################################XXX XXXXX, Oct 14 2019 ################################## 

library(grid)
library(scales)
library(TeachingDemos)
library(ape)
library(phangorn)
library(Hmisc,)
library(brms)
library(lme4)
library(lmerTest)
library(MuMIn)
library(bayestestR)
library(tidyverse)
library(egg)
library(broom.mixed)
library(tidybayes)

# ----  Adjust all file paths to read in data and save plots
gdrive_path <- file.path('/Users/jgradym/Documents/GitHub/Foodweb_thermal_asymmetries/Data')


################################## Plotting ###########################33

# Colors
invert_col <- "deepskyblue2"
fish_col <- "navy"
mamm_col <- "red"
bird_col <- "darkred"
ecto_col <- "dodgerblue"
endo_col <- "red2"

theme_single<- theme(panel.grid = element_blank(), 
                          aspect.ratio = .70,
                          axis.text = element_text(size = 18, color = "black"), 
                          axis.ticks.length=unit(0.2,"cm"),
                          axis.title = element_text(size = 18),
                          axis.title.y.right = element_text(margin = margin(r = 0)),
                          axis.title.y.left = element_text(margin = margin(l = 0)),
                          axis.title.x = element_text(margin = margin(t = 0)),
                          axis.title.x.top = element_text(margin = margin(b = 0)),
                          plot.title = element_text(size = 18, hjust = 10),
                          panel.border = element_rect(colour = "black", fill=NA, size=1),
                          panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),
                          panel.background = element_rect(fill = "transparent",colour = NA),
                          plot.background = element_rect(fill = "transparent",colour = NA),
                          legend.background = element_rect(color = NA),
                          legend.key=element_blank(), 
                          text = element_text(family = 'Helvetica')) 

theme_no_y <- theme(axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank())

theme_no_x <- theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())

theme_ticks <- theme(axis.ticks.length=unit(-0.25, "cm"),
                     axis.text.x = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")), 
                     axis.text.x.top = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")), 
                     axis.text.y.right = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")),
                     axis.text.y.left = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")))


######################################################
############## Mortality Data first
######################################################


###################################  Add endo temp data ################################## 

# Create 1/kT column; 
mortality0 <- read_csv(file.path(gdrive_path,'McCoy_mortality_updated.csv')) 
mortality0$one_kT <- 1/(8.617e-5*(mortality0$Temp_C + 273.15)) # Inverse Temperature

mortality1 <- mortality0 %>%
  rename(species = Species) %>%
  arrange(Group, species) %>%
  mutate(log_mort = log(Mortality_yr), log_mass = log(dry_mass_g),
         Genus = word(species, sep = "_"), species2 = species)
mortality1
 
##########################################################################
#--------------- Analyze Groups ----------

#------------------- Ectotherms -----------------
model_ecto_mort <- read_rds("/Users/jgradym/Google Drive/Gibert Paper/model_ecto_mort.rds")

# Data
mort_ecto <- mortality1 %>%
  filter(Group == "Fish" | Group == "Invertebrate") 

# Bayesian framework with species as random effect
priors_ecto <- c(set_prior("normal(-0.65, 1)", class = "b", coef = "one_kT"), #slope = 0 for endotherms
             set_prior("normal(-0.25, 1)", class = "b", coef = "log_mass"),
             set_prior("normal(0, 1)", class = "b")
             
)

# Model, takes a minute
system.time(model_ecto_mort <- brm(
  log_mort ~ log_mass + one_kT + (log_mass + one_kT|Group/species), 
  data = mort_ecto, 
  family = gaussian(), 
  prior = priors_ecto,
  cores = parallel::detectCores() -1) 
) 
saveRDS(model_ecto_mort, "/Users/jgradym/Google Drive/Gibert Paper/model_ecto_mort.rds")
#Results
print(summary(model_ecto_mort), digits = 5) 
coef(model_ecto_mort)$Group
bayes_R2(model_ecto_mort) #full model
bayes_R2(model_ecto_mort, re.form = NA) # exclude random effects
pd_ecto_mort <- p_direction(model_ecto_mort) # p direction
pd_to_p(pd_ecto_mort$pd[2], "two-sided") # bayesian equivalent to p value for mass scaling
pd_to_p(pd_ecto_mort$pd[3], "two-sided") # bayesian equivalent to p value for temp scaling



c_mort <- coef(model_ecto_mort)
c_mort$Group



system.time(model_attack<- brm(
  log_attack ~  log_mass + one_kT +
    (log_mass|predator.ana.group/predator.species + one_kT|predator.ana.group/predator.species),
  data = attack, 
  family = gaussian(), 
  prior = priors3,
  cores = parallel::detectCores() -1) 
) 
c_attack <- coef(model_attack)
c_attack$predator.ana.group
#--------------Mammals: comprehensive phylogeny available-----------------
# Data
mort_mammal <- mortality1 %>%
  filter(Group == "Mammal") 

#-- Add phylogeny ---
mam_tree_mcc <- read.tree('/Users/jgradym/Google Drive/Phylo_trees/mammal_mcc.tree')

#------------ can skip below
# Read in 100 phylogenetic trees - download at https://data.vertlife.org/  mammaltree > Completed_5911sp_topoCons_FBDasZhouEtAl.zip

mamm_tree_files <- list.files('/Users/jgradym/Google Drive/Phylo_trees/Upham_Phylos/Mammals/Completed_5911sp_topoCons_FBDasZhouEtAl', full.names = T)[1:100]
mam_trees <- lapply(mamm_tree_files, FUN = ape::read.tree)
class(mam_trees) <- "multiPhylo"

# get maximum credible clade
mam_tree_mcc0 <- maxCladeCred(mam_trees) 

#prune down tree to our species
mam_spmatch <- match(mam_tree_mcc0$tip.label, #get matched species
                     mort_mammal_0$species)
mam_tree_mcc <- drop.tip(mam_tree_mcc0, 
                         mam_tree_mcc0$tip.label[is.na(mam_spmatch)]) #prune unused spp
#write.tree(mam_tree_mcc, "/Users/jgradym/Google Drive/Phylo_trees/mammal_mcc.tree")


mort_mammal <- mort_mammal_0[mort_mammal_0$species %in% mam_tree_mcc$tip.label,]
#------------ begin again
#----------------- bayesian model for mammals using brms, species as random effect -----------------

model_mort_mammal <- read_rds('/Users/jgradym/Google Drive/Gibert Paper/Mammal_fit.rds')
# variance and covariance from phylogeny
A <- ape::vcv.phylo(mam_tree_mcc)
# Priors
get_prior(formula = log_mort ~ log_mass + one_kT +
            (1|gr(species, cov = A)) + (1|species2), data = mort_mammal)

priors_endo1 <- c(set_prior("normal(0, 1)", class = "b", coef = "one_kT"), #slope = 0 for endotherms
             set_prior("normal(-0.25, 1)", class = "b", coef = "log_mass"),
             set_prior("normal(0, 1)", class = "b")
             
)


system.time(model_mort_mammal <- brm(
    log_mort ~ log_mass + one_kT + (1|gr(species, cov = A)) + (1|species2), 
    data = mort_mammal, 
    family = gaussian(), 
    data2 = list(A = A),
    prior = priors_endo1,
    cores = parallel::detectCores() -1) 
  )


# Model Results
print(summary(model_mort_mammal), digits = 5) 
bayes_R2(model_mort_mammal, re.form = NA) # exclude random effects
bayes_R2(model_mort_mammal) # full model
pd_mamm <- p_direction(model_mort_mammal) # p direction
pd_to_p(pd_mamm$pd[2], "two-sided") # bayesian equivalent to p value for mass scaling
pd_to_p(pd_mamm$pd[3], "two-sided") # bayesian equivalent to p value for temp scaling

pdf('~/Desktop/mammal_bayes_plot.pdf')
plot(model_mort_mammal, N = 6)
dev.off()
saveRDS(model_mort_mammal, file = "/Users/jgradym/Google Drive/Gibert Paper/Mammal_fit.rds")
check <- readRDS("/Users/jgradym/Google Drive/Gibert Paper/Mammal_fit.rds")

mammal_fitted <- fitted(model_mort_mammal, newdata = NULL, re_formula = NA, summary = T)

# ----------- Birds ------------------

mort_bird <- mortality1 %>%
  filter(Group == "Bird") #%>%

bird_tree_mcc  <- read.tree('/Users/jgradym/Google Drive/Phylo_trees/bird_tree_mcc.tree')

model_mort_bird <- read_rds("/Users/jgradym/Google Drive/Gibert Paper/bird_fit.rds")

#-------can skip below 
# Read in 100 phylogenetic trees - download at https://data.vertlife.org/  birdtree > PatchClade > Stage 1 > EricsonStage1_0001_1000.zip
bird_tree <- read.tree('/Users/jgradym/Google Drive/Phylo_trees/Upham_Phylos/Birds/mnt/data/projects/birdphylo/Tree_sets/Stage1_full_data/CombinedTrees/EricsonStage1Full_1.tre')[1:100]

# get maximum credible tree
bird_tree_mcc0 <- maxCladeCred(bird_tree) 

#prune down tree to species in dataset
bird_spmatch <- match(bird_tree_mcc0$tip.label, #get matched species
                     mort_bird0$species)

bird_tree_mcc <- drop.tip(bird_tree_mcc0 , bird_tree_mcc0$tip.label[is.na(bird_spmatch)]) #prune unused spp
#write.tree(bird_tree_mcc , "/Users/jgradym/Google Drive/Phylo_trees/bird_tree_mcc.tree")

# read in model results
model_mort_bird <- read_rds('/Users/jgradym/Google Drive/Gibert Paper/bird_fit.rds')
# generate model results (can skip)-----------------------

# variance and covariance from phylogeny
A <- ape::vcv.phylo(bird_tree_mcc)

# Run model
system.time(model_mort_bird <- brm(
    log_mort ~ log_mass + one_kT + (1|gr(species, cov = A)) + (1|species2), 
    data = mort_bird, 
    family = gaussian(), 
    data2 = list(A = A),
    prior = priors_endo1,
    cores = parallel::detectCores() -1) 
  ) 
saveRDS(model_mort_bird, file = "/Users/jgradym/Google Drive/Gibert Paper/bird_fit.rds")
#----------------------------------------
# Model Results
print(summary(model_mort_bird), digits = 6) 
bayes_R2(model_mort_bird) 
bayes_R2(model_mort_bird, re.form = NA)  # fixed effects only

pd_bird <- p_direction(model_mort_bird) # p direction
pd_to_p(pd_bird$pd[2], "two-sided") # Bayesian equivalent to p value for mass scaling
pd_to_p(pd_bird$pd[3], "two-sided") # Bayesian equivalent to p value for temp scaling

pdf('~/Desktop/bird_bayes_plot.pdf')
plot(model_mort_bird, N = 6)
dev.off()
# Get intercept of mass-corrected scaling for plotting
mort_bird <- mort_bird  %>%
  mutate(log_mort_corr = log(exp(mort_bird$log_mort)/
                               exp(mort_bird$log_mass)^ -0.1642))
bird_fitted <- fitted(model_mort_bird, newdata = NULL, re_formula = NA, summary = T)
head(bird_fitted)
#------- Add mass-corrected mortality for plotting ----

mortality2 <- mortality1 %>% 
  mutate(mass_corr_mortality = NA)

bird <- mortality2[mortality2$Group == "Bird",]
mortality2$mass_corr_mortality[mortality2$Group == "Bird"] <- bird$Mortality_yr/bird$dry_mass_g^-0.169266

mammal <- mortality2[mortality2$Group == "Mammal",]
mortality2$mass_corr_mortality[mortality2$Group == "Mammal"] <- mammal$Mortality_yr/mammal$dry_mass_g^-0.178

invert <- mortality2[mortality2$Group == "Invertebrate",]
mortality2$mass_corr_mortality[mortality2$Group == "Invertebrate"] <- invert$Mortality_yr/invert$dry_mass_g^-0.249

fish <- mortality2[mortality2$Group == "Fish",]
mortality2$mass_corr_mortality[mortality2$Group == "Fish"] <- fish$Mortality_yr/fish$dry_mass_g^-0.291


#get credible bands

#-------Adjust and columns for genus level analysis ------------

# Get number measurements per species or per genus, if you want to restrict by individuals per species or genus:, 
mortality3 <- mortality2 %>% #firs
  mutate(Genus = word(species,1, sep = "_")) %>%
  group_by(Genus) %>% 
  add_count(Genus, name = "n_genus") %>%# individuals per genus
  mutate(temp_range_genus = max(Temp_C) - min(Temp_C)) %>%
  ungroup()

#Add a temperature range for individuals in a species (only relevant for ectos)
mortality4 <- mortality3 %>%
  group_by(species) %>% 
  add_count(species, name = "n_species") %>%
  mutate(temp_range_spp = max(Temp_C) - min(Temp_C)) %>%
  ungroup()

# Add thermy
mortality5 <- mortality4 %>%
  filter(Group == "Bird" | Group == "Mammal"|  Group == "Invertebrate" | Group == "Fish") %>%
  mutate(Thermy = if_else(Group == "Bird" | Group == "Mammal", "Endotherm", "Ectotherm"))

# Final dataset
mortality <- mortality5 # final dataset
mortality$log_mort_corr <- log(mortality$mass_corr_mortality)
#--------------- with mass-corrected mortality ---------------
bird_mort <- mortality[mortality$Group == "Bird",]
mammal_mort <- mortality[mortality$Group == "Mammal",]
ecto_mort <- mortality[mortality$Group == "Invertebrate" | mortality$Group == "Fish",]

# Birds Run model
priors_endo2 <- c(set_prior("normal(0, 1)", class = "b", coef = "one_kT"), #slope = 0 for endotherms
                  set_prior("normal(0, 1)", class = "b")
                  
)

                  
A <- ape::vcv.phylo(bird_tree_mcc)
system.time(model_mort_bird_corr <- brm(
  log_mort_corr ~ one_kT + (1|gr(species, cov = A)) + (1|species2), 
  data = bird_mort, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors_endo2,
  cores = parallel::detectCores() -1) 
) 
saveRDS(model_mort_bird_corr, file = "/Users/jgradym/Google Drive/Gibert Paper/model_mort_bird_corr.rds")
plot(model_mort_bird_corr)
print(summary(model_mort_bird_corr), digits = 5)
get_prior(formula = log_mort_corr ~  one_kT + (one_kT|Group/species), data = ecto_mort)


priors_ecto2 <- c(set_prior("normal(0, 1)", class = "b"),
                  set_prior("normal(-0.65, 1)", class = "b", coef = "one_kT") )


# Model, takes a minute
system.time(model_ecto_mort_corr <- brm(
  log_mort_corr  ~ one_kT + (one_kT|Group/species), 
  data = ecto_mort, 
  family = gaussian(), 
  prior = priors_ecto2,
  cores = parallel::detectCores() -1) 
) 
saveRDS(model_ecto_mort_corr, "/Users/jgradym/Google Drive/Gibert Paper/model_ecto_mort_corr.rds")
print(summary(model_ecto_mort_corr), digits = 5)
##########################################
 #### Panel 2a
##########################################

bird_lab  <- expression(paste('y = -0.022x + 0.20')) #from model summary
mamm_lab  <- expression(paste('y = 0.0057x - 0.48'))

# for manual plotting, results from Bayesian mixed effect analysis with phylogeny
bird_fit <- data.frame(x = c(min(mort_bird$one_kT),  max(mort_bird$one_kT)), 
                     y = c(exp((0.203 + (-0.0220 * min(mort_bird$one_kT)))), 
                     exp(0.203 + (-0.0220 * max(mort_bird$one_kT)))),
                     Group = c("Bird", "Bird"))

mammal_fit <- data.frame(x = c(min(mort_mammal$one_kT),  max(mort_mammal$one_kT)), 
                         y = c( exp((-0.478 + (0.00572 * min(mort_mammal$one_kT)))), 
                         exp((-0.478 + (0.00572 *  max(mort_mammal$one_kT))))),
                         Group = c("Mammal", "Mammal"))

endo_fit <- data.frame(rbind(bird_fit, mammal_fit))



#------------- Plot Fig 2A  -----------------

# add fitted credible bands

# birds
bird_fitted0 <- as_tibble(fitted(model_mort_bird_corr, scale = "response", newdata = NULL, re_formula = NA, summary = T))
bird_fitted <- as_tibble(cbind(bird_mort_simp, bird_fitted0))

ecto_fitted0 <- as_tibble(fitted(model_ecto_mort_corr,  newdata = NULL, re_formula = NA, summary = T))
ecto_simp <- ecto_mort %>% select(one_kT)
ecto_fitted <- as_tibble(cbind(ecto_simp , ecto_fitted0 ))

# mammals

#ectos




mortal_plot_endo <- ggplot(mortality %>% filter(Group == "Mammal" | Group == "Bird") %>%
                             arrange(Group),
                           aes(x = one_kT, y = mass_corr_mortality, shape = Group)) + 
  theme_single + theme(legend.position = "none") +  theme_ticks +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{-alpha},')')), 
                     trans = "log", breaks = c(0.1, 1, 10), labels = c("0.1", "1", "10"), 
                     limits = c(0.02, 40)) +  
  scale_x_reverse(limits = c(45.1, 38), name = NULL, breaks = c(46,44, 42,40, 38),
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, breaks = c(-20, -10, 0, 10, 20, 30),
                                      name = "Ambient Temperature (ºC)")) +
  geom_point(size = 2.5, shape = 21, color = "black",  stroke = .5, aes(fill = Group)) + 
  scale_fill_manual(values =c("Mammal" = mamm_col, "Bird" = "darkred")) +
  scale_color_manual(values =c("Mammal" = mamm_col, "Bird" = "brown4")) +
  geom_line(data = endo_fit, aes(x = x, y = y, color = Group), size = 1) + 
  annotate("text", x = 45, size = 6, y = 17, hjust = 0,  fontface = "bold",label = "Bird", color = bird_col) +
  annotate("text", x = 45, size = 6, y = 35,  hjust = 0, fontface = "bold",label = "Mammal", color = mamm_col) +
  annotate("text", x = 41, size = 6, y = .045,  hjust = 0, label = mamm_lab, color = "red1")+
  annotate("text", x = 41, size = 6, y = .025,  hjust = 0, label = bird_lab, color = "brown4") +
  #stat_lineribbon(data = bird_fits, aes(y = .value, x = one_kT)) +
  geom_ribbon(data = bird_fitted, # %>% filter(fit == 'light per area', dbh < 156), 
              aes(x = one_kT, y = exp(Estimate), ymin = exp(Q2.5), ymax = exp(Q97.5)), alpha = 0.4) 

grid.newpage()
grid.draw(mortal_plot_endo)
model_mort_bird_corr

mortality %>%
  filter(Group == "Fish" | Group == "Invertebrate") %>%
  group_by(Group) %>%
  #data_grid(hp = seq_range(hp, n = 51)) %>%
  add_fitted_draws(model_ecto_mort_corr, re_formula = NA) %>%
  #add_predicted_draws(model_ecto_mort_corr, re_formula = NA) %>%
  ggplot(aes(x = one_kT, y = log_mort_corr)) +
  geom_point(data = mortality %>% filter(Group == "Fish" | Group == "Invertebrate"), 
             aes(x = one_kT, y = log_mort_corr, color = Group)) +
  theme_single +
  #geom_ribbon(data = ecto_fitted, 
   #           aes(x = one_kT, y = Estimate, ymin = Q2.5, ymax = Q97.5), alpha = 0.4) 
  geom_smooth(data = mortality %>%
                filter(Group == "Fish" | Group == "Invertebrate"),
              aes(y = log_mort_corr, x = one_kT), method = "lm", se = T) +
  geom_smooth(data = ecto_fitted ,
              aes(y = Estimate), color = "red") +
  #geom_ribbon(data = bird_fitted, 
   #           aes(x = one_kT, y = Estimate, ymin = Q2.5, ymax = Q97.5), alpha = 0.4, fill = "red") +
  #stat_lineribbon( alpha = .1) +
  stat_lineribbon(aes(y = Estimate), alpha = .1)  +
  
  
  

ggplot() + 
  geom_point(data = mortality %>% filter(Group == "Bird") %>%
               arrange(Group),
             aes(x = one_kT, y = mass_corr_mortality, color = Group, fill = Group, shape = Group)) +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{-alpha},')')), 
                     trans = "log", breaks = c(0.1, 1, 10), labels = c("0.1", "1", "10"), 
                     limits = c(0.02, 40)) +  
  scale_x_reverse(limits = c(45.1, 38), name = NULL, breaks = c(46,44, 42,40, 38),
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, breaks = c(-20, -10, 0, 10, 20, 30),
                                      name = "Ambient Temperature (ºC)")) +
  theme_single + theme(legend.position = "none") +  theme_ticks +
  scale_fill_manual(values =c("Mammal" = mamm_col, "Bird" = "darkred")) +
  scale_color_manual(values =c("Mammal" = mamm_col, "Bird" = "brown4")) +
  add_fitted_draws(model_mort_bird_corr) %>%
  stat_lineribbon(aes(y = .value)) +
  geom_ribbon(data = bird_fitted, 
              aes(x = one_kT, y = exp(Estimate), ymin = exp(Q2.5), ymax = exp(Q97.5)), alpha = 0.4) +
  geom_ribbon(data = ecto_fitted, 
              aes(x = one_kT, y = exp(Estimate), ymin = exp(Q2.5), ymax = exp(Q97.5)), alpha = 0.4) 
  

ggsave(mortal_plot_endo, height = 4.6, width = 6.7, filename = '~/Desktop/Fig_2/Fig2A.pdf')

# Save plot
ggsave(mortal_plot_endo, height = 4.6, width = 6.7, filename = file.path(gdrive_path, "NAME YOUR FILE.pdf"))


###########################################
#### Panel 2b - species as random effect
##########################################

# Prepare plot
invert_lab  <- expression(paste('y = -0.55x + 22.5'))
fish_lab  <- expression(paste('y = -0.52x + 21.8'))

# from models summary
invert_fit <- data.frame(x = c(min(mort_invert$one_kT),  max(mort_invert$one_kT)), 
                       y = c(exp((22.47 + (-0.548 * min(mort_invert$one_kT)))), 
                             exp((22.47 + (-0.548 *  max(mort_invert$one_kT))))),
                       Group = c("Invertebrate", "Invertebrate"))


fish_fit <- data.frame(x = c(min(mort_fish$one_kT), max(mort_fish$one_kT)), 
                         y = c( exp((21.80 + (-0.523 * min(mort_fish$one_kT)))), 
                                exp((21.80 + (-0.523 *  max(mort_fish$one_kT))))),
                         Group = c("Fish", "Fish"))

ecto_fit <- data.frame(rbind(invert_fit, fish_fit))
ecto_fit 

#--------------  Plot Fig 2b ------------
mortal_plot_ecto <- ggplot(mortality %>% 
                             filter(Group == "Fish" | Group == "Invertebrate") %>%
                             arrange(Group),
                           aes(x = one_kT, y = mass_corr_mortality,fill = Group)) + 
  scale_y_continuous(trans = "log10",  name = expression(paste('Mortality (yr'^-1,' g'^{-alpha},')')), 
                    breaks = c(0.1,  1, 10), labels = c("0.1", "1",  "10"), 
                    limits = c(0.03, 30), 
                    position = "right")+ 
  scale_x_reverse(name = expression("1/kT"), limits = c(43, 38.2), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, 
                                      name = "Ambient Temperature (ºC)")) +
  theme(legend.position = "none") +theme_single + 
  geom_point(size = 2.5, shape = 21, color = "black", 
             stroke = .5, aes(fill = Group)) +
  scale_fill_manual(values =c("Fish" = fish_col, "Invertebrate" =  invert_col)) +
  scale_color_manual(values =c("Fish" = fish_col, "Invertebrate" =  invert_col)) +
  geom_line(data = ecto_fit, aes(x = x, y = y, color = Group), size = 1) + 
   annotate("text", x = 43, size = 6, y = 14, hjust = 0, 
                          fontface = "bold",label = "Fish", color = fish_col) +
  annotate("text", x = 43, size = 6, y = 25, hjust = 0, 
           fontface = "bold",label = "Invertebrate", color = invert_col) +
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text.x = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")), 
        axis.text.x.top = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")), 
        axis.text.y.right = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")),
        axis.text.y.left = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")))
grid.newpage()
grid.draw(mortal_plot_ecto) 


# Save plot
ggsave(mortal_plot_ecto, height = 4.6, width = 6.7, filename = '~/Desktop/Fig_2/Fig2b.pdf')

####################################################################################
#### Panel 2c - genus level regressions, min 5 individual per genus, temp range 5 C 
####################################################################################

# Thermy
endo_genus <- mortality %>%
  filter(Thermy == "Endotherm") %>%
  pull(Genus) %>%
  unique()

# Linear regressions for genus level plotting 
genus_mort_lm <- mortality %>%
  filter(n_genus >= 5, temp_range_genus >= 5) %>%
  nest(-Genus) %>% 
  mutate(
    fit = map(data, ~ lm(log(mass_corr_mortality) ~ one_kT, data = .x)),
    tidied = purrr::map(fit, tidy, conf.int = T),
    pred = map(fit, broom::augment)
  ) 

genus_mort_lm <- unnest(genus_mort_lm, pred) %>%
  rename(log_mass_corr_mortality = `log(mass_corr_mortality)`)

# add column info back for plotting
genus_mort_lm$Thermy<-if_else(genus_mort_lm$Genus %in% endo_genus, "Endotherm", "Ectotherm" )
genus_mort_group <- map(genus_mort_lm $data, ~ .x[["Group"]][[1]])
genus_mort_lm$Group <- unlist(genus_mort_group) 
unnest(genus_mort_lm, tidied) %>% filter(term == "one_kT")

#output results
genus_mort_lm2 <- mortality %>%
  filter(n_genus >= 5, temp_range_genus >= 5) %>%
  nest(-Genus) %>% 
  mutate(
    fit = map(data, ~ lm(log(mass_corr_mortality) ~ one_kT, data = .x)),
    tidied = purrr::map(fit, tidy, conf.int = T),
  ) %>%
  unnest(tidied)
genus_mort_lm2$Thermy<-if_else(genus_mort_lm2$Genus %in% endo_genus, "Endotherm", "Ectotherm" )
genus_mort_group <- map(genus_mort_lm2$data, ~ .x[["Group"]][[1]])
genus_mort_lm2$Group <- unlist(genus_mort_group) 

genus_mort_lm2a <- genus_mort_lm2 %>%
  filter(term == "one_kT") %>%
  dplyr::select(-data, -fit) %>%
  dplyr::select(Thermy, Group, Genus, everything())
genus_mort_lm2a 
write_csv(genus_mort_lm2a, "~/Desktop/Fig_2/genus_regressions.csv")


# get credible bands 
bird_simp <- bird %>% 
  dplyr::select(one_kT) %>%
  mutate(log_mass = mean(bird$log_mass))
bird_fitted <- fitted(model_mort_bird, newdata = bird_simp, re_formula = NA, summary = T, resp = log_mort)
length(bird_fitted)
head(bird_fitted)
mammal_simp <- mammal %>% 
  dplyr::select(one_kT) %>%
  mutate(log_mass = mean(mammal$log_mass)) #use average mass
mammal_fitted <- fitted(model_mort_mammal, newdata = mammal_simp, re_formula = NA, summary = T, resp = log_mort)
length(mammal_fitted)
head(mammal_fitted)
#--- plot Fig 2C ------

mortal_genus_regr <- ggplot(data = mortality %>% 
                              filter(n_genus >= 5, temp_range_genus >= 5),
                            aes(y = mass_corr_mortality, x = one_kT)) +
    geom_smooth(method = "lm", se =- F, size = 0.3, aes(group = Genus, color = Group)) + 
  geom_line(data = ecto_fit, aes(x = x, y = y, color = Group), size = 1.5) +
  geom_line(data = endo_fit, aes(x = x, y = y, color = Group), size = 1.5) + 
  scale_color_manual(values =c("Invertebrate" = invert_col, "Fish" = fish_col, 
                                 "Mammal" = mamm_col, "Bird" = bird_col)) +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{-alpha},')')), 
                     trans = "log", breaks = c(0.3, 1, 3,10), labels = c("0.3", "1", "3", "10"), 
                     limits = c(.15, 10)) +  
  scale_x_reverse(name =NULL, limits = c(45.1, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15,  name = NULL)) + 
  theme_single + theme(legend.position = "none") +theme_ticks 
grid.newpage()
grid.draw(mortal_genus_regr)

# Save plot
ggsave(mortal_genus_regr, height = 4.3, width = 6.3, filename = '~/Desktop/Fig2C.pdf')


##########################################
#### Panel 2d, Mortality violin plot
 #first get regressions
##########################################

# get parameter coefficients
genus_mort_lm2 <- mortality %>%
  filter(n_genus >= 5, temp_range_genus >= 5) %>%
  nest(-Genus) %>% 
  mutate(
    fit = map(data, ~ lm(log(mass_corr_mortality) ~ one_kT, data = .x)),
    tidied = purrr::map(fit, tidy, conf.int = T)
  ) %>% 
  unnest(tidied) 
genus_mort_lm2$Thermy<-if_else(genus_mort_lm2$Genus %in% endo_genus, "Endotherm", "Ectotherm" )
genus_mort_lm2

# Plot - remove y limits to see outliers

mort_vio_Ea <- ggplot(genus_mort_lm2 %>% filter(term == "one_kT"), 
                      aes(x = Thermy, y = -1*estimate))+
  geom_violin(size = 1, aes(color = Thermy)) + 
  geom_jitter(aes(fill = Thermy),shape = 21, size = 3, color = "black", width = 0.08, stroke = .5 ) +
  scale_y_continuous(limits = c(-0.6, 1.1),
                     position = "right",
                     name = expression(paste('Thermal Sensitivity (E'[a],')')))+ 
  scale_x_discrete(name = NULL) +
  scale_color_manual(values =c("Ectotherm" = ecto_col, "Endotherm" = endo_col)) +
  scale_fill_manual(values =c("Ectotherm" = ecto_col, "Endotherm" = endo_col)) +
  theme_single + theme(legend.position = "none") +  
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + theme_ticks +
  ggtitle("Mortality")  
 
mort_vio_Ea
ggsave(mort_vio_Ea , height = 4.3, width = 6.3, filename = '~/Desktop/Fig2D.pdf')


####################################################
##v  Fig 2E - Mortality & Attack regressions (ectos)
####################################################
##---------------------------------
## Plotting Fig 2e, Part I (mortality)


mortal_scaling <- ggplot(mortality %>% 
                                filter(Group == "Fish" | Group == "Invertebrate")  %>%
                                filter(n_genus >= 5) %>%
                                filter(temp_range_genus >= 5),
                              aes(x = one_kT, y = mass_corr_mortality)) + 
  geom_smooth(method = 'lm', size = 1,  alpha = 0, aes(group = Genus, color = Group)) +
  scale_color_manual(values =c("Invertebrate" = invert_col, "Fish" = fish_col, 
                               "Mammal" = mamm_col, "Bird" = bird_col)) +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{-alpha},')')), position = "right",
                     trans = "log10", breaks = c(0.1, 1, 10), limits = c(0.1, 10), 
                     labels = c("0.1", "1", "10")) +  
  scale_x_reverse(name =NULL, limits = c(42.5, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, name = NULL)) + 
  theme_single + theme_ticks +theme(legend.position = "none")
grid.draw(mortal_scaling )
# Note: range of y limits is equal to attack rates below - can observe difference in spread of intercepts

# Save plot
ggsave(mortal_scaling, height = 4.7, width = 6.7, filename = '~/Desktop/mortal_scaling.pdf')

ggsave(mortal_regress, height = 4.7, width = 6.7, filename = file.path(gdrive_path,'NAME YOUR FILE.pdf'))

##---------------------------------
## PART 2 (attack)

# Read in data
attack <- read_csv(file.path(gdrive_path, 'Lietal_oikos_2017_data.csv'))
attack$one_kT <- 1/(8.617e-5*(attack$temperature.degree.celcius +273.15))

# Add columns
attack2 <- attack %>%
  mutate(Pred_Genus = word(predator.species)) %>%
  add_count(predator.species, name = "n_species") %>%# exclude taxa with few members
  add_count(Pred_Genus, name = "n_genus")#adjust n by changing filter(n_species/n_genus = )

attack3  <- attack2  %>%
  group_by(predator.species ) %>% 
  mutate(temp_range_spp = max(temperature.degree.celcius) - min(temperature.degree.celcius))

attack4  <- attack3  %>%
  group_by(Pred_Genus) %>% 
  mutate(temp_range_genus = max(temperature.degree.celcius) - min(temperature.degree.celcius))

attack <- attack4
attack$log_attack <- log(attack$attack.rate)
attack$log_mass <- log(attack$predator.mass.mg)

# Mixed model test - Group and species as random effect  (species nested in group)
ecto_attack_lmer <-  lmer(data = attack, log_attack  ~ one_kT + log_mass +
                            (one_kT|predator.ana.group:predator.species) +
                            (log_mass|predator.ana.group:predator.species))
str(coef(ecto_attack_lmer))
invert_mass_scaling <- coef(ecto_attack_lmer)
summary(ecto_attack_lmer)
confint.merMod(ecto_attack_lmer, method = "Wald")
r.squaredGLMM(ecto_attack_lmer)

# Mixed model for inverts separately
# Data


invert_attack <- attack %>% filter(predator.ana.group == "invertebrate")
invert_attack_lmer <- lmer(data = invert_attack, log_attack ~ one_kT + log(predator.mass.mg) +
                             (log(predator.mass.mg)|Pred_Genus) +
                             (one_kT|Pred_Genus) + (1|predator.species) + (1|Pred_Genus:predator.species))
summary(invert_attack_lmer)  #0.56031 mass scaling slope
confint.merMod(invert_attack_lmer, method = "Wald")
r.squaredGLMM(invert_attack_lmer)
invert_attack$attack_mass_corr <- invert_attack$attack.rate/invert_attack$predator.mass.mg^0.538

# Bayesian framework with species as random effect
priors2 <- c(set_prior("normal(0, 1)", class = "sd", coef = "Intercept", group = "predator.species"),
             set_prior("normal(-0.65, 1)", class = "b", coef = "one_kT"), #slope = 0 for endotherms
             set_prior("normal(1, 1)", class = "b", coef = "log_mass"),
             set_prior("normal(0, 1)", class = "b")
             
)
get_prior(formula = log_attack ~  log_mass + one_kT + (predator.ana.group|log_mass/predator.species) +
            (one_kT|predator.ana.group/predator.species), data = attack)
lmer_attack <- lmer(log_attack ~  log_mass + one_kT + (log_mass|predator.ana.group/predator.species) +
           (one_kT|predator.ana.group/predator.species), data = attack)
summary(lmer_attack)
str(coef(lmer_attack))
coef(lmer_attack)$predator.ana.group
#(Intercept)  log_mass     one_kT
#invertebrate   -3.287899 0.5075178 -0.3610663
#vertebrate     -3.147275 0.7868092 -0.3628188

lmer_attack_invert <- lmer(log_attack ~  log_mass + one_kT + (log_mass|predator.species) +
                      (one_kT|predator.species), data = invert_attack)
summary(lmer_attack_invert)
lmer_attack_vert <- lmer(log_attack ~  log_mass + one_kT + (log_mass|predator.species) +
                             (one_kT|predator.species), data = vert_attack)
summary(lmer_attack_vert)

priors3 <- c(set_prior("normal(-0.65, 1)", class = "b", coef = "one_kT"), #slope = 0 for endotherms
             set_prior("normal(1, 1)", class = "b", coef = "log_mass"),
             set_prior("normal(0, 1)", class = "b")
            )
lm <- lm(log_attack ~  log_mass + one_kT * predator.ana.group, data = attack)
summary(lm)
system.time(model_attack<- brm(
  log_attack ~  log_mass + one_kT + (predator.ana.group|log_mass/predator.species) +
    (one_kT|predator.ana.group/predator.species),
  data = attack, 
  family = gaussian(), 
  prior = priors3,
  cores = parallel::detectCores() -1) 
) 

# Results
print(summary(model_attack, waic = T), digits = 5) 
str(coef(model_attack_invert))
coef(model_attack_invert)$Estimate
head(coef(model_attack_invert))
print(summary(model_attack_invert), digits = 5) 
print(summary(model_attack_vert), digits = 5) 
bayes_R2(model_attack_invert  ) 
bayes_R2(model_attack_invert , re.form = NA) # exclude random effects
pd_invert_attack <- p_direction(model_attack_invert ) # p direction
pd_to_p(pd_invert_attack $pd[2], "two-sided") # bayesian equivalent to p value for mass scaling
pd_to_p(pd_invert_attack $pd[3], "two-sided") # bayesian equivalent to p value for temp scaling



vert_attack <- attack %>% filter(predator.ana.group == "vertebrate")

system.time(model_attack_vert <- brm(
  log_attack ~  log_mass + one_kT + (1|predator.species), 
  data = vert_attack, 
  family = gaussian(), 
  prior = priors2,
  cores = parallel::detectCores() -1) 
) 

# Results
print(summary(model_attack_vert ), digits = 5) 
bayes_R2(model_attack_vert  ) 
bayes_R2(model_attack_vert , re.form = NA) # exclude random effects
pd_vert_attack <- p_direction(model_attack_vert ) # p direction
pd_to_p(pd_vert_attack $pd[2], "two-sided") # bayesian equivalent to p value for mass scaling
pd_to_p(pd_vert_attack $pd[3], "two-sided") # bayesian equivalent to p value for temp scaling


invert_attack <- attack %>% filter(predator.ana.group == "invertebrate")
invert_attack_lmer <- lmer(data = invert_attack, log_attack ~ one_kT + log(predator.mass.mg) +
                             (log(predator.mass.mg)|Pred_Genus) +
                             (one_kT|Pred_Genus) + (1|predator.species)) 
invert_attack_lmer_sum <- summary(invert_attack_lmer) 
coef_invert <- invert_attack_lmer_sum$coefficients
invert_slope <- coef_invert[3]
summary(invert_attack_lmer)  #0.56031 mass scaling slope
confint.merMod(invert_attack_lmer, method = "Wald")
r.squaredGLMM(invert_attack_lmer)
invert_attack$attack_mass_corr <- invert_attack$attack.rate/invert_attack$predator.mass.mg^invert_slope




# Mixed model for vertebrates separately
vert_attack <- attack %>% filter(predator.ana.group == "vertebrate")
vert_attack_lmer <- lmer(data = vert_attack, log_attack ~ one_kT + log(predator.mass.mg) +
                           (1 + one_kT|Pred_Genus) + (1|predator.species))

vert_attack_lmer_sum <- summary(vert_attack_lmer) 
coef <- vert_attack_lmer_sum$coefficients
vert_slope <- coef[3]
confint.merMod(vert_attack_lmer, method = "Wald")
r.squaredGLMM(vert_attack_lmer)
str(vert_attack_lmer)
# add mass-corrected attack rate as column
vert_attack$attack_mass_corr <- vert_attack$attack.rate/vert_attack$predator.mass.mg^vert_slope

# combine vert and invert attack
attack <- as_tibble(rbind(vert_attack, invert_attack))

# get attack regressions with temperature
attack_genus_lm0 <- attack %>%
  filter(n_genus >= 5, temp_range_genus >= 5) %>%
  nest(-Pred_Genus) %>% 
  mutate(
    fit = map(data, ~ lm(log(attack_mass_corr) ~ one_kT, data = .x)),
    tidied = purrr::map(fit, tidy, conf.int = T)
  ) %>% unnest(tidied) %>%
  filter(term == "one_kT") %>% 
  select(-data, -fit)


attack_genus_lm0 
attack_genus_lm0
attack_genus_lm <- attack %>%
  filter(n_genus >= 5, temp_range_genus >= 5) %>%
  nest(-Pred_Genus) %>% 
  mutate(
    fit = map(data, ~ lm(log(attack_mass_corr) ~ one_kT, data = .x)),
    tidied = purrr::map(fit, tidy, conf.int = T),
    pred = map(fit, broom.mixed::augment)
  ) 
attack_genus_lm
#Add group names
attack_genus_lm$Group <- unlist(map(attack_genus_lm$data, ~ .x[["predator.ana.group"]][[1]])) 
attack_genus_lm <- unnest(attack_genus_lm,  pred)%>% 
  rename(log_attack_mass_corr = `log(attack_mass_corr)`)
attack_genus_lm

##---------------  Plotting Fig 2e, Part II-----------
## Attack and mortality plots were manually superimposed

attack_regr_genus <- ggplot(data = attack %>% filter(n_genus >= 5, temp_range_genus >= 5),
            aes(x = one_kT, y =attack_mass_corr, group = Pred_Genus)) + 
  geom_smooth(method = "lm", 
              aes(color = predator.ana.group), size = 1,  se = F,linetype = "dashed") +
  scale_color_manual(values =c("invertebrate" = invert_col, "vertebrate" = "navy")) +
  scale_y_continuous(name = expression(paste('Attack Rate (m'^2,' s'^-1,')')),
                     limits = c(10^-8.5, 10^-5.5), #note: same total range as mortality
                     labels = trans_format("log10", math_format(10^.x)),
                     breaks = c(10^-8, 10^-7, 10^-6), trans = "log10") +  
  scale_x_reverse(name = expression("1/kT"), limits = c(42.5, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, name = NULL)) +
  annotate("text", y = 10^-9.5, x = 40, hjust = 0, size = 6, label = NA, color = invert_col) +
  theme_single + theme_ticks + theme(legend.position = "none")
  
grid.newpage()
grid.draw(attack_regr_genus )

ggsave(attack_regr_genus , height = 4.7, width = 6.7, filename = '~/Desktop/attack_scaling.pdf')

ggsave(attack_regr_genus , height = 4.7, width = 6.7, filename = file.path(gdrive_path,'NAME YOUR FILE.pdf'))


############################################
##  Fig 2F - violin plot
# Plot genus-level slopes of thermal sensitivity
############################################

# attack slopes
attack_lm <- attack %>%
  filter(temp_range_genus >= 5, n_genus >= 5) %>%
  nest(-Pred_Genus) %>% # the group variable
  mutate(
    fit = map(data, ~ lm(log(attack_mass_corr) ~ one_kT, data = .x)), #this is a regression (y = attack.rate, x = one_kT), but you could adapt to whatever
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied) 
attack_lm

attack_lm <- filter(attack_lm, term != '(Intercept)') %>%
  mutate(Type = "Attack") %>%
  rename(Genus = Pred_Genus)
attack_lm

#mortality slopes
mort_lm_ecto <- mortality %>%
  filter(temp_range_genus >= 5, n_genus >= 5, Thermy == "Ectotherm") %>%
  nest(-Genus) %>% # the group variable
  mutate(
    fit = map(data, ~ lm(log(mass_corr_mortality) ~ one_kT, data = .x)), #this is a regression (y = attack.rate, x = one_kT), but you could adapt to whatever
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied) 

mort_lm_ecto <- filter(mort_lm_ecto, term != '(Intercept)') %>%
  mutate(Type = "Mortality")
mort_lm_ecto

#combine
mort_attack <- as_tibble(rbind(mort_lm_ecto, attack_lm)) 

# Plot

## 2) Plots the violin plot
vio_mort_attack_ea <- ggplot(mort_attack, aes(x = Type, y = -1 * estimate))+
  geom_violin(size = 1, color = "dodgerblue", aes(linetype = Type)) + 
  scale_linetype_manual(aes(Group = Type), values=c("Attack" = "dashed","Mortality" =  "solid")) +
  geom_jitter(shape = 21, color = "black",fill = "dodgerblue", width = 0.06, stroke = .6, size = 3) +
  scale_y_continuous(
    limits = c(-0.5, 1.1), #remove to see outliers
    position = "right", name=expression(paste('Thermal Sensitivity (E'[a],')')))+ 
  scale_x_discrete(name = NULL) +
  theme_single + theme(legend.position = "none")  + theme_ticks + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
vio_mort_attack_ea
ggsave(vio_mort_attack_ea,  height = 4.22, width = 5.9,  filename = '~/Desktop/Fig2F.pdf')

ggsave(vio_mort_attack_ea,  height = 4.22, width = 5.9,  filename = file.path(gdrive_path,'NAME YOUR FIGURE.pdf'))


###################################################################
############# Compare Mortality Variance  between Thermies #########
###################################################################

#Stats

mort_lm_endo <- mortality %>%
  filter(temp_range_genus >= 5, n_genus >= 5, Thermy == "Endotherm") %>%
  nest(-Genus) %>% # the group variable
  mutate(
    fit = map(data, ~ lm(log(mass_corr_mortality) ~ one_kT , data = .x)), #this is a regression (y = attack.rate, x = one_kT), but you could adapt to whatever
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied) 

mort_lm_endo <- filter(mort_lm_endo, term != '(Intercept)') %>%
  mutate(Type = "Mortality")
mort_lm_endo

# ecto vs endo range
range(mort_lm_ecto$estimate)
range(mort_lm_ecto$estimate)[2] - range(mort_lm_ecto$estimate)[1] #1.33 magnitude of ecto E range

range(mort_lm_endo$estimate)
range(mort_lm_endo$estimate)[2] - range(mort_lm_endo$estimate)[1] #1.37 magnitude of Endo E range - about the same

#standard dev  
sd(mort_lm_ecto$estimate) 
sd(mort_lm_endo$estimate) 

var_ecto <-var(mort_lm_ecto$estimate)
var_ecto
var_endo <- var(mort_lm_endo$estimate)
var_endo

#compare variance - not significantly different
sigma.test(mort_lm_ecto$estimate, sigmasq = var_endo)
