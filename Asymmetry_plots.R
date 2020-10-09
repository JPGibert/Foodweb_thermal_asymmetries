################################XXX XXXXX, Oct 14 2019 ################################## 

library(grid)
library(scales)
library(TeachingDemos)
library(ape)
library(phangorn)
library(Hmisc, pos = 100)
library(brms)
library(lme4)
library(lmerTest)
library(MuMIn)
library(bayestestR)
library(tidyverse)
library(egg)
library(broom.mixed)

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
                          plot.title = element_text(size = 18, face = "plain", hjust = 10),
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
  mutate(log_mort = log(Mortality_yr), log_mass = log(dry_mass_g))

 
##########################################################################
#--------------- Analyze Groups ----------

#------------------- Fish -----------------
# Data
mort_fish <- mortality1 %>%
  filter(Group == "Fish") %>%
  select(species, log_mass, log_mort, one_kT)

# Model, species as random effect
fish_mort_lmer <- lmer(log_mort ~ log_mass + one_kT + (1 |species), data = mort_fish)

# Results
summary(fish_mort_lmer) # mass slope = -0.29283
confint.merMod(fish_mort_lmer)
r.squaredGLMM(fish_mort_lmer)

# Mass-corrected  for plotting
mort_fish$log_mort_corr <- log(exp(mort_fish$log_mort)/
                                   exp(mort_fish$log_mass)^-0.29283)
#--------------Invertebrates -----------------
# Data
mort_invert <- mortality1 %>%
  filter(Group == "Invertebrate")

# Model, species as random effect
invert_mort_lmer <- lmer(log_mort ~ log_mass + one_kT + (1|species), data = mort_invert )

# Results
summary(invert_mort_lmer)  # mass slope = -0.24932
confint.merMod(invert_mort_lmer)
r.squaredGLMM(invert_mort_lmer)

# Or evaluate ectos together 
mort_ecto <- mortality1 %>%
  filter(Group == "Fish" | Group == "Invertebrate") 
mort_ecto_lmer <- lmer(log_mort ~ log_mass + one_kT + (1 + one_kT|Group) +(1|Group:species) +
                         (1 + log_mass|Group), data = mort_ecto)
summary(mort_ecto_lmer)
confint.merMod(mort_ecto_lmer, method = "Wald") # Wald handles multiple terms better
r.squaredGLMM(mort_ecto_lmer)

#--------------Mammals: comprehensive phylogeny available-----------------
# Data
mort_mammal_0 <- mortality1 %>%
  filter(Group == "Mammal") %>%
  select(species, log_mass, log_mort, one_kT)

#-- Add phylogeny ---

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
mort_mammal <- mort_mammal_0[mort_mammal_0$species %in% mam_tree_mcc$tip.label,]

#----------------- bayesian model for mammals using brms, species as random effect -----------------


# variance and covariance from phylogeny
A <- ape::vcv.phylo(mam_tree_mcc)
# Priors
get_prior(formula = log_mort ~ log_mass + one_kT + (1|gr(species, cov = A)), data = mort_mammal)

priors1 <- c(set_prior("normal(0, 1)", class = "Intercept"),
             set_prior("normal(0, 1)", class = "b", coef = "one_kT"), #slope = 0 for endotherms
             set_prior("normal(-0.25, 1)", class = "b", coef = "log_mass")
)


system.time(model_mort_mammal <- brm(
    log_mort ~ log_mass + one_kT + (1|gr(species, cov = A)), 
    data = mort_mammal, 
    family = gaussian(), 
    data2 = list(A = A),
    prior = priors1,
    cores = parallel::detectCores() -1) 
  )


# Model Results
print(summary(model_mort_mammal), digits = 5) 
bayes_R2(model_mort_mammal) 
bayes_R2(model_mort_mammal, re.form = NA) # exclude random effects
pd_mamm <- p_direction(model_mort_mammal) # p direction
pd_to_p(pd_mamm$pd[2], "two-sided") # bayesian equivalent to p value for mass scaling
pd_to_p(pd_mamm$pd[3], "two-sided") # bayesian equivalent to p value for temp scaling


# Get intercept of mass-corrected scaling for plotting
mort_mammal$log_mort_corr <- log(exp(mort_mammal$log_mort)/
                                       exp(mort_mammal$log_mass)^ -0.179)

priors2 <- c(set_prior("normal(0, 1)", class = "Intercept"),
             set_prior("normal(0, 1)", class = "b", coef = "one_kT")
)

model_mort_mammal_corr <- brm(
    log_mort_corr ~ one_kT + (1|gr(species, cov = A)), 
    data = mort_mammal, 
    family = gaussian(), 
    data2 = list(A = A),
    prior = priors2,
    cores = parallel::detectCores() -1) 
   

print(summary(model_mort_mammal_corr), digits = 5) 


# ----------- Birds ------------------

mort_bird0 <- mortality1 %>%
  filter(Group == "Bird") #%>%
 
# Read in 100 phylogenetic trees - download at https://data.vertlife.org/  birdtree > PatchClade > Stage 1 > EricsonStage1_0001_1000.zip
bird_tree <- read.tree('/Users/jgradym/Google Drive/Phylo_trees/Upham_Phylos/Birds/mnt/data/projects/birdphylo/Tree_sets/Stage1_full_data/CombinedTrees/EricsonStage1Full_1.tre')[1:100]

# get maximum credible tree
bird_tree_mcc0 <- maxCladeCred(bird_tree) 

#prune down tree to species in dataset
bird_spmatch <- match(bird_tree_mcc0$tip.label, #get matched species
                     mort_bird0$species)

bird_tree_mcc <- drop.tip(bird_tree_mcc0 , bird_tree_mcc0$tip.label[is.na(bird_spmatch)]) #prune unused spp

mort_bird <- mort_bird0[mort_bird0$species %in% bird_tree_mcc$tip.label,]


# variance and covariance from phylogeny
A <- ape::vcv.phylo(bird_tree_mcc)

# Run model
system.time(model_mort_bird <- brm(
    log_mort ~ log_mass + one_kT + (1|gr(species, cov = A)), 
    data = mort_bird, 
    family = gaussian(), 
    data2 = list(A = A),
    prior = priors1,
    cores = parallel::detectCores() -1) 
  ) 

# Model Results
print(summary(model_mort_bird), digits = 4) 
bayes_R2(model_mort_bird) 
bayes_R2(model_mort_bird, re.form = NA)  # fixed effects only

pd_bird <- p_direction(model_mort_bird) # p direction
pd_to_p(pd_bird$pd[2], "two-sided") # Bayesian equivalent to p value for mass scaling
pd_to_p(pd_bird$pd[3], "two-sided") # Bayesian equivalent to p value for temp scaling

# Get intercept of mass-corrected scaling for plotting
mort_bird <- mort_bird  %>%
  mutate(log_mort_corr = log(exp(mort_bird$log_mort)/
                               exp(mort_bird$log_mass)^ -0.1642))

# Run Model to look at mass-corrected mortality with temperature
model_mort_bird_corr <- brm(
    log_mort_corr ~ one_kT + (1|gr(species, cov = A)), 
    data = mort_bird, 
    family = gaussian(), 
    data2 = list(A = A),
    prior = priors2,
    cores = parallel::detectCores() -1,
    control = list(adapt_delta = 0.9) 
  ) 
print(summary(model_mort_bird_corr), digits = 5) 


#------- Add mass-corrected mortality to mortality dataset for plotting ----

mortality1 <- mortality1 %>% 
  mutate(mass_corr_mortality = NA)
bird <- mortality1[mortality1$Group == "Bird",]
mortality1$mass_corr_mortality[mortality1$Group == "Bird"] <- bird$Mortality_yr/bird$dry_mass_g^-0.164

mammal <- mortality1[mortality1$Group == "Mammal",]
mortality1$mass_corr_mortality[mortality1$Group == "Mammal"] <- mammal$Mortality_yr/mammal$dry_mass_g^-0.179

invert <- mortality1[mortality1$Group == "Invertebrate",]
mortality1$mass_corr_mortality[mortality1$Group == "Invertebrate"] <- invert$Mortality_yr/invert$dry_mass_g^-0.249

fish <- mortality1[mortality1$Group == "Fish",]
mortality1$mass_corr_mortality[mortality1$Group == "Fish"] <- fish$Mortality_yr/fish$dry_mass_g^-0.293


#-------Adjust and columns for genus level analysis ------------

# Get number measurements per species or per genus, if you want to restrict by individuals per species or genus:, 
mortality2 <- mortality1 %>% #firs
  mutate(Genus = word(species,1, sep = "_")) %>%
  group_by(Genus) %>% 
  add_count(Genus, name = "n_genus") %>%# individuals per genus
  mutate(temp_range_genus = max(Temp_C) - min(Temp_C)) %>%
  ungroup()

#Add a temperature range for individuals in a species (only relevant for ectos)
mortality3 <- mortality2 %>%
  group_by(species) %>% 
  add_count(species, name = "n_species") %>%
  mutate(temp_range_spp = max(Temp_C) - min(Temp_C)) %>%
  ungroup()

# Add thermy
mortality4 <- mortality3 %>%
  filter(Group == "Bird" | Group == "Mammal"|  Group == "Invertebrate" | Group == "Fish") %>%
  mutate(Thermy = if_else(Group == "Bird" | Group == "Mammal", "Endotherm", "Ectotherm"))

# Final dataset
mortality <- mortality4 # final dataset


##########################################
 #### Panel 2a
##########################################

bird_lab  <- expression(paste('y = -0.026x + 0.39')) #from model summary
mamm_lab  <- expression(paste('y = 0.010x - 0.46'))

# for manual plotting, results from Bayesian mixed effect analysis with phylogeny
bird_fit <- data.frame(x = c(min(mort_bird$one_kT),  max(mort_bird$one_kT)), 
                     y = c(exp((0.386 + (-0.0258 * min(mort_bird$one_kT)))), 
                     exp(0.386 + (-0.0258 *  max(mort_bird$one_kT)))),
                     Group = c("Bird", "Bird"))

mammal_fit <- data.frame(x = c(min(mort_mammal$one_kT),  max(mort_mammal$one_kT)), 
                         y = c( exp((-0.464 + (0.00980 * min(mort_mammal$one_kT)))), 
                         exp((-0.464 + (0.00980 *  max(mort_mammal$one_kT))))),
                         Group = c("Mammal", "Mammal"))

endo_fit <- data.frame(rbind(bird_fit, mammal_fit))



#------------- Plot Fig 2A  -----------------
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
  annotate("text", x = 41, size = 6, y = .025,  hjust = 0, label = bird_lab, color = "brown4") 
grid.newpage()
grid.draw(mortal_plot_endo)
ggsave(mortal_plot_endo, height = 4.6, width = 6.7, filename = '~/Desktop/Fig2A.pdf')

# Save plot
ggsave(mortal_plot_endo, height = 4.6, width = 6.7, filename = file.path(gdrive_path, "NAME YOUR FILE.pdf"))


###########################################
#### Panel 2b - species as random effect
##########################################

# Prepare plot
invert_lab  <- expression(paste('y = -0.58x + 23.8'))
fish_lab  <- expression(paste('y = -0.50x + 21.1'))

# from models summary
invert_fit <- data.frame(x = c(min(mort_invert$one_kT),  max(mort_invert$one_kT)), 
                       y = c(exp((23.84 + (-0.581 * min(mort_invert$one_kT)))), 
                             exp((23.84 + (-0.581 *  max(mort_invert$one_kT))))),
                       Group = c("Invertebrate", "Invertebrate"))


fish_fit <- data.frame(x = c(min(mort_fish$one_kT), max(mort_fish$one_kT)), 
                         y = c( exp((21.05 + (-0.504 * min(mort_fish$one_kT)))), 
                                exp((21.05 + (-0.504 *  max(mort_fish$one_kT))))),
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
ggsave(mortal_plot_ecto, height = 4.6, width = 6.7, filename = '~/Desktop/ecto_file.pdf')

####################################################################################
#### Panel 2c - genus level regressions, min 5 individual per genus, temp range 5 C 
####################################################################################

# Thermy
endo_genus <- mortality %>%
  filter(Thermy == "Endotherm") %>%
  pull(Genus) %>%
  unique()

# Regressions for genus level plotting 
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

# add column info back
genus_mort_lm$Thermy<-if_else(genus_mort_lm$Genus %in% endo_genus, "Endotherm", "Ectotherm" )
genus_mort_group <- map(genus_mort_lm $data, ~ .x[["Group"]][[1]])
genus_mort_lm$Group <- unlist(genus_mort_group) 
unnest(genus_mort_lm, tidied) %>% filter(term == "one_kT")

# Group level, species as random effect (intercept)
group_mort_lmer <- mortality %>%
  filter(n_genus >= 5, temp_range_genus >= 5) %>%
  nest(-Group) %>% 
  mutate(
    fit = purrr::map(data, ~ lmer(
      log(mass_corr_mortality) ~ one_kT + (1|species), data = .x)),
    pred = map(fit, re.form = NA, broom.mixed::augment)
  ) 
group_mort_lmer<- unnest(group_mort_lmer ,pred)
group_mort_lmer$Thermy<-if_else(word(group_mort_lmer$species, sep = "_") %in% endo_genus, "Endotherm", "Ectotherm" )
group_mort_lmer

#--- plot Fig 2C ------

mortal_genus_regr <- ggplot(data = mortality %>% 
                              filter(n_genus >= 5, temp_range_genus >= 5),
                            aes(y = mass_corr_mortality, x = one_kT)) +
    geom_smooth(method = "lm", se =- F, size = 0.3, aes(group = Genus, color = Group)) + 
  geom_line(data = group_mort_lmer, size = 1.5,
            aes(x = one_kT, y = exp(.fitted), group = Group, color = Group)) +
  scale_color_manual(values =c("Invertebrate" = invert_col, "Fish" = fish_col, 
                                 "Mammal" = mamm_col, "Bird" = bird_col)) +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{-alpha},')')), 
                     trans = "log", breaks = c(0.3, 1, 3,10), labels = c("0.3", "1", "3", "10"), 
                     limits = c(.15, 10)) +  
  scale_x_reverse(name =NULL, limits = c(44.4, 38), 
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
                     position = "right",name = expression(paste('Thermal Sensitivity (E'[a],')')))+ 
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


# Mixed model test - Group and species as random effect  (species nested in group)
ecto_attack_lmer <-  lmer(data = attack, log(attack.rate) ~ one_kT + log(predator.mass.mg) +
                            (1 + one_kT|predator.ana.group) + (1|predator.ana.group:predator.species) +
                            (1 + log(predator.mass.mg)|predator.ana.group)
                          )
summary(ecto_attack_lmer)
confint.merMod(ecto_attack_lmer, method = "Wald")
r.squaredGLMM(ecto_attack_lmer)

# Mixed model for inverts separately
invert_attack <- attack %>% filter( predator.ana.group == "invertebrate")
invert_attack_lmer <- lmer(data = invert_attack, log(attack.rate) ~ one_kT + log(predator.mass.mg) +
                             (1|predator.species))
summary(invert_attack_lmer)  #0.56031 mass scaling slope
confint.merMod(invert_attack_lmer)
r.squaredGLMM(invert_attack_lmer)
invert_attack$attack_mass_corr <- invert_attack$attack.rate/invert_attack$predator.mass.mg^0.560

#vMixed model for vertebrates separately
vert_attack <- attack %>% filter(predator.ana.group == "vertebrate")
vert_attack_lmer <- lmer(data = vert_attack, log(attack.rate) ~ one_kT + log(predator.mass.mg) +
                           (1|predator.species))
summary(vert_attack_lmer) 
confint.merMod(vert_attack_lmer)
r.squaredGLMM(vert_attack_lmer)

# add mass-corrected attack rate as column
vert_attack$attack_mass_corr <- vert_attack$attack.rate/vert_attack$predator.mass.mg^0.757

# combine vert and invert attack
attack <- as_tibble(rbind(vert_attack, invert_attack))

# get attack regressions with temperature
attack_genus_lm <- attack %>%
  filter(n_genus >= 5, temp_range_genus >= 5) %>%
  nest(-Pred_Genus) %>% 
  mutate(
    fit = map(data, ~ lm(log(attack_mass_corr) ~ one_kT, data = .x)),
    tidied = purrr::map(fit, tidy, conf.int = T),
    pred = map(fit, broom.mixed::augment)
  ) 

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
