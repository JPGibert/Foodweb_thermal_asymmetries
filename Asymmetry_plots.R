################################XXX XXXXX, Oct 14 2019 ################################## 

library(grid)
library(scales)
library(TeachingDemos)
library(ape)
library(phangorn)
library(Hmisc)
library(brms)
library(lme4)
library(lmerTest)
library(MuMIn)
library(bayestestR)
library(tidyverse)
library(egg)
library(broom.mixed)
#library(tidybayes)

# ----  Adjust all file paths to read in data and save plots
github_path <- file.path('/Users/jgradym/Documents/GitHub/Foodweb_thermal_asymmetries/Data')


################################## Plotting ###########################33

# Colors

ecto_col <- "deepskyblue2"
ecto_col2 <- "dodgerblue"
endo_col <- "red2"
endo_col2 <- "red3"
ecto_mort_col <- "deepskyblue2"
ecto_attack_col <- "navy"

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
mortality0 <- read_csv(file.path(github_path, 'McCoy_mortality_updated.csv')) 
mortality0$one_kT <- 1/(8.617e-5*(mortality0$Temp_C + 273.15)) # Inverse Temperature

mortality1 <- mortality0 %>%
  rename(species = Species) %>%
  arrange(Group, species) %>%
  mutate(log_mort = log(Mortality_yr), log_mass = log(dry_mass_g),
         Genus = word(species, sep = "_"))

# Get number measurements per species or per genus, if you want to restrict by individuals per species or genus:, 
mortality2 <- mortality1 %>% #firs
  mutate(Genus = word(species, 1, sep = "_")) %>%
  group_by(Genus) %>% 
  add_count(Genus, name = "n_genus") %>%# individuals per genus
  mutate(temp_range_genus = max(Temp_C) - min(Temp_C)) %>%
  mutate(min_one_kT = min(one_kT), max_one_kT = max(one_kT)) %>%
  ungroup() %>%
  filter(Group == "Bird" | Group == "Mammal"|  Group == "Invertebrate" | Group == "Fish") %>%
  mutate(Thermy = if_else(Group == "Bird" | Group == "Mammal", "Endotherm", "Ectotherm")) %>%
  arrange(Thermy, Genus) 


##########################################################################
#--------------- Analyze Groups ----------

#########################################################################3
# Mixed Model Ectotherms
mort_ecto <- mortality2 %>% filter(Thermy == "Ectotherm")

# Model
lmer_mort_ecto <- lmer(log_mort ~ one_kT + log_mass +  (one_kT|Genus),  data = mort_ecto)
lmer_mort_ecto_full <- lmer(log_mort ~ one_kT + log_mass + (one_kT|Group/Genus),  data = mort_ecto)

# Results
rand(lmer_mort_ecto) # AIC
rand(lmer_mort_ecto_full) # AIC
summary(lmer_mort_ecto) #summary
summary_mort_ecto <- as.data.frame(summary(lmer_mort_ecto)$coefficients) 
confint.merMod(lmer_mort_ecto, method = "Wald") #confidence intervals
r.squaredGLMM(lmer_mort_ecto) #r2

# dataframe for plotting
mort_ecto_results_all <- tibble::rownames_to_column(summary_mort_ecto) %>% 
  rename(parameter = rowname) %>%
  select(parameter, Estimate) %>%
  pivot_wider(names_from = parameter, values_from = Estimate) %>%
  rename(intercept = `(Intercept)`, temp_slope = one_kT, mass_slope = log_mass)
mort_ecto_results_all$min_one_kT <- min(mort_ecto$min_one_kT)
mort_ecto_results_all$max_one_kT <- max(mort_ecto$max_one_kT)
mort_ecto_results_all$min_mort <- exp(mort_ecto_results_all$intercept + mort_ecto_results_all$temp_slope * mort_ecto_results_all$min_one_kT)
mort_ecto_results_all$max_mort <- exp(mort_ecto_results_all$intercept + mort_ecto_results_all$temp_slope * mort_ecto_results_all$max_one_kT)
mort_ecto_results_all$Thermy <- "Ectotherm"
mort_ecto_results_all



# Mixed Model Ectotherms
attack <- mortality2 %>% filter(Thermy == "Endotherm")

# Model
lmer_mort_endo <-  lmer(log_mort ~ log_mass + one_kT + (one_kT|Genus), data = mort_endo )
lmer_mort_endo_full <-  lmer(log_mort ~ log_mass+ one_kT + (one_kT|Group/Genus), data = mort_endo )

# Results
rand(lmer_mort_endo) # AIC
rand(lmer_mort_endo_full) # AIC
print(summary(lmer_mort_endo), digits = 5)
summary_mort_endo <- as.data.frame(summary(lmer_mort_endo)$coefficients)
confint.merMod(lmer_mort_endo, method = "Wald") #confidence intervals
r.squaredGLMM(lmer_mort_endo) #r2
# dataframe for plotting

mort_endo_results_all <- tibble::rownames_to_column(summary_mort_endo) %>% 
  rename(parameter = rowname) %>%
  select(parameter, Estimate) %>%
  pivot_wider(names_from = parameter, values_from = Estimate) %>%
  rename(intercept = `(Intercept)`, temp_slope = one_kT, mass_slope = log_mass)
mort_endo_results_all$min_one_kT <- min(mort_endo$min_one_kT)
mort_endo_results_all$max_one_kT <- max(mort_endo$max_one_kT)
mort_endo_results_all$min_mort <- exp(mort_endo_results_all$intercept + mort_endo_results_all$temp_slope * mort_endo_results_all$min_one_kT)
mort_endo_results_all$max_mort <- exp(mort_endo_results_all$intercept + mort_endo_results_all$temp_slope * mort_endo_results_all$max_one_kT)
mort_endo_results_all$Thermy <- "Endotherm"
mort_endo_results_all



# add mass correction for plotting

mort_ecto$mort_mass_corr <- mort_ecto$Mortality_yr/mort_ecto$dry_mass_g^mort_ecto_results_all$mass_slope
mortality2$mort_mass_corr[mortality2$Thermy == "Ectotherm"] <- mort_ecto$mort_mass_corr 

mort_endo$mort_mass_corr <- mort_endo$Mortality_yr/mort_endo$dry_mass_g^mort_endo_results_all$mass_slope
mortality2$mort_mass_corr[mortality2$Thermy == "Endotherm"] <- mort_endo$mort_mass_corr 


###########################################
#### Panel 2a- Endotherms
##########################################
#------------- Plot Fig 2A  -----------------
mort_endo_results_all
endo_lab <- expression(paste('y = -0.011x + 0.093'))
mortal_plot_endo <- ggplot(mortality2 %>% filter(Thermy == "Endotherm"),
                           aes(x = one_kT, y = mort_mass_corr)) + 
  theme_single + theme(legend.position = "none") +  theme_ticks +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{-alpha},')')), 
                     trans = "log", breaks = c(0.1, 1, 10), labels = c("0.1", "1", "10"), 
                     #limits = c(0.02, 10)
  ) +  
  scale_x_reverse( name = NULL, breaks = c(46,44, 42,40, 38),
                   limits = c(45.1, 38.1),
                   sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, breaks = c(-20, -10, 0, 10, 20, 30),
                                       name = "Ambient Temperature (ºC)")) +
  geom_point(size = 2.5, shape = 21, color = "black",  stroke = .5,fill = "red3") + 
  geom_segment(data = mort_endo_results_all, size = 1,  color = "#4A120A",
            aes(x = min_one_kT, xend = max_one_kT, y = min_mort, yend = max_mort)) +
  annotate("text", x = 42, size = 6, y = 0.04, hjust = 0,  fontface = "bold",label = endo_lab , color = "red3") 


grid.newpage()
grid.draw(mortal_plot_endo)
mortal_plot_endo <- set_panel_size(mortal_plot_endo, width = unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(mortal_plot_endo )

pdf('~/Desktop/Fig_2/Fig2A.pdf')
grid.draw(mortal_plot_endo )
dev.off()
system2(command = "pdfcrop", 
        args  = c('~/Desktop/Fig_2/Fig2A.pdf', 
                  '~/Desktop/Fig_2/Fig2A.pdf'))
# Save plot


###########################################
#### Panel 2b - Ectotherms
##########################################


mort_ecto_results_all
ecto_lab <- expression(paste('y = -0.53x + 25'))

mortal_plot_ecto <- ggplot(mortality2 %>% filter(Thermy == "Ectotherm"),
                           aes(x = one_kT, y = mort_mass_corr)) + 
  scale_y_continuous(trans = "log10",  name = expression(paste('Mortality (yr'^-1,' g'^{-alpha},')')), 
                     breaks = c(0.1,  1, 10), labels = c("0.1", "1",  "10"), 
                     limits = c(0.02, 30), 
                     position = "right")+ 
  scale_x_reverse(name = expression("1/kT"),
                  limits = c(43, 38.1), 
                  breaks = c(42, 40, 38),
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, 
                                      name = "Ambient Temperature (ºC)")) +
  theme(legend.position = "none") +theme_single + 
  geom_point(size = 2.5, shape = 21, color = "black",  stroke = .5, fill  = ecto_mort_col ) +
  geom_segment(data = mort_ecto_results_all, size = 1,  color = "#23598F",
               aes(x = min_one_kT, xend = max_one_kT, y = min_mort, yend = max_mort)) +
  annotate("text", x = 40.3, size = 6, y = 0.026, hjust = 0,  fontface = "bold",label = ecto_lab , color = ecto_mort_col ) +
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text.x = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")), 
        axis.text.x.top = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")), 
        axis.text.y.right = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")),
        axis.text.y.left = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")))
grid.newpage()
grid.draw(mortal_plot_ecto) 

mortal_plot_ecto <- set_panel_size(mortal_plot_ecto, width = unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(mortal_plot_ecto)

pdf('~/Desktop/Fig_2/Fig2B.pdf')
grid.draw(mortal_plot_ecto)
dev.off()

system2(command = "pdfcrop", 
        args  = c('~/Desktop/Fig_2/Fig2B.pdf', 
                  '~/Desktop/Fig_2/Fig2B.pdf'))

# Save plot


####################################################################################
#### Panel 2c - genus level regressions, min 5 individual per genus, temp range 5 C 
####################################################################################

mort_results_all <- as_tibble(rbind(mort_endo_results_all, mort_ecto_results_all))

# Genus analysis
mortality3 <- mortality2 %>%
  filter(temp_range_genus >= 5, n_genus >= 5)

# Plotting data, by Genus
mort_ecto_genus <- mortality3 %>% 
  filter(Thermy == "Ectotherm") %>%
  group_by(Genus) %>%
  slice_head() %>%
  arrange(Genus) %>%
  select(Thermy, Group, Genus, min_one_kT,  max_one_kT)
mort_ecto_genus 


mort_ecto3 <- mortality3 %>% filter(Thermy == "Ectotherm")

# Model
lmer_mort_gen_ecto <- lmer(log_mort ~ one_kT + log_mass +  (one_kT|Genus),  data = mort_ecto3)
lmer_mort_gen_full_ecto <- lmer(log_mort ~ one_kT + log_mass + (one_kT|Group/Genus),  data = mort_ecto3)

# Results
rand(lmer_mort_gen_ecto ) # AIC
rand(lmer_mort_gen_full_ecto) # AIC
summary(lmer_mort_gen_ecto )
confint.merMod(lmer_mort_gen_ecto, method = "Wald") #confidence intervals
coef(lmer_mort_gen_ecto)
r.squaredGLMM(lmer_mort_gen_ecto) #r2


#-------- Regression fits for plotting ---------
# Ectotherm Genera
mort_ecto_results <- as_tibble(tibble::rownames_to_column(coef(lmer_mort_gen_ecto)$Genus, "Genus") %>%
                                 rename(intercept = `(Intercept)`, temp_slope = one_kT, mass_slope = log_mass)) %>%
  arrange(Genus)
mort_ecto_results$Group <- mort_ecto_genus$Group[mort_ecto_genus$Genus %in% mort_ecto_results$Genus]

mort_ecto_results$min_one_kT <- mort_ecto_genus$min_one_kT
mort_ecto_results$max_one_kT <- mort_ecto_genus$max_one_kT
mort_ecto_results$min_y <- exp(mort_ecto_results$intercept + mort_ecto_results$temp_slope * mort_ecto_results$min_one_kT)
mort_ecto_results$max_y <- exp(mort_ecto_results$intercept + mort_ecto_results$temp_slope * mort_ecto_results$max_one_kT)
mort_ecto_results$Thermy <- "Ectotherm"
mort_ecto_results$Group <- mort_ecto_genus$Group[mort_ecto_genus$Genus %in% mort_ecto_results$Genus]
mort_ecto_results


mortality_genus_ecto <- mortality_genus %>% filter(Thermy == "Ectotherm")


# Endotherm Genera


# Model
mort_endo_genus <- mortality3 %>% 
  filter(Thermy == "Endotherm") %>%
  group_by(Genus) %>%
  slice_head() %>%
  arrange(Genus) %>%
  select(Thermy, Group, Genus, min_one_kT,  max_one_kT)

mort_endo_genus 

# data
mort_endo3 <- mortality3 %>% filter(Thermy == "Endotherm")

# model
lmer_mort_gen_endo <- lmer(log_mort ~ one_kT + log_mass +  (one_kT|Genus),  data = mort_endo3 )
lmer_mort_gen_full_endo <- lmer(log_mort ~ one_kT + log_mass + (one_kT|Group/Genus),  data = mort_endo3 )

# Results
rand(lmer_mort_gen_endo) # AIC
rand(lmer_mort_gen_full_endo) # AIC
summary(lmer_mort_gen_endo )
confint.merMod(lmer_mort_gen_endo , method = "Wald") #confidence intervals
coef(lmer_mort_gen_endo )
r.squaredGLMM(lmer_mort_gen_endo ) #r2


mort_endo_results <- as_tibble(tibble::rownames_to_column(coef(lmer_mort_gen_endo)$Genus, "Genus") %>%
                                 rename(intercept = `(Intercept)`, temp_slope = one_kT, mass_slope = log_mass)) %>%
  arrange(Genus)

mort_endo_results$Group <- mort_endo_genus$Group[mort_endo_genus$Genus %in% mort_endo_results$Genus]
mort_endo_results$min_one_kT <- mort_endo_genus$min_one_kT
mort_endo_results$max_one_kT <- mort_endo_genus$max_one_kT
mort_endo_results$min_y <- exp(mort_endo_results$intercept + mort_endo_results$temp_slope * mort_endo_results$min_one_kT)
mort_endo_results$max_y <- exp(mort_endo_results$intercept + mort_endo_results$temp_slope * mort_endo_results$max_one_kT)
mort_endo_results$Thermy <- "Endotherm"
mort_endo_results$Group <- mort_endo_genus$Group[mort_endo_genus$Genus %in% mort_endo_results$Genus]
mort_endo_results

# Genus regressions
lmer_mort_gen_ecto <- coef(lmer_mort_ecto)$`Group:Genus`

lmer_mort_gen_ecto <- as_tibble(rownames_to_column(lmer_mort_gen_ecto)) %>%
  rename(intercept = `(Intercept)`, temp_slope = one_kT, mass_slope = log_mass) %>%
  mutate(Group = word(rowname, 1, sep = ":"), Genus = word(rowname, 2, sep = ":")) %>%
  select(-rowname) %>%
  select(Group, Genus, everything())
library(data.table)


# Combine 
mort_results <- as_tibble(rbind(mort_endo_results, mort_ecto_results))
mort_results$Type <- "mortality" 
#--- plot Fig 2C ------

mortal_genus_regr <- ggplot(data =mort_results) +
  geom_segment(aes(x = min_one_kT, xend = max_one_kT, y = min_y, yend = max_y, color = Thermy), size = 0.3) +
  geom_segment(data =mort_results_all,
               aes(x = min_one_kT, xend = max_one_kT, y = min_mort, yend = max_mort, color = Thermy), size = 2) +
  scale_color_manual(values =c("Endotherm" = endo_col2, "Ectotherm" = ecto_mort_col )) +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{-alpha},')')), 
                     trans = "log", 
                     limits = c(.15, 10),
                     breaks = c(0.3, 1, 3,10), labels = c("0.3", "1", "3", "10")) + 
  scale_x_reverse(name =NULL, 
                  limits = c(45, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15,  name = NULL)) + 
  theme_single + theme(legend.position = "none") +theme_ticks 
grid.newpage()
grid.draw(mortal_genus_regr)

# Save plot
mortal_genus_regr <- set_panel_size(mortal_genus_regr , width = unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(mortal_genus_regr)

pdf('~/Desktop/Fig_2/Fig2C.pdf')
grid.draw(mortal_genus_regr)
dev.off()

system2(command = "pdfcrop", 
        args  = c('~/Desktop/Fig_2/Fig2C.pdf', 
                  '~/Desktop/Fig_2/Fig2C.pdf'))


# supplmemental, genus color coded with data
mort_colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')

# Ectotherms
mortality3$mort_mass_corr <- mortality3$Mortality_yr/mortality3$dry_mass_g^mort_ecto_results_all$mass_slope

mortal_ecto_regr <- ggplot(data = mort_results %>% filter(Thermy == "Ectotherm")) +
  geom_point(data = mortality3 %>% filter (Thermy == "Ectotherm"), shape = 21, stroke = 0.75, size =2,
             aes(x = one_kT, y = mort_mass_corr, color = Genus)) +
  geom_segment(aes(x = min_one_kT, xend = max_one_kT, y = min_y, yend = max_y, color = Genus), size = 0.7) +
  scale_color_manual(values = mort_colors ) +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{-alpha},')')), 
                     trans = "log", 
                    # limits = c(.3, 10),
                     breaks = c(0.3, 1, 3,10), labels = c("0.3", "1", "3", "10")) + 
  scale_x_reverse(name =NULL, 
                  #limits = c(42, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15,  name = NULL)) + 
  theme_single + theme(legend.position = "bottom") +theme_ticks 
grid.newpage()
grid.draw(mortal_ecto_regr)

mortal_ecto_regr  <- set_panel_size(mortal_ecto_regr  , width = unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(mortal_ecto_regr )

pdf('~/Desktop/Fig_2/sup_ecto_mort_regr.pdf')
grid.draw(mortal_ecto_regr )
dev.off()

system2(command = "pdfcrop", 
        args  = c('~/Desktop/Fig_2/sup_ecto_mort_regr.pdf', 
                  '~/Desktop/Fig_2/sup_ecto_mort_regr.pdf'))


# Endotherms
mortality3$mort_mass_corr <- mortality3$Mortality_yr/mortality3$dry_mass_g^mort_endo_results_all$mass_slope

mortal_endo_regr <- ggplot(data = mort_results %>% filter(Thermy == "Endotherm")) +
  geom_point(data = mortality3 %>% filter (Thermy == "Endotherm"), shape = 21, stroke = 0.75, size =2,
             aes(x = one_kT, y = mort_mass_corr, color = Genus)) +
  geom_segment(aes(x = min_one_kT, xend = max_one_kT, y = min_y, yend = max_y, color = Genus), size = 0.7) +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{-alpha},')')), 
                     trans = "log", 
                     #limits = c(.3, 10),
                     breaks = c(0.3, 1, 3,10), labels = c("0.3", "1", "3", "10")) + 
  scale_x_reverse(name =NULL, 
                  #limits = c(42, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15,  name = NULL)) + 
  theme_single + theme(legend.position = "bottom") +theme_ticks 
grid.newpage()
grid.draw(mortal_endo_regr)

mortal_endo_regr  <- set_panel_size(mortal_endo_regr  , width = unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(mortal_endo_regr )

pdf('~/Desktop/Fig_2/sup_endo_mort_regr.pdf')
grid.draw(mortal_endo_regr )
dev.off()

system2(command = "pdfcrop", 
        args  = c('~/Desktop/Fig_2/sup_endo_mort_regr.pdf', 
                  '~/Desktop/Fig_2/sup_endo_mort_regr.pdf'))

##########################################
#### Panel 2d, Mortality violin plot

##########################################


vio_mort_attack_ea <- ggplot(mort_results, aes(x = Thermy, y = -1 * temp_slope)) +
  geom_violin(size = 1, aes(color = Thermy)) + 
  geom_jitter(aes(fill = Thermy), shape = 21, color = "black", width = 0.06, stroke = .6, size = 3) +
  scale_fill_manual(values =c("Ectotherm" = ecto_mort_col, "Endotherm" = endo_col2)) +
  scale_color_manual(values =c("Ectotherm" = ecto_mort_col , "Endotherm" = endo_col2)) +
  scale_y_continuous(position = "right", limits = c(-0.15, 0.85), breaks = c(0, 0.2, 0.4, 0.6, 0.8),
                     name=expression(paste('Thermal Sensitivity (E'[a],')')))+ 
  scale_x_discrete(name = NULL) +
  theme_single + theme(legend.position = "none")  + theme_ticks + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
vio_mort_attack_ea

# Save plot
# Save plot
vio_mort_attack_ea<- set_panel_size(vio_mort_attack_ea, width = unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(vio_mort_attack_ea)

pdf('~/Desktop/Fig_2/Fig2D.pdf')
grid.draw(vio_mort_attack_ea)
dev.off()

system2(command = "pdfcrop", 
        args  = c('~/Desktop/Fig_2/Fig2D.pdf', 
                  '~/Desktop/Fig_2/Fig2D.pdf'))




####################################################
##v  Fig 2E - Mortality & Attack regressions (ectos)
####################################################
##---------------------------------
## Plotting Fig 2e, Part I (mortality)

mortal_genus_regr <- ggplot(data =mort_results %>% filter(Thermy == "Ectotherm")) +
  geom_segment(aes(x = min_one_kT, xend = max_one_kT, y = min_y, yend = max_y, color = Thermy), size = .5) +

  scale_color_manual(values =c("Endotherm" = endo_col, "Ectotherm" = ecto_mort_col)) +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{-alpha},')')), 
                     trans = "log10", position = "right", limits = c(0.3, 30),
                     breaks = c(0.3, 3, 30), labels = c("0.3", "3", "30")) + 
  scale_x_reverse(name =NULL, 
                  limits = c(42.5, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15,  name = NULL)) + 
  theme_single + theme(legend.position = "none") +theme_ticks 
grid.newpage()
grid.draw(mortal_genus_regr)

mortal_genus_regr<- set_panel_size(mortal_genus_regr, width = unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(mortal_genus_regr)

pdf('~/Desktop/Fig_2/mortal_genus_2E.pdf')
grid.draw(mortal_genus_regr)
dev.off()

system2(command = "pdfcrop", 
        args  = c('~/Desktop/Fig_2/mortal_genus_2E.pdf', 
                  '~/Desktop/Fig_2/mortal_genus_2E.pdf'))


##---------------------------------
## PART 2 (attack)

# Read in data
attack <- read_csv(file.path(github_path, 'Lietal_oikos_2017_data.csv'))
attack$one_kT <- 1/(8.617e-5*(attack$temperature.degree.celcius +273.15))

# Add columns
attack1 <- attack %>%
  mutate(Genus = word(predator.species)) %>%
  add_count(predator.species, name = "n_species") %>%# exclude taxa with few members
  add_count(Genus, name = "n_genus")#adjust n by changing filter(n_species/n_genus = )

attack2  <- attack1  %>%
  group_by(Genus) %>% 
  mutate(temp_range_genus = max(temperature.degree.celcius) - min(temperature.degree.celcius)) %>%
  mutate(log_attack = log(attack.rate), log_mass = log(predator.mass.mg)) %>%
  mutate(min_one_kT = min(one_kT), max_one_kT = max(one_kT)) %>%
  rename(Group = predator.ana.group ) %>%
  ungroup() %>%
  arrange(Genus)



# Model
lmer_attack<-  lmer(log_attack ~ log_mass + one_kT + (one_kT|Genus), data = attack2 )
lmer_attack_full <-  lmer(log_attack ~ log_mass+ one_kT + (one_kT|Group/Genus), data = attack2)

# Results
rand(lmer_attack) # AIC
rand(lmer_attack_full ) # AIC
summary(lmer_attack)
attack_results_all <- summary(lmer_attack)$coefficients
confint.merMod(lmer_attack, method = "Wald") #confidence intervals
r.squaredGLMM(lmer_attack) #r2


# Genus analysis
attack3 <-  attack2 %>%
  filter(n_genus >= 5, temp_range_genus >= 5)

lmer_attack_gen <-  lmer(log_attack ~ one_kT + log_mass + (one_kT|Genus), data = attack3 )
lmer_attack_gen_full <-  lmer(log_attack ~ one_kT + log_mass + (one_kT|Group/Genus), data = attack3 )

attack_genus <- attack3 %>%
  group_by(Genus) %>%
  slice_head()
# Results
rand(lmer_attack_gen) # AIC
rand(lmer_attack_gen_full ) # AIC
summary(lmer_attack_gen)
confint.merMod(lmer_attack_gen, method = "Wald") #confidence intervals
coef(lmer_attack_gen)
# dataframe for genus fits


attack_results <- as_tibble(tibble::rownames_to_column(coef(lmer_attack_gen)$Genus, "Genus") %>%
                                 rename(intercept = `(Intercept)`, temp_slope = one_kT, mass_slope = log_mass)) %>%
  arrange(Genus)

attack_results$Group <- attack_genus$Group[attack_genus$Genus %in% attack_results$Genus]
attack_results$min_one_kT <- attack_genus$min_one_kT
attack_results$max_one_kT <- attack_genus$max_one_kT
attack_results$min_y<- exp(attack_results$intercept + attack_results$temp_slope * attack_results$min_one_kT)
attack_results$max_y <- exp(attack_results$intercept + attack_results$temp_slope * attack_results$max_one_kT)
attack_results$Thermy<- "Ectotherm"
attack_results$Type<- "attack"



##---------------  Plotting Fig 2e, Part II-----------
## Attack and mortality plots were manually superimposed
attack_genus_regr <- ggplot(data = attack_results) +
  geom_segment(aes(x = min_one_kT, xend = max_one_kT, y = min_y, yend = max_y), 
               col = "navy",size = .5) +
  scale_y_continuous(name = expression(paste('Attack Rate (m'^2,' s'^-1,' g'^{-alpha},')')),
                     limits = c(10^-8.5, 10^-6.5), #note: same total range as mortality
                     labels = trans_format("log10", math_format(10^.x)),
                     breaks = c(10^-8, 10^-7, 10^-6), trans = "log10") +  
  scale_x_reverse(name = expression("1/kT"), limits = c(42.5, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, name = NULL)) +
  theme_single + theme(legend.position = "none") +theme_ticks 
grid.newpage()
grid.draw(attack_genus_regr )


attack_genus_regr <- set_panel_size(attack_genus_regr , width = unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(attack_genus_regr )

pdf('~/Desktop/Fig_2/attack_scaling_2E.pdf')
grid.draw(attack_genus_regr)
dev.off()

system2(command = "pdfcrop", 
        args  = c('~/Desktop/Fig_2/attack_scaling_2E.pdf', 
                  '~/Desktop/Fig_2/attack_scaling_2E.pdf'))


#showing data
attack_color <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', "black", "gray", "darkslategray1")
attack_mass_slope <- attack_results_all[2]
attack3$attack_mass_corr <- attack3$attack.rate/attack3$predator.mass.mg^attack_mass_slope
attack_data_regr <- ggplot() +
  geom_point(data = attack3, shape = 21, size = 2, stroke = 0.75,
             aes(x = one_kT, y = attack_mass_corr , color = Genus)) +
  geom_segment(data = attack_results, aes(x = min_one_kT, xend = max_one_kT, y = min_y, yend = max_y, color = Genus), 
               size = .5) +
  scale_color_manual(values = attack_color) +
  scale_y_continuous(name = expression(paste('Attack Rate (m'^2,' s'^-1,' g'^{-alpha},')')),
                    # limits = c(10^-8.5, 10^-6.5), #note: same total range as mortality
                     labels = trans_format("log10", math_format(10^.x)),
                     breaks = c(10^-8, 10^-7, 10^-6), trans = "log10") +  
  scale_x_reverse(name = expression("1/kT"), #limits = c(42.5, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, name = NULL)) +
  theme_single + theme(legend.position = "bottom") +theme_ticks 
grid.newpage()
grid.draw(attack_data_regr )


attack_genus_regr <- set_panel_size(attack_data_regr , width = unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(attack_data_regr)

pdf('~/Desktop/Fig_2/supp_attack_data.pdf')
grid.draw(attack_data_regr)
dev.off()

system2(command = "pdfcrop", 
        args  = c('~/Desktop/Fig_2/supp_attack_data.pdf', 
                  '~/Desktop/Fig_2/supp_attack_data.pdf'))

############################################
##  Fig 2F - violin plot
# Plot genus-level slopes of thermal sensitivity
############################################
mort_ecto_results$Type = "mortality"
mort_results

ecto_results <- as_tibble(rbind(attack_results, mort_ecto_results))
all_results <- rbind(attack_results, mort_results) %>%
  arrange(Type, Group, Genus) %>%
  select(Type, Group, Genus, everything())
write_csv(all_results, '~/Desktop/Fig_2/Ea_results.csv')

vio_mort_attack_ea <- ggplot(ecto_results , aes(x = Type, y = -1 * temp_slope)) +
  geom_violin(size = 1, aes(color = Type)) + 
  geom_jitter(aes(fill = Type), shape = 21, color = "black", width = 0.06, stroke = .6, size = 3) +
  scale_fill_manual(values =c("mortality" = ecto_mort_col, "attack" = "navy")) +
  scale_color_manual(values =c("mortality" = ecto_mort_col, "attack" = "navy")) +
  scale_y_continuous(position = "right", limits = c(0, 0.8), 
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8),
                     name=expression(paste('Thermal Sensitivity (E'[a],')')))+ 
  scale_x_discrete(name = NULL) +
  theme_single + theme(legend.position = "none")  + theme_ticks + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
vio_mort_attack_ea

# Save plot

vio_mort_attack_ea <- set_panel_size(vio_mort_attack_ea, width = unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(vio_mort_attack_ea)

pdf('~/Desktop/Fig_2/Fig2.pdf')
grid.draw(vio_mort_attack_ea)
dev.off()

system2(command = "pdfcrop", 
        args  = c('~/Desktop/Fig_2/Fig2.pdf', 
                  '~/Desktop/Fig_2/Fig2.pdf'))


###################################################################
############# Compare Mortality Variance  between Thermies #########
###################################################################

#Stats


# ecto vs endo mortality Ea range
range(mort_ecto_results$temp_slope)
range(mort_ecto_results$temp_slope)[1]- range(mort_ecto_results$temp_slope)[2] #0.705 magnitude of ecto E range

range(mort_endo_results$temp_slope)
range(mort_endo_results$temp_slope)[1]- range(mort_endo_results$temp_slope)[2] #0.180 magnitude of Endo E range - about the same

#standard dev  
sd(mort_ecto_results$temp_slope) 
sd(mort_endo_results$temp_slope) 

var_ecto <-var(mort_ecto_results$temp_slope)
var_ecto
var_endo <- var(mort_endo_results$temp_slope)
var_endo

#compare variance - not significantly different
sigma.test(mort_ecto_results$temp_slope, sigmasq = var_endo, alternative = "greater") #significantly greater variance in ectotherms




########################################################################
#------------ Extra analysis: Phylogenetic covariate ------------------
#######################################################################





#------- load output and view results
model_mort_mammal <- read_rds('/Users/jgradym/Google Drive/Gibert Paper/best_model/Mammal_fit.rds')

# Results
print(summary(model_mort_mammal), digits = 5) 
plot(model_mort_mammal, N = 6)

# r2
bayes_R2(model_mort_mammal, re.form = NA) # exclude random effects
bayes_R2(model_mort_mammal) # full model

# p value equivalent
pd_mamm <- p_direction(model_mort_mammal) # p direction
pd_to_p(pd_mamm$pd[2], "two-sided") # bayesian equivalent to p value for mass scaling
pd_to_p(pd_mamm$pd[3], "two-sided") # bayesian equivalent to p value for temp scaling


#--------------Generate output----------------
mam_tree_mcc <- read.tree('/Users/jgradym/Google Drive/Phylo_trees/mammal_mcc.tree') # to skip making
# Data
mort_mammal <- mortality1 %>%
  filter(Group == "Mammal") 

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


mort_mammal$species2 <- mort_mammal$species #random effect for multiple individuals per species

#----------------- bayesian model for mammals using brms, species and phylogeny as random effect -----------------

# variance and covariance from phylogeny
A <- ape::vcv.phylo(mam_tree_mcc)

# Priors
get_prior(formula = log_mort ~ log_mass + one_kT +
            (1|gr(species, cov = A)) + (1|species2), data = mort_mammal)


priors_endo1 <- c(set_prior("normal(0, 1)", class = "b", coef = "one_kT"), #slope = 0 for endotherms
                  set_prior("normal(-0.25, 1)", class = "b", coef = "log_mass")
                  
)

# run model
system.time(model_mort_mammal <- brm(
  log_mort ~ log_mass + one_kT + (1|gr(species, cov = A)) + (1|species2), 
  data = mort_mammal, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors_endo1,
  control = list(adapt_delta = 0.99),
  iter = 5000,
  cores = parallel::detectCores() -1) 
)

#
# Model Results
# r2
bayes_R2(model_mort_mammal, re.form = NA) # exclude random effects
bayes_R2(model_mort_mammal) # full model

# p value equivalent
pd_mamm <- p_direction(model_mort_mammal) # p direction
pd_to_p(pd_mamm$pd[2], "two-sided") # bayesian equivalent to p value for mass scaling
pd_to_p(pd_mamm$pd[3], "two-sided") # bayesian equivalent to p value for temp scaling

#saveRDS(model_mort_mammal, file = "/Users/jgradym/Google Drive/Gibert Paper/best_model/Mammal_fit.rds")



#------------------- Birds ---------------------
# ----------- Birds ------------------


model_mort_bird <- read_rds("/Users/jgradym/Google Drive/Gibert Paper/best_model/bird_fit.rds")

# Results
print(summary(model_mort_bird), digits = 5) 
plot(model_mort_bird, N = 6)

# r2
bayes_R2(model_mort_bird, re.form = NA) # exclude random effects
bayes_R2(model_mort_bird) # full model

# p value equivalent
pd_bird <- p_direction(model_mort_bird) # p direction
pd_to_p(pd_bird$pd[2], "two-sided") # bayesian equivalent to p value for mass scaling
pd_to_p(pd_bird$pd[3], "two-sided") # bayesian equivalent to p value for temp scaling

plot(model_mort_bird, N = 6)

#------ Generate bird output ----------------
mort_bird <- mortality1 %>%
  filter(Group == "Bird") #%>%

bird_tree_mcc  <- read.tree('/Users/jgradym/Google Drive/Phylo_trees/bird_tree_mcc.tree') # to skip making

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

priors_endo1 <- c(set_prior("normal(0, 1)", class = "b", coef = "one_kT"), #slope = 0 for endotherms
                  set_prior("normal(-0.25, 1)", class = "b", coef = "log_mass")
                  
)
# Run model
system.time(model_mort_bird <- brm(
  log_mort ~ log_mass + one_kT + (1|gr(species, cov = A)) + (1|species2), 
  data = mort_bird, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors_endo1,
  iter = 5000,
  cores = parallel::detectCores() -1) 
) 
#saveRDS(model_mort_bird, file = "/Users/jgradym/Google Drive/Gibert Paper/best_model/bird_fit.rds")
#----------------------------------------
# Model Results
print(summary(model_mort_bird), digits = 6) 
bayes_R2(model_mort_bird, re.form = NA) 
bayes_R2(model_mort_bird) 
 # fixed effects only

pd_bird <- p_direction(model_mort_bird) # p direction
pd_to_p(pd_bird$pd[2], "two-sided") # Bayesian equivalent to p value for mass scaling
pd_to_p(pd_bird$pd[3], "two-sided") # Bayesian equivalent to p value for temp scaling

pdf('~/Desktop/bird_bayes_plot.pdf')
plot(model_mort_bird, N = 6)
dev.off()


 #----------- Ectotherm mortality 

# ------------ to load results ----------
model_mort_ecto <- read_rds("/Users/jgradym/Google Drive/Gibert Paper/best_model/model_mort_ecto.rds")
print(summary(model_mort_ecto), digits = 5) #main effects haverhat near 1
coef(model_mort_ecto)$Group 

bayes_R2(model_mort_ecto, re.form = NA) # exclude random effects
bayes_R2(model_mort_ecto) # full model

# p value equivalent
pd_mort_ecto <- p_direction(pd_mort_ecto ) # p direction
pd_to_p(pd_mort_ecto $pd[2], "two-sided") # bayesian equivalent to p value for mass scaling

#----------- generate results -----------------
mort_ecto <- mortality1 %>%
  filter(Group == "Fish" | Group == "Invertebrate") 

# Model
#------------------- Ectotherm Mortality -----------------
priors_ecto <- c(set_prior("normal(-0.65, 1)", class = "b", coef = "one_kT"), #slope = 0 for endotherms
                 set_prior("normal(-0.25, 1)", class = "b", coef = "log_mass")
                 )
priors_ecto2 <- c(set_prior("normal(-0.65, 1)", class = "b", coef = "one_kT")) #slope = 0 for endotherms



mort_ecto <- mortality1 %>%
  filter(Group == "Fish" | Group == "Invertebrate") 

#
system.time(model_mort_ecto <- brm(
  log_mort ~ log_mass + one_kT + (log_mass + one_kT|Group) + (1|species), 
  data = mort_ecto, 
  family = gaussian(),
  prior = priors_ecto,
  iter = 50000,
  control = list(adapt_delta = 0.99999),
  cores = parallel::detectCores() -1) 
)  

#Results
print(summary(model_mort_ecto), digits = 5) #fixed coefficents have low rhat, more uncertainty in 1|species
coef(model_mort_ecto)$Group
plot(model_mort_ectol, N = 3)
bayes_R2(model_mort_ecto) #full model
bayes_R2(model_mort_ecto, re.form = NA) # exclude random effects
pd_mort_ecto <- p_direction(model_mort_ecto) # p direction
pd_to_p(pd_mort_ecto$pd[2], "two-sided") # bayesian equivalent to p value for mass scaling
pd_to_p(pd_mort_ecto$pd[3], "two-sided") # bayesian equivalent to p value for temp scaling




#--------------- ---------------
# ecto genus 
priors_ecto2 <- c(set_prior("normal(-0.65, 1)", class = "b", coef = "one_kT")) #slope = 0 for endotherms



mort_ecto_genus <- mortality %>%
  filter(n_genus >= 5, temp_range_genus >= 5, Group == "Fish" | Group == "Invertebrate" ) 
mort_ecto_genus$log_mass_corr_mort <- log(mort_ecto_genus$mass_corr_mortality)
#
system.time(model_mort_ecto_genus  <- brm(
  log_mass_corr_mort ~  one_kT + (one_kT|Group/Genus), 
  data = mort_ecto_genus , 
  family = gaussian(),
  prior = priors_ecto,
  iter = 100000,
  control = list(adapt_delta = 0.99999))
  #cores = parallel::detectCores() -1) 
)  
##########################################
 #### Panel 2a
##########################################

bird_lab  <- expression(paste('y = -0.022x + 0.19')) #from model summary
mamm_lab  <- expression(paste('y = 0.0050x - 0.45'))

# for manual plotting, results from Bayesian mixed effect analysis with phylogeny
bird_fit <- data.frame(x = c(min(mort_bird$one_kT),  max(mort_bird$one_kT)), 
                     y = c(exp((bird_intercept  + (bird_temp_slope * min(mort_bird$one_kT)))), 
                     exp(bird_intercept  + (bird_temp_slope * max(mort_bird$one_kT)))),
                     Group = c("Bird", "Bird"))

mammal_fit <- data.frame(x = c(min(mort_mammal$one_kT),  max(mort_mammal$one_kT)), 
                         y = c(exp((mammal_intercept + (mammal_temp_slope * min(mort_mammal$one_kT)))), 
                         exp(mammal_intercept + (mammal_temp_slope*  max(mort_mammal$one_kT)))),
                         Group = c("Mammal", "Mammal"))

endo_fit <- data.frame(rbind(bird_fit, mammal_fit))
endo_fit 


#------------- Plot Fig 2A  -----------------

mortal_plot_endo <- ggplot(mortality %>% filter(Group == "Mammal" | Group == "Bird") %>%
                             arrange(Group),
                           aes(x = one_kT, y = mass_corr_mortality, shape = Group)) + 
  theme_single + theme(legend.position = "none") +  theme_ticks +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{-alpha},')')), 
                     trans = "log", breaks = c(0.1, 1, 10), labels = c("0.1", "1", "10"), 
                    limits = c(0.02, 40)
                     ) +  
  scale_x_reverse( name = NULL, breaks = c(46,44, 42,40, 38),
                   limits = c(45.1, 38),
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


ggsave(mortal_plot_endo, height = 4.6, width = 6.7, filename = '~/Desktop/Fig_2/Fig2A.pdf')

# Save plot
ggsave(mortal_plot_endo, height = 4.6, width = 6.7, filename = file.path(github_path, "NAME YOUR FILE.pdf"))


###########################################
#### Panel 2b - species as random effect
##########################################

# Prepare plot
invert_lab  <- expression(paste('y = -0.55x + 22.6'))
fish_lab  <- expression(paste('y = -0.53x + 22.0'))

# from models summary - fit from model
mort_invert <-  mortality[mortality$Group == "Invertebrate",]
invert_fit <- data.frame(x = c(min(mort_invert$one_kT),  max(mort_invert$one_kT)), 
                       y = c(exp((invert_intercept+ (invert_temp_slope * min(mort_invert$one_kT)))), 
                             exp((invert_intercept + (invert_temp_slope  *  max(mort_invert$one_kT))))),
                       Group = c("Invertebrate", "Invertebrate"))

mort_fish <-  mortality[mortality$Group == "Fish",]
fish_fit <- data.frame(x = c(min(mort_fish$one_kT), max(mort_fish$one_kT)), 
                         y = c( exp((fish_intercept + (fish_temp_slope * min(mort_fish$one_kT)))), 
                                exp((fish_intercept + (fish_temp_slope *  max(mort_fish$one_kT))))),
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
  scale_x_reverse(name = expression("1/kT"),
                  limits = c(43, 38.1), 
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
                     trans = "log", 
                     #limits = c(.15, 10),
                     breaks = c(0.3, 1, 3,10), labels = c("0.3", "1", "3", "10")) + 
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
# Thermy
endo_genus <- mortality %>%
  filter(Thermy == "Endotherm") %>%
  pull(Genus) %>%
  unique()

# get parameter coefficients for plotting
genus_mort_lm <- mortality %>%
  filter(n_genus >= 5, temp_range_genus >= 5) %>%
  nest(-Genus) %>% 
  mutate(
    fit = map(data, ~ lm(log(mass_corr_mortality) ~ one_kT, data = .x)),
    tidied = purrr::map(fit, tidy, conf.int = T)
  ) %>% 
  unnest(tidied) 
genus_mort_lm$Thermy<-if_else(genus_mort_lm$Genus %in% endo_genus, "Endotherm", "Ectotherm" )
genus_mort_lm

# output results
genus_mort_lm2 <- genus_mort_lm %>%
  filter(term == "one_kT") %>%
  dplyr::select(-data, -fit) %>%
  dplyr::select(Thermy,Genus, everything())
genus_mort_lm2
write_csv(genus_mort_lm2, "~/Desktop/Fig_2/genus_regressions.csv")


# compare approaches 
mortality_5C <- mortality %>% filter(n_genus >= 5, temp_range_genus >= 5)

# linear model, genera calculated separately
mort_lm1a <- mortality_5C  %>%
  nest(-Genus) %>% 
  mutate(
    fit = map(data, ~ lm(log_mort ~ one_kT +log_mass  , data = .x)),
    tidied = purrr::map(fit, tidy, conf.int = T)
  ) %>% unnest(tidied) %>%
  filter(term == "one_kT") %>% 
  select(-data, -fit) 
mort_lm1a

mort_group <- mortality_5C %>% 
  group_by(Genus) %>%
  slice_head() %>%
  select(Genus, Group) 

mort_lm1 <- mort_lm1a %>% left_join(mort_group, by = "Genus") %>%
  arrange(Group, Genus) %>%
  select(Group, Genus, estimate) %>%
  rename(lm_estimate_separate = estimate)

accipter <- mortality_5C %>% filter(Genus == "Accipiter")
circus <- mortality_5C %>% filter(Genus == "Circus")
lm_acc <- lm(log_mort ~ log_mass + one_kT, data = accipter)
print(summary(lm_acc), digits = 8)
lm_circ <- lm(log_mort ~   one_kT + log_mass, data = circus)
summary(lm_circ )
### mixed model
mort_lmer <-  lmer(log_mort~ one_kT +log_mass  +
                     (log_mass + one_kT|Group:Genus),
                   data = mortality_5C )

mort_lmer_results0 <- tibble::rownames_to_column(coef(mort_lmer)$`Group:Genus`,"Genus")
mort_lmer_results <- tidyr::separate(mort_lmer_results0,Genus, sep = ":", into = c("Group", "Genus")) %>%
  select(Group, Genus, one_kT) %>%
  arrange(Group, Genus) %>%
  rename(lmer_estimate = one_kT)
mort_lmer_results  

mort_lmer2 <-  lmer(log_mort~ one_kT +log_mass  +
                     (one_kT|Group:Genus),
                   data = mortality_5C )
mort_lmer_results20 <- tibble::rownames_to_column(coef(mort_lmer2)$`Group:Genus`,"Genus")
mort_lmer_results2 <- tidyr::separate(mort_lmer_results20,Genus, sep = ":", into = c("Group", "Genus")) %>%
  select(Group, Genus, one_kT) %>%
  arrange(Group, Genus) %>%
  rename(lmer_estimate = one_kT)
mort_lmer_results2  

mort_lmer3 <-  lmer(log_mort~ one_kT +log_mass  +
                      (one_kT|Genus),
                    data = mortality)
summary(mort_lmer3)
coef <- coef(mort_lmer3)
str(coef)
coef(mort_lmer3)$Group
coef(mort_lmer3)$Genus
mort_lmer_results30 <- tibble::rownames_to_column(coef(mort_lmer)$`Group:Genus`,"Genus")
mort_lmer_results3 <- tidyr::separate(mort_lmer_results30,Genus, sep = ":", into = c("Group", "Genus")) %>%
  select(Group, Genus, one_kT) %>%
  arrange(Group, Genus) %>%
  rename(lmer_estimate = one_kT)
mort_lmer_results3  


write_csv(mort_lmer_results, "~/Desktop/Fig_2/mort_lmer_results.csv")

# combine
combine <- as.data.frame(cbind(mort_lmer_results, mort_lmer_results2,mort_lm1[,3]))
write_csv(combine, "~/Desktop/Fig_2/mort_genus_results.csv")



mort_lm2a <- lm(log_mort ~ one_kT * Genus + log_mass, data = mortality_5C)
mort_lm2b <- as.data.frame(coef(mort_lm2a))
mort_lm2c <- tibble::rownames_to_column(mort_lm2b, "Coefficient")
mort_lm2 <- separate(mort_lm2c, Coefficient, sep = ":", into = c("Genus", "term" )) %>%
  filter(term == "one_kT") %>%
  rename(lm_estimate_interaction = "coef(mort_lm2a)")
mort_lm2$Genus <- str_remove_all(mort_lm2$Genus, "Genus")
mort_lm2 <- mort_lm2 %>% 
  left_join(mort_group, by = "Genus") 
mort_lm2  
mort_lmer_full <-  lmer(log_mort~ log_mass + one_kT +
                     (log_mass + one_kT|Group) + (log_mass + one_kT|Group:Genus),
                   data = mortality_5C )

mort_lmer <-  lmer(log_mort~ log_mass + one_kT +
                    (log_mass + one_kT|Group:Genus),
                   data = mortality_5C )

mort_lm <-  lm(log_mort~ log_mass + one_kT +
                     (log_mass + one_kT|Group:Genus),
                   data = mortality_5C )






mort_vio_Ea <- ggplot(genus_mort_lm %>% filter(term == "one_kT"), 
                      aes(x = Thermy, y = -1*estimate))+
  geom_violin(size = 1, aes(color = Thermy)) + 
  geom_jitter(aes(fill = Thermy),shape = 21, size = 3, color = "black", width = 0.08, stroke = .5 ) +
  scale_y_continuous(position = "right",
                     limits = c(-0.6, 1.1),
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
                     trans = "log10", breaks = c(0.1, 1, 10), 
                     limits = c(0.1, 10), 
                     labels = c("0.1", "1", "10")) +  
  scale_x_reverse(name =NULL, limits = c(42.5, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, name = NULL)) + 
  theme_single + theme_ticks +theme(legend.position = "none")
grid.draw(mortal_scaling )
# Note: range of y limits is equal to attack rates below - can observe difference in spread of intercepts

# Save plot
ggsave(mortal_scaling, height = 4.7, width = 6.7, filename = '~/Desktop/mortal_scaling.pdf')

ggsave(mortal_regress, height = 4.7, width = 6.7, filename = file.path(github_path,'NAME YOUR FILE.pdf'))

##---------------------------------
## PART 2 (attack)

# Read in data
attack <- read_csv(file.path(github_path, 'Lietal_oikos_2017_data.csv'))
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


# linear and mixed model
# compare approaches 
attack_5C <- attack %>% filter(n_genus >= 5, temp_range_genus >= 5)

# linear model, genera calculated separately
attack_lm1a <- attack_5C  %>%
  nest(-Pred_Genus) %>% 
  mutate(
    fit = map(data, ~ lm(log_attack ~ log_mass + one_kT, data = .x)),
    tidied = purrr::map(fit, tidy, conf.int = T)
  ) %>% unnest(tidied) %>%
  filter(term == "one_kT") %>% 
  select(-data, -fit) 
attack_lm1a

attack_group <- attack_5C %>% 
  group_by(Pred_Genus) %>%
  slice_head() %>%
  select(Pred_Genus, predator.ana.group )  

attack_lm1 <- attack_lm1a %>% left_join(attack_group, by = "Pred_Genus") %>%
  arrange(predator.ana.group, Pred_Genus) %>%
  select(predator.ana.group, Pred_Genus, estimate) %>%
  rename(lm_estimate_separate = estimate)
attack_lm1

### mixed model
library(lmerTest)
attack_lmer <-  lmer(log_attack~ one_kT + log_mass + ( one_kT|predator.ana.group:Pred_Genus),
                   data = attack_5C )
attack_lmer2 <-  lmer(log_attack~ one_kT + log_mass + (one_kT|predator.ana.group/Pred_Genus),
                     data = attack_5C )
attack_lmer3 <-  lmer(log_attack~ one_kT + log_mass + (one_kT|Pred_Genus),
                      data = attack_5C )
attack_lmer4<-  lmer(log_attack~ one_kT + log_mass + (one_kT|predator.ana.group) +
                        ( one_kT|predator.ana.group:Pred_Genus),data = attack_5C )
attack_lmer5<-  lmer(log_attack~ one_kT + log_mass + (one_kT|predator.ana.group) +
                       ( one_kT|Pred_Genus),data = attack_5C )
attack_lm<-  lm(log_attack~ one_kT + log_mass + predator.ana.group +
                        one_kT*Pred_Genus,data = attack_5C )
coef(attack_lm)
coef(attack_lmer3)
coef(attack_lmer)
coef(attack_lmer2)
coef(attack_lmer4)
coef(attack_lmer5)
library(jtools)
summ(attack_lmer,confint = T, digits = 4)
rand(attack_lmer)
coef(attack_lmer2)
plotREsim(REsim(attack_lmer, n.sims = 500))
library(parameters)
model_parameters(attack_lmer), robust = TRUE, vcov_estimation = "CL", vcov_type = "HC1")

library(merTools)
randomSims <- REsim(attack_lmer, n.sims = 500)
FEsim(attack_lmer, n.sims = 200, oddsRatio = FALSE, seed = NULL)
confint.merMod(attack_lmer , method = "Wald", digits = 4)
r.squaredGLMM(attack_lmer )
attack_lmer_results0 <- tibble::rownames_to_column(coef(attack_lmer)$`predator.ana.group:Pred_Genus`,"Genus")
attack_lmer_results <- tidyr::separate(attack_lmer_results0, Genus, sep = ":", into = c("Group", "Genus")) %>%
  select(Group, Genus, one_kT) %>%
  arrange(Group, Genus) %>%
  rename(lmer_estimate = one_kT)
attack_lmer_results  
write_csv(attack_lmer_results, "~/Desktop/Fig_2/attack_lmer_results.csv")

# combine
attack_combine<- as.data.frame(cbind(attack_lmer_results, attack_lm1[,3]))
attack_combine
write_csv(attack_combine, "~/Desktop/Fig_2/attack_genus_results.csv")


# together 
model_ecto_attack <- read_rds('/Users/jgradym/Google Drive/Gibert Paper/best_model/model_ecto_attack.rds')

#------- can skip -------
system.time(model_ecto_attack <- brm(
  log_attack  ~ log_mass + one_kT + (log_mass + one_kT|predator.ana.group) +(1|predator.species), 
  data = attack, 
  family = gaussian(), 
  prior = priors_ecto,
  iter = 100000,
  control = list(adapt_delta = 0.99999),
  cores = parallel::detectCores() -1) 
) 
attack2$log_mass_corr_attack <- log(attack2$mass_corr_attack)
attack2 <- attack2 %>% filter(n_genus >= 5, temp_range_genus >= 5)

attack_lmer <- lmer(log_attack ~ log_mass + one_kT +
                      (log_mass + one_kT|predator.ana.group:Pred_Genus),
                    data = attack2)
attack_genus_results0 <- tibble::rownames_to_column(coef(attack_lmer)$`predator.ana.group:Pred_Genus`,"Genus")
attack_genus_results <- tidyr::separate(attack_genus_results,Genus, sep = ":", into = c("Group", "Genus"))
attack_genus_results
str(attack_genus_results0 )
write_csv(attack_genus_results, "~/Desktop/Fig_2/attack_mixed.csv")


system.time(model_ecto_attack_genus <- brm(
  mass_corr_attack ~ one_kT + (one_kT|predator.ana.group/Pred_Genus), 
  data = attack2, 
  family = gaussian(), 
  prior = priors_ecto2,
  iter = 100000,
  control = list(adapt_delta = 0.99999),
  cores = parallel::detectCores() -1) 
) 
saveRDS(model_ecto_attack_genus , "/Users/jgradym/Google Drive/Gibert Paper/model_ecto_attack_genus.rds")
#------- can skip -------

print(summary(model_ecto_attack ), digits = 5)
saveRDS(model_ecto_attack, "/Users/jgradym/Google Drive/Gibert Paper/model_ecto_attack.rds")
coef(model_ecto_attack)$predator.ana.group
bayes_R2(model_vert_attack, re.form = NA) # exclude random effects
bayes_R2(model_vert_attack) #full model
pd_vert_attack <- p_direction(model_vert_attack) # p direction
pd_to_p(pd_vert_attack$pd[2], "two-sided") # bayesian equivalent to p value for mass scaling
pd_to_p(pd_vert_attack$pd[3], "two-sided") # bayesian equivalent to p value for temp scaling

pdf('~/Desktop/ecto_attack_bayes_plot.pdf')
plot(model_ecto_attack, N = 6)
dev.off()


######33
# mass slopes, intercept
attack_coef <- coef(model_ecto_attack )$predator.ana.group
invert_attack_intercept <- attack_coef[1]
invert_attack_mass_slope <- attack_coef[9]
invert_attack_temp_slope <- attack_coef[17]

fish_attack_intercept <-  attack_coef[2]
fish_attack_mass_slope <- attack_coef[10]
fish_attack_temp_slope <- attack_coef[18]


#------- Add mass-corrected mortality for plotting ----

attack2 <- attack %>% 
  mutate(mass_corr_attack= NA)

invert <- attack2[attack2$predator.ana.group == "invertebrate",]
attack2$mass_corr_attack[attack2$predator.ana.group == "invertebrate"] <- invert$attack.rate/invert$predator.mass.mg^invert_attack_mass_slope  

fish <- attack2[attack2$predator.ana.group == "vertebrate",]
attack2$mass_corr_attack[attack2$predator.ana.group == "vertebrate"] <- fish$attack.rate/fish$predator.mass.mg^fish_attack_mass_slope




# get attack regressions with temperature
attack_genus_lm0 <- attack2 %>%
  filter(n_genus >= 5, temp_range_genus >= 5) %>%
  nest(-Pred_Genus) %>% 
  mutate(
    fit = map(data, ~ lm(log(mass_corr_attack) ~ one_kT, data = .x)),
    tidied = purrr::map(fit, tidy, conf.int = T)
  ) %>% unnest(tidied) %>%
  filter(term == "one_kT") %>% 
  select(-data, -fit)
attack_genus_lm0 
write_csv(attack_genus_lm0 , "~/Desktop/Fig_2/attack_slopes.csv")

attack_genus_lm0 
attack_genus_lm0
attack_genus_lm <- attack2 %>%
  filter(n_genus >= 5, temp_range_genus >= 5) %>%
  nest(-Pred_Genus) %>% 
  mutate(
    fit = map(data, ~ lm(log(mass_corr_attack) ~ one_kT, data = .x)),
    tidied = purrr::map(fit, tidy, conf.int = T),
    pred = map(fit, broom.mixed::augment)
  ) 
attack_genus_lm
#Add group names
attack_genus_lm$Group <- unlist(map(attack_genus_lm$data, ~ .x[["predator.ana.group"]][[1]])) 
attack_genus_lm <- unnest(attack_genus_lm,  pred)%>% 
  rename(log_attack_mass_corr = `log(attack_mass_corr)`)
attack_genus_lm
