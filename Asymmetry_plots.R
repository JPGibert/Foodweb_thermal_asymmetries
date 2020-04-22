################################XXX XXXXX, Oct 14 2019 ################################## 

library(tidyverse)
library(egg)
library(grid)
library(gridExtra)
library(gtable)
library(ggpubr)
library(broom)
library(scales)
library(cvcqv)
library(TeachingDemos)

# Add your own google drive path to the shared folder to read in data and save plots
gdrive_path <- file.path('/Users/XXXX/Google Drive/XXXXXX Paper') 

invert_col <- "deepskyblue2"
fish_col <- "navy"
mamm_col <- "red"
bird_col <- "darkred"
ecto_col <- "dodgerblue"
endo_col <- "red2"

################################## Plotting ###########################33


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
# Read in data

attack <- read_csv(file.path(gdrive_path, 'Data/Lietal_oikos_2017_data.csv'))
attack$one_kT <- 1/(8.617e-5*(attack$temperature.degree.celcius +273.15))

mortality1 <- read_csv(file.path(gdrive_path,'Data/McCoy_Gillooly_Data.csv')) #2,117 values
mortality1[is.na(mortality1$Species),]

######################################################
############## Mortality Data first
######################################################


#
# All slopes, calculated separately
##########################################################################

###################################  Add endo temp data ################################## 
  # This temp determined from matching lat/lon centroid of species range to climate data from ...
bird_temp <- read_csv(file.path(gdrive_path,"Data/bird_temp.csv"))
mammal_temp <- read_csv(file.path(gdrive_path,"Data/mammal_temp.csv"))
mammal_temp2 <- mammal_temp %>%
  dplyr::select(Species, Mean_Temp_C) 
 # use new endo data and old ecto data in column called 'Temp_C_good'
mortality2 <- mortality1 %>%
  left_join(mammal_temp, by = "Species") %>%
  left_join(bird_temp, by = "Species") %>%
  mutate(new_temps = coalesce(Mean_Temp_C.x,  Mean_Temp_C.y)) %>%
  mutate(Temp_C_good =  if_else(Group == "Mammal" | Group == "Bird", new_temps, Temp_C))
mortality2 #2,103 spp
unique(mortality2$Group)


# Create 1/kT column; 
mortality2$one_kT <- 1/(8.617e-5*(mortality2$Temp_C_good +273.15)) # my default temp
mortality2$mortality_corr <- mortality2$Mortality_yr*mortality2$dry_mass_g^0.25


##################### remove NA's, NAN, Inf's for analysis ################################## 
mortality2.1 <- mortality2[is.na(log(mortality2$mortality_corr)) == F,]
mortality3 <- mortality2[is.nan(log(mortality2$mortality_corr)) == F,]
mortality4 <- mortality3[is.infinite(log(mortality3$mortality_corr)) == F,]
mortality4 <- mortality4 %>%
  dplyr::select(Group, Species, dry_mass_g,  one_kT, mortality_corr, Temp_C_good,  Mortality_yr, one_kT_mccoy, Temp_C)

# Get number measurements per species or per genus, if you want to restrict by individuals per species or genus:, 
mortality5 <- mortality4 %>% #firs
  mutate(Genus = word(Species,1)) %>%
  dplyr::select(1:10) %>%
  group_by(Genus) %>% 
  add_count(Genus, name = "n_genus") %>%# individuals per genus
  mutate(temp_range_genus = max(Temp_C_good) - min(Temp_C_good)) %>%
  ungroup()

#Add a temperature range for individuals in a species (only relevant for ectos)
mortality6 <- mortality5 %>%
  group_by(Species) %>% 
  add_count(Species, name = "n_species") %>%
  mutate(temp_range_spp = max(Temp_C_good) - min(Temp_C_good)) %>%
  ungroup()

mortality7 <- mortality6%>%
  filter(Group == "Bird" | Group == "Mammal"|  Group == "Invertebrate" | Group == "Fish") %>%
  mutate(Thermy = if_else(Group == "Bird" | Group == "Mammal", "Endotherm", "Ectotherm"))
# Give a name for final dataset
mortality_final <- mortality7
mortality_final_tax <- mortality7 %>%
  dplyr::select(Genus, Thermy)
mortality <- mortality_final
# get regression line info 

mamm <- mortality[mortality$Group == "Mammal",]
bird <- mortality[mortality$Group == "Bird",]
lm_mamm <- lm(log(mortality_corr) ~ one_kT, data = mamm)
summary(lm_mamm) #0.08935, Adjusted R-squared:  0.03642
lm_bird <- lm(log(mortality_corr) ~ one_kT, data = bird)
summary(lm_bird) #-0.0243, Multiple R-squared:  0.002222


bird_lab  <- expression(paste('y = -0.089x, r'^2,' = 0.036'))
mamm_lab  <- expression(paste('y = 0.024x, r'^2,' = 0.0022'))

##########################################
 #### Panel 2a
##########################################
mortal_plot_endo <- ggplot(mortality %>% 
                             filter(Group == "Mammal" | Group == "Bird") %>%
                             arrange(Group),
                           aes(x = one_kT, y = mortality_corr, shape = Group)) + 
  geom_point(size = 2.5, shape = 21, color = "black", 
             stroke = .5, aes(fill = Group)) + 
  geom_smooth(method = 'lm', size = 1.25, alpha = 0.2,
              aes(fill = Group,color = Group)) +
  scale_fill_manual(values =c("Mammal" = mamm_col, "Bird" = "darkred")) +
  scale_color_manual(values =c("Mammal" = mamm_col, "Bird" = "brown4")) +
  
  #scale_color_manual(values =c("Mammal" = "red1", "Bird" = "brown4")) +
  #scale_fill_manual(values =c("Mammal" = "red", "Bird" = "brown3")) +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{1/4},')')), 
                     trans = "log10", breaks = c(0.1, 1, 10), 
                     labels = c("0.1", "1", "10"), limits = c(0.05, 20)) +  
  scale_x_reverse(limits = c(44.8, 38), name = NULL, breaks = c(46,44, 42,40, 38),
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, breaks = c(-20, -10, 0, 10, 20, 30),
                                      name = "Ambient Temperature (ºC)")) +
  guides(color=guide_legend(override.aes=list(fill=NA))) + theme(legend.position = "none")+
  theme_single + annotate("text", x = 44.8, size = 6, y = 10, hjust = 0, 
                          fontface = "bold",label = "Bird", color = bird_col) +
  annotate("text", x = 44.8, size = 6, y = 16, 
           hjust = 0, fontface = "bold",label = "Mammal", color = mamm_col) +
  annotate("text", x = 41.8, size = 6, y = .1, 
           hjust = 0, label = mamm_lab, color = "red1")+
  annotate("text", x = 41.8, size = 6, y = .06, 
           hjust = 0, label = bird_lab, color = "brown4") +
  theme_ticks
grid.newpage()
grid.draw(mortal_plot_endo)

# Save plot
ggsave(mortal_plot_endo, height = 4.6, width = 6.7, filename = file.path(gdrive_path, "Plots/Fig_2/Fig2a_endo_mortality.pdf"))


# all ectos combined (note have to decide on color, navy (fish) and darkred (birds) don't look as good but increase contrast, maybe worth it)
# Ecto Mortality
invert <- mortality[mortality$Group == "Invertebrate",]
fish <- mortality[mortality$Group == "Fish",]
lm_invert <- lm(log(mortality_corr) ~ one_kT, data = invert)
summary(lm_invert)#0.60, r2 =0.43
lm_fish<- lm(log(mortality_corr) ~ one_kT, data = fish)
summary(lm_fish) #0.53, r2 =0.47

invert_lab  <- expression(paste('y = -0.60x,   r'^2,' = 0.43'))
fish_lab  <- expression(paste('y = -0.53x, r'^2,' = 0.47'))

##########################################
#### Panel 2b
##########################################
mortal_plot_ecto <- ggplot(mortality %>% 
                             filter(Group == "Fish" | Group == "Invertebrate") %>%
                             arrange(Group),
                           aes(x = one_kT, y = mortality_corr,fill = Group)) + 
  geom_point(size = 2.5, shape = 21, color = "black", 
             stroke = .5, aes(fill = Group)) +
  geom_smooth(method = 'lm', size = 1.25, alpha = 0.25,
              aes(fill = Group, color = Group)) +
  #scale_fill_manual(values =c("Fish" = "royalblue2", "Invertebrate" = "deepskyblue")) +
  scale_fill_manual(values =c("Fish" = fish_col, "Invertebrate" =  invert_col)) +
  scale_color_manual(values =c("Fish" = fish_col, "Invertebrate" =  invert_col)) +
  
  #scale_color_manual(values =c("Fish" = "royalblue3", "Invertebrate" = "lightskyblue")) +
  scale_y_continuous(trans = "log10",  name = expression(paste('Mortality (yr'^-1,' g'^{1/4},')')),
                     breaks = c(0.1,  1, 10), labels = c("0.1", "1",  "10"), 
                     limits = c(0.03, 30), position = "right")+ 
  scale_x_reverse(name = expression("1/kT"), limits = c(43, 38.2), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, 
                                      #name = NULL))+
                                      name = "Ambient Temperature (ºC)")) +
  theme(legend.position = "none") +
  theme_single + annotate("text", x = 43, size = 6, y = 14, hjust = 0, 
                          fontface = "bold",label = "Fish", color = fish_col) +
  annotate("text", x = 43, size = 6, y = 25, hjust = 0, 
           fontface = "bold",label = "Invertebrate", color = invert_col) +#+
  annotate("text", x = 40.7, size = 6, y = .08, hjust = 0, label = invert_lab, color = invert_col)+
  annotate("text", x = 40.7, size = 6, y = .04, hjust = 0, label = fish_lab, color = fish_col) +
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text.x = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")), 
        axis.text.x.top = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")), 
        axis.text.y.right = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")),
        axis.text.y.left = element_text(margin=unit(c(0.4,0.4,0.4,0.4), "cm")))
#theme_no_y
grid.newpage()
grid.draw(mortal_plot_ecto) 

# Save plot
ggsave(mortal_plot_ecto, height = 5, width = 7, filename = file.path(gdrive_path, 'Plots/Fig_2/Fig2b_ecto_mortality2.pdf'))


# Restrict to 5 individuals per genus, 5 temp C range
mort_5 <- mortality_final %>% 
  filter(n_genus >= 5) %>%
  filter(temp_range_genus >= 5)

##########################################
#### Panel 2c, Mortality regressions
##########################################
# 5 indiv per genus, temp range of 5, showing group regression line too
##########################################

mortal_plot_5C_5ind <- ggplot(mort_5 %>% 
                                filter(Group == "Fish" | Group == "Invertebrate" |
                                         Group == "Mammal" | Group == "Bird") %>%
                                filter(n_genus >= 5) %>%
                                filter(temp_range_genus >= 5),
                              aes(x = one_kT, y = mortality_corr)) + 
  geom_smooth(method = 'lm', size = 0.3,  alpha = 0,
              aes(group = Genus, color = Group)) +
  geom_smooth(method = 'lm', size = 1.5,  alpha = 0,
              aes(group = Group, color = Group)) +
  scale_color_manual(values =c("Invertebrate" = invert_col, "Fish" = fish_col, 
                               "Mammal" = mamm_col, "Bird" = bird_col)) +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{1/4},')')),
                     trans = "log10", breaks = c(0.3, 1, 3,10), labels = c("0.3", "1", "3", "10"), limits = c(0.3, 10)) +  
  scale_x_reverse(name =NULL, limits = c(44.4, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, 
                                      name = NULL)) + 
  #annotate("text", label = "Bird", hjust = 0, color = bird_col, x = 44.8, y = 6.9, size = 6, fontface = "bold") +
  #annotate("text", label = "Mammal", hjust = 0,color = mamm_col, x = 44.8, y = 9, size = 6, fontface = "bold") +
  #annotate("text", label = "Invertebrate", hjust = 0,color =  invert_col, x = 44.8, y = 5.3, size = 6, fontface = "bold") +
  #annotate("text", label = "Fish", hjust = 0,color = fish_col, x = 44.8, y = 4, size = 6, fontface = "bold") +
  theme_single + theme(legend.position = "none") +theme_ticks 
grid.newpage()
grid.draw(mortal_plot_5C_5ind)

# Save plot
ggsave(mortal_plot_5C_5ind, height = 4.3, width = 6.3, filename = file.path(gdrive_path,'Plots/Fig_2/Fig2C_mortality_genus_5_indiv_5C_and_group.pdf'))
ggsave(mortal_plot_5C_5ind, height = 4.3, width = 6.3, filename = 'figfig.pdf')


##########################################
#### Panel 2d, Mortality violin plot
 #first get regressions
##########################################

#Restricted to 5 indiv per genus, 5 C range

mort_ecto <- mort_5 %>% 
  filter(Group == "Fish" | Group == "Invertebrate")  

# All slopes, calculated separately
mort_lm0 <- mort_5 %>%
  nest(-Group) %>% 
  mutate(
    fit = map(data, ~ lm(log(mortality_corr) ~ one_kT, data = .x)),
    tidied = map(fit, tidy, conf.int = T)
    #glanced = map(fit, glance) %>% filter(r.squared),
    #augmented = map(fit, augment)
  ) %>% 
  unnest(tidied)
mort_lm0

#get r2
mort_lm1 <- mort_5 %>%
  nest(-Group) %>% 
  mutate(
    fit = map(data, ~ lm(log(mortality_corr) ~ one_kT, data = .x)),
    #tidied = map(fit, tidy, conf.int = T),
    glanced = map(fit, glance)
    #augmented = map(fit, augment)
  ) %>% 
  unnest(glanced)
mort_lm1 

#combine
mort_lm <- mort_lm0 %>%
  filter(term != '(Intercept)') %>%
  mutate(r2 = mort_lm1$r.squared)# note: 'estimate' = slope
mort_lm


#### look within genus 
mort_lm2 <- mort_5  %>%
  nest(-Genus) %>% 
  mutate(
    fit = map(data, ~ lm(log(mortality_corr) ~ one_kT, data = .x)),
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied) 
mort_lm2 
# remove intercept
mort_lm2 <- filter(mort_lm2, term != '(Intercept)') # note: 'estimate' = slope
mort_lm2

# Add therym
mort_lm2 <- mort_lm2 %>%
  left_join(mortality_final_tax, by = "Genus") %>%
  distinct()
mort_lm2 

#B violin plot
mort_vio_Ea <- ggplot(mort_lm2, aes(x = Thermy, y = -1*estimate))+
  geom_violin(size = 1, aes(color = Thermy)) + geom_jitter(aes(fill = Thermy),shape = 21, size = 3, color = "black", width = 0.05, stroke = .5 ) +
  scale_y_continuous(position = "right",name=expression(paste('Thermal Sensitivity (E'[a],')')))+ 
  scale_x_discrete(name = NULL) +
  scale_color_manual(values =c("Ectotherm" = ecto_col, "Endotherm" = endo_col)) +
  scale_fill_manual(values =c("Ectotherm" = ecto_col, "Endotherm" = endo_col)) +
  theme_single + theme(legend.position = "none") + ggtitle("Mortality") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + theme_ticks
mort_vio_Ea

# Save plot
ggsave(mort_vio_Ea,  height = 4.22, width = 6.05, filename = file.path(gdrive_path,'Plots/Fig_2/Fig2d_vio_mort.pdf'))


####### Mortality regressions
mortal_plot_5C_5ind <- ggplot(mort_5 %>% 
                                filter(Group == "Fish" | Group == "Invertebrate")  %>%
                                filter(n_genus >= 5) %>%
                                filter(temp_range_genus >= 5),
                              aes(x = one_kT, y = mortality_corr)) + 
  geom_smooth(method = 'lm', size = 1,  alpha = 0, 
              aes(group = Genus, color = Group)) +

  scale_color_manual(values =c("Invertebrate" = invert_col, "Fish" = fish_col, 
                               "Mammal" = mamm_col, "Bird" = bird_col)) +
  #scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{1/4},')')), position = "right",
  #                   trans = "log10", breaks = c(0.3, 1, 3,10), labels = c("0.3", "1", "3", "10"), limits = c(0.5, 10)) +  
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{1/4},')')), position = "right",
                     trans = "log10", breaks = c(0.1, 10,  1000), limits = c(0.03, 3000), labels = c("0.1", "10", "1000")) +  
  scale_x_reverse(name =NULL, limits = c(42, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, 
                                      name = NULL)) + 

  theme_single + theme_ticks +theme(legend.position = "none")
grid.draw(mortal_plot_5C_5ind)

# Save plot
ggsave(mortal_plot_5C_5ind, height = 4.7, width = 6.7, filename = file.path(gdrive_path,'Plots/Fig_2/Fig2e_2v2_mort_by_genus_5_indiv_5C_overlay.pdf'))



############################################
##v  Fig 2F - violin plot
############################################

#################################  Fig 2e Attack Rates & Mortality Regressions  #################################

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

attack_5<- attack %>%
  filter(n_genus >= 5, temp_range_genus >= 5)

lm_attack_all <- lm(log(attack.rate) ~ one_kT, data = attack )
summary(lm_attack_all) #-0.77

lm_attack_all_spp <- lm(log(attack.rate) ~ one_kT + predator.species, data = attack )
summary(lm_attack_all_spp) #-0.43
lm_attack_all_genus <- lm(log(attack.rate) ~ one_kT + Pred_Genus, data = attack )
summary(lm_attack_all_genus) #-0.46

lm_attack_5<- lm(log(attack.rate) ~ one_kT + Pred_Genus, data =attack_5 )
summary(lm_attack_5) #y = -0.47x, r2 = 0.87


attack_lab  <- expression(paste('y = -0.4x, r'^2,' = 0.87'))

attack_lm <- attack_5 %>%
  nest(-Pred_Genus) %>% # the group variable
  mutate(
    fit = map(data, ~ lm(log(attack.rate) ~ one_kT, data = .x)), #this is a regression (y = attack.rate, x = one_kT), but you could adapt to whatever
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied) 
attack_lm


attack_lm <- filter(attack_lm, term != '(Intercept)') %>%
  mutate(Thermy = "Ectotherm") %>%
  mutate(Type = "Attack") %>%
  rename(Taxa = Pred_Genus)
attack_lm

mort_lm <- mort_lm2 %>%
  filter(Thermy == "Ectotherm") %>%
  #dplyr::select(-Thermy) %>%
  mutate(Type = "Mortality") %>%
  rename(Taxa = Genus)

mort_attack <-bind_rows(mort_lm, attack_lm)

vio_mort_attack_ea <- ggplot(mort_attack , aes(x = Type, y = -1*estimate))+
  geom_violin(size = 1, color = "dodgerblue", aes(linetype = Type)) + 
  scale_linetype_manual(values=c("dashed", "solid"))+
  geom_jitter(shape = 21, color = "black",fill = "dodgerblue", width = 0.06, stroke = .6, size = 3) +
  scale_y_continuous(limits = c(-0.5, 1), position = "right", name=expression(paste('Thermal Sensitivity (E'[a],')')))+ 
  scale_x_discrete(name = NULL) +
  #scale_color_manual(values =c("Ectotherm" = "skyblue3", "Endotherm" = "red2")) +
  theme_single + theme(legend.position = "none")  + theme_ticks + #+ggtitle("Ectotherm Thermal Sensitivity") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
vio_mort_attack_ea
ggsave(vio_mort_attack_ea,  height = 4.22, width = 5.9,  filename = file.path(gdrive_path,'Plots/Fig_2/Fig2f_vio_mort_attack_ea_dashed.pdf'))


##########################################
# Attack genus regressions
##########################################
#### Regressions within species, first add n per species and temp range


# Plot by species, at least 3 spp



attack_regr_genus <- ggplot(attack_5 %>% 
                              arrange(Pred_Genus),
                            aes(x = one_kT, y = attack.rate)) + 
  geom_smooth(method = 'lm', size = 1,  alpha = 0, linetype = "dashed",
              aes(group = Pred_Genus, color = predator.met.group)) +
  scale_color_manual(values =c("invert" = invert_col, "ectovert" = "navy"))+
  scale_y_continuous(name = expression(paste('Attack Rate (m'^2,' s'^-1,')')),
                     limits = c(10^-9, 10^-4), 
                     labels = trans_format("log10", math_format(10^.x)),
                     breaks = c(10^-11, 10^-9, 10^-7, 10^-5, 10^-3), trans = "log10") +  
  scale_x_reverse(name = expression("1/kT"), limits = c(42, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, 
                                      name = NULL)) +
  
  annotate("text", y = 10^-9.5, x = 40, hjust = 0, size = 6, label = NA, color = invert_col) +
  theme_single + theme_ticks + theme(legend.position = "none")

grid.newpage()
grid.draw(attack_regr_genus)
ggsave(attack_regr_genus , height = 4.7, width = 6.7, filename = file.path(gdrive_path,'Plots/Fig_2/Fig2e_1_attack_by_genus_5_indiv_5C.pdf'))

###################################################################
#############Compare Mortality Variance  between Thermies #########
###################################################################

#Stats
mort_lm_ecto <- mort_lm %>% filter(Thermy == "Ectotherm")
mort_lm_endo <- mort_lm2 %>% filter(Thermy == "Endotherm")

# ecto vs endo range
range(mort_lm_ecto$estimate)
range(mort_lm_ecto$estimate)[2] - range(mort_lm_ecto$estimate)[1] #1.33 magnitude of ecto E range

range(mort_lm_endo$estimate)
range(mort_lm_endo$estimate)[2] - range(mort_lm_endo$estimate)[1] #1.37 magnitude of Endo E range - about the same

#standard dev
sd(mort_lm_ecto$estimate) #0.3550812
sd(mort_lm_endo$estimate) #0.2647

var_ecto <-var(mort_lm_ecto$estimate)
var_ecto
var_endo <- var(mort_lm_endo$estimate)
var_endo
#CI of variance
sigma.test(mort_lm_ecto$estimate, sigmasq = var_endo)
# ecto vs endo range

#interquartile distance
IQR(mort_lm_ecto$estimate) #0.4150781 Interquartile range
IQR(mort_lm_endo$estimate) #0.225 Interquartile range

# CI's for interquartile range - found this on the web, doesn't seem to be very sensitive
# get IQR CI's with bootstrapping
library(boot)
iqr.fun <- function(data, indices) {
  
  d <- data[indices]  
  iqr <- IQR(d)  
  return(iqr)
}

set.seed(1000) # for reproducibility

ecto_iqr_boot <- boot(data = mort_lm_ecto$estimate, statistic = iqr.fun, R = 100000)

ecto_iqr_boot 
#Bootstrap Statistics :
# original      bias    std. error
#t1* 0.4150781 -0.01578602   0.1718928

ecto_mort_IQR_CIs <- boot.ci(ecto_iqr_boot )
ecto_mort_IQR_CIs
#Intervals : 
# Level      Normal              Basic         
#95%   ( 0.0940,  0.7678 )   ( 0.1299,  0.7574 )  

#Level     Percentile            BCa          
#95%   ( 0.0727,  0.7003 )   ( 0.0769,  0.8530 ) 

endo_iqr_boot <- boot(data = mort_lm_endo$estimate, statistic = iqr.fun, R = 100000)
#Bootstrap Statistics :
# original      bias    std. error
#t1* 0.4150781 -0.01578602   0.1718928
endo_iqr_boot 
#Bootstrap Statistics :
#  original      bias    std. error
#t1* 0.2258414 -0.01793671  0.08412654
endo_mort_IQR_CIs <- boot.ci(endo_iqr_boot )
endo_mort_IQR_CIs 
#Intervals : 
#  Level      Normal              Basic         
#95%   ( 0.0789,  0.4087 )   ( 0.0739,  0.3787 )  
# Confidence intervals



#coefficience of variation
cv_versatile(mort_lm_ecto$estimate, method = "All")
cv_versatile(mort_lm_endo$estimate, method = "All")
#kelley    -93.5 -304.8  -59.5   cv with Kelley 95% CI

#coefficience of quartile variation
cqv_versatile(mort_lm_ecto$estimate, method = "all")
cqv_versatile(mort_lm_endo$estimate, method = "all")

range(mort_lm_endo$estimate)
quantile(mort_lm_endo$estimate)



##########################################################################################
##########################################################################################
##########################################################################################
###################################  Mccoy & gillooly's Figure #1 ################################## 
# note: some species have multiple data points
mortal_plot_orig <- ggplot(mortality %>% 
                         filter(Group == "Fish" | Group == "Invertebrate") %>%
                         arrange(Group),
                       aes(x = one_kT_mccoy, y = mortality_corr)) + 
  geom_point(size = 2.5, stroke = .5, aes(shape = Group)) + 
  scale_shape_manual(values = c("Invertebrate" = 2, "Fish" = 19)) +
  geom_smooth(method = 'lm', size = 1,  alpha = 0.3, color = "black") +
  scale_y_continuous(trans = "log10",
                     name = expression(paste('Mortality (yr'^-1,' g'^{1/4},')'))) +  
  scale_x_continuous(name = expression("1/kT (1/T - 1/T20 ºC)"), limits = c(-2,4))+
  theme_single
grid.newpage()
grid.draw(mortal_plot_orig)

# The way I would do it; Ectos

mortal_plot_ecto <- ggplot(mortality %>% 
                       filter(Group == "Fish" | Group == "Invertebrate") %>%
                       arrange(Group),
                    aes(x = one_kT, y = mortality_corr, shape = Group)) + 
  geom_point(size = 2.5, stroke = .5, aes(shape = Group), color = "blue2") + 
  scale_shape_manual(values = c("Invertebrate" = 2, "Fish" = 19)) +
  geom_smooth(method = 'lm', size = 1,  alpha = 0.3, color = "darkblue") +
  scale_y_continuous(trans = "log10",  name = expression(paste('Mortality (yr'^-1,' g'^{1/4},')')),
                     limits = c(0.05, 20), breaks = c(0.1,  1, 10), labels = c("0.1", "1",  "10")) +  
  scale_x_reverse(name = expression("1/kT"), limits = c(43.2, 39.5), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, 
                                      name = "Temperature (ºC)")) +

  theme_single
grid.newpage()
grid.draw(mortal_plot_ecto)

# Endos
mortal_plot_endo <- ggplot(mortality %>% 
                             filter(Group == "Bird" |
                                      Group == "Mammal") %>%
                             arrange(Group),
                           aes(x = one_kT, y = mortality_corr, shape = Group)) + 
  geom_point(size = 2.5, stroke = .5, aes(shape = Group), color = "brown2") + 
  scale_shape_manual(values = c("Mammal" = 11, "Bird" = 18)) +
  geom_smooth(method = 'lm', size = 1,  alpha = 0.3, color = "brown3") +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{1/4},')')), 
                     trans = "log10", breaks = c(0.1, 1, 10), labels = c("0.1", "1", "10")) +  
  scale_x_reverse(name = expression("1/kT"), limits = c(),#43.2, 39.5), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, 
                                      name = "Temperature (ºC)")) +
  
  theme_single
grid.newpage()
grid.draw(mortal_plot_endo)

#### All
mortal_plot_all <- ggplot(mortality %>% 
                         filter(Group == "Fish" | Group == "Invertebrate" | Group == "Bird" |
                                  Group == "Mammal") %>%
                         arrange(Group),
                       aes(x = one_kT, y = mortality_corr, shape = Group)) + 
  geom_point(size = 2.5, stroke = .5, aes(shape = Group, color = Group)) + 
  scale_shape_manual(values = c("Invertebrate" = 2, "Fish" = 19, "Mammal" = 11, "Bird" = 18)) +
  geom_smooth(method = 'lm', size = 1,  alpha = 0.3, aes(color = Group)) +
  scale_color_manual(values =c("Invertebrate" = invert_col, "Fish" = fish_col, 
                               "Mammal" = mamm_col, "Bird" = bird_col)) +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{1/4},')')),
                     trans = "log10", breaks = c(0.1, 1, 10), labels = c("0.1", "1", "10")) + 
  scale_x_reverse(name = expression("1/kT"), limits = c(43.2, 39.5), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, 
                                      name = "Temperature (ºC)")) +
  theme_single
grid.newpage()
grid.draw(mortal_plot_all)
# A mess!


# Compare with McCoy results
    # regression of fish and inverts together
mortality5 <- mortality4 %>%
  filter(Group == "Invertebrate" | Group == "Fish")

lm1 <- lm(log(mortality_corr) ~ one_kT_mccoy , data = mortality5)
summary(lm1) # slope = -0.56971 (matches -0.57!) 

lm2a <- lm(log(mortality_corr) ~ one_kT_mccoy + Group, data = mortality5)
summary(lm2a) #-0.56190 (slightly different than -0.57!); 

#one_kT vs more complicated setup (one_kt_mccoy), mnakes no difference
lm2b <- lm(log(mortality_corr) ~ one_kT + Group, data = mortality5)
summary(lm2b)  #-0.56190 



#ectos
mean(mort_lm_ecto$estimate,na.rm=T)
#[1] -0.379574
median(mort_lm_ecto$estimate,na.rm=T)
#-0.5248579
t.test(mort_lm_ecto$estimate )$"conf.int"
#[1] -0.5941474 -0.1650006

#endos
mean(mort_lm_endo$estimate,na.rm=T)
#[1] 0.02355298
median(mort_lm_endo$estimate,na.rm=T)
#0.03015333
t.test(mort_lm_endo$estimate )$"conf.int"
#[1] -0.1003453  0.1474513

t.test(mort_lms2c$estimate ~ mort_lms2c$Thermy)
#t = -3.5084, df = 20.544, p-value = 0.002144

#b Birds
mean(mort_lm_endo$estimate[1:10])
#[1] 0.0623127
median(mort_lm_endo$estimate[1:10],na.rm=T)
#[1] 0.07359674
t.test(mort_lm_endo$estimate[1:10] )$"conf.int"
#[1] -0.0440225  0.1686479

# Mammals
mean(mort_lm_endo$estimate[11:20])
#[1] -0.01520674
median(mort_lm_endo$estimate[11:20],na.rm=T)
#[1] 0.02241846
t.test(mort_lm_endo$estimate[11:20] )$"conf.int"
#[1] -0.2656006  0.2351872


# Fish
mean(mort_lm_ecto$estimate[1:5])
#[1] -0.2095318
median(mort_lm_ecto$estimate[1:5],na.rm=T)
#[1] -0.151948
t.test(mort_lm_ecto$estimate[1:5] )$"conf.int"
#[1] -0.51757396  0.09851029


# Inverts
mean(mort_lm_ecto$estimate[6:13])
#[1] -0.4858504
median(mort_lm_ecto$estimate[6:13],na.rm=T)
#[1] -0.5480302
t.test(mort_lm_ecto$estimate[6:13] )$"conf.int"
#-0.8068326 -0.1648681

t.test(mort_lms2c$estimate ~ mort_lms2c$Thermy)
#t = -3.5084, df = 20.544, p-value = 0.002144


#kelley  1124    222.1  -234.6     cv with Kelley 95% CI
# Restricted to 5 indiv per genus, 5 C range


#one_kT vs more complicated setup (one_kt_mccoy), mnakes no difference

lm3 <- lm(log(mortality_final$mortality_corr) ~ mortality_final$one_kT + mortality_final$Group + mortality_final$Species )
summary(lm3)  #-0.56190 
aov(lm3)
library(rsq)
rsq.partial(lm3)
########################################################################
################################# Plots ##############################
########################################################################


#Phytoplankton don't look great, only two "species"
unique(mortality_final$Species[mortality_final$Group == "Phytoplankton"])
mortal_plot_phyto <- ggplot(mortality_final %>% 
                              filter(temp_range_spp >= 2) %>%
                         filter(Group == "Phytoplankton" | Group == "Multicellular plant") %>%
                         arrange(Species),
                       aes(x = one_kT, y = mortality_corr)) + 
  geom_smooth(method = 'lm', size = 1,  alpha = 0, color = "forestgreen",
              aes(group = Species)) +
  scale_y_log10(name = expression(paste('Mortality (yr'^-1,' g'^{1/4},')')))+#,
                     #limits = c(-4, 0.5)) +  
  scale_x_reverse(name = expression("1/kT"), limits = c(43, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, 
                                      name = "Temperature (ºC)")) +
  theme_single
grid.newpage()
grid.draw(mortal_plot_phyto)

# Ectotherms
# restricted to temp range of 5 c
mortal_plot_ectos <- ggplot(mortality_final %>% 
                              filter(Group == "Fish" | Group == "Invertebrate") %>%
                              #filter(n_species >= 2) %>%
                              filter(temp_range_spp >= 5) %>%
                              arrange(Species),
                            aes(x = one_kT, y = mortality_corr)) + 
  geom_smooth(method = 'lm', size = 0.5,  alpha = 0,
              aes(group = Genus, color = Group)) +
  scale_color_manual(values =c("Invertebrate" = invert_col, "Fish" = fish_col, 
                               "Mammal" = mamm_col, "Bird" = bird_col)) +
  scale_y_continuous(trans = "log10", name = expression(paste('Mortality (yr'^-1,' g'^{1/4},')')),
                     limits = c(.1, 10), breaks = c(0.1, 1, 10), labels = c(0.1, 1, 10)) +  
  scale_x_reverse(name = expression("1/kT"), limits = c(42.5, 39.5), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, 
                                      name = "Temperature (ºC)")) +
  theme_single
grid.newpage()
grid.draw(mortal_plot_ectos)

### All, by species for ectos (no phytoplankton), by genus for endos - looks  a little crazy
mortal_plot6 <- ggplot() + 
  geom_smooth(data = mortality_final %>%
                filter(Group == "Fish" | Group == "Invertebrate"),
                      # temp_range_spp >=5, n_species >= 3), 
              method = 'lm', size = 0.5,  alpha = 0,
              aes(group = Species, color = Group,
                  x = one_kT, y = mortality_corr)) +
  geom_smooth(data = mortality_final %>%
                filter(Group == "Mammal" | Group == "Bird"),
                      # temp_range_genus >=5, n_genus >= 3),
              method = 'lm', size = 0.5,  alpha = 0,
              aes(group = Genus, color = Group,
                  x = one_kT, y = mortality_corr)) +
  #geom_smooth(method = 'lm', size = 1.5,  alpha = 0,
  #           aes(group = Group, color = Group)) +
  scale_color_manual(values =c("Invertebrate" = invert_col, "Fish" = fish_col, 
                               "Mammal" = mamm_col, "Bird" = bird_col)) +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{1/4},')')),
                     trans = "log10", breaks = c(0.3, 1, 3,10), labels = c("0.3", "1", "3", "10"), limits = c(0.18, 10)) +  
  scale_x_reverse(name =  expression("1/kT"), limits = c(44.8, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, 
                                      name = "Temperature (ºC)")) +
  annotate("text", label = "Bird", hjust = 0, color = "brown4", x = 44.8, y = 6.9, size = 6, fontface = "bold") +
  annotate("text", label = "Mammal", hjust = 0,color = "red", x = 44.8, y = 9, size = 6, fontface = "bold") +
  annotate("text", label = "Invertebrate", hjust = 0,color = "deepskyblue2", x = 44.8, y = 5.3, size = 6, fontface = "bold") +
  annotate("text", label = "Fish", hjust = 0,color = "navy", x = 44.8, y = 4, size = 6, fontface = "bold") +
  theme_single +theme(legend.position = "none")
grid.newpage()
grid.draw(mortal_plot6)

ggsave(mortal_plot6, height = 5, width = 7, filename = 'figfig.pdf')

# Best compromise??
# 5 indiv per genus, temp range of 5
mortal_plot_all <- ggplot(mortality_final %>% 
                         filter(Group == "Fish" | Group == "Invertebrate" |
                                  Group == "Mammal" | Group == "Bird") %>%
                         #filter(n_genus >= 5) %>%
                         filter(temp_range_genus >= 5),
                       aes(x = one_kT, y = mortality_corr)) + 
  geom_smooth(method = 'lm', size = 0.5,  alpha = 0,
              aes(group = Genus, color = Group)) +
  scale_color_manual(values =c("Invertebrate" = invert_col, "Fish" = fish_col, 
                               "Mammal" = mamm_col, "Bird" = bird_col)) +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{1/4},')')),
                     trans = "log10", breaks = c(0.3, 1, 3,10), labels = c("0.3", "1", "3", "10"), limits = c(0.3, 10)) +  
  scale_x_reverse(name = expression("1/kT"), limits = c(44.8, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, 
                                      name = "Temperature (ºC)")) +
  #annotate("text", label = "Bird", hjust = 0, color = bird_col, x = 44.8, y = 6.9, size = 6, fontface = "bold") +
  #annotate("text", label = "Mammal", hjust = 0,color = mamm_col, x = 44.8, y = 9, size = 6, fontface = "bold") +
  #annotate("text", label = "Invertebrate", hjust = 0, color = invert_col, x = 44.8, y = 5.3, size = 6, fontface = "bold") +
  #annotate("text", label = "Fish", hjust = 0,color = fish_col, x = 44.8, y = 4, size = 6, fontface = "bold") +
  theme_single +theme(legend.position = "none")
grid.newpage()
grid.draw(mortal_plot_all)
ggsave(mortal_plot_all, height = 5, width = 7, filename = file.path(gdrive_path, 'Plots/mortality_by_genus_5_indiv_5C.pdf'))




###############Supplemental Mortality: regression for all groups, showing the data points

#combine supplements
g_ecto <- ggplotGrob(mortal_plot_ecto)
g_endo <- ggplotGrob(mortal_plot_endo)
g3 <- rbind(g_endo, g_ecto , size = "first")
g3$widths <- unit.pmax(g_ecto$widths, g_ecto$widths)
grid.newpage()
grid.draw(g3)
ggsave(g3, height = 7.6, width = 6, filename = file.path(gdrive_path,'Plots/comb_mortality.pdf'))


##################################################################  
######################  Attack Rates #############################
##################################################################  

unique(attack$predator.species)
# slope and r2 depends on what variables are included in addition to temp
lm_attack <- lm( data= attack, log(attack.rate) ~ one_kT)
summary(lm_attack) # E = -0.75, r2 = 0.03
lm_attack_dim <- lm(data= attack, log(attack.rate) ~ one_kT +dimensionality)
summary(lm_attack_dim) # E = -0.75, r2 = 0.05

lm_attack_dim_sp <- lm(data= attack, log(attack.rate) ~ one_kT +dimensionality + predator.species)
summary(lm_attack_dim_sp) # E = -0.429, int = 4.005, r2 = 0.93

lm_attack_sp <- lm(data= attack, log(attack.rate) ~ one_kT + predator.species)
summary(lm_attack_sp) # E = -0.425, int = 43.856
attack_lab  <- expression(paste('y = -0.425x, r'^2,' = 0.94'))

## All data points
attack_all2 <- ggplot(attack,
                       aes(x = one_kT, y = attack.rate)) + #, shape = Group)) + 
  geom_point(size = 2.5, stroke = .5,  shape =21, fill = "deepskyblue3", color = "black") +
  geom_smooth(method = 'lm', size = 1,  alpha = 0.3, color = "deepskyblue3") +
  scale_y_log10(breaks = trans_breaks(n = 4,"log10", function(x) 10^x),
                                   labels = trans_format("log10", math_format(10^.x)),
                                   name = expression(paste('Attack Rate (m'^2,' s'^-1,')'))) +
             
  scale_x_reverse(name = expression("1/kT"), limits = c(42.2, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, 
                                      name = "Temperature (ºC)")) +
  annotate("text", label = "y = -0.76x, r") +
  theme_single
grid.newpage()
grid.draw(attack_all2)


#combined with mortality
attack_regr_genus <- ggplot(attack %>% 
         #filter(dimensionality == "2D") %>% 
         filter(n_genus >= 5) %>%
         filter(temp_range_genus >= 5) %>%
         arrange(Pred_Genus),
       aes(x = one_kT, y = attack.rate)) + 
  geom_smooth(method = 'lm', size = 1,  alpha = 0,
              aes(group = Pred_Genus, color = predator.met.group), linetype = "dashed") +
  scale_color_manual(values =c("invert" = invert_col, "ectovert" = "navy"))+
  scale_y_continuous(name = expression(paste('Attack Rate (m'^2,' s'^-1,')')),
                     limits = c(10^-10, 10^-3), 
                     labels = trans_format("log10", math_format(10^.x)),
                     breaks = c(10^-11, 10^-9, 10^-7, 10^-5, 10^-3), trans = "log10") +  
  scale_x_reverse(name = expression("1/kT"), limits = c(42, 38), 
                  sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, 
                                      name = NULL)) +
  #annotate("text", y = 10^-9.5, x = 40, hjust = 0, size = 6, label = attack_lab, color = invert_col) +
  theme_single + theme(legend.position = "none")

grid.newpage()
grid.draw(attack_regr_genus)
ggsave(attack_regr_genus , height = 4.7, width = 6.7, filename = file.path(gdrive_path,'Plots/attack_by_genus_5_indiv_5C_over.pdf'))


# All slopes, calculated separately


box_mort_attack_ea <- ggplot(mort_attack , aes(x = Type, y = -1*estimate, fill = Type))+
  geom_boxplot(fill = "dodgerblue") + 
  scale_y_continuous(position = "right", name=expression(paste('Thermal Sensitivity (E'[a],')')))+ 
  scale_x_discrete(name = NULL) +
  # scale_fill_manual(values =c("Ectotherm" = "skyblue3", "Endotherm" = "red2")) +
  theme_single + theme(legend.position = "none") 
box_mort_attack_ea
ggsave(box_mort_attack_ea,  height = 4.22, width = 5.9,  filename = file.path(gdrive_path,'Plots/box_mort_attack_ea.pdf'))


# combine plots
#combine supplements
g_ecto <- ggplotGrob(mortal_plot_ecto)
g_endo <- ggplotGrob(mortal_plot_endo)
g3 <- rbind(g_endo, g_ecto , size = "first")
g3$widths <- unit.pmax(g_ecto$widths, g_ecto$widths)
grid.newpage()
grid.draw(g3)
ggsave(g3, height = 7.6, width = 6, filename = file.path(gdrive_path,'Plots/comb_mortality.pdf'))

mortal_plot_endo
mortal_plot_ecto
mortal_plot_5C_5ind
box_mort_Ea 
attack_regr
box_mort_attack_ea 


g_mortal_plot_endo<- ggplotGrob(mortal_plot_endo)
g_mortal_plot_ecto <- ggplotGrob(mortal_plot_ecto)
g_mortal_plot_5C_5ind <- ggplotGrob(mortal_plot_5C_5ind)
g_box_mort_Ea<- ggplotGrob(box_mort_Ea)
g_attack_regr <- ggplotGrob(attack_regr)
g_box_mort_attack_ea <- ggplotGrob(box_mort_attack_ea)



g3a <- rbind(g_mortal_plot_endo, g_mortal_plot_5C_5ind , g_attack_regr, size = "first")
g3a$widths <- unit.pmax(g_mortal_plot_endo$widths, g_mortal_plot_5C_5ind$widths,g_attack_regr$widths)
grid.newpage()
grid.draw(g3a)
g3b <- rbind(g_mortal_plot_ecto, g_box_mort_Ea, g_box_mort_attack_ea , size = "first")
g3b$widths <- unit.pmax(g_mortal_plot_ecto$widths, g_box_mort_Ea$widths,g_box_mort_attack_ea$widths)
grid.newpage()
grid.draw(g3b)

g4 <- cbind(g3a, g3b)
grid.newpage()
grid.draw(g4)

ggsave(g4, height = 7.6, width = 6, filename = file.path(gdrive_path,'Plots/all_figs.pdf'))
