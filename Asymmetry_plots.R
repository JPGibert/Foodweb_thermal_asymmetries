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
gdrive_path <- file.path('YOUR PATH HERE')

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

attack <- read_csv(file.path(gdrive_path, 'Lietal_oikos_2017_data.csv'))
attack$one_kT <- 1/(8.617e-5*(attack$temperature.degree.celcius +273.15))

mortality1 <- read_csv(file.path(gdrive_path,'McCoy_Gillooly_Data.csv')) #2,117 values
mortality1[is.na(mortality1$Species),]

######################################################
############## Mortality Data first
######################################################


#
# All slopes, calculated separately
##########################################################################

###################################  Add endo temp data ################################## 
  # This temp determined from matching lat/lon centroid of species range to climate data from ...
bird_temp <- read_csv(file.path(gdrive_path,"bird_temp.csv"))
mammal_temp <- read_csv(file.path(gdrive_path,"mammal_temp.csv"))
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
ggsave(mortal_plot_endo, height = 4.6, width = 6.7, filename = file.path(gdrive_path, "NAME YOUR FILE.pdf"))


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
ggsave(mortal_plot_ecto, height = 5, width = 7, filename = file.path(gdrive_path, 'NAME YOUR FILE.pdf'))


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
ggsave(mortal_plot_5C_5ind, height = 4.3, width = 6.3, filename = file.path(gdrive_path,'NAME YOUR FILE.pdf'))
ggsave(mortal_plot_5C_5ind, height = 4.3, width = 6.3, filename = 'NAME YOUR FILE.pdf')


##########################################
#### Panel 2d, Mortality violin plot
 #first get regressions
##########################################

#Restricted to 5 indiv per genus, 5 C range
# IN prep for the plot------------------------------------
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

##-----------------
#ACTUALLY MAKES THE VIOLIN PLOT:
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
ggsave(mort_vio_Ea,  height = 4.22, width = 6.05, filename = file.path(gdrive_path,'NAME YOUR FILE.pdf'))


####################################################
##v  Fig 2E - Mortality & Attack regressions (ectos)
####################################################

##---------------------------------
## PART 1 (mortality)
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
ggsave(mortal_plot_5C_5ind, height = 4.7, width = 6.7, filename = file.path(gdrive_path,'NAME YOUR FILE.pdf'))

##---------------------------------
## PART 2 (attack)
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
ggsave(attack_regr_genus , height = 4.7, width = 6.7, filename = file.path(gdrive_path,'NAME YOUR FILE.pdf'))
## Attack and mortality plots were manually superimposed

############################################
##  Fig 2F - violin plot
############################################
## 1) Gets the data for Violin plot
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

## 2) Plots the violin plot
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
ggsave(vio_mort_attack_ea,  height = 4.22, width = 5.9,  filename = file.path(gdrive_path,'NAME YOUR FIGURE.pdf'))


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
