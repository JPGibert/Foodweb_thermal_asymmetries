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
#library(egg)
library(broom.mixed)
#library(tidybayes)

# ----  Adjust all file paths to read in data and save plots
github_path <- file.path('/Users/jgradym/Documents/GitHub/Foodweb_thermal_asymmetries/Data')


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
mortality0 <- read_csv(file.path(github_path, 'McCoy_mortality_updated.csv')) 
mortality0$one_kT <- 1/(8.617e-5*(mortality0$Temp_C + 273.15)) # Inverse Temperature

mortality1 <- mortality0 %>%
  rename(species = Species) %>%
  arrange(Group, species) %>%
  mutate(log_mort = log(Mortality_yr), log_mass = log(dry_mass_g),
         Genus = word(species, sep = "_"), species2 = species)
mortality1

# Get number measurements per species or per genus, if you want to restrict by individuals per species or genus:, 
mortality3 <- mortality1 %>% #firs
  mutate(Genus = word(species, 1, sep = "_")) %>%
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
  mutate(Thermy = if_else(Group == "Bird" | Group == "Mammal", "Endotherm", "Ectotherm")) %>%
  arrange(Thermy, Group, Genus)

mortality <- mortality5

# Plotting data, by Group
mortality_group <- mortality %>% 
  group_by(Group) %>%
  mutate(min_one_kT = min(one_kT), max_one_kT = max(one_kT)) %>%
  slice_head() %>%
  arrange(Thermy, Group, Genus)

mortality_group_ecto <- mortality_group %>% filter(Thermy == "Ectotherm")

mortality_group_endo <- mortality_group %>% filter(Thermy == "Endotherm")

# Plotting data, by Genus
mortality_genus <- mortality %>% 
  group_by(Genus) %>%
  mutate(min_one_kT = min(one_kT), max_one_kT = max(one_kT)) %>%
  slice_head() %>%
  arrange(Thermy, Group, Genus)

##########################################################################
#--------------- Analyze Groups ----------

#########################################################################3
#Ectotherms

mort_ecto <- mortality %>% filter(Thermy == "Ectotherm")
lmer_mort_ecto3 <- lmer(log_mort ~ one_kT + log_mass +
                         (one_kT|Genus),  data = mort_ecto)
lmer_mort_ecto <- lmer(log_mort ~ one_kT + log_mass + (one_kT|Group) +
                       (one_kT|Group:Genus),  data = mort_ecto)
lmer_mort_ecto <- lmer(log_mort ~ one_kT + log_mass + (one_kT|Group) +
                         (one_kT|Group:Genus),  data = mort_ecto)
lmer_mort_ecto1 <- lmer(log_mort ~ one_kT + log_mass +
                         (one_kT|Group:Genus), data = mort_ecto)

lmer_mort_ecto2 <- lmer(log_mort ~ log_mass + one_kT + (one_kT|Group/Genus),
                       data = mort_ecto) #genus slopes fine but barely vary

lmer_mort_ecto2 <- lmer(log_mort ~ log_mass + one_kT + (one_kT|Group/Genus),
                        data = mort_ecto) #genus slopes fine but barely vary


lmer_mort_ecto2 <- lmer(log_mort ~ log_mass + (one_kT|Group/Genus),
                        data = mort_ecto) #genus slopes close to zero

lmer_mort_ecto2 <- lmer(log_mort ~ log_mass + (one_kT|Group) + (one_kT|Group:Genus),
                        data = mort_ecto) 
(x|site/block) = (x|site)+(x|site:block) = (1 + x|site)+(1+x|site:block)
lmer_mort_ecto2 <- lmer(log_mort ~ log_mass + one_kT + (1|Group) + (0+ one_kT|Group/Genus),
                        data = mort_ecto)

lmer_mort_ecto2 <- lmer(log_mort ~ log_mass + one_kT + (one_kT|Group) + (one_kT|Group/Genus),
                        data = mort_ecto)

lmer_mort_ecto2 <- lmer(log_mort ~ log_mass + one_kT + (one_kT|Group) + (one_kT|Genus),
                        data = mort_ecto)

lmer_mort_ecto2 <- lmer(log_mort ~ log_mass + one_kT + (one_kT|Group) + (one_kT|Genus),
                        data = mort_ecto)

(x|site)+(x|site:block) = (1 + x|site)+(1+x|site:block)

lmer_mort_ecto2 <- lmer(log_mort ~ log_mass + one_kT + (1|Group) + (0+ one_kT|Group/Genus),
                        data = mort_ecto)

lmer_mort_all<- lmer(log_mort ~ log_mass + one_kT + (one_kT|Group/Genus),
                        data = mortality)

lmer_mort_all<- lmer(log_mort ~ log_mass + one_kT + (one_kT|Group) + (one_kT|Genus),
                     data = mortality)

lmer_mort_ecto2a <- lmer(log_mort ~ log_mass + one_kT + (1|Group) + (0+ one_kT|Group/Genus),
                        data = mort_ecto)
x*site+(x|site:block)

summary(lmer_mort_all)

mort_genus <-  coef(lmer_mort_all)$Genus
coef(lmer_mort_all)$Genus
coef(lmer_mort_all)$Group

coef(lmer_mort_all)$`Genus:Group`
coef <- coef(lmer_mort_all)
str(coef)
coef$`Genus:Group`
coef(lmer_mort_ecto2)$`Genus:Group`


lmer_mort_group_ecto <- coef(lmer_mort_ecto)
lmer_mort_group_ecto$`Group:Genus`


#Endotherms
mort_endo <- mortality1 %>% filter(Group == "Mammal" | Group == "Bird")


lmer_mort_endo <- lmer(log_mort ~ log_mass + (one_kT|Group/Genus),
                       data = mort_endo)

(x|site/block) 
# Results
summary(lmer_mort_endo)
coef(lmer_mort_endo)$Group
r.squaredLR(lmer_mort_endo)
r.squaredGLMM(lmer_mort_endo)


confint.merMod(lmer_mort_endo, method = "Wald")

lmer_mort_group_endo <- coef(lmer_mort_endo)$Group

coef_endo <- coef(lmer_mort_endo2)
coef_endo$`Genus:Group`
# Group regressions
     #Endotherms 

lmer_mort_group_endo <- as_tibble(rownames_to_column(lmer_mort_group_endo)) %>%
  rename(intercept = `(Intercept)`, temp_slope = one_kT, mass_slope = log_mass, Group = rowname) 
lmer_mort_group_endo

lmer_mort_group_endo$min_one_kT <-  mortality_group_endo$min_one_kT
lmer_mort_group_endo$min_mort <-  (lmer_mort_group_endo$intercept + lmer_mort_group_endo$min_one_kT *lmer_mort_group_endo$temp_slope)

lmer_mort_group_endo$max_one_kT <-  mortality_group_endo$max_one_kT
lmer_mort_group_endo$max_mort <-  lmer_mort_group_endo$intercept + lmer_mort_group_endo$max_one_kT *lmer_mort_group_endo$temp_slope
lmer_mort_group_endo

lmer_mort_group_endo2 <- lmer_mort_group_endo %>%
  pivot_longer(cols = c(min_one_kT, max_one_kT),
               names_to = "names") %>%
  rename(one_kT = value) %>%
  mutate(mort = c(min_mort[1], max_mort[1], min_mort[3], max_mort[3])) %>%
  select(-c(names, min_mort, max_mort))
lmer_mort_group_endo2



# add plotting information

# add plotting information
lmer_mort_gen_endo$min <-  mortality_genus_endo$min_one_kT
lmer_mort_gen_endo$max <-  mortality_genus_endo$max_one_kT
lmer_mort_gen_endo$n_genus <-  mortality_genus_endo$n_genus
lmer_mort_gen_endo$temp_range <-  mortality_genus_endo$temp_range_genus


# Parameters for plotting
fish_intercept <-  coef(lmer_mort_ecto)$Group[1, 1]
fish_temp_slope <- coef(lmer_mort_ecto)$Group[1, 2]
ecto_mass_slope <- coef(lmer_mort_ecto)$Group[1, 3]

invert_intercept <- coef(lmer_mort_ecto)$Group[2, 1]
invert_temp_slope <- coef(lmer_mort_ecto)$Group[2, 2]




# mortality slopes, intercept for plotting

bird_intercept <- coef(lmer_mort_endo)$Group[1, 1]
bird_temp_slope <- coef(lmer_mort_endo)$Group[1, 2]
endo_mass_slope <- coef(lmer_mort_endo)$Group[1, 3]

mammal_intercept <- coef(lmer_mort_endo)$Group[2, 1]
mammal_temp_slope <- coef(lmer_mort_endo)$Group[2, 2]

# add mass correction for plotting
mort_ecto$mort_mass_corr <- mort_ecto$Mortality_yr/mort_ecto$dry_mass_g^ecto_mass_slope
mortality$mort_mass_corr[mortality$Thermy == "Ectotherm"] <- mort_ecto$mort_mass_corr 

mort_endo$mort_mass_corr <- mort_endo$Mortality_yr/mort_endo$dry_mass_g^endo_mass_slope 
mortality$mort_mass_corr[mortality$Thermy == "Endotherm"] <- mort_endo$mort_mass_corr 


# plotting
#------------- Plot Fig 2A  -----------------

mortal_plot_endo <- ggplot(mortality %>% filter(Thermy == "Endotherm") %>%
                             arrange(Group),
                           aes(x = one_kT, y = mort_mass_corr, shape = Group)) + 
  theme_single + theme(legend.position = "none") +  theme_ticks +
  scale_y_continuous(name = expression(paste('Mortality (yr'^-1,' g'^{-alpha},')')), 
                     trans = "log", breaks = c(0.1, 1, 10), labels = c("0.1", "1", "10"), 
                     #limits = c(0.02, 40)
  ) +  
  scale_x_reverse( name = NULL, breaks = c(46,44, 42,40, 38),
                   #limits = c(45.1, 38),
                   sec.axis = sec_axis(~ 1/(.*0.00008617) - 273.15, breaks = c(-20, -10, 0, 10, 20, 30),
                                       name = "Ambient Temperature (ºC)")) +
  geom_point(size = 2.5, shape = 21, color = "black",  stroke = .5, aes(fill = Group)) + 
  scale_fill_manual(values =c("Mammal" = mamm_col, "Bird" = "darkred")) +
  scale_color_manual(values =c("Mammal" = mamm_col, "Bird" = "brown4")) +
  geom_line(data = lmer_mort_group_endo2, size = 1,
            aes(x = one_kT, y = exp(mort),  color = Group))
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
                         y = c(exp((invert_intercept + (invert_temp_slope * min(mort_invert$one_kT)))), 
                               exp((invert_intercept + (invert_temp_slope  *  max(mort_invert$one_kT))))),
                         Group = c("Invertebrate", "Invertebrate"))

mort_fish <-  mortality[mortality$Group == "Fish",]
fish_fit <- data.frame(x = c(min(mort_fish$one_kT), max(mort_fish$one_kT)), 
                       y = c( exp((fish_intercept + (fish_temp_slope * min(mort_fish$one_kT)))), 
                              exp((fish_intercept + (fish_temp_slope *  max(mort_fish$one_kT))))),
                       Group = c("Fish", "Fish"))

ecto_fit <- data.frame(rbind(invert_fit, fish_fit))
ecto_fit 


# Group regressions
#Ectotherms
# Group regressions
lmer_mort_group_ecto <- coef(lmer_mort_ecto)$Group


lmer_mort_group_ecto <- as_tibble(rownames_to_column(lmer_mort_group_ecto)) %>%
  rename(intercept = `(Intercept)`, temp_slope = one_kT, mass_slope = log_mass, Group = rowname) 
lmer_mort_group_ecto

lmer_mort_group_ecto$min_one_kT <-  mortality_group_ecto$min_one_kT
lmer_mort_group_ecto$min_mort <-  (lmer_mort_group_ecto$intercept + lmer_mort_group_ecto$min_one_kT *lmer_mort_group_ecto$temp_slope)

lmer_mort_group_ecto$max_one_kT <-  mortality_group_ecto$max_one_kT
lmer_mort_group_ecto$max_mort <-  lmer_mort_group_ecto$intercept + lmer_mort_group_ecto$max_one_kT *lmer_mort_group_ecto$temp_slope
lmer_mort_group_ecto

lmer_mort_group_ecto2 <- lmer_mort_group_ecto %>%
  pivot_longer(cols = c(min_one_kT, max_one_kT),
               names_to = "names") %>%
  rename(one_kT = value) %>%
  mutate(mort = c(min_mort[1], max_mort[1], min_mort[3], max_mort[3])) %>%
  select(-c(names, min_mort, max_mort))
lmer_mort_group_ecto2



#--------------  Plot Fig 2b ------------
mortal_plot_ecto <- ggplot(mortality %>% 
                             filter(Group == "Fish" | Group == "Invertebrate") %>%
                             arrange(Group),
                           aes(x = one_kT, y = mort_mass_corr,fill = Group)) + 
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
  geom_line(data = lmer_mort_group_ecto2, aes(x = one_kT, y = exp(mort), color = Group), size = 1) + 
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

mortality_genus_ecto <- mortality_genus %>% filter(Thermy == "Ectotherm")

# Genus regressions
lmer_mort_gen_ecto <- coef(lmer_mort_ecto)$`Group:Genus`

lmer_mort_gen_ecto <- as_tibble(rownames_to_column(lmer_mort_gen_ecto)) %>%
  rename(intercept = `(Intercept)`, temp_slope = one_kT, mass_slope = log_mass) %>%
  mutate(Group = word(rowname, 1, sep = ":"), Genus = word(rowname, 2, sep = ":")) %>%
  select(-rowname) %>%
  select(Group, Genus, everything())
library(data.table)


lmer_mort_gen_ecto$min_one_kT <-  mortality_genus_ecto$min_one_kT
lmer_mort_gen_ecto$min_mort <-  (lmer_mort_gen_ecto$intercept + lmer_mort_gen_ecto$min_one_kT *lmer_mort_gen_ecto$temp_slope)

lmer_mort_gen_ecto$max_one_kT <-  mortality_genus_ecto$max_one_kT
lmer_mort_gen_ecto$max_mort <-  lmer_mort_gen_ecto$intercept + lmer_mort_gen_ecto$max_one_kT *lmer_mort_gen_ecto$temp_slope

lmer_mort_gen_ecto$n_genus <-  mortality_genus_ecto$n_genus
lmer_mort_gen_ecto$temp_range <-  mortality_genus_ecto$temp_range_genus


lmer_mort_gen_ecto2 <- lmer_mort_gen_ecto %>%
  pivot_longer(cols = c(min_one_kT, max_one_kT),
               names_to = "names") %>%
  rename(one_kT = value) %>%
  select(-names)
lmer_mort_gen_ecto2
lmer_mort_gen_ecto2$index <- 1:384

mort1 <- lmer_mort_gen_ecto2 %>% 
  filter(row_number() %% 2 == 1) %>%
  select(min_mort, index) %>%
  rename(mort = min_mort)
mort2 <- lmer_mort_gen_ecto2 %>% 
  filter(row_number() %% 2 == 0) %>%
  select(max_mort, index)  %>%
  rename(mort = max_mort)
mort3 <- as_tibble(rbind(mort1, mort2 )) %>%
  arrange(index)
mort4 <- mort3 %>%pull(mort) 

lmer_mort_gen_ecto2 <-  lmer_mort_gen_ecto2 %>%
  mutate(mort = mort4) %>%
  select(-min_mort, -max_mort)




# Genus
mortality_genus_endo <- mortality_genus %>% filter(Thermy == "Endotherm")

genus_mort_endo <- coef(lmer_mort_endo)$`Group:Genus`
genus_mort_endo 

lmer_mort_gen_endo <- coef(lmer_mort_endo)$`Group:Genus`

lmer_mort_gen_endo <- as_tibble(rownames_to_column(lmer_mort_gen_endo)) %>%
  rename(intercept = `(Intercept)`, temp_slope = one_kT, mass_slope = log_mass) %>%
  mutate(Group = word(rowname, 1, sep = ":"), Genus = word(rowname, 2, sep = ":")) %>%
  select(-rowname) %>%
  select(Group, Genus, everything())

# add plotting information
lmer_mort_gen_endo$min_one_kT <-  mortality_genus_endo$min_one_kT
lmer_mort_gen_endo$min_mort <-  lmer_mort_gen_endo$intercept + lmer_mort_gen_endo$min_one_kT *lmer_mort_gen_endo$temp_slope

lmer_mort_gen_endo$max_one_kT <-  mortality_genus_endo$max_one_kT
lmer_mort_gen_endo$max_mort <-  lmer_mort_gen_endo$intercept + lmer_mort_gen_endo$max_one_kT *lmer_mort_gen_endo$temp_slope

lmer_mort_gen_endo$n_genus <-  mortality_genus_endo$n_genus
lmer_mort_gen_endo$temp_range <-  mortality_genus_endo$temp_range_genus


lmer_mort_gen_endo2 <- lmer_mort_gen_endo %>%
  pivot_longer(cols = c(min_one_kT, max_one_kT),
               names_to = "names") %>%
  rename(one_kT = value) %>%
  select(-names)
lmer_mort_gen_endo2
lmer_mort_gen_endo2$index <- 1:992

mort1 <- lmer_mort_gen_endo2 %>% 
  filter(row_number() %% 2 == 1) %>%
  select(min_mort, index) %>%
  rename(mort = min_mort)
mort2 <- lmer_mort_gen_endo2 %>% 
  filter(row_number() %% 2 == 0) %>%
  select(max_mort, index)  %>%
  rename(mort = max_mort)
mort3 <- as_tibble(rbind(mort1, mort2 )) %>%
  arrange(index)
mort4 <- mort3 %>%pull(mort) 

lmer_mort_gen_endo2 <-  lmer_mort_gen_endo2 %>%
  mutate(mort = mort4) %>%
  select(-min_mort, -max_mort)

lmer_mort_gen_endo2 
lmer_mort_gen_ecto2 

lmer_mort_genus <- as_tibble(rbind(lmer_mort_gen_endo2, lmer_mort_gen_ecto2 ))
#--- plot Fig 2C ------

mortal_genus_regr <- ggplot(data = lmer_mort_genus %>% 
                              filter(n_genus >= 5, temp_range>= 5)) +
  geom_smooth(method = "lm", se =- F, size = 0.3, aes(x = one_kT, y = exp(mort),
                                                      group = Genus, color = Group)) + 
  #geom_line(data = lmer_mort_group, aes(x = one_kT, y = exp(mort), color = Group), size = 1.5) +
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

mortality$log_mort_mass <- log(mortality$mort_mass_corr)
mort_ecto$log_mort_mass <- log(mort_ecto$mort_mass_corr)

genus_lmer <- mort_ecto %>%
  filter(n_genus >= 5, temp_range_genus >= 5) %>%
  nest(-Genus) %>% 
  mutate(
    fit = purrr::map(data, ~ lm(
      log_mort ~ log_mass + one_kT, data = .x )),#, 
      #control = lmerControl(check.nlev.gtr.1 = "ignore",
       #                     check.conv.singular = "ignore", 
        #                    check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE= "ignore"))
     
      tidied = map(fit, tidy, conf.int = T)
  ) 
genus_lmer 
genus_lmer <- unnest(genus_lmer, tidied)

genus_lmer <- mort_ecto %>%
  filter(n_genus >= 5, temp_range_genus >= 5) %>%
  nest(-Group) %>% 
  mutate(
    fit = purrr::map(data, ~ lmer(
      log_mort ~ log_mass + one_kT + (one_kT|Genus), data = .x ),#, 
    control = lmerControl(check.nlev.gtr.1 = "ignore",
                         check.conv.singular = "ignore", 
                        check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE= "ignore")),
    tidied = map(fit, tidy, conf.int = T)) 

    

genus_lmer 
genus_lmer <- unnest(genus_lmer, tidied)  %>%  filter(term == "one_kT")

#------- Add mass-corrected mortality for plotting ----

mortality2 <- mortality1 %>% 
  mutate(mass_corr_mortality = NA)

bird <- mortality2[mortality2$Group == "Bird",]
mortality2$mass_corr_mortality[mortality2$Group == "Bird"] <- bird$Mortality_yr/bird$dry_mass_g^bird_mass_slope 

mammal <- mortality2[mortality2$Group == "Mammal",]
mortality2$mass_corr_mortality[mortality2$Group == "Mammal"] <- mammal$Mortality_yr/mammal$dry_mass_g^mammal_mass_slope

invert <- mortality2[mortality2$Group == "Invertebrate",]
mortality2$mass_corr_mortality[mortality2$Group == "Invertebrate"] <- invert$Mortality_yr/invert$dry_mass_g^invert_mass_slope 

fish <- mortality2[mortality2$Group == "Fish",]
mortality2$mass_corr_mortality[mortality2$Group == "Fish"] <- fish$Mortality_yr/fish$dry_mass_g^fish_mass_slope

# 



#-------Add columns for genus level analysis ------------

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
mort_bird <- mortality[mortality$Group == "Bird",]
mort_mammal <- mortality[mortality$Group == "Mammal",]

 #----------------Mammals ------------------------------------------

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




# Plot - remove y limits to see outliers

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

ggsave(attack_regr_genus , height = 4.7, width = 6.7, filename = file.path(github_path,'NAME YOUR FILE.pdf'))


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

ggsave(vio_mort_attack_ea,  height = 4.22, width = 5.9,  filename = file.path(github_path,'NAME YOUR FIGURE.pdf'))


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
