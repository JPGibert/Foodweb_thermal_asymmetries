library(tidyverse)
#mammal mortality data

#mammal diet data
mammal_diet0 <- read_tsv('/Users/jgradym/Google Drive/Gibert Paper/Data/Diet/Mammals/doi_10.5061_dryad.6cd0v__v1/MammalDIET_v1.0.txt')
mammal_diet1 <- mammal_diet0 %>%
  mutate(Species = paste(Genus, Species, sep = "_")) %>%
  select(Species, TrophicLevel)
unique(mammal_diet1$TrophicLevel)
mammal_diet1$trophic_level <- NA
mammal_diet1$trophic_level[mammal_diet1$TrophicLevel == "Herbivore"] <- 2
mammal_diet1$trophic_level[mammal_diet1$TrophicLevel == "Carnivore"] <- 3
mammal_diet1$trophic_level[mammal_diet1$TrophicLevel == "Omnivore"] <- 2.5
mammal_diet1$trophic_level[mammal_diet1$TrophicLevel == "NotAssigned"] <- NA
mammal_diet1

#combine with mortality
mortality0 <- read_csv(file.path(gdrive_path,'McCoy_mortality_updated.csv')) 
mammal_mort0 <- mortality0 %>% filter(Group == "Mammal")
mammal_mort <- left_join(mammal_mort0, mammal_diet1, by = "Species")
mammal_mort$Genus <- word(mammal_mort$Species, 1, sep = "_")
missing_mamm <- mammal_mort[is.na(mammal_mort$trophic_level),]
missing_mamm_gen <- word(missing_mamm$Species, 1, sep = "_")
unique(missing_mamm_gen )

#manual fix
mammal_mort$trophic_level[mammal_mort$Species == "Vulpes_lagopus"] <- 3
mammal_mort$trophic_level[mammal_mort$Genus == "Arctocephalus"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Balaena"] <- 3
mammal_mort$trophic_level[mammal_mort$Genus == "Balaenoptera"] <- 3
mammal_mort$trophic_level[mammal_mort$Genus == "Berardius"] <- 4.5
mammal_mort$trophic_level[mammal_mort$Genus == "Sapajus"] <- 2.5
mammal_mort$trophic_level[mammal_mort$Genus == "Callorhinus"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Cephalorhynchus"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Cystophora"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Delphinapterus"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Delphinus"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Dugong"] <- 2
mammal_mort$trophic_level[mammal_mort$Genus == "Erignathus"] <- 3.5
mammal_mort$trophic_level[mammal_mort$Genus == "Eschrichtius"] <- 3
mammal_mort$trophic_level[mammal_mort$Genus == "Eumetopias"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Halichoerus"] <- 3
mammal_mort$trophic_level[mammal_mort$Genus == "Hydrurga"] <- 4.5
mammal_mort$trophic_level[mammal_mort$Genus == "Hyperoodon"] <- 4.5
mammal_mort$trophic_level[mammal_mort$Genus == "Lagenorhynchus"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Lobodon"] <- 3
mammal_mort$trophic_level[mammal_mort$Genus == "Megaptera"] <- 3
mammal_mort$trophic_level[mammal_mort$Genus == "Mirounga"] <- 4.5
mammal_mort$trophic_level[mammal_mort$Genus == "Monodon"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Neophocaena"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Odobenus"] <- 3
mammal_mort$trophic_level[mammal_mort$Genus == "Orcinus"] <- 4.5
mammal_mort$trophic_level[mammal_mort$Genus == "Otaria"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Pusa"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Peponocephala"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Histriophoca"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Pagophilus"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Phoca"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Phocoena"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Phocoenoides"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Physeter"] <- 5
mammal_mort$trophic_level[mammal_mort$Genus == "Pontoporia"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Pseudorca"] <- 4.5
mammal_mort$trophic_level[mammal_mort$Genus == "Leontocebus"] <- 2.5
mammal_mort$trophic_level[mammal_mort$Genus == "Stenella"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Trichechus"] <- 2
mammal_mort$trophic_level[mammal_mort$Genus == "Tursiops"] <- 4
mammal_mort$trophic_level[mammal_mort$Genus == "Zalophus"] <- 4

mammal_mort[is.na(mammal_mort$trophic_level),]
#---- birds

bird_diet0 <- read_csv('/Users/jgradym/Google Drive/Gibert Paper/Data/Diet/Birds/41559_2019_1070_MOESM3_ESM.csv')
unique(bird_diet0$TrophicLevel)

bird_diet1 <- bird_diet0 %>%
  rename(Species = Binomial) %>%
  select(Species, TrophicLevel)
bird_diet1$trophic_level <- NA
bird_diet1$trophic_level[bird_diet1$TrophicLevel == "Herbivore"] <- 2
bird_diet1$trophic_level[bird_diet1$TrophicLevel == "Carnivore"] <- 3
bird_diet1$trophic_level[bird_diet1$TrophicLevel == "Omnivore"] <- 2.5
bird_diet1$trophic_level[bird_diet1$TrophicLevel == "Scavenger"] <- 2

#combine
bird_mort0 <- mortality0 %>% filter(Group == "Bird")
bird_mort <- left_join(bird_mort0 , bird_diet1, by = "Species")
bird_mort
bird_mort[is.na(bird_mort$trophic_level),]

#-------- invertebrates

invert_mort0 <- mortality0 %>% filter(Group == "Invertebrate")
invert_genera <- as_tibble(unique(word(invert_mort0$Species), 1, sep = "_"))
invert_genera <- invert_genera %>%
  rename(genera = value)
write_csv(invert_genera, '~/Desktop/invert_genera.csv')

invert_diet <- read_csv('/Users/jgradym/Google Drive/Gibert Paper/Data/Diet/invert_diet.csv')

invert_mort0$genus <- word(invert_mort0$Species, 1, sep = " ")

invert_mort <- left_join(invert_mort0, invert_diet, by= "genus")

#fish 
fish_mort0 <- mortality0 %>% filter(Group == "Fish") #234

fish_genera <- unique(word(fish_mort0$Species, 1, sep = " "))
fish_spp <- fish_mort0$Species 

write.csv(fish_genera, "~/Desktop/fish_genera.csv")
write.csv(fish_spp, "~/Desktop/fish_spp.csv")

library(rfishbase)
library(taxize)
fish0 <- species_list(Class = "Actinopterygii")
fish1 <- fish0[fish0 %in% fish_spp]
missing_fish <- fish_mort0 %>%
  filter(Species %nin% fish0)
#write_csv(missing_fish, "~/Desktop/missing_fish.csv")

#manual fixes
fish_mort1 <- fish_mort0
fish_mort1$Species[fish_mort1$Species == "Clupea pallassii"] <- "Clupea pallasii" 
fish_mort1$Species[fish_mort1$Species == "Restrelliger kanagurta"] <- "Rastrelliger kanagurta"
fish_mort1$Species[fish_mort1$Species == "Restrelliger neglectus"] <- "Rastrelliger neglectus"
fish_mort1$Species[fish_mort1$Species == "Lethrinops longispinis"] <- "Lethrinops longipinnis"
fish_mort1$Species[fish_mort1$Species == "Pseudopeneus maculatus"] <- "Pseudupeneus maculatus"
fish_mort1$Species[fish_mort1$Species == "Sebastes paucispinus"] <- "Sebastes paucispinis"
fish_mort1$Species[fish_mort1$Species == "Bathylagus milleri"] <- "Pseudobathylagus milleri"
fish_mort1$Species[fish_mort1$Species == "Blennius pholis"] <- "Lipophrys pholis"
fish_mort1$Species[fish_mort1$Species == "Clupea pallasii"] <- "Clupea pallasii pallasii"
fish_mort1$Species[fish_mort1$Species == "Cynoglossus macrolepidus"] <- "Cynoglossus arel"
fish_mort1$Species[fish_mort1$Species == "Cynolebias adloffi"] <- "Austrolebias adloffi"
fish_mort1$Species[fish_mort1$Species == "Engraulis encrasicholus"] <- "Engraulis encrasicolus"
fish_mort1$Species[fish_mort1$Species == "Gadus minimus"] <- "Raniceps raninus"
fish_mort1$Species[fish_mort1$Species == "Gadus minutus"] <- "Trisopterus minutus"
fish_mort1$Species[fish_mort1$Species == "Gadus minitus"] <- "Trisopterus minutus"
fish_mort1$Species[fish_mort1$Species == "Haemulon plumieri"] <- "Haemulon plumierii"
fish_mort1$Species[fish_mort1$Species == "Haplochromis anaphyrmus"] <- "Mylochromis anaphyrmus"
fish_mort1$Species[fish_mort1$Species == "Haplochromis mloto"] <- "Copadichromis mloto"
fish_mort1$Species[fish_mort1$Species == "Lampanyctus regalis"] <- "Nannobrachium regale"
fish_mort1$Species[fish_mort1$Species == "Leiognathus splendens"] <- "Eubleekeria splendens"
fish_mort1$Species[fish_mort1$Species == "Leucichthys artedi"] <- "Coregonus artedi"
fish_mort1$Species[fish_mort1$Species == "Leucichthys sardinella"] <- "Coregonus sardinella"
fish_mort1$Species[fish_mort1$Species == "Lithrinus enigmaticus"] <- "Lethrinus enigmaticus"
fish_mort1$Species[fish_mort1$Species == "Merluccius gayi"] <- "Merluccius gayi gayi"
fish_mort1$Species[fish_mort1$Species == "Nemipterus bleekeri"] <- "Nemipterus bipunctatus"
fish_mort1$Species[fish_mort1$Species == "Nemipterus delagoe"] <- "Nemipterus bipunctatus"
fish_mort1$Species[fish_mort1$Species == "Nemipterus tolu"] <- "Nemipterus peronii"
fish_mort1$Species[fish_mort1$Species == "Pneumatophorus japonicus"] <- "Scomber japonicus"
fish_mort1$Species[fish_mort1$Species == "Pseudosciaena diacanthus"] <- "Protonibea diacanthus"
fish_mort1$Species[fish_mort1$Species == "Rastrelliger neglectus"] <- "Rastrelliger brachysoma"
fish_mort1$Species[fish_mort1$Species == "Sardinops caerrula"] <- "Sardinops sagax"
fish_mort1$Species[fish_mort1$Species == "Sardinops melanosticta"] <- "Sardinops sagax "
fish_mort1$Species[fish_mort1$Species == "Sebastes dalli"] <- "Sebastes dallii"
fish_mort1$Species[fish_mort1$Species == "Sebastes jorani"] <- "Sebastes jordani"
fish_mort1$Species[fish_mort1$Species == "Sebastes paucipinis"] <- "Sebastes paucispinis"
fish_mort1$Species[fish_mort1$Species == "Sebastes ruberrinus"] <- "Sebastes ruberrimus"
fish_mort1$Species[fish_mort1$Species == "Soela vulgaris"] <- "Solea solea"
fish_mort1$Species[fish_mort1$Species == "Stizostedion canadensis"] <- "Sander canadensis"
fish_mort1$Species[fish_mort1$Species == "Thunnus germo"] <- "Thunnus alalunga"
fish_mort1$Species[fish_mort1$Species == "Thunnus alaunga"] <- "Thunnus alalunga"
fish_mort1$Species[fish_mort1$Species == "Thunnus macoyi"] <- "Thunnus maccoyii"
fish_mort1$Species[fish_mort1$Species == "Tracharus japonicus"] <- "Trachurus japonicus"
fish_mort1$Species[fish_mort1$Species == "Tilapia esculenta"] <- "Oreochromis esculentus"
fish_mort1$Species[fish_mort1$Species == "Cheilodactylus macropterus"] <- "Nemadactylus macropterus"
fish_mort1$Species[fish_mort1$Species == "Cynolebias bellottii"] <- "Austrolebias bellottii"
fish_mort1$Species[fish_mort1$Species == "Cynoscion macdonaldi"] <- "Totoaba macdonaldi"
fish_mort1$Species[fish_mort1$Species == "Cynoscion nobilis"] <- "Atractoscion nobilis"
fish_mort1$Species[fish_mort1$Species == "Cynolebias wolterstarfii"] <- "Austrolebias wolterstorffi"
fish_mort1$Species[fish_mort1$Species == "Sardinops sagax "] <- "Sardinops sagax"
fish_mort1$Species[fish_mort1$Species == "Acipsnser fulvescens"] <- "Acipenser fulvescens"
fish_mort1$Species[fish_mort1$Species == "Aphinius fasciatus"] <- "Aphanius fasciatus"
fish_mort1$Species[fish_mort1$Species == "Centengraulis mysticetus"] <- "Cetengraulis mysticetus"
fish_mort1$Species[fish_mort1$Species == "Cololabis aira"] <- "Cololabis saira"
fish_mort1$Species[fish_mort1$Species == "Coryphaennoides acrolepis"] <- "Coryphaenoides acrolepis"
fish_mort1$Species[fish_mort1$Species == "Chelodactylus macropterus"] <- "Nemadactylus macropterus"
fish_mort1$Species[fish_mort1$Species == "Pseudoupeneus macularus"] <- "Pseudupeneus maculatus"

fish_mort2 <- fish_mort1
fish_mort2$TL <- estimate(fish_mort2$Species)$Troph
fish_mort2$Species[is.na(fish_mort2$TL)]

#combine mortality
fish_mort_fin <- fish_mort2 %>%
  rename(trophic_level = TL)
fish_mort_fin$type <- NA
fish_mort_fin$TrophicLevel <- NA
invert_mort$genus <- NULL
invert_mort$TrophicLevel <- NA
invert_mort$TrophicLevel <- as.character(invert_mort$TrophicLevel )
mammal_mort$Genus <- NULL
mammal_mort <- mammal_mort %>%
  rename(TrophicLevel = trophic_level)
endo_mort <- rbind(mammal_mort, bird_mort)
endo_mort$type <- NA
mort_TL <- rbind(endo_mort, invert_mort, fish_mort_fin)
write_csv(mort_TL, "~/Desktop/mortality_TL.csv" )
  
  
#attack rates
attack <- read_csv(file.path(gdrive_path, 'Lietal_oikos_2017_data.csv')) %>%
  select(predator.ana.group, predator.species, predator.mass.mg, temperature.degree.celcius, attack.rate, everything() ) %>%
  arrange(predator.ana.group, predator.species, attack.rate) %>%
  mutate(index = 1:451)

attack_verts <- attack %>% 
  filter(predator.ana.group == "vertebrate") %>%
  select(predator.species) %>%
  rename(Species = predator.species) 
attack_verts$TL <- estimate(attack_verts$Species)$Troph
missing_attack_vert <- unique(attack_verts$Species[is.na( attack_verts$TL)])
missing_attack_vert 
attack_verts$Species[attack_verts$Species == "Brachydanio rerio"] <- "Danio rerio"
attack_verts$Species[attack_verts$Species == "Perca ﬂuviatilis"] <- "Perca fluviatilis"
attack_verts$TL <- estimate(attack_verts$Species)$Troph 
attack_verts$pred_type <- NA
attack_verts$pred_type <- as.character(attack_verts$pred_type)
attack_verts_TL <- attack_verts %>%
  rename(predator.species = Species)

#Attack invertebrates 
attack_inverts <- attack %>% 
  filter(predator.ana.group == "invertebrate") %>%
  select(predator.species) 
attack_invert_spp <-  unique(attack_inverts$Species)
write.csv(unique(attack_inverts$Species), "~/Desktop/attack_invert_spp_new.csv")

attack_invert_TL0 <- read_csv('~/Desktop/attack_invert_spp.csv') %>%
  rename(predator.species = Species) 

attack_invert_TL <- left_join(attack_inverts, attack_invert_TL0 , by = "predator.species") 

attack_TL0 <- as_tibble(rbind(attack_invert_TL, attack_verts_TL)) %>%
  mutate(index = 1:451)

# Together
attack_TL <- left_join(attack, attack_TL0, by = "index") %>%
  select(TL, predator.ana.group, pred_type, everything())
write_csv(attack_TL, "~/Desktop/attack_TL.csv" )


mort_Ea <- mortality %>%
  filter(temp_range_genus >= 5, n_genus >= 5) %>%
  nest(-Genus) %>% # the group variable
  mutate(
    fit = map(data, ~ lm(log(mass_corr_mortality) ~ one_kT , data = .x)), #this is a regression (y = attack.rate, x = one_kT), but you could adapt to whatever
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied) 

attack_Ea <- mortality %>%
  filter(temp_range_genus >= 5, n_genus >= 5) %>%
  nest(-Genus) %>% # the group variable
  mutate(
    fit = map(data, ~ lm(log(mass_corr_mortality) ~ one_kT , data = .x)), #this is a regression (y = attack.rate, x = one_kT), but you could adapt to whatever
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied) 

