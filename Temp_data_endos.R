gdrive_path <- file.path('/Users/jgradym/Google Drive/Gibert Paper') 


library(tidyverse)
library(sf)
library(raster)
library(lwgeom)
library(taxize)
library(fasterize)
library(RColorBrewer)
library(pryr)
library(rangeBuilder)
library(Hmisc)
library(compare)


ea_tf <- st_crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs") #equal area projection
lat_lon <- st_crs("+init=epsg:4326 +proj=longlat +ellps=WGS84") #conventional WGS84 lat long projection

raster110by110 <- raster(file.path(gdrive_path, "Data/raster110by110.tif")) #Load the raster used in analyses, ratio of endos to ectos
coord <- coordinates(raster110by110)
extent1 <- extent(raster110by110)
res1 = res(raster110by110)


##### Read in Gillooly data and Species
  # correct bird names
mortality1 <- read_csv(file.path(gdrive_path,'Data/McCoy_Gillooly_Data.csv'))
mortality1[duplicated(mortality1),] #interesting to look at duplicate rows
length(mortality1$Species) #2117
unique(mortality1$Group ) #6
length(mortality1$Species[mortality1$Group == "Multicellular plant"]) #348
sum(!is.na(mortality1$Species[mortality1$Group == "Multicellular plant"])) #348
  
#check numbers
length(mortality1$Species[mortality1$Group == "Bird"]) #779
length(unique(mortality1$Species[mortality1$Group == "Bird"])) #563 unique bird species

length(mortality1$Species[mortality1$Group == "Mammal"]) #524
length(unique(mortality1$Species[mortality1$Group == "Mammal"])) #323 unique mammal species

length(mortality1$Species) #2117
sum(is.na(mortality1$Species)) #0
mort_birds <- mortality1[mortality1$Group == "Bird",]

# Fix names using Birdlife International (id = 175); source of range maps
result_bird <- mort_birds$Species %>% #using Taxize
  gnr_resolve(data_source_ids = 175, 
              with_canonical_ranks=T) #557
better_bird <- result_bird  %>%
  rename(Species = user_supplied_name) %>%
  dplyr::select(Species, matched_name2)
mort_birds2 <- mort_birds %>%
  left_join(better_bird, by = "Species")

#Another round of name error checks using library rangeBuilder
mort_birds2$new_name <- NA 
for(i in 1: length(mort_birds2$Species)) {
  try(mort_birds2$new_name[i] <- getAcceptedFromSynonym(mort_birds2$matched_name2[i], db = "birds")) 
  try(print(mort_birds2$new_name[i]))
}
mort_birds2$new_name <- gsub("_", " ",mort_birds2$new_name)
mort_birds2$new_name <- as.character(mort_birds2$new_name)
mort_birds2$new_name2 <- if_else(is.na(mort_birds2$new_name), mort_birds2$matched_name2, mort_birds2$new_name)
mort_birds2


#check
mort_birds2<- mort_birds2 %>% 
  mutate(row = row_number())
mort_birds2$Species[mort_birds2$Species != mort_birds2$new_name2]
mort_birds2$row[mort_birds2$Species != mort_birds2$new_name2]
mort_birds2$new_name2[mort_birds2$new_name2 != mort_birds2$Species]
bird_name_check <- as.data.frame(cbind(row = mort_birds2$row[mort_birds2$Species != mort_birds2$new_name2], old_name = mort_birds2$Species[mort_birds2$Species != mort_birds2$new_name2],
           new_name = mort_birds2$new_name2[mort_birds2$new_name2 != mort_birds2$Species]
           ))
bird_name_check
length(bird_name_check$old_name) #146
which(is.na(mort_birds2$new_name2)) #no NA's

# Read in bird shapefiles
birds <- st_read(file.path(gdrive_path, 'Data/gibert_birds.gpkg')) #shortened to match Gillooly species
mort_birds$Species[mort_birds$Species %nin% birds$binomial] #214 spp  down not in both mortality and shapefiles - mostly spelling differences
#improvement?
mort_birds2$new_name2[mort_birds2$new_name2 %nin% birds$binomial] # 86 spp don't match - less but not great

#check

# Some manual fixes
which(is.na(mort_birds2$new_name2))
#[1]  11 230 231 257 493 494 650 709 710
mort_birds2$new_name2[mort_birds2$Species == "Accipiter n. nisus"] <- "Accipiter nisus"
mort_birds2$new_name2[11] 
mort_birds2[230,] 
mort_birds2$new_name2[mort_birds2$Species == "Columba fasciata"] <- "Patagioenas fasciata"
mort_birds2[257 ,] 
mort_birds2$new_name2[mort_birds2$Species == "Cygnus c. columbianus"] <- "Cygnus columbianus"
mort_birds2[493 ,] 
mort_birds2$new_name2[mort_birds2$Species == "Nyctea scandiaca"] <- "Bubo scandiacus"
mort_birds2[650 ,] 
mort_birds2$new_name2[mort_birds2$Species == "Seiurus novaeboracensis"] <- "Seiurus noveboracensis"
mort_birds2[ 710,]
mort_birds2$new_name2[mort_birds2$Species == "Sula (= Morus) bassanus"] <- "Morus bassanus"
mort_birds2$new_name2[mort_birds2$Species == "Ammodramus caudacutus"] <- "Ammospiza caudacuta"
mort_birds2$new_name2[mort_birds2$Species == "Ammodramus maritimus"] <- "Ammospiza maritima"
mort_birds2$new_name2[mort_birds2$Species == "Butorides striatus"] <- "Butorides striata"
mort_birds2$new_name2[mort_birds2$Species == "Carduelis flammea"] <- "Acanthis flammea"
mort_birds2$new_name2[mort_birds2$Species == "Acanthis hornemanni"] <- "Carduelis hornemanni"
mort_birds2$new_name2[mort_birds2$new_name2 == "Acanthis hornemanni"] <- "Carduelis hornemanni"
mort_birds2$new_name2[mort_birds2$Species == "Carduelis pinus"] <- "Spinus pinus"
mort_birds2$new_name2[mort_birds2$Species == "Carduelis psaltria"] <- "Spinus psaltria"
mort_birds2$new_name2[mort_birds2$Species == "Carduelis tristis"] <- "Spinus tristis"
mort_birds2$new_name2[mort_birds2$Species == "Carpodacus cassinii"] <- "Haemorhous cassinii"
mort_birds2$new_name2[mort_birds2$Species == "Carpodacus mexicanus"] <- "Haemorhous mexicanus"
mort_birds2$new_name2[mort_birds2$Species == "Carpodacus purpureus"] <- "Haemorhous purpureus"
mort_birds2$new_name2[mort_birds2$Species == "Dendragapus canadensis"] <- "Falcipennis canadensis"
mort_birds2$new_name2[mort_birds2$Species == "Diomedea chlororhynchos"] <- "Thalassarche chlororhynchos"
mort_birds2$new_name2[mort_birds2$Species == "Diomedea melanophris"] <- "Thalassarche melanophris"
mort_birds2$new_name2 <- gsub("Dendroica", "Setophaga", mort_birds2$new_name2)
#mort_birds2$new_name2[mort_birds2$Species == "Haematopus bachmani"] <- "Haematopus bachmani" real species apparently but not in Birdlife Int
#mort_birds2$new_name2[mort_birds2$Species == "Himantopus mexicanus"] <- 
mort_birds2$new_name2[mort_birds2$Species == "Oporornis tolmiei"] <- "Geothlypis tolmiei"
mort_birds2$new_name2[mort_birds2$new_name2== "Oporornis tolmiei"] <- "Geothlypis tolmiei"
mort_birds2$new_name2[mort_birds2$Species == "Oporornis formosus"] <- "Geothlypis formosa"
mort_birds2$new_name2[mort_birds2$Species == "Oporornis philadelphia"] <- "Geothlypis philadelphia"
mort_birds2$new_name2[mort_birds2$Species == "Parula americana"] <- "Setophaga americana"
mort_birds2$new_name2[mort_birds2$Species == "Parus atricapillus"] <- "Poecile atricapillus"
mort_birds2$new_name2[mort_birds2$Species == "Parus carolinensis"] <- "Poecile carolinensis"
mort_birds2$new_name2[mort_birds2$Species == "Parus cinctus"] <- "Poecile cinctus"
mort_birds2$new_name2[mort_birds2$Species == "Parus gambeli"] <- "Poecile gambeli"
mort_birds2$new_name2[mort_birds2$Species == "Parus inornatus"] <- "Baeolophus inornatus"
mort_birds2$new_name2[mort_birds2$Species == "Parus rufescens"] <- "Poecile rufescens"
mort_birds2$new_name2[mort_birds2$Species == "Parus sclateri"] <- "Poecile sclateri"
mort_birds2$new_name2[mort_birds2$Species == "Parus hudsonicus"] <- "Poecile hudsonicus"
mort_birds2$new_name2[mort_birds2$Species == "Pica nuttalli"] <- "Pica nutalli"
mort_birds2$new_name2[mort_birds2$Species == "Melozone fuscus"] <- "Melozone fusca"
mort_birds2$new_name2[mort_birds2$new_name2 == "Melozone fuscus"] <- "Melozone fusca"
mort_birds2$new_name2[mort_birds2$Species == "Polyborus plancus"] <- "Caracara plancus"
mort_birds2$new_name2[mort_birds2$Species == "Seiurus noveboracensis"] <- "Parkesia noveboracensis"
mort_birds2$new_name2[mort_birds2$new_name2 == "Seiurus noveboracensis"] <- "Parkesia noveboracensis"
mort_birds2$new_name2[mort_birds2$Species == "Spizella arborea"] <- "Spizelloides arborea"
mort_birds2$new_name2[mort_birds2$Species == "Spizelloides arborea"] <- "Passerella arborea"
mort_birds2$new_name2[mort_birds2$new_name2 == "Spizelloides arborea"] <- "Passerella arborea"
mort_birds2$new_name2[mort_birds2$Species == "Sterna nilotica"] <- "Gelochelidon nilotica"
mort_birds2$new_name2[mort_birds2$Species == "Stigmatopelia chinensis"] <- "Spilopelia chinensis"
mort_birds2$new_name2[mort_birds2$Species == "Vermivora celata"] <- "Leiothlypis celata"
mort_birds2$new_name2[mort_birds2$Species == "Vermivora luciae"] <- "Leiothlypis luciae"
mort_birds2$new_name2[mort_birds2$Species == "Vermivora ruficapilla"] <- "Leiothlypis ruficapilla"
mort_birds2$new_name2[mort_birds2$Species == "Vermivora virginiae"] <- "Leiothlypis virginiae"
mort_birds2$new_name2[mort_birds2$Species == "Wilsonia canadensis"] <- "Cardellina canadensis"
mort_birds2$new_name2[mort_birds2$Species == "Wilsonia citrina"] <- "Setophaga citrina"
mort_birds2$new_name2[mort_birds2$Species == "Wilsonia pusilla"] <- "Cardellina pusilla"
mort_birds2$new_name2[mort_birds2$Species == "Turdoides squamiceps"] <- "Argya squamiceps"
mort_birds2$new_name2[mort_birds2$Species == "Stigmatopelia chinensis"] <- "Spilopelia chinensis"
mort_birds2$new_name2[mort_birds2$new_name2 == "Stigmatopelia chinensis"] <- "Spilopelia chinensis"
mort_birds2$new_name2[mort_birds2$Species == "Vermivora peregrina"] <- "Leiothlypis peregrina"
which(mort_birds2[mort_birds2$new_name2 == "Pagophila eburnea",])
which(mort_birds2$new_name2 == "Pagophila eburnea") #515, 516
mort_birds2[which(mort_birds2$new_name2 == "Pagophila eburnea"),] #bad centroid
mort_birds2 <- mort_birds2[-which(mort_birds2$new_name2 == "Pagophila eburnea"),]
bird_spp <- unique(mort_birds2$new_name2)

bird_spp[bird_spp %nin% birds$binomial] #only 3 species don't match now
length(unique(bird_spp)) #545 birds now accounted for

# try bird stats again
mort_birds3 <- birds %>%
  filter(binomial %in% bird_spp)
length(unique(mort_birds3$binomial)) #543 of 546  species from Gillooly
length(mort_birds3$binomial) #2069 - lots of extra shapefiles

mort_birds4 <- mort_birds3 %>%
  group_by(binomial) %>%
  filter(Shape_Area == max(Shape_Area)) # only keep largest polygon per species
length(unique(mort_birds4$binomial)) #543
length(mort_birds4$binomial) #543


mort_birds_shp <- mort_birds4
length(mort_birds_shp$binomial )
#mort_birds_shp <- mort_birds_shp[-524,] # bad row

bird_centroid <- mort_birds_shp %>%
  group_by(binomial) %>%
  st_centroid() %>%
  dplyr::select(binomial, geom)

#select largest polygon where multiples exist - eg islands

mort_birds_shp_centr <- cbind(mort_birds_shp , st_coordinates(st_centroid(mort_birds_shp$geom))) # get centroid coordinates
length(mort_birds_shp_centr$binomial)
bird_spp[bird_spp %nin% mort_birds_shp_centr$binomial]

mort_birds_shp_centr_simp <- st_simplify(mort_birds_shp_centr, dTolerance = 10000) # simplifying countours to reduce file size and plot in reasonable time

#Plotting example - check that centroids look right
ggplot(data = filter(mort_birds_shp_centr_simp,binomial == "Turdus migratorius" )) + #wolf
  geom_sf(aes(fill = binomial), fill = NA) + theme_bw() +
  geom_point(aes(x = X, y = Y), shape = 21,  fill = "blue", size = 4) +
  ggtitle("American Robin")

#transform back to lat long
lat_lon <- crs("+init=epsg:4326 +proj=longlat +ellps=WGS84")
raster_wgs <- raster(file.path(gdrive_path, 'Data/raster_wgs.tif'))
bird_shp_lat_long <- st_transform(mort_birds_shp_centr_simp, lat_lon)
length(bird_shp_lat_long$binomial)

bird_centroid_lat_long0 <- st_transform(bird_centroid, lat_lon) #transform equal area centroid to get equivalent lat long
#bird_centroid_lat_long0 [526,] #Pagophila eburnea       EMPTY
length(bird_centroid_lat_long0$binomial)
st_coordinates(st_centroid(bird_centroid$geom)) # in  meters, equal area projection
st_coordinates(st_centroid(bird_centroid_lat_long0$geom)) # in degrees, WGS84

mort_birds_shp_lat_lon <- cbind(bird_shp_lat_long,st_coordinates(st_centroid(bird_centroid_lat_long0$geom))) # get centroid coordinates
bird_centroid_lat_long_simp <- st_simplify(mort_birds_shp_lat_lon, dTolerance = 0)

# Equal Area - note latitude is about 45º
ggplot(data = filter(mort_birds_shp_centr_simp,binomial == "Turdus migratorius" )) + #coyote
  geom_sf(aes(fill = binomial), fill = NA) + theme_bw() +
  geom_point(aes(x = X, y = Y), shape = 21,  fill = "blue", size = 4) +
  ggtitle("American Robin")


ggplot(data = filter(mort_birds_shp_lat_lon ,binomial == "Turdus migratorius" )) + #coyote
  geom_sf(aes(fill = binomial), fill = NA) + theme_bw() +
  geom_point(aes(x = X.1, y = Y.1), shape = 21,  fill = "blue", size = 4) +
  ggtitle("American Robin")


# Add Temp ffrom World Clim for land, and from Predator paper for ocean
lat_lon <- crs("+init=epsg:4326 +proj=longlat +ellps=WGS84")

ocean_temp <- raster(file.path(gdrive_path, 'Data/sstC_2006_2015.tif')) #Data compiled in 2019 Grady Predator paper
temp_ocean_latlon <- projectRaster(ocean_temp, raster_wgs)
temps_raster <- raster::extract(temp_ocean_latlon, data.frame(x = mort_birds_shp_lat_lon$X.1, y = mort_birds_shp_lat_lon$Y.1))
temps_raster_df <- as.data.frame(temps_raster, xy = T)

# from 
#from world clim
r <- getData("worldclim", var="bio", res=10)
str(r)
coords <- data.frame(x= mort_birds_shp_lat_lon$X.1, y= mort_birds_shp_lat_lon$Y.1)
length(coords$x)
points <- SpatialPoints(coords, proj4string = r@crs)
length(points)
values0 <- extract(r, points)
length(values0)
values1 <- values0[,1]/10
values2 <- temps_raster_df$temps_raster
plot(values1 ~values2)
abline(a = 0, b = 1, col = "red")
df_bird <- as_tibble(cbind.data.frame(coordinates(points),values1, values2, mort_birds_shp_lat_lon$binomial))

# swimming bird species that forage in the ocean - Alcids, penguins, some grebes, loons, and cormorants
sw_bird_spp <-read_csv(file.path(gdrive_path,'Data/sw_bird_names.csv'))
bird_temp <- df_bird %>%
  rename(Species = `mort_birds_shp_lat_lon$binomial`, 
         World_Clim_Temp_C = values1,
         Ocean_Temp_C = values2) %>%
  mutate(Mean_Temp_C = if_else(Species %in% sw_bird_spp$x, Ocean_Temp_C, World_Clim_Temp_C)) %>%
  dplyr::select(Mean_Temp_C, World_Clim_Temp_C,Ocean_Temp_C, Species)


write.csv(bird_temp, file.path(gdrive_path, "Data/bird_temp.csv"))


#st_write(mort_birds4, "/Users/jgradym/Google Drive/Gibert Paper/Data/gibert_birds.gpkg")

# Now Resolve Mammals
mort_mammals <- mortality1[mortality1$Group == "Mammal",]
unique(mort_mammals$Species) #323, some domestic
result_mammals <- mort_mammals$Species%>%
  gnr_resolve(data_source_ids = 174, 
              with_canonical_ranks=T) #557
better_mammals <- result_mammals  %>%
  rename(Species = user_supplied_name) %>%
  dplyr::select(Species, matched_name2)
mort_mammals2 <- mort_mammals %>%
  left_join(better_mammals, by = "Species")
length(mort_mammals2$Species) #525
length(mort_mammals$Species) #524
mort_mammals2$new_name <- NA
for(i in 1: length(mort_mammals2$Species)) {
  try(mort_mammals2$new_name[i] <- getAcceptedFromSynonym(mort_mammals2$matched_name2[i], db = "mammals")) 
  try(print(mort_mammals2$new_name[i]))
}
mort_mammals2$new_name <- gsub("_", " ",mort_mammals2$new_name)
mort_mammals2$new_name <- as.character(mort_mammals2$new_name)
mort_mammals2$new_name2 <- if_else(is.na(mort_mammals2$new_name), mort_mammals2$matched_name2, mort_mammals2$new_name)
mort_mammals2

########################
#check old vs new names
mort_mammals2$Species[mort_mammals2$Species != mort_mammals2$new_name2]
mort_mammals2$row[mort_mammals2$Species != mort_mammals2$new_name2]
mort_mammals2$new_name2[mort_mammals2$new_name2 != mort_mammals2$Species]
mammal_name_check <- as.data.frame(cbind(row = mort_mammals2$row[mort_mammals2$Species != mort_mammals2$new_name2],
                           old_name = mort_mammals2$Species[mort_mammals2$Species != mort_mammals2$new_name2],
                           new_name = mort_mammals2$new_name2[mort_mammals2$new_name2 != mort_mammals2$Species]))
mammal_name_check
which(is.na(mort_mammals2$new_name2)) #Alopex lagopus
mort_mammals2[which(is.na(mort_mammals2$new_name2)),] #no NA's

#[1] 13
mort_mammals2[13,] 
mort_mammals2$new_name2[383:400] <- mort_mammals2$Species[383:400]
mort_mammals2$new_name2[mort_mammals2$Species == "Alopex lagopus"] <- "Vulpes lagopus"
mort_mammals2$new_name2[mort_mammals2$Species == "Equus burchelli"] <- "Equus quagga"
mort_mammals2$new_name2[mort_mammals2$Species == "Hylobates syndactylus"] <- "Symphalangus syndactylus"
mort_mammals2$new_name2[mort_mammals2$Species == "Lobodon carcinophagus"] <- "Lobodon carcinophaga"
mort_mammals2$new_name2[mort_mammals2$Species == "Erethizon dorsata"] <- "Erethizon dorsatum"
mort_mammals2$new_name2[mort_mammals2$new_name2== "Erethizon dorsata"] <- "Erethizon dorsatum"
mort_mammals2$new_name2[mort_mammals2$Species == "Ochotona dalli"] <- "Ovis dalli"
mort_mammals2$new_name2[mort_mammals2$Species == "Saguinus caffer"] <- "Syncerus caffer"
mort_mammals2$new_name2[mort_mammals2$Species == "Parantechinus bilarni"] <- "Pseudantechinus bilarni"
mort_mammals2$new_name2[mort_mammals2$Species == "Paguma hispida"] <- "Phoca hispida"
mort_mammals2$new_name2[mort_mammals2$Species == "Phoca caspica"] <- "Pusa caspica"
mort_mammals2$new_name2[mort_mammals2$Species == "Phoca hispida"] <- "Pusa hispida"
mort_mammals2$new_name2[mort_mammals2$Species == "Cephalophus monticola"] <- "Philantomba monticola"
mort_mammals2$new_name2[mort_mammals2$new_name2== "Cephalophus monticola"] <- "Philantomba monticola"
mort_mammals2$new_name2[mort_mammals2$Species == "Cebus apella"] <- "Sapajus apella"
mort_mammals2$new_name2[mort_mammals2$Species == "Equus asinus"] <- "Equus africanus"
mort_mammals2$new_name2[mort_mammals2$Species == "Equus caballus"] <- "Equus ferus"
mort_mammals2$new_name2[mort_mammals2$Species == "Liomys adspersus"] <- "Heteromys adspersus"
mort_mammals2$new_name2[mort_mammals2$Species == "Saguinus nigricollis"] <- "Leontocebus nigricollis"
mort_mammals2$new_name2[mort_mammals2$Species == "Damaliscus korrigum"] <- "Damaliscus lunatus"


mort_mammals2$new_name2[mort_mammals2$new_name2 != mort_mammals2$Species]
mort_mammals2$new_name2[mort_mammals2$Species == "Cephalophus monticola"]

mammal_spp <- unique(mort_mammals2$new_name2)[order(unique(mort_mammals2$new_name2))]

##### Get centroids and associated temp # this a big file and I didn't share the path - got it from IUCN
  #Import data
sf_import <- function(x) { 
  raw = st_read(x)
  names(raw) = tolower(names(raw)) #standardize spelling
  raw2 = filter(raw, origin != 4) #get rid of vagrants
  ready = st_transform(raw2, ea_tf) #equal area projection
}
mammals_shp <- sf_import(file.path(gdrive_path, 'Data/MAMMALS/MAMMALS.shp'))

  #Limit to species in our data
#all_mamm_rich <- fastest_rich(terr_mammals)

mort_mammals2$new_name2[mort_mammals2$new_name2== "Phoca hispida"] <- "Pusa hispida"
mort_mammals2$new_name2[mort_mammals2$new_name2== "Phoca sibirica"] <- "Pusa sibirica"
mort_mammals2$new_name2[mort_mammals2$new_name2== "Physeter catodon"] <- "Physeter macrocephalus"
mort_mammals2$new_name2[mort_mammals2$new_name2== "Phoca groenlandica"] <- "Pagophilus groenlandicus"
mort_mammals2$new_name2[mort_mammals2$new_name2== "Otaria flavescens"] <- "Otaria byronia"
mort_mammals2$new_name2[mort_mammals2$new_name2== "Phoca fasciata"] <- "Histriophoca fasciata"
mammal_spp <- unique(mort_mammals2$new_name2)

mort_mammals_shp <- mammals_shp  %>%
  filter(binomial %in% mammal_spp)%>%
  #st_simplify( dTolerance = 10000) %>% # reduce border complexity for faster plotting, error checking
  group_by(binomial) %>%
  filter(shape_area == max(shape_area)) %>%
  ungroup() %>%
  mutate(validity = st_is_valid(geometry,reason = TRUE)) %>%
  filter(validity == "Valid Geometry") #173.8 
#rm(mammals_shp)
length(unique(mort_mammals_shp$binomial)) #314
mamm_names_shp <- unique(mort_mammals_shp$binomial)
mammal_spp[mammal_spp %nin% mamm_names_shp]  # species that don't match - check names

#rm(mammals_shp)
mamm_centroid <- mort_mammals_shp %>%
  group_by(binomial) %>%
  st_centroid() %>%
  dplyr::select(binomial, geometry)
#select largest polygon where multiples exist - eg islands
mort_mammals_shp_centr <- cbind(mort_mammals_shp , st_coordinates(st_centroid(mort_mammals_shp$geometry))) # get centroid coordinates
mammal_spp[mammal_spp %nin% mort_mammals_shp_centr$binomial] # only 4 missing

mort_mammals_shp_centr_simp <- st_simplify(mort_mammals_shp_centr, dTolerance = 10000) # for plotting check
  
#Plotting example of centroids
ggplot(data = filter(mort_mammals_shp_centr_simp,binomial == "Canis lupus" )) + #wolf
  geom_sf(aes(fill = binomial), fill = NA) + theme_bw() +
  geom_point(aes(x = X, y = Y), shape = 21,  fill = "blue", size = 4) +
  ggtitle("Canis lupus - Equal area plot")

ggplot(data = filter(mort_mammals_shp_centr_simp,binomial == "Canis latrans" )) + #coyote
  geom_sf(aes(fill = binomial), fill = NA) + theme_bw() +
  geom_point(aes(x = X, y = Y), shape = 21,  fill = "blue", size = 4) +
  ggtitle("Canis latrans - Equal area plot")

#transform back to lat long
mamm_centroid_lat_long0 <- st_transform(mamm_centroid, lat_lon) #transform equal area centroid to get equivalent lat long
  st_coordinates(st_centroid(mamm_centroid$geometry)) 

mort_mammals_shp_lat_lon <- cbind(mamm_shp_lat_long,st_coordinates(st_centroid(mamm_centroid_lat_long0$geometry))) # get centroid coordinates

# Equal Area - centroid is about 45º
ggplot(data = filter(mort_mammals_shp_centr_simp,binomial == "Canis latrans" )) + #coyote
  geom_sf(aes(fill = binomial), fill = NA) + theme_bw() +
  geom_point(aes(x = X, y = Y), shape = 21,  fill = "blue", size = 4) +
  ggtitle("Canis latrans - Equal Area")


ggplot(data = filter(mort_mammals_shp_lat_lon ,binomial == "Canis latrans" )) + #coyote
  geom_sf(aes(fill = binomial), fill = NA) + theme_bw() +
  geom_point(aes(x = X.1, y = Y.1), shape = 21,  fill = "blue", size = 4) +
  ggtitle("Canis latrans - Degrees")


temps_raster <- raster::extract(temp_ocean_latlon, data.frame(x = mort_mammals_shp_lat_lon$X.1, y = mort_mammals_shp_lat_lon$Y.1))
temps_raster_df <- as.data.frame(temps_raster, xy = T)

# from 
#from world clim
r <- getData("worldclim", var="bio", res=10)

coords <- data.frame(x= mort_mammals_shp_lat_lon$X.1, y= mort_mammals_shp_lat_lon$Y.1)
points <- SpatialPoints(coords, proj4string = r@crs)

values0 <- extract(r,points)
values1 <-values0[,1]/10
values2 <- temps_raster_df$temps_raster
plot(values1 ~values2) #some interesting differences
abline(a = 0, b = 1, col = "red")
mammal_temp0 <- as_tibble(cbind.data.frame(coordinates(points),values1,values2, mort_mammals_shp_lat_lon$binomial,mort_mammals_shp_lat_lon$family, 
                                  mort_mammals_shp_lat_lon$order_, mort_mammals_shp_lat_lon$freshwater,mort_mammals_shp_lat_lon$terrestial))

capitalize(tolower(mammal_temp0$`mort_mammals_shp_lat_lon$family`))
Marine_fam <- c("Phocidae", "Otariidae","	Odobenidae", #marine mammal families
                "Delphinidae", "Monodontidae", "Phocoenidae", "Pontoporiidae", "Kogiidae", "Physeteridae", "Ziphidae",
                "Eschrichtiidae", "Balaenidae", "Cetotheriidae", "Balaenopteridae","Neobalaenidae")
mar_spp <- c("Lontra felina", "Enhydra lutris") #marine otters
mammal_temp0$`mort_mammals_shp_lat_lon$family`<-  capitalize(tolower(mammal_temp0$`mort_mammals_shp_lat_lon$family`))
Marine_fam2 <- capitalize(Marine_fam)
mammal_temp <- mammal_temp0 %>%
  rename(Species = `mort_mammals_shp_lat_lon$binomial`, 
         Family = `mort_mammals_shp_lat_lon$family`,
         Order = `mort_mammals_shp_lat_lon$order_`, 
         Terrestrial =`mort_mammals_shp_lat_lon$terrestial`,
         Freshwater = `mort_mammals_shp_lat_lon$freshwater`,
         World_Clim_Temp_C = values1,
         Ocean_Temp_C = values2) %>%
  mutate(Mean_Temp_C = if_else(Family %in% Marine_fam | Species %in% mar_spp, Ocean_Temp_C, World_Clim_Temp_C)) %>%
  arrange(Order, Family, Species) %>%
  dplyr::select(Mean_Temp_C, World_Clim_Temp_C,Ocean_Temp_C, Order, Family, Species, Terrestrial, Freshwater)
  
head(mammal_temp)
mammal_temp <- mammal_temp %>% dplyr::select(Mean_Temp_C, Order, Family, Species)
write.csv(mammal_temp, file.path(gdrive_path,"Data/mammal_temp.csv"))
