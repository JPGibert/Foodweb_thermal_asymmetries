path <- file.path('YOUR FOLDER HERE') 

library(tidyverse)
library(sf)
library(raster)
library(taxize)
library(rangeBuilder)
library(Hmisc, pos =100)

##### Read in Mortality data from MccCoy & Gillooly 
# correct bird names
mortality1 <- read_csv(file.path(path,'McCoy_Gillooly_Data.csv'))

#----------------------------------------------------------------------------------------------------
#---------Check Species Names and Update to match names for Birdlife International Range Maps--------
#----------------------------------------------------------------------------------------------------

#--------- Fix names using Birdlife International (id = 175); source of range maps

mort_birds <- mortality1[mortality1$Group == "Bird",]

result_bird <- mort_birds$Species %>% #using Taxize
  gnr_resolve(data_source_ids = 175, 
              with_canonical_ranks=T) #557

better_bird <- result_bird  %>%
  rename(Species = user_supplied_name) %>%
  dplyr::select(Species, matched_name2)

mort_birds2 <- mort_birds %>%
  left_join(better_bird, by = "Species")

#------- Another round of name error checks using library rangeBuilder
mort_birds2$new_name <- NA 
for(i in 1: length(mort_birds2$Species)) {
  try(mort_birds2$new_name[i] <- getAcceptedFromSynonym(mort_birds2$matched_name2[i], db = "birds")) 
  try(print(mort_birds2$new_name[i]))
}

mort_birds2$new_name <- gsub("_", " ",mort_birds2$new_name)
mort_birds2$new_name <- as.character(mort_birds2$new_name)
mort_birds2$new_name2 <- if_else(is.na(mort_birds2$new_name), mort_birds2$matched_name2, mort_birds2$new_name)
mort_birds2

#------ check
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
which(is.na(mort_birds2$new_name2))
#[1]  11 230 231 257 493 494 650 709 710


# -------- Some manual fixes, comparing to Birdlife international names

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
mort_birds2 <- mort_birds2[-which(mort_birds2$new_name2 == "Pagophila eburnea"),] #bad polygon data

#updated Species List
bird_spp <- unique(mort_birds2$new_name2)


#-----------------------------------------------------------------------------------------------------
#--------------------------------- Spatial Analysis - Get Centroid -----------------------------------
#-----------------------------------------------------------------------------------------------------

#---  Spatial projections

#Behrmann Equal Area Projection
ea_tf <- st_crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs") #equal area projection

# Read in bird shapefiles
birds <- st_read(file.path(path, 'bird_range_sample.gpkg')) #shortened to match Gillooly species

#---- Limit Birdlife Int species to those in mortality file
mort_birds3 <- birds %>%
  filter(SCINAME  %in% bird_spp) 

#--Convert to equal area projection, standardize column headers, exclude vagrants
reshape_birds <- function(x) {
  raw1 = st_transform(x, ea_tf) 
  raw2 = filter(raw1, ORIGIN != 4)
  raw3 = rename(raw2,binomial = SCINAME)
  #raw4 = rename(raw3, geom = SHAPE)
}


mort_birds4 <- reshape_birds(mort_birds3)
mort_birds4 


#----- for species with multiple maps, limit to largest polygon
mort_birds5 <- mort_birds4 %>%
  group_by(binomial) %>%
  filter(Shape_Area == max(Shape_Area)) # only keep largest polygon per species

mort_birds_shp <- mort_birds5

#----- Get Centroid
bird_centroid <- mort_birds_shp %>%
  group_by(binomial) %>%
  st_centroid() %>%
  dplyr::select(binomial, geom)

#---- Bind centroid to datafile
mort_birds_shp_centr <- cbind(mort_birds_shp , st_coordinates(st_centroid(mort_birds_shp$geom))) # get centroid coordinates


#---- Simplify contours
mort_birds_shp_centr_simp <- st_simplify(mort_birds_shp_centr, dTolerance = 10000) # simplifying countours to reduce file size and plot in reasonable time

#---- Plotting example 
ggplot(data = filter(mort_birds_shp_centr_simp,binomial == "Accipiter nisus" )) + #wolf
  geom_sf(aes(fill = binomial), fill = NA) + theme_bw() +
  geom_point(aes(x = X, y = Y), shape = 21,  fill = "blue", size = 4) +
  ggtitle("Eurasian Sparrowhawk")

#---- Transform back to lat long
lat_lon <- crs("+init=epsg:4326 +proj=longlat +ellps=WGS84")

#----- Raster lat long
raster_latlon <- raster(ncol = 360, nrow = 180, xmn = -180, xmx = 180, ymn = -90, ymx = 90,
                        crs = lat_lon)

bird_shp_lat_long <- st_transform(mort_birds_shp_centr_simp, lat_lon)
length(bird_shp_lat_long$binomial)

bird_centroid_lat_long0 <- st_transform(bird_centroid, lat_lon) #transform equal area centroid to get equivalent lat long


st_coordinates(st_centroid(bird_centroid$geom)) # in  meters, equal area projection
st_coordinates(st_centroid(bird_centroid_lat_long0$geom)) # in degrees, WGS84

mort_birds_shp_lat_lon <- cbind(bird_shp_lat_long,st_coordinates(st_centroid(bird_centroid_lat_long0$geom))) # get centroid coordinates
bird_centroid_lat_long_simp <- st_simplify(mort_birds_shp_lat_lon, dTolerance = 0)

# Equal Area - note latitude is about 45º
ggplot(data = filter(mort_birds_shp_centr_simp,binomial == "Accipiter nisus" )) + #coyote
  geom_sf(aes(fill = binomial), fill = NA) + theme_bw() +
  geom_point(aes(x = X, y = Y), shape = 21,  fill = "blue", size = 4) +
  ggtitle("Eurasian Sparrowhawk")

# Lat Long
ggplot(data = filter(mort_birds_shp_lat_lon ,binomial == "Accipiter nisus" )) + #coyote
  geom_sf(aes(fill = binomial), fill = NA) + theme_bw() +
  geom_point(aes(x = X.1, y = Y.1), shape = 21,  fill = "blue", size = 4) +
  ggtitle("Eurasian Sparrowhawk")


# Add Temp from World Clim for land, and from XXXX paper for ocean for 
#swimming birds (Alcidae, Shpheniscidae and marine going Phalacrocoracidae, Gaviidae, Podicipedidae)

# ocean 
ocean_temp <- raster(file.path(path, 'sstC_2006_2015.tif')) #Data compiled in 2019 XXXX paper
temp_ocean_latlon <- projectRaster(ocean_temp, raster_latlon)
temps_raster <- raster::extract(temp_ocean_latlon, data.frame(x = mort_birds_shp_lat_lon$X.1, y = mort_birds_shp_lat_lon$Y.1))
temps_raster_df <- as.data.frame(temps_raster, xy = T)

#--------------- Land Temperature from WorldClim -------------
#---- Get mean Temperature, latitude, longitude 
r <- raster::getData("worldclim", var = "bio", res = 10) 
coords <- data.frame(x= mort_birds_shp_lat_lon$X.1, y= mort_birds_shp_lat_lon$Y.1)
points <- SpatialPoints(coords, proj4string = r@crs)

values0 <- extract(r, points)
values1 <- values0[,1]/10 # Land Temp
values2 <- temps_raster_df$temps_raster #sea surface temp

#----- Combine Temps
df_bird <- as_tibble(cbind.data.frame(coordinates(points),values1, values2, mort_birds_shp_lat_lon$binomial))

# Species list of marine swimming birds from above families
sw_bird_spp <-read_csv(file.path(path,'sw_bird_names.csv'))

# Mean Temp = If bird is marine, use sea surface temp, otherwise WorldClim 
bird_temp <- df_bird %>%
  rename(Species = `mort_birds_shp_lat_lon$binomial`, 
         World_Clim_Temp_C = values1,
         Ocean_Temp_C = values2) %>%
  mutate(Mean_Temp_C = if_else(Species %in% sw_bird_spp$x, Ocean_Temp_C, World_Clim_Temp_C)) %>%
  dplyr::select(Mean_Temp_C, World_Clim_Temp_C,Ocean_Temp_C, Species)


bird_temp

#----------------------------------------------------------------------------------------------------
#---------Check Mammal Species Names and Update to match names for IUCN Range Maps--------
#----------------------------------------------------------------------------------------------------

#--------------------------- Mammal Temps -------------------------

#-------- Resolve Mammals Species names - standardize so match IUCN
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

# get accepted synonyms
mort_mammals2$new_name <- NA
for(i in 1: length(mort_mammals2$Species)) {
  try(mort_mammals2$new_name[i] <- getAcceptedFromSynonym(mort_mammals2$matched_name2[i], db = "mammals")) 
  try(print(mort_mammals2$new_name[i]))
}

#------ Remove "_" in names, keep correct name
mort_mammals2$new_name <- gsub("_", " ",mort_mammals2$new_name)
mort_mammals2$new_name <- as.character(mort_mammals2$new_name)
mort_mammals2$new_name2 <- if_else(is.na(mort_mammals2$new_name), mort_mammals2$matched_name2, mort_mammals2$new_name)
mort_mammals2


#-------check old vs new names

mammal_name_check <- as.data.frame(cbind(old_name = mort_mammals2$Species[mort_mammals2$Species != mort_mammals2$new_name2],
                                         new_name = mort_mammals2$new_name2[mort_mammals2$new_name2 != mort_mammals2$Species]))

# Manual Correction of names

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
mort_mammals2$new_name2[mort_mammals2$new_name2== "Phoca sibirica"] <- "Pusa sibirica"
mort_mammals2$new_name2[mort_mammals2$new_name2== "Physeter catodon"] <- "Physeter macrocephalus"
mort_mammals2$new_name2[mort_mammals2$new_name2== "Phoca groenlandica"] <- "Pagophilus groenlandicus"
mort_mammals2$new_name2[mort_mammals2$new_name2== "Otaria flavescens"] <- "Otaria byronia"
mort_mammals2$new_name2[mort_mammals2$new_name2== "Phoca fasciata"] <- "Histriophoca fasciata"
mammal_spp <- unique(mort_mammals2$new_name2)

mort_mammals2$new_name2[mort_mammals2$new_name2 != mort_mammals2$Species]
mort_mammals2$new_name2[mort_mammals2$Species == "Cephalophus monticola"]

mammal_spp <- unique(mort_mammals2$new_name2)[order(unique(mort_mammals2$new_name2))]

#-------------  Get centroids and associated temp 

#---- read in range shapefiles
#mammals_shp0 <- st_read(file.path(path, 'Data/MAMMALS/MAMMALS.shp'))
mammals_shp0 <- st_read(file.path(path, 'mammal_range_sample.shp'))

#------- Limit to species in mortality dataset
mort_mammals_shp1 <- mammals_shp0 %>%
  filter(binomial %in% mammal_spp)

#----- Format shapefiles
reshape_mamm <- function(x) { 
  raw = filter(x, origin != 4) #get rid of vagrants
  ready = st_transform(raw, ea_tf) #equal area projection
}

mort_mammals_shp2 <- reshape_mamm(mort_mammals_shp1)

#----  Simplify boundaries, remove invalid polygons
mort_mammals_shp <- mort_mammals_shp2   %>%
  st_simplify( dTolerance = 10000) %>% 
  group_by(binomial) %>%
  filter(SHAPE_Area == max(SHAPE_Area)) %>%
  ungroup() %>%
  mutate(validity = st_is_valid(geometry, reason = TRUE)) %>%
  filter(validity == "Valid Geometry") #173.8 

#----- Get centroid in equal area plot, add to dataset
mamm_centroid <- mort_mammals_shp %>%
  group_by(binomial) %>%
  st_centroid() %>%
  dplyr::select(binomial, geometry)

mort_mammals_shp_centr <- cbind(mort_mammals_shp , st_coordinates(st_centroid(mort_mammals_shp$geometry))) # get centroid coordinates

mort_mammals_shp_centr_simp <- st_simplify(mort_mammals_shp_centr, dTolerance = 10000) # for plotting check

#------- Plotting example of centroids

ggplot(data = filter(mort_mammals_shp_centr_simp,binomial == "Canis latrans" )) + #coyote
  geom_sf(aes(fill = binomial), fill = NA) + theme_bw() +
  geom_point(aes(x = X, y = Y), shape = 21,  fill = "blue", size = 4) +
  ggtitle("Canis latrans - Equal area plot")

#---- Transform back to lat long to match with World Clim

mamm_shp_lat_long <-  st_transform(mort_mammals_shp_centr_simp, lat_lon)
mamm_centroid_lat_long0 <- st_transform(mamm_centroid, lat_lon) #transform equal area centroid to get equivalent lat long

#---- Combine and remove NA's
mort_mammals_shp_lat_lon <- cbind(mamm_shp_lat_long, st_coordinates(st_centroid(mamm_centroid_lat_long0$geometry))) # get centroid coordinates
mort_mammals_shp_lat_lon <- mort_mammals_shp_lat_lon[!is.na(mort_mammals_shp_lat_lon$X.1),]


#------- Plotting example of centroids

   # Equal Area - centroid is about 45º
ggplot(data = filter(mort_mammals_shp_centr_simp, binomial == "Canis latrans" )) + #coyote
  geom_sf(aes(fill = binomial), fill = NA) + theme_bw() +
  geom_point(aes(x = X, y = Y), shape = 21,  fill = "blue", size = 4) +
  ggtitle("Canis latrans - Equal Area")

# Projected back to Lat lon

ggplot(data = filter(mort_mammals_shp_lat_lon ,binomial == "Canis latrans" )) + #coyote
  geom_sf(aes(fill = binomial), fill = NA) + theme_bw() +
  geom_point(aes(x = X.1, y = Y.1), shape = 21,  fill = "blue", size = 4) +
  ggtitle("Canis latrans - Degrees")

#----------- Get ocean temps associated with marine mammal centroids
temps_raster <- raster::extract(temp_ocean_latlon, data.frame(x = mort_mammals_shp_lat_lon$X.1, y = mort_mammals_shp_lat_lon$Y.1))
temps_raster_df <- as.data.frame(temps_raster, xy = T)
temps_raster_df  <- na.omit(temps_raster_df )
# from 

#---- from world clim -  extract temp from centroid
coords <- data.frame(x = mort_mammals_shp_lat_lon$X.1, y= mort_mammals_shp_lat_lon$Y.1)
points <- SpatialPoints(coords, proj4string = r@crs)

values0 <- extract(r, points)
values1 <-values0[,1]/10 # World Clim temps
values2 <- temps_raster_df$temps_raster # Ocean temps

#--- Combine rows
mammal_temp0 <- as_tibble(cbind.data.frame(coordinates(points), values1, values2, 
                                           mort_mammals_shp_lat_lon$binomial,
                                           mort_mammals_shp_lat_lon$family, 
                                           mort_mammals_shp_lat_lon$order_, 
                                           mort_mammals_shp_lat_lon$freshwater,
                                           mort_mammals_shp_lat_lon$terrestial))

#------- Some formatting, define marine species
capitalize(tolower(mammal_temp0$`mort_mammals_shp_lat_lon$family`))
Marine_fam <- c("Phocidae", "Otariidae","	Odobenidae", #marine mammal families
                "Delphinidae", "Monodontidae", "Phocoenidae", "Pontoporiidae", "Kogiidae", "Physeteridae", "Ziphidae",
                "Eschrichtiidae", "Balaenidae", "Cetotheriidae", "Balaenopteridae","Neobalaenidae")
mar_spp <- c("Lontra felina", "Enhydra lutris") #marine otters

mammal_temp0$`mort_mammals_shp_lat_lon$family`<-  capitalize(tolower(mammal_temp0$`mort_mammals_shp_lat_lon$family`))
Marine_fam2 <- capitalize(Marine_fam)
#-------- If terrestrial use Worldclime, otherwise marine temps
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


