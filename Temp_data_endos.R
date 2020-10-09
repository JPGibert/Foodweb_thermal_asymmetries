# libraries 
library(Hmisc)
library(sf)
library(raster)
library(taxize)
library(rangeBuilder)
library(phangorn)
library(ape)
library(tidyverse)



# Note: adjust file paths for reading data, saving plots to relevant locations
    # github data path
path <- file.path('/Users/jgradym/Documents/GitHub/Foodweb_thermal_asymmetries/Data') #Enter file path to data here


##### Read in Mortality data from MccCoy & Gillooly 
# correct bird names
mortality_mccoy <- read_csv(file.path(path,'McCoy_Gillooly_Data.csv'))

#----------------------------------------------------------------------------------------------------
#---------Check Species Names and Update to match names for Birdlife International Range Maps--------
#--------- Fix names using Birdlife International (id = 175); source of range maps

mortality1 <- mortality_mccoy 
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


# -------- Manual fixes, comparing to Birdlife international names

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


#-----------------------------------------------------------------------------------------------------
#--------------------------------- Spatial Analysis - Get Centroid -----------------------------------
#-----------------------------------------------------------------------------------------------------

#---  Spatial projections

#Behrmann Equal Area Projection
ea_tf <- st_crs("+proj=cea +lat_ts=30")

# Read in bird shapefiles
birds <- st_read('/Users/jgradym/Google Drive/Gibert Paper/Data/gibert_birds.gpkg') #shortened to match Gillooly species
#birds <- st_read(file.path(path, 'FILE_NAME.gpkg', layer="All_Species")) #if from Birdlife Int

length(birds$binomial)
mort_birds2$new_name2[mort_birds2$new_name2 %nin% birds$binomial] 

mort_bird_spp <- unique(mort_birds2$new_name2)
#---- Limit Birdlife Int species to those in mortality file

mort_birds3 <- birds %>%
  filter(binomial  %in% mort_bird_spp) 
length(mort_birds3$binomial)

#--Convert to equal area projection, standardize column headers, exclude vagrants
## For Full dataset from birdlife, remove pound sign
reshape_birds <- function(x) {
  raw1 = st_transform(x, ea_tf) 
  raw2 = filter(raw1, ORIGIN != 4)
  #raw3 = rename(raw2,binomial = SCINAME)
  #raw4 = rename(raw3, geom = SHAPE) 
}

mort_birds4 <- reshape_birds(mort_birds3)
length(mort_birds4$binomial) 

#----- for species with multiple maps - ie non-contiguous habitat -  limit to largest polygon
mort_birds5 <- mort_birds4 %>%
  group_by(binomial) %>%
  filter(Shape_Area == max(Shape_Area)) # only keep largest polygon per species

mort_birds_shp <- mort_birds5
length(mort_birds_shp$binomial)

#----- Get Centroid
bird_centroid <- mort_birds_shp %>%
  group_by(binomial) %>%
  st_centroid() %>%
  dplyr::select(binomial, geom)

#---- Bind centroid to datafile
mort_birds_shp_centr <- cbind(mort_birds_shp , st_coordinates(st_centroid(mort_birds_shp$geom))) # get centroid coordinates

#---- Simplify contours
mort_birds_shp_centr_simp <- st_simplify(mort_birds_shp_centr, dTolerance = 10000) # simplifying countours to reduce file size and plot in reasonable time
length(mort_birds_shp_centr_simp$binomial )

#---- Transform back to lat long
lat_lon <- crs("+init=epsg:4326 +proj=longlat +ellps=WGS84")

#----- Get latitude and longitude to match with wordlclim
raster_latlon <- raster(ncol = 360, nrow = 180, xmn = -180, xmx = 180, ymn = -90, ymx = 90,
                        crs = lat_lon)

bird_shp_lat_long <- st_transform(mort_birds_shp_centr_simp, lat_lon)
length(bird_shp_lat_long$binomial)

bird_centroid_lat_long0 <- st_transform(bird_centroid, lat_lon) #transform equal area centroid to get equivalent lat long

# get centroid in equal area projection (meters) 
st_coordinates(st_centroid(bird_centroid$geom)) # 

# convert centroid back to latitude/longitude 
st_coordinates(st_centroid(bird_centroid_lat_long0$geom)) # convert back in degrees, WGS84

mort_birds_shp_lat_lon <- cbind(bird_shp_lat_long,st_coordinates(st_centroid(bird_centroid_lat_long0$geom))) # get centroid coordinates
bird_centroid_lat_long_simp <- st_simplify(mort_birds_shp_lat_lon, dTolerance = 0.1)

# Equal Area - note latitude is about 45º
ggplot(data = filter(mort_birds_shp_centr_simp, binomial == "Accipiter nisus" )) + #coyote
  geom_sf(aes(fill = binomial), fill = NA) + theme_bw() +
  geom_point(aes(x = X, y = Y), shape = 21,  fill = "blue", size = 4) +
  ggtitle("Eurasian Sparrowhawk")

# Lat Long
ggplot(data = filter(mort_birds_shp_lat_lon, binomial == "Accipiter nisus" )) + #coyote
  geom_sf(aes(fill = binomial), fill = NA) + theme_bw() +
  geom_point(aes(x = X.1, y = Y.1), shape = 21,  fill = "blue", size = 4) +
  ggtitle("Eurasian Sparrowhawk")

# Add Temp from World Clim for land, and from XXXX paper for ocean for 
#swimming birds (Alcidae, Shpheniscidae and marine going Phalacrocoracidae, Gaviidae, Podicipedidae)

# ocean 
ocean_temp <- raster(file.path(path, 'sstC_2006_2015.tif')) #Satellite Data, 10 year mean, compiled in 2019 Grady et al Science
plot(ocean_temp)
temp_ocean_latlon <- projectRaster(ocean_temp, raster_latlon)
temps_raster <- raster::extract(temp_ocean_latlon, data.frame(x = mort_birds_shp_lat_lon$X.1, y = mort_birds_shp_lat_lon$Y.1))
temps_raster_df <- as.data.frame(temps_raster, xy = T)

#--------------- Land Temperature from WorldClim -------------
#---- Get mean Temperature, latitude, longitude 
r <- raster::getData("worldclim", var = "bio", res = 10) #also avavilable at https://www.worldclim.org/data/worldclim21.html
coords <- data.frame(x= mort_birds_shp_lat_lon$X.1, y= mort_birds_shp_lat_lon$Y.1)
points <- SpatialPoints(coords, proj4string = r@crs)

values0 <- raster::extract(r, points)
values1 <- values0[,1]/10 # Land Temp
values2 <- temps_raster_df$temps_raster #sea surface temp

#----- Combine Temps
df_bird <- as_tibble(cbind.data.frame(coordinates(points),values1, values2, mort_birds_shp_lat_lon$binomial))
df_bird 

# Mean Temp = If bird is marine, use sea surface temp, otherwise WorldClim 
bird_temp <- df_bird %>%
  rename(Species = `mort_birds_shp_lat_lon$binomial`, 
         World_Clim_Temp_C = values1,
         Ocean_Temp_C = values2) %>%
  mutate(Mean_Temp_C = if_else(is.na(World_Clim_Temp_C), Ocean_Temp_C, World_Clim_Temp_C)) %>%
  dplyr::select(Mean_Temp_C, World_Clim_Temp_C,Ocean_Temp_C, Species)

bird_temp_simp <- bird_temp %>%
  dplyr::select(Species, Mean_Temp_C)

mort_birds2b <- mort_birds2 %>%
  rename(old_name = Species, Species = new_name2)

# remove old columns
mort_birds2b$new_name <- NULL
mort_birds2b$matched_name2 <- NULL
mort_birds2b$row <- NULL
mort_birds2b$Species[mort_birds2b$Species %nin%  bird_temp$Species]

mort_bird3 <- mort_birds2b %>%
  left_join(bird_temp_simp, by = "Species") 

# match bird names to phylogeny

#update birds 
bird_tree <- read.tree('/Users/jgradym/Google Drive/Phylo_trees/Upham_Phylos/Birds/mnt/data/projects/birdphylo/Tree_sets/Stage1_full_data/CombinedTrees/EricsonStage1Full_1.tre')[[1]] #read in first tree to get names
str(bird_tree)

mort_bird3$Species <- gsub(" ", "_", mort_bird3$Species)
mort_bird3$Species[mort_bird3$Species %nin% bird_tree$tip.label]

#manually update bird names
mort_bird3$Species[mort_bird3$Species == "Acanthis_flammea"] <- "Carduelis_flammea"
mort_bird3$Species[mort_bird3$Species == "Ammospiza_caudacuta"] <- "Ammodramus_caudacutus"
mort_bird3$Species[mort_bird3$Species == "Ammospiza_maritima" ] <- "Ammodramus_maritimus"
mort_bird3$Species[mort_bird3$Species == "Amphispiza_quinquestriata"] <- "Aimophila_quinquestriata"
mort_bird3$Species[mort_bird3$Species == "Anser_caerulescens"] <- "Chen_caerulescens"
mort_bird3$Species[mort_bird3$Species == "Anser_canagicus"] <- "Chen_canagica"
mort_bird3$Species[mort_bird3$Species == "Anser_rossii"] <- "Chen_rossii"
mort_bird3$Species[mort_bird3$Species == "Antigone_canadensis"] <- "Grus_canadensis"
mort_bird3$Species[mort_bird3$Species == "Antrostomus_vociferus"] <- "Caprimulgus_vociferus"
mort_bird3$Species[mort_bird3$Species == "Ardenna_gravis"] <- "Puffinus_gravis"
mort_bird3$Species[mort_bird3$Species == "Ardenna_pacifica"] <- "Puffinus_pacificus"
mort_bird3$Species[mort_bird3$Species == "Ardenna_tenuirostris"] <- "Puffinus_tenuirostris"
#mort_bird3$Species[mort_bird3$Species == "Argya_squamiceps"] <- "" # not available
mort_bird3$Species[mort_bird3$Species == "Bubo_scandiacus"] <- "Bubo_scandiaca"
#mort_bird3$Species[mort_bird3$Species == "Calothorax_lucifer"] <- "" # not available
mort_bird3$Species[mort_bird3$Species == "Cardellina_canadensis"] <- "Wilsonia_canadensis"
mort_bird3$Species[mort_bird3$Species == "Cardellina_pusilla"] <- "Wilsonia_pusilla"
#mort_bird3$Species[mort_bird3$Species == "Circus_macrourus"] <- "" # not available
mort_bird3$Species[mort_bird3$Species == "Clanga_pomarina"] <- "Aquila_pomarina"
mort_bird3$Species[mort_bird3$Species == "Dryobates_pubescens"] <- "Picoides_pubescens"
mort_bird3$Species[mort_bird3$Species == "Dryobates_scalaris"] <- "Picoides_scalaris"
mort_bird3$Species[mort_bird3$Species == "Falcipennis_canadensis"] <- "Dendragapus_canadensis"
mort_bird3$Species[mort_bird3$Species == "Gelochelidon_nilotica"] <- "Sterna_nilotica"
mort_bird3$Species[mort_bird3$Species == "Geothlypis_formosa"] <- "Oporornis_formosus"
mort_bird3$Species[mort_bird3$Species == "Geothlypis_philadelphia"] <- "Oporornis_philadelphia"
mort_bird3$Species[mort_bird3$Species == "Geothlypis_tolmiei"] <- "Oporornis_tolmiei"
mort_bird3$Species[mort_bird3$Species == "Haemorhous_purpureus"] <- "Carpodacus_purpureus"
mort_bird3$Species[mort_bird3$Species == "Hieraaetus_wahlbergi"] <- "Aquila_wahlbergi"
mort_bird3$Species[mort_bird3$Species == "Hydrobates_furcatus"] <- "Oceanodroma_furcata"
mort_bird3$Species[mort_bird3$Species == "Hydrobates_homochroa"] <- "Oceanodroma_homochroa"
mort_bird3$Species[mort_bird3$Species == "Hydrobates_leucorhous"] <- "Oceanodroma_leucorhoa"
mort_bird3$Species[mort_bird3$Species == "Hydrobates_tristrami"] <- "Oceanodroma_tristrami"
mort_bird3$Species[mort_bird3$Species == "Hydroprogne_caspia"] <- "Sterna_caspia"
mort_bird3$Species[mort_bird3$Species == "Hylatomus_pileatus"] <- "Dryocopus_pileatus"
mort_bird3$Species[mort_bird3$Species == "Leiothlypis_celata"] <- "Vermivora_celata"
mort_bird3$Species[mort_bird3$Species == "Leiothlypis_luciae"] <- "Vermivora_luciae"
mort_bird3$Species[mort_bird3$Species == "Leiothlypis_peregrina"] <- "Vermivora_peregrina"
mort_bird3$Species[mort_bird3$Species == "Leiothlypis_ruficapilla"] <- "Vermivora_ruficapilla"
mort_bird3$Species[mort_bird3$Species == "Leiothlypis_virginiae"] <- "Vermivora_virginiae"
mort_bird3$Species[mort_bird3$Species == "Leuconotopicus_borealis"] <- "Picoides_borealis"
mort_bird3$Species[mort_bird3$Species == "Leuconotopicus_villosus"] <- "Picoides_villosus"
mort_bird3$Species[mort_bird3$Species == "Mareca_americana"] <- "Anas_americana"
mort_bird3$Species[mort_bird3$Species == "Mareca_penelope"] <- "Anas_penelope"
mort_bird3$Species[mort_bird3$Species == "Mareca_strepera"] <- "Anas_strepera"
mort_bird3$Species[mort_bird3$Species == "Melanitta_deglandi"] <- "Melanitta_fusca"
mort_bird3$Species[mort_bird3$Species == "Melanitta_stejnegeri"] <- "Melanitta_fusca"
mort_bird3$Species[mort_bird3$Species == "Melozone_aberti"] <- "Pipilo_aberti"
mort_bird3$Species[mort_bird3$Species == "Melozone_fusca"] <- "Pipilo_fuscus"
mort_bird3$Species[mort_bird3$Species == "Onychoprion_fuscatus"] <- "Sterna_fuscata"
mort_bird3$Species[mort_bird3$Species == "Onychoprion_lunatus"] <- "Sterna_lunata"
mort_bird3$Species[mort_bird3$Species == "Parkesia_noveboracensis"] <- "Coturnicops_noveboracensis"
mort_bird3$Species[mort_bird3$Species == "Passerella_arborea"] <- "Spizella_arborea"
mort_bird3$Species[mort_bird3$Species == "Peucaea_aestivalis"] <- "Aimophila_aestivalis"
mort_bird3$Species[mort_bird3$Species == "Pica_nutalli"] <- "Pica_nuttalli"
mort_bird3$Species[mort_bird3$Species == "Poecile_atricapillus"] <- "Parus_atricapillus"
mort_bird3$Species[mort_bird3$Species == "Poecile_cinctus"] <- "Parus_cinctus"
mort_bird3$Species[mort_bird3$Species == "Poecile_gambeli"] <- "Parus_gambeli"
mort_bird3$Species[mort_bird3$Species == "Poecile_hudsonicus"] <- "Parus_hudsonicus"
mort_bird3$Species[mort_bird3$Species == "Poecile_rufescens"] <- "Parus_rufescens"
mort_bird3$Species[mort_bird3$Species == "Poecile_sclateri"] <- "Parus_sclateri"
#mort_bird3$Species[mort_bird3$Species == "Procelsterna_cerulea"] <- "" # not available in phylogeny (checked Anous cerulea)
mort_bird3$Species[mort_bird3$Species == "Psiloscops_flammeolus"] <- "Otus_flammeolus"
mort_bird3$Species[mort_bird3$Species == "Poecile_carolinensis"] <- "Parus_carolinensis"
#mort_bird3$Species[mort_bird3$Species == "Pterodroma_leucoptera"] <- "" # not available
#mort_bird3$Species[mort_bird3$Species == "Puffinus_auricularis"] <- "" # not available
mort_bird3$Species[mort_bird3$Species == "Selasphorus_calliope"] <- "Stellula_calliope"
mort_bird3$Species[mort_bird3$Species == "Setophaga_americana"] <- "Parula_americana"
mort_bird3$Species[mort_bird3$Species == "Setophaga_caerulescens"] <- "Dendroica_caerulescens"
mort_bird3$Species[mort_bird3$Species == "Setophaga_castanea"] <- "Dendroica_castanea"
mort_bird3$Species[mort_bird3$Species == "Setophaga_coronata"] <- "Dendroica_coronata"
mort_bird3$Species[mort_bird3$Species == "Setophaga_discolor"] <- "Dendroica_discolor"
mort_bird3$Species[mort_bird3$Species == "Setophaga_citrina"] <- "Wilsonia_citrina"
mort_bird3$Species[mort_bird3$Species == "Setophaga_flavescens"] <- "Dendroica_dominica"
mort_bird3$Species[mort_bird3$Species == "Setophaga_fusca"] <- "Dendroica_fusca"
mort_bird3$Species[mort_bird3$Species == "Setophaga_kirtlandii"] <- "Dendroica_kirtlandii"
mort_bird3$Species[mort_bird3$Species == "Setophaga_magnolia"] <- "Dendroica_magnolia"
mort_bird3$Species[mort_bird3$Species == "Setophaga_palmarum"] <- "Dendroica_palmarum"
mort_bird3$Species[mort_bird3$Species == "Setophaga_pensylvanica"] <- "Dendroica_pensylvanica"
mort_bird3$Species[mort_bird3$Species == "Setophaga_petechia"] <- "Dendroica_petechia"
mort_bird3$Species[mort_bird3$Species == "Setophaga_pinus"] <- "Dendroica_pinus"
mort_bird3$Species[mort_bird3$Species == "Setophaga_striata"] <- "Dendroica_striata"
mort_bird3$Species[mort_bird3$Species == "Setophaga_tigrina"] <- "Dendroica_tigrina"
mort_bird3$Species[mort_bird3$Species == "Setophaga_townsendi"] <- "Dendroica_townsendi"
mort_bird3$Species[mort_bird3$Species == "Setophaga_virens"] <- "Dendroica_virens"
mort_bird3$Species[mort_bird3$Species == "Spatula_clypeata"] <- "Anas_clypeata"
mort_bird3$Species[mort_bird3$Species == "Spatula_cyanoptera"] <- "Anas_cyanoptera"
mort_bird3$Species[mort_bird3$Species == "Spatula_discors"] <- "Anas_discors"
mort_bird3$Species[mort_bird3$Species == "Spilopelia_chinensis"] <- "Stigmatopelia_chinensis"
mort_bird3$Species[mort_bird3$Species == "Spinus_pinus"] <- "Carduelis_pinus"
mort_bird3$Species[mort_bird3$Species == "Spinus_psaltria"] <- "Carduelis_psaltria"
mort_bird3$Species[mort_bird3$Species == "Spinus_tristis"] <- "Carduelis_tristis"
mort_bird3$Species[mort_bird3$Species == "Sternula_albifrons"] <- "Sterna_albifrons"
mort_bird3$Species[mort_bird3$Species == "Sternula_antillarum"] <- "Sterna_antillarum"
#mort_bird3$Species[mort_bird3$Species == "Thalassarche_melanophris"] <-  not available
mort_bird3$Species[mort_bird3$Species == "Thalasseus_elegans"] <- "Sterna_elegans"
mort_bird3$Species[mort_bird3$Species == "Thalasseus_maximus"] <- "Sterna_maxima"
mort_bird3$Species[mort_bird3$Species == "Thalasseus_sandvicensis"] <- "Sterna_sandvicensis"
mort_bird3$Species[mort_bird3$Species == "Tringa_semipalmata"] <- "Catoptrophorus_semipalmatus"
mort_bird3$Species[mort_bird3$Species == "Vermivora_cyanoptera"] <- "Vermivora_pinus"
mort_bird3$Species[mort_bird3$Species == "Haemorhous_cassinii"] <- "Carpodacus_cassinii"
mort_bird3$Species[mort_bird3$Species == "Haemorhous_mexicanus"] <- "Carpodacus_mexicanus"

mort_bird3$Species[mort_bird3$Species %nin% bird_tree$tip.label] # a few could not be included

mort_bird_fin <- mort_bird3 %>%
  filter(Species %in% bird_tree$tip.label) %>%
  dplyr::select(Group, Species, dry_mass_g, Mean_Temp_C, Mortality_yr) %>%
     rename(Temp_C = Mean_Temp_C)
mort_bird_fin



#----------------------------------------------------------------------------------------------------
#---------Check Mammal Species Names and Update to match names for IUCN Range Maps--------
#----------------------------------------------------------------------------------------------------

#--------------------------- Mammal Temps -------------------------

#-------- Resolve Mammals Species names - standardize so match IUCN range maps
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

mort_mammals3 <- mort_mammals2
#-------check old vs new names

mammal_name_check <- as.data.frame(cbind(old_name = mort_mammals3$Species[mort_mammals3$Species != mort_mammals3$new_name2],
                                         new_name = mort_mammals3$new_name2[mort_mammals3$new_name2 != mort_mammals3$Species]))
mammal_name_check 

# Manual Correction of names

mort_mammals3$new_name2[383:400] <- mort_mammals3$Species[383:400]
mort_mammals3$new_name2[mort_mammals3$Species == "Alopex lagopus"] <- "Vulpes lagopus"
mort_mammals3$new_name2[mort_mammals3$Species == "Equus burchelli"] <- "Equus quagga"
mort_mammals3$new_name2[mort_mammals3$Species == "Hylobates syndactylus"] <- "Symphalangus syndactylus"
mort_mammals3$new_name2[mort_mammals3$Species == "Lobodon carcinophagus"] <- "Lobodon carcinophaga"
mort_mammals3$new_name2[mort_mammals3$Species == "Erethizon dorsata"] <- "Erethizon dorsatum"
mort_mammals3$new_name2[mort_mammals3$new_name2== "Erethizon dorsata"] <- "Erethizon dorsatum"
mort_mammals3$new_name2[mort_mammals3$Species == "Ochotona dalli"] <- "Ovis dalli"
mort_mammals3$new_name2[mort_mammals3$Species == "Saguinus caffer"] <- "Syncerus caffer"
mort_mammals3$new_name2[mort_mammals3$Species == "Parantechinus bilarni"] <- "Pseudantechinus bilarni"
mort_mammals3$new_name2[mort_mammals3$Species == "Paguma hispida"] <- "Phoca hispida"
mort_mammals3$new_name2[mort_mammals3$Species == "Phoca caspica"] <- "Pusa caspica"
mort_mammals3$new_name2[mort_mammals3$new_name2 == "Phoca hispida"] <- "Pusa hispida"
mort_mammals3$new_name2[mort_mammals3$Species == "Cephalophus monticola"] <- "Philantomba monticola"
mort_mammals3$new_name2[mort_mammals3$new_name2== "Cephalophus monticola"] <- "Philantomba monticola"
mort_mammals3$new_name2[mort_mammals3$Species == "Cebus apella"] <- "Sapajus apella"
mort_mammals3$new_name2[mort_mammals3$Species == "Equus asinus"] <- "Equus africanus"
mort_mammals3$new_name2[mort_mammals3$Species == "Equus caballus"] <- "Equus ferus"
mort_mammals3$new_name2[mort_mammals3$Species == "Liomys adspersus"] <- "Heteromys adspersus"
mort_mammals3$new_name2[mort_mammals3$Species == "Saguinus nigricollis"] <- "Leontocebus nigricollis"
mort_mammals3$new_name2[mort_mammals3$Species == "Damaliscus korrigum"] <- "Damaliscus lunatus"
mort_mammals3$new_name2[mort_mammals3$new_name2== "Phoca sibirica"] <- "Pusa sibirica"
mort_mammals3$new_name2[mort_mammals3$new_name2== "Physeter catodon"] <- "Physeter macrocephalus"
mort_mammals3$new_name2[mort_mammals3$new_name2== "Phoca groenlandica"] <- "Pagophilus groenlandicus"
mort_mammals3$new_name2[mort_mammals3$new_name2== "Otaria flavescens"] <- "Otaria byronia"
mort_mammals3$new_name2[mort_mammals3$new_name2== "Phoca fasciata"] <- "Histriophoca fasciata"
mammal_spp <- unique(mort_mammals3$new_name2[order(mort_mammals3$new_name2)])

mort_mammals3$new_name2[mort_mammals3$new_name2 == "Phoca hispida"] 
mort_mammals3$new_name2[mort_mammals3$new_name2 != mort_mammals3$Species] #changed names
mort_mammals3$new_name2[mort_mammals3$Species == "Cephalophus monticola"]

mort_mammals3 <-  mort_mammals3  %>% 
  rename(old_name = Species) %>%
  rename(Species = new_name2)
mammal_spp <- unique(mort_mammals3$Species)[order(unique(mort_mammals3$Species))]
mammal_spp #320

#
#-------------  Get centroids and associated temp 

#---- read in range shapefiles
mammals_shp0 <- st_read('/Users/jgradym/Google Drive/Gibert Paper/Data/Big_data/MAMMALS/MAMMALS.shp')
#mammals_shp0 <- st_read(file.path(path, 'mammal_range_sample.shp'))

#------- Limit to species in mortality dataset
mort_mammals_shp1 <- mammals_shp0 %>%
  filter(binomial %in% mammal_spp)
unique(mort_mammals_shp1$binomial) #316
rm(mammals_shp0)

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
  filter(validity == "Valid Geometry") 

#----- Get centroid in equal area plot, add to dataset
mamm_centroid <- mort_mammals_shp %>%
  group_by(binomial) %>%
  st_centroid() %>%
  dplyr::select(binomial, geometry)

mort_mammals_shp_centr <- cbind(mort_mammals_shp , st_coordinates(st_centroid(mort_mammals_shp$geometry))) # get centroid coordinates

mort_mammals_shp_centr_simp <- st_simplify(mort_mammals_shp_centr, dTolerance = 10000) # for plotting check

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

values0 <- raster::extract(r, points)
values1 <-values0[,1]/10 # World Clim mean temps
values2 <- temps_raster_df$temps_raster # Ocean temps

#--- Combine data
mammal_temp0 <- as_tibble(cbind.data.frame(coordinates(points), values1, values2, 
                                           mort_mammals_shp_lat_lon$binomial,
                                           mort_mammals_shp_lat_lon$family, 
                                           mort_mammals_shp_lat_lon$order_))

#-------- If terrestrial use Worldclim, otherwise marine temps
mammal_temp <- mammal_temp0 %>%
  rename(Species = `mort_mammals_shp_lat_lon$binomial`, 
         Family = `mort_mammals_shp_lat_lon$family`,
         Order = `mort_mammals_shp_lat_lon$order_`, 
         World_Clim_Temp_C = values1, 
         Ocean_Temp_C = values2) %>%
  mutate(Mean_Temp_C = if_else(is.na(World_Clim_Temp_C), Ocean_Temp_C, World_Clim_Temp_C)) %>%
  arrange(Order, Family, Species) %>%
  dplyr::select(Mean_Temp_C, World_Clim_Temp_C,Ocean_Temp_C, Order, Family, Species)

head(mammal_temp )

mammal_temp <- mammal_temp %>% dplyr::select(Mean_Temp_C, Order, Family, Species)

unique(mammal_temp$Species) #314
unique(mort_mammals$Species)

missing_mamm <- mort_mammals3$Species[mort_mammals3$Species %nin% mammal_temp$Species] 
missing_mamm #domesticated animals not included

mammal_temp_simp <- mammal_temp %>%
  dplyr::select(Species, Mean_Temp_C) %>%
  arrange(Species)

mort_mammals_4 <- mort_mammals3 %>%
  filter(Species %nin% missing_mamm) %>%
  left_join(mammal_temp_simp, by = "Species")
  

#####################################################
#------------ match species names to phylogeny --------
#######################################################


# read in mammal tree 
mamm_tree_files <- list.files('/Users/jgradym/Google Drive/Phylo_trees/Upham_Phylos/Mammals/Completed_5911sp_topoCons_FBDasZhouEtAl', full.names = T)[1] #100 trees, access from data.vertlife.org
mam_tree <- read.tree(mamm_tree_files)
mort_mammals_4$Species <- gsub(" ", "_", mort_mammals_4$Species)
mort_mammals_4$Species[mort_mammals_4$Species %nin% mam_tree$tip.label]
#"Heteromys_adspersus" "Otaria_byronia" 

# change names to match phylogeny
mort_mammals_4$Species[mort_mammals_4$Species == "Heteromys_adspersus"] <- "Liomys_adspersus"
mort_mammals_4$Species[mort_mammals_4$Species == "Otaria_byronia"] <- "Otaria_bryonia"
mort_mammals_fin <- mort_mammals_4 %>%
  filter(Species %in% mam_tree$tip.label) %>%
  dplyr::select(Group, Species, dry_mass_g, Mean_Temp_C, Mortality_yr) %>%
  rename(Temp_C = Mean_Temp_C)

mort_mammals_fin$matched_name2 <- NULL
mort_mammals_fin$new_name <- NULL
mort_mammals_fin
#combine mammals and birds and update names in McCoy and Gillooly

#combine endo data
mortality_endo <- as_tibble(rbind(mort_mammals_fin , mort_bird_fin)) 

mortality_ecto <- mortality_mccoy %>%
  filter(Group == "Invertebrate" | Group == "Fish") %>%
  dplyr::select(-Ref)
mortality_updated <- as_tibble(rbind(mortality_endo, mortality_ecto))


write_csv(mortality_updated,'/Users/jgradym/Documents/GitHub/Foodweb_thermal_asymmetries/Data/McCoy_mortality_updated.csv')
