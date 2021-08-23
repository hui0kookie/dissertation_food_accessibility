library(sf)
library(tidyverse)
library(tmap)
library(tidycensus)
library(tigris)
library(rmapshaper)
library(matrixStats)
library(SpatialAcc)

#Reading population file
pop <- read.csv("casa2020/i2p/dissertation/dissertation202106/001/ward_pop.csv")

#Reading shapefile
ward <- st_read("casa2020/i2p/dissertation/dissertation202106/001/shp01/London_Ward_CityMerged.shp") %>%
  st_transform(., 27700)
qtm(ward)


#Merge data
joined_df <- merge(ward, pop,
                   by.x='GSS_CODE',
                   by.y='New_code') %>%
  distinct(.keep_all=TRUE)

stores <- read_csv("casa2020/i2p/dissertation/dissertation202106/001/geo201509.csv") %>%
  st_as_sf(., coords = c("LongWGS84", "LatWGS84"), 
           crs = 4326) %>%
  st_transform(., 27700)
qtm(stores)

#clip points data
stores_subset <- stores[ward, ]
qtm(stores_subset)

stores_subset <- stores_subset %>%
  select(GLUID,StoreName,size_band, geometry)


green <- st_read('R_dissertation/from_qgis/green_grocery_Lon.geojson') %>%
  st_transform(., 27700)

stores_subset <- rbind(stores_subset, green) %>%
  distinct()

#plot
tm_shape(joined_df) +
  tm_polygons() +
  tm_shape(stores_subset) +
  tm_dots(col = "blue") +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(main.title = "Grocery store locations in London",
            main.title.size = 1.2, frame = FALSE,
            legend.outside = TRUE, legend.outside.position = "right")

#count store numbers in each ward
lon_food <- stores_subset  %>% 
  st_join(joined_df) %>%
  group_by(GSS_CODE) %>% 
  summarize(stores_counts = n())

lon_food <- st_drop_geometry(lon_food)

joined_df <- joined_df %>%
  left_join(lon_food, by = "GSS_CODE") %>%
  mutate(stores_counts = replace_na(stores_counts, 0))

joined_df <- joined_df %>%
  mutate(storesperpop = (stores_counts/Population_2015)*1000)

summary(joined_df$storesperpop)

#choropleth map of stores per 1000 population
tm_shape(joined_df) +
  tm_polygons(col = "storesperpop", style = "jenks",palette = "YlGnBu", 
              border.alpha = 0, title = "Stores per\n1k population") +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(main.title = "Grocery stores per 1000 people in London",
            main.title.size = 1.2, frame = FALSE,
            legend.outside = TRUE, legend.outside.position = "right")

# population density

tm_shape(joined_df) +
  tm_polygons(col = "Population_density", style = "jenks",palette = "YlGnBu", 
              border.alpha = 0, title = "population density") +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(main.title = "Population density in London",
            main.title.size = 1.2, frame = FALSE,
            legend.outside = TRUE, legend.outside.position = "right")



#Point proximity buffers
tract.centroids <- st_centroid(joined_df$geometry)

tm_shape(joined_df) +
  tm_polygons(col = "blue") +
  tm_shape(tract.centroids) +
  tm_dots(col = "red")

tract.centroids = tract.centroids %>%
  st_sf

tract.centroids <- tract.centroids %>%
  st_join(joined_df)

tract.buff  <-st_buffer(tract.centroids, dist = 1609)
class(tract.buff)
qtm(tract.buff)

tract.buff = tract.buff %>%
  st_sf %>%
  st_cast

buff.food <- stores_subset  %>% 
  st_join(tract.buff) 
  

buff.food <- buff.food %>% 
  group_by(GSS_CODE) %>%
  summarize(stores1m = n()) 

buff.food <- st_drop_geometry(buff.food)

joined_df <- joined_df %>%
  left_join(buff.food, by = "GSS_CODE") %>%
  mutate(stores1m = replace_na(stores1m, 0)) %>%
  mutate(storesbuff1m = (stores1m/joined_df$Population_2015)*1000)


tmap_mode("plot")

#Distance to nearest bank

#Euclidean distance
food.dist<-st_distance(tract.centroids, stores_subset)

joined_df <- joined_df %>%
  mutate(storemin = rowMins(food.dist))

tm_shape(joined_df) +
  tm_polygons(col = "storemin", style = "jenks",palette = "inferno", 
              border.alpha = 0, title = "Distance to nearest \nstore (m)") +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(main.title = "Distance to nearest store",
            main.title.size = 1.25, frame = FALSE,
            legend.outside = TRUE, legend.outside.position = "right")

# FCA

buff.tracts <- tract.centroids  %>% 
  select(Population_2015) %>% 
  st_join(tract.buff) %>%
  group_by(GSS_CODE) %>% 
  summarize(buffpop = sum(Population_2015.x)) %>%
  ungroup()

buff.tracts <- st_drop_geometry(buff.tracts)

joined_df <- joined_df %>%
  left_join(buff.tracts, by = "GSS_CODE") %>%
  mutate(buffpop = replace_na(buffpop, 0)) %>%
  mutate(fca = (stores1m/buffpop)*10000)

tm_shape(joined_df, unit = "mi") +
  tm_polygons(col = "fca", style = "jenks",palette = "YlGn", 
              border.alpha = 0, title = "FCA") +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(main.title = "Store spatial accessibility in London",
            main.title.size = 0.95, frame = FALSE,
            legend.outside = TRUE, legend.outside.position = "right")

# 2SFCA

class(food.dist)

centroid.coords <- st_coordinates(tract.centroids)
store.coords <- st_coordinates(stores_subset)

dist.matrix <- distance(centroid.coords, store.coords, type = "euclidean")
class(dist.matrix)

TSFCA <- ac(p = joined_df$Population_2015, 
            n = stores_subset$size_band, 
            D = dist.matrix, d0 = 1609, family = "2SFCA")

joined_df <- joined_df %>%
  mutate(TSFCA = TSFCA)

normalize <- function(x) {
  return (((x - min(x)) / (max(x) - min(x))) * 10)
}

joined_df$TSFCA <- normalize(joined_df$TSFCA)

tm_shape(joined_df) +
  tm_polygons(col = "TSFCA", style = "jenks",palette = "plasma", 
              border.alpha = 0, title = "2SFCA Score") +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(main.title = "healthy food spatial accessibility in London",
            main.title.size = 1.2, frame = FALSE,
            legend.outside = TRUE, legend.outside.position = "right")

#bar plots
library(ggplot2)

joined_df <- joined_df %>%
  mutate(tsfca_level = case_when(TSFCA > 5.54 ~ "Very High Accessibility",
                                 TSFCA > 3.72 ~ "High Accessibility",
                                 TSFCA > 2.57 ~ "Moderate Accessibility",
                                 TSFCA > 1.48 ~ "Low  Accessibility",
                                 TRUE ~ "Very Low  Accessibility"))

ggplot(joined_df, aes(fill=tsfca_level, y=TSFCA, x=BOROUGH)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  scale_fill_brewer(palette = 'Accent') +
  labs(title = "Accessibility Levels for Each of the Boroughs", 
       y = "2SFCA ", x = "Borough", fill = "Accessibility Level")


# summaries
joined_dfz <- joined_df %>%
  group_by(tsfca_level)%>%
  summarise(count=n(), Average=mean(stores_counts), ratio=mean(storesperpop))






























