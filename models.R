library(tidyverse)
library(tmap)
library(geojsonio)
library(plotly)
library(rgdal)
library(broom)
library(mapview)
library(crosstalk)
library(sf)
library(sp)
library(spdep)
library(car)
library(fs)
library(janitor)

#SES data
ses <- read.csv('R_dissertation/London_Ward_Atlas.csv')

#Ward with population and 2sfca score
ward <- st_read('R_dissertation/tsfca0.geojson') %>%
  st_transform(crs=27700)
qtm(ward)

#normalize data
ward$TSFCA[is.na(ward$TSFCA)] <- median(ward$TSFCA, na.rm = TRUE)

#function min-max scaler (1,10)
normalize <- function(x) {
  return (((x - min(x)) / (max(x) - min(x))) * 10)
}

ward$TSFCA <- normalize(ward$TSFCA)

ses <- ses %>%
  mutate(white_rate = White /
           Population_2011)


#Store plot

#descriptive statistics
histplot <- ggplot(data=ward, aes(x=TSFCA))
histplot +geom_histogram(fill = "white", colour = "black")

histplot <- ggplot(data=ses, aes(x=Unemployed_rate))
histplot +geom_histogram(fill = "white", colour = "black")

histplot <- ggplot(data=ses, aes(x=Median_House_Price_2014))
histplot +geom_histogram(fill = "white", colour = "black")

histplot <- ggplot(data=ses, aes(x=No_cars_household_rate))
histplot +geom_histogram(fill = "white", colour = "black")

histplot <- ggplot(data=ses, aes(x=ALevel_2014))
histplot +geom_histogram(fill = "white", colour = "black")

histplot <- ggplot(data=ses, aes(x=Housing_Benefit_rates_2015))
histplot +geom_histogram(fill = "white", colour = "black")

histplot <- ggplot(data=ses, aes(x=white_rate))
histplot +geom_histogram(fill = "white", colour = "black")

#join two data sets
LonWardProfiles <- ward %>%
  left_join(., 
            ses, 
            by = c("GSS_CODE" = "New_ward_code"))
  
#now OLS model
model1 <- LonWardProfiles %>%
  lm((TSFCA)^0.5 ~
       Deprivation_Score,
     data=.)

summary(model1)

# normal distribution: Tukey¡¯s ladder of transformations

symbox(~TSFCA, 
       LonWardProfiles, 
       na.rm=T,
       powers=seq(-3,3,by=.5))

symbox(~Median_House_Price_2014, 
       LonWardProfiles, 
       na.rm=T,
       powers=seq(-3,3,by=.5))

model2 <- lm((TSFCA)^0.5 ~ 
                 ALevel_2014 + white_rate  + log(Median_House_Price_2014) +
               No_cars_household_rate +
                 Unemployed_rate + Housing_Benefit_rates_2015, 
             data = LonWardProfiles)

summary(model2)

#the residuals
model_data <- model2 %>%
  augment(., LonWardProfiles)

LonWardProfiles <- LonWardProfiles %>%
  mutate(model2resids = residuals(model2))

#plot residuals
model_data%>%
  dplyr::select(.resid)%>%
  pull()%>%
  qplot()+ 
  geom_histogram() 

#model diagnositcs
par(mfrow=c(2,2))
plot(model2)

#Standard Autocorrelation:durbin-watson test
DW <- durbinWatsonTest(model2)
tidy(DW)

#check Multicolinearity
library(corrr)

Correlation <- LonWardProfiles %>%
  st_drop_geometry()%>%
  dplyr::select(No_cars_household_rate,
                ALevel_2014,
                white_rate,
                Median_House_Price_2014,
                Housing_Benefit_rates_2015,
                Unemployed_rate) %>%
  correlate()
 
#visualise the correlation matrix
rplot(Correlation)

vif(model2)


coordsW <- LonWardProfiles%>%
  st_centroid()%>%
  st_geometry()


plot(coordsW)

LWard_nb <- LonWardProfiles %>%
  poly2nb(., queen=T)

#or nearest neighbours
knn_wards <-coordsW %>%
  knearneigh(., k=4)

LWard_knn <- knn_wards %>%
  knn2nb()

#plot them
plot(LWard_nb, st_geometry(coordsW), col="red")

plot(LWard_knn, st_geometry(coordsW), col="blue")

#create a spatial weights matrix object from these weights

Lward.queens_weight <- LWard_nb %>%
  nb2listw(., style="C")

Lward.knn_4_weight <- LWard_knn %>%
  nb2listw(., style="C")

Queen <- LonWardProfiles %>%
  st_drop_geometry()%>%
  dplyr::select(model2resids)%>%
  pull()%>%
  moran.test(., Lward.queens_weight)%>%
  tidy()

Nearest_neighbour <- LonWardProfiles %>%
  st_drop_geometry()%>%
  dplyr::select(model2resids)%>%
  pull()%>%
  moran.test(., Lward.knn_4_weight)%>%
  tidy()

Queen

Nearest_neighbour

#Spatial Regression Models

#Spatial Lag model
library(spatialreg)

slag_dv_model2_queen <- lagsarlm((TSFCA)^0.5 ~ 
                                   ALevel_2014 + white_rate  + log(Median_House_Price_2014) +
                                   No_cars_household_rate +
                                   Unemployed_rate + Housing_Benefit_rates_2015,
                                 data = LonWardProfiles, 
                                 nb2listw(LWard_nb, style="C"), 
                                 method = "eigen")
tidy(slag_dv_model2_queen)
glance(slag_dv_model2_queen)


#run a spatially-lagged regression model
slag_dv_model2_knn4 <- lagsarlm((TSFCA)^0.5 ~ 
                                  ALevel_2014 + white_rate  + log(Median_House_Price_2014) +
                                  No_cars_household_rate +
                                  Unemployed_rate + Housing_Benefit_rates_2015, 
                                data = LonWardProfiles, 
                                nb2listw(LWard_knn, 
                                         style="C"), 
                                method = "eigen")
glance(slag_dv_model2_knn4)

LonWardProfiles <- LonWardProfiles %>%
  mutate(slag_dv_model2_knn_resids = residuals(slag_dv_model2_knn4))


KNN4Moran <- LonWardProfiles %>%
  st_drop_geometry()%>%
  dplyr::select(slag_dv_model2_knn_resids)%>%
  pull()%>%
  moran.test(., Lward.knn_4_weight)%>%
  tidy()

KNN4Moran

#The Spatial Error Model

sem_model1 <- errorsarlm((TSFCA)^0.5 ~ 
                           ALevel_2014 + white_rate  + log(Median_House_Price_2014) +
                           No_cars_household_rate +
                           Unemployed_rate + Housing_Benefit_rates_2015, 
                         data = LonWardProfiles,
                         nb2listw(LWard_knn, style="C"), 
                         method = "eigen")

tidy(sem_model1)

#GWR model

myvars <- LonWardProfiles %>%
  dplyr::select(TSFCA,
                No_cars_household_rate,
                ALevel_2014,
                white_rate,
                Median_House_Price_2014,
                Housing_Benefit_rates_2015,
                Unemployed_rate)

model_final <- lm((TSFCA)^0.5 ~ 
                    ALevel_2014 + white_rate  + log(Median_House_Price_2014) +
                    No_cars_household_rate +
                    Unemployed_rate + Housing_Benefit_rates_2015, 
                  data = myvars)
tidy(model_final)
glance(model_final)

LonWardProfiles <- LonWardProfiles %>%
  mutate(model_final_res = residuals(model_final))

qtm(LonWardProfiles, fill = "model_final_res")

final_model_Moran <- LonWardProfiles %>%
  st_drop_geometry()%>%
  dplyr::select(model_final_res)%>%
  pull()%>%
  moran.test(., Lward.knn_4_weight)%>%
  tidy()

final_model_Moran

library(spgwr)

st_crs(LonWardProfiles) = 27700

LonWardProfilesSP <- LonWardProfiles %>%
  as(., "Spatial")

st_crs(coordsW) = 27700

GWRbandwidth <- gwr.sel((TSFCA)^0.5 ~ 
                          ALevel_2014 + white_rate  + log(Median_House_Price_2014) +
                          No_cars_household_rate +
                          Unemployed_rate + Housing_Benefit_rates_2015, 
                        data = LonWardProfilesSP, 
                        coords=coordsWSP,
                        adapt=T)

gwr.model = gwr((TSFCA)^0.5 ~ 
                  ALevel_2014 + white_rate  + log(Median_House_Price_2014) +
                  No_cars_household_rate +
                  Unemployed_rate + Housing_Benefit_rates_2015, 
                data = LonWardProfilesSP, 
                coords=coordsWSP, 
                adapt=GWRbandwidth, 
                hatmatrix=TRUE, 
                se.fit=TRUE)

gwr.model

results <- as.data.frame(gwr.model$SDF)
names(results)

tmap_mode("plot")
tm_shape(LonWardProfiles) +
  tm_polygons(col = "model_final_res",palette = "RdYlBu", 
              border.alpha = 0, title = "residuals") +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(main.title = "OLS model residuals",
            main.title.size = 0.95, frame = FALSE,
            legend.outside = TRUE, legend.outside.position = "right")


results <- as.data.frame(gwr.model$SDF)
names(results)

LonWardProfiles2 <- LonWardProfiles %>%
  mutate(coefAlevel = results$ALevel_2014,
         coefHousePrice = results$log.Median_House_Price_2014.,
         coefWrite = results$white_rate,
         coefNocar = results$No_cars_household_rate,
         coefUnemployed = results$Unemployed_rate,
         coefBenefit = results$Housing_Benefit_rates_2015,
         localR2 = results$localR2)

tmap_mode("plot")
tm1 <- tm_shape(LonWardProfiles2) + 
  tm_polygons(col = "coefAlevel",
              border.alpha = 0,
              style="jenks", 
              palette='RdYlBu') + 
  tm_legend(show=TRUE,scale=0.6)+
  tm_layout(frame=FALSE) +
  tm_credits("(a)", position=c(0.1,0.1), size=1.2)

tm2 <- tm_shape(LonWardProfiles2)+ 
  tm_polygons(col = "coefHousePrice",
              border.alpha = 0,
              style="jenks",
              palette="RdYlBu")+
  tm_legend(show=TRUE,scale=0.6) +
  tm_compass(north=0, position=c(0.8,0.65)) +
  tm_layout(frame=FALSE)+
  tm_credits("(b)", position=c(0.1,0.1), size=1.2)

tm3 <- tm_shape(LonWardProfiles2)+ 
  tm_polygons(col = "coefWrite",
              border.alpha = 0,
              style="jenks",
              midpoint = NA,
              palette="RdYlBu")+
  tm_legend(show=TRUE,scale=0.6)+
  tm_layout(frame=FALSE)+
  tm_credits("(c)", position=c(0.1,0.1), size=1.2)


tm4 <- tm_shape(LonWardProfiles2) +
  tm_polygons(col = 'coefNocar',
              border.alpha = 0,
              style="jenks",
              palette="RdYlBu") +
  tm_legend(show=TRUE,scale=0.6)+
  tm_layout(frame=FALSE)+
  tm_credits("(d)", position=c(0.1,0.1), size = 1.2)

tm5 <- tm_shape(LonWardProfiles2) +
  tm_polygons(col = 'coefUnemployed',
              border.alpha = 0,
              style="jenks",
              palette="RdYlBu") +
  tm_legend(show=TRUE,scale=0.6)+
  tm_scale_bar(position=c(0.10,0.02)) +
  tm_layout(frame=FALSE)+
  tm_credits("(e)", position=c(0.1,0.1), size = 1.2)

tm6 <- tm_shape(LonWardProfiles2) +
  tm_polygons(col = 'coefBenefit',
              border.alpha = 0,
              style="jenks",
              palette="RdYlBu") +
  tm_legend(show=TRUE,scale=0.6)+
  tm_layout(frame=FALSE)+
  tm_credits("(f)", position=c(0.1,0.1), size = 1.2)

t0=tmap_arrange(tm1, tm2, tm3, tm4, tm5,tm6, ncol=2)

t0

#local r2
tm_shape(LonWardProfiles2) +
  tm_polygons(col = "localR2",palette = "RdYlGn", 
              border.alpha = 0.05, title = "R Square") +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(main.title = "GWR model local R square",
            main.title.size = 0.95, frame = FALSE,
            legend.outside = TRUE, legend.outside.position = "right")































ses$white_rate <- normalize(ses$white_rate)

ses$Public_Transport_Accessibility <- normalize(ses$Public_Transport_Accessibility)

ses$Median_House_Price_2014 <- normalize(ses$Median_House_Price_2014)

ses$Unemployed_rate <- normalize(ses$Unemployed_rate)

ses$Housing_Benefit_rates_2015 <- normalize(ses$Housing_Benefit_rates_2015)

ses$children_Claimant_Households <- normalize(ses$children_Claimant_Households)

ses$GCSE_2014 <- normalize(ses$GCSE_2014)

ses$ALevel_2014 <- normalize(ses$ALevel_2014)

ses$No_cars_household_rate <- normalize(ses$No_cars_household_rate)

ses$Deprivation_Score <- normalize(ses$Deprivation_Score)

