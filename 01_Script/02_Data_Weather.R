################################################################
################################################################
#####                        R-Script                      #####
#####      Sebastian Mader & Tobias Rüttenauer (2022):     #####
#####          The effects of non-pharmaceutical           #####
#####         interventions on COVID-19 mortality          ##### 
#####                      28.10.2021				               #####
################################################################
################################################################

rm(list = ls())

library(WDI)

library(plm)
library(splm)
library(feisr)
library(texreg)
library(extrafont)
loadfonts()
library(ggplot2)
library(cowplot)

library(ecmwfr)
library(ncdf4)
library(mapview)

setwd("Path/02_Data")



#########################
#### Load Covid data ####
#########################

load("Corona_all_weekly.RData")





#############################
#### Copernicus API data ####
#############################

user <- "XX" # specify user
key <- "YY" # Specify key

### Set up API access
wf_set_key(user = user, 
           key = key, 
           service = "cds")


### Request ERA5 monthly data of single levels
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=overview
request <- list(
  format = "netcdf",
  variable = c("2m_temperature", "total_cloud_cover", "total_precipitation"),
  product_type = "monthly_averaged_reanalysis",
  year = c("2017", "2018", "2019", "2020", "2021"),
  month = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"),
  time = "00:00",
  dataset_short_name = "reanalysis-era5-single-levels-monthly-means",
  target = "era5.nc"
)

ncfile <- wf_request(user = user,
                     request = request,   
                     transfer = TRUE,  
                     path = ".",
                     verbose = FALSE)

# Read the nc data
era5.nc <- nc_open(ncfile)

lon <- ncvar_get(era5.nc, "longitude")
lat <- ncvar_get(era5.nc, "latitude", verbose = F)
t <- ncvar_get(era5.nc, "time")
fillvalue <- era5.nc$var[["t2m"]]$missval

data.array1 <- ncvar_get(era5.nc, c("t2m"))
data1 <- data.frame(matrix(data.array1, ncol = dim(data.array1)[4], byrow = FALSE))
#data1 <- data1[1:(nrow(data1)/2), ]
names(data1) <- paste0("t2m_", t)
rm(data.array1)

data.array2 <- ncvar_get(era5.nc, c("tcc"))
data2 <- data.frame(matrix(data.array2, ncol = dim(data.array2)[4], byrow = FALSE))
#data2 <- data2[1:(nrow(data2)/2), ]
names(data2) <- paste0("tcc_", t)
rm(data.array2)

data.array3 <- ncvar_get(era5.nc, c("tp"))
data3 <- data.frame(matrix(data.array3, ncol = dim(data.array3)[4], byrow = FALSE))
#data3 <- data3[1:(nrow(data3)/2), ]
names(data3) <- paste0("tp_", t)
rm(data.array3)
Sys.sleep(3); gc(); Sys.sleep(3)

nc_close(era5.nc)


### Request ERA5 monthly data of pressure levels
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=overview
request <- list(
  format = "netcdf",
  product_type = "monthly_averaged_reanalysis",
  variable = "specific_humidity",
  year = c("2017", "2018", "2019", "2020", "2021"),
  month = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"),
  time = "00:00",
  pressure_level = "1000",
  dataset_short_name = "reanalysis-era5-pressure-levels-monthly-means",
  target = "era5_2.nc"
)

ncfile <- wf_request(user = user,
                     request = request,   
                     transfer = TRUE,  
                     path = ".",
                     verbose = FALSE)


# Read the nc data
era5.nc <- nc_open(ncfile)

fillvalue2 <- era5.nc$var[["q"]]$missval

data.array4 <- ncvar_get(era5.nc, c("q"))
data4 <- data.frame(matrix(data.array4, ncol = dim(data.array4)[4], byrow = FALSE))
# data4 <- data4[1:(nrow(data4)/2), ]
names(data4) <- paste0("q_", t)
rm(data.array4)

nc_close(era5.nc)




### Combined coordinates
londf <- rep(lon, length(lat))
latdf <- rep(lat, each = length(lon))

coords <- data.frame(lon = londf, lat = latdf)



###############################
#### Compute country means ####
###############################

### Reproject the longitudes to WGS84 (latitudes are already in WGS84)

# Easy (potentially imprecise): just bring lon on -180 to 180 scale
coords_new <- coords
coords_new$lon <- ifelse(coords_new$lon > 180, coords_new$lon - 360, coords_new$lon)

grid <- cbind(coords_new, id = c(1:nrow(coords)))
grid <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)


### Country shapes
data(wrld_simpl, package = "maptools")
wrld_simpl <- st_as_sf(wrld_simpl)
wrld_simpl <- st_transform(wrld_simpl, crs = st_crs(grid))
wrld_simpl <- st_make_valid(wrld_simpl)

### Country codes

# Use Sudan for south sudan
rep <- wrld_simpl[wrld_simpl$ISO3 == "SDN", ]
rep$ISO3 <- "SSD"
bb <- c(xmin = 23.8869795809, xmax = 35.2980071182, ymax = 12.2480077571, ymin = 3.50917)
bb <- st_bbox(bb, crs = 4326)
rep <- st_crop(rep, bb)

wrld_simpl <- rbind(wrld_simpl, rep)

# Kosovo
rep <- wrld_simpl[wrld_simpl$ISO3 == "SRB", ]
rep$ISO3 <- "RKS"
bb <- c(xmin = 20.002913, xmax = 21.790395, ymax = 43.268039, ymin = 41.857474)
bb <- st_bbox(bb, crs = 4326)
rep <- st_crop(rep, bb)

wrld_simpl <- rbind(wrld_simpl, rep)


### For small countries, add 0.2 degree buffer
sc <- c("ATG", 
"BRB", 
"BMU", 
"CYM", 
"COK", 
"DMA", 
"FSM", 
"GRD", 
"LIE", 
"MSR", 
"MLT", 
"MDV", 
"NIU", 
"HKG", 
"MNP", 
"AND", 
"GIB", 
"MAC", 
"MCO", 
"MYT", 
"NFK", 
"CCK", 
"BVT", 
"IOT", 
"CXR", 
"UMI", 
"NRU", 
"REU", 
"KNA", 
"SGP", 
"TKL", 
"TON", 
"TUV", 
"VGB", 
"WLF", 
"ANT", 
"PCN", 
"SPM", 
"SHN", 
"SMR", 
"TCA", 
"VAT", 
"MAF", 
"BLM", 
"GGY", 
"JEY")

wrld_small <- wrld_simpl[which(wrld_simpl$ISO3 %in% sc), ]
wrld_small <- st_buffer(wrld_small, 0.2)

# Append
wrld_simpl <- rbind(wrld_simpl[which(!wrld_simpl$ISO3 %in% sc),], wrld_small)


### Reduce data to land cover
land <- st_intersects(grid, wrld_simpl)
oo <- which(lapply(land, length) > 0)

grid <- grid[oo, ]
data1 <- data1[oo, ]
data2 <- data2[oo, ]
data3 <- data3[oo, ]
data4 <- data4[oo, ]


### Intersection between climate grid and countries
ints.spdf <- st_intersection(grid, wrld_simpl)
save(ints.spdf, file = "Intersection_grid_countries.RData")

# load("Intersection_grid_countries.RData")

# Check any country not merged 
(ni <- setdiff(wrld_simpl$ISO3, ints.spdf$ISO3))



### Merge variables
data1$id <- rownames(data1)
data2$id <- rownames(data2)
data3$id <- rownames(data3)
data4$id <- rownames(data4)

weather.spdf <- merge(ints.spdf, data1, by = "id")
weather.spdf <- merge(weather.spdf, data2, by = "id")
weather.spdf <- merge(weather.spdf, data3, by = "id")
weather.spdf <- merge(weather.spdf, data4, by = "id")


### Aggregate data
weather.df <- st_drop_geometry(weather.spdf)
weather.df <- aggregate(weather.df[, which(grepl("_", names(weather.df)))], 
                        by = list(iso2 = weather.df$ISO2, 
                                  iso3 = weather.df$ISO3,
                                  name = weather.df$NAME,
                                  region = weather.df$REGION),
                   FUN = function(x) mean(x, na.rm = TRUE))


# Save wide
save(weather.df, file = "Weather_wide.RData")



#################
#### Reshape ####
#################

# Reshape
weather_long.df <- reshape(weather.df, varying = names(weather.df)[5:ncol(weather.df)],
                           sep = "_", direction = "long")

# time as date (minutes since 1990)
weather_long.df$time <- as.Date(weather_long.df$time/24, origin = "1900-01-01")
weather_long.df$year <- as.numeric(format(weather_long.df$time, "%Y"))
weather_long.df$month <- as.numeric(format(weather_long.df$time, "%m"))

# Make temperature to Celsius (is Kelvin)
weather_long.df$t2m <- weather_long.df$t2m - 273.15
summary(weather_long.df$t2m)

# Create monthly averages averages (to use for recent 3 months) 
vars <- c("t2m", "tcc", "tp", "q")
for(v in vars){
  n <- paste0("avg_", v)
  weather_long.df[, n] <- ave(weather_long.df[, v],
                              weather_long.df$iso3, weather_long.df$month,
                              FUN = function(x) mean(x, na.rm = TRUE))
}


save(weather_long.df, file = "Weather_long.RData")






################################
#### Merge with weekly data ####
################################

load("Corona_all_weekly.RData")



weather_long.df$iso_code <- weather_long.df$iso3

### Add monthly data
corona_all.df$year <- as.numeric(format(as.Date(corona_all.df$date, format = "%Y%m%d"), "%Y"))
corona_all.df$month <- as.numeric(format(as.Date(corona_all.df$date, format = "%Y%m%d"), "%m"))

corona_all.df <- merge(corona_all.df, 
                       weather_long.df[, c("iso_code", "month", "year", "t2m", "tcc", "tp", "q")],
                       by = c("iso_code", "year", "month"), all.x = TRUE)

### Add the average
corona_all.df <- merge(corona_all.df, 
                       unique(weather_long.df[, c("iso_code", "month","avg_t2m", "avg_tcc", "avg_tp", "avg_q")]),
                       by = c("iso_code", "month"), all.x = TRUE)

# Impute the averages of last three years if missing
vars <- c("t2m", "tcc", "tp", "q")
for(v in vars){
  n <- paste0("avg_", v)
  oo <- which(is.na(corona_all.df[, v]))
  corona_all.df[oo, v] <- corona_all.df[oo, n]
}


save(corona_all.df, file = "Corona_all_weekly_2.RData")







###############################
#### Merge with daily data ####
###############################

load("Corona_all_daily.RData")



weather_long.df$iso_code <- weather_long.df$iso3

### Add monthly data
corona_all.df$year <- as.numeric(format(as.Date(corona_all.df$date, format = "%Y%m%d"), "%Y"))
corona_all.df$month <- as.numeric(format(as.Date(corona_all.df$date, format = "%Y%m%d"), "%m"))

corona_all.df <- merge(corona_all.df, 
                       weather_long.df[, c("iso_code", "month", "year", "t2m", "tcc", "tp", "q")],
                       by = c("iso_code", "year", "month"), all.x = TRUE)

### Add the average
corona_all.df <- merge(corona_all.df, 
                       unique(weather_long.df[, c("iso_code", "month","avg_t2m", "avg_tcc", "avg_tp", "avg_q")]),
                       by = c("iso_code", "month"), all.x = TRUE)

# Impute the averages of last three years if missing
vars <- c("t2m", "tcc", "tp", "q")
for(v in vars){
  n <- paste0("avg_", v)
  oo <- which(is.na(corona_all.df[, v]))
  corona_all.df[oo, v] <- corona_all.df[oo, n]
}


save(corona_all.df, file = "Corona_all_daily_2.RData")




