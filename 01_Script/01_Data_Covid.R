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

setwd("Path/02_Data")



##############################
#### Get World Covid data ####
##############################

### Use Our World in Data
#https://github.com/owid/covid-19-data/tree/master/public/data

covid.df <- read.csv("https://covid.ourworldindata.org/data/owid-covid-data.csv",
                     header = TRUE, stringsAsFactors = FALSE)

names(covid.df)[which(names(covid.df) == "location")] <- "countryname"

# Drop international
covid.df <- covid.df[which(covid.df$countryname != "International"), ]
covid.df <- covid.df[which(covid.df$countryname != "World"), ]


# recode date
covid.df$date <- gsub("-", "", covid.df$date)


########################################
#### Get World Bank data and merge  ####
########################################

wd.df <- WDI(country = "all", 
             indicator = c('stat_cap_all' = "iq.sci.ovrl", 
                           'stat_cap_meth' = "iq.sci.mthd", 
                           'stat_cap_source' = "iq.sci.srce", 
                           'stat_cap_time' = "iq.sci.prdc"),
             start = 2019, end = 2019)
wd.df[] <- lapply(wd.df, function(x) { attributes(x) <- NULL; x })
names(wd.df)[which(names(wd.df) == "country")] <- "countryname"
wd.df$year <- NULL

### Test country names
cn <- unique(covid.df$countryname)
oo <- which(!cn %in% wd.df$countryname)
cn[oo]

# Get unmatched in wd
cn2 <- unique(wd.df$countryname)
oo2 <- which(!cn2 %in% cn)
cn2[oo2][order(cn2[oo2])]

# Do approx merge
apmatch <- lapply(cn[oo], function(x) agrep(x, cn2[oo2], ignore.case = TRUE, 
                                          max.distance = 0.4, value = TRUE))
apmatch <- lapply(seq_along(cn[oo]), function(x) c(cn[oo][x], apmatch[[x]]))
apmatch
  
# Rename countryname
# ren <- c("Bahamas"                      =   "Bahamas, The",
#          "Brunei"                       =        "Brunei Darussalam",
#          "Cape Verde"                   = "Cabo Verde",
#          "Congo"                        =       "Congo, Rep.",
#          "Czechia"                      =     "Czech Republic",
#          "Democratic Republic of Congo" = "Congo, Dem. Rep.",
#          "Egypt"                        =     "Egypt, Arab Rep.",
#          "Gambia"                       ="Gambia, The",
#          "Hong Kong"                    =         "Hong Kong SAR, China",
#          "Iran"                         = "Iran, Islamic Rep.",
#          "Kyrgyzstan"                   = "Kyrgyz Republic",
#          "Laos"                         =                "Lao PDR",
#          "Palestine"                    = "Palestine", ## no match
#          "Russia"                       = "Russian Federation",
#          "Saint Kitts and Nevis"        = "St. Kitts and Nevis",
#          "Saint Vincent and the Grenadines" = "St. Vincent and the Grenadines",
#          "Slovakia"                    = "Slovak Republic",
#          "South Korea"                 = "Korea, Rep.",
#          "Syria"                       = "Syrian Arab Republic",
#          "Timor"                       = "Timor-Leste",
#          "Vatican"                     = "Vatican",
#          "Venezuela"                   = "Venezuela, RB",
#          "Yemen"                       = "Yemen, Rep.")

ren <- c("Africa"                           = "Africa", 
         "Anguilla"                         = "Anguilla",
         "Asia"                             = "Asia",
         "Bahamas"                          = "Bahamas, The",
         "Bonaire Sint Eustatius and Saba"  = "Bonaire Sint Eustatius and Saba",
         "Brunei"                           = "Brunei Darussalam",
         "Cape Verde"                       = "Cabo Verde",
         "Congo"                            = "Congo, Rep.",
         "Cook Islands"                     = "Cook Islands",
         "Czechia"                          = "Czech Republic",
         "Democratic Republic of Congo"     = "Congo, Dem. Rep.",
         "Egypt"                            = "Egypt, Arab Rep.",
         "Europe"                           = "Europe", 
         "Faeroe Islands"                   = "Faroe Islands", 
         "Falkland Islands"                 = "Falkland Islands",
         "Gambia"                           = "Gambia, The",
         "Guernsey"                         = "Guernsey",
         "Hong Kong"                        = "Hong Kong SAR, China",
         "Iran"                             = "Iran, Islamic Rep.",
         "Jersey"                           = "Jersey",
         "Kyrgyzstan"                       = "Kyrgyz Republic",
         "Laos"                             =  "Lao PDR",
         "Macao"                            = "Macao SAR, China", 
         "Micronesia (country)"             = "Micronesia (country)", 
         "Montserrat"                       = "Montserrat", 
         "Northern Cyprus"                  = "Northern Cyprus",
         "Oceania"                          = "Oceania",
         "Palestine"                        = "Palestine", ## no match
         "Pitcairn"                         = "Pitcairn",
         "Russia"                           = "Russian Federation",
         "Saint Helena"                     = "Saint Helena",
         "Saint Kitts and Nevis"            = "St. Kitts and Nevis",
         "Saint Lucia"                      = "St. Lucia", 
         "Saint Vincent and the Grenadines" = "St. Vincent and the Grenadines",
         "Slovakia"                         = "Slovak Republic",
         "South America"                    = "South America",
         "South Korea"                      = "Korea, Rep.",
         "Syria"                            = "Syrian Arab Republic",
         "Taiwan"                           = "Taiwan", 
         "Timor"                            = "Timor-Leste",
         "Vatican"                          = "Vatican",
         "Venezuela"                        = "Venezuela, RB",
         "Wallis and Futuna"                = "Wallis and Futuna",
         "Yemen"                            =   "Yemen, Rep.")

for(i in ren){
  wd.df$countryname[which(wd.df$countryname == i)] <- names(ren)[which(ren == i)]
}


### Merge, keep only if present in both data sets
covid.df <- merge(covid.df, wd.df, by = "countryname",
                  all.x = FALSE, all.y = FALSE)





######################################
#### Merge covid response tracker ####
######################################

### Covid response tracker
# https://www.bsg.ox.ac.uk/research/research-projects/coronavirus-government-response-tracker#data 

response.df <- read.csv("https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv",
                        header = TRUE, stringsAsFactors = FALSE)

names(response.df) <- tolower(names(response.df))
names(response.df) <- gsub(".", "", names(response.df), fixed = TRUE)


# ### Add vaccination data full
# vaccine.df <- read.csv("https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_vaccines_full.csv",
#                         header = TRUE, stringsAsFactors = FALSE)
# 
# names(vaccine.df) <- tolower(names(vaccine.df))
# names(vaccine.df) <- gsub(".", "", names(vaccine.df), fixed = TRUE)
# 
# # Merge
# response.df <- merge(response.df, vaccine.df, by = c("countryname", "countrycode", "date"))


### Keep only national totals
table(response.df$jurisdiction)
response.df <- response.df[which(response.df$jurisdiction == "NAT_TOTAL"), ]


### Test codes
names(response.df)[which(names(response.df) == "countrycode")] <- "iso_code"

oo <- which(!covid.df$iso_code %in% response.df$iso_code)
table(covid.df$iso_code[oo])
oo2 <- which(!response.df$iso_code %in% covid.df$iso_code)
table(response.df$iso_code[oo2])

View(covid.df[oo,])
View(response.df[oo2,])

# Replace Kosovo
covid.df$iso_code[which(covid.df$iso_code == "OWID_KOS")] <- "RKS"


### Merge to covid data
response.df$countryname <- NULL; response.df$regionname <- NULL; response.df$regioncode <- NULL



corona_all.df <- merge(covid.df, response.df, by = c("iso_code", "date"),
                       all.x = TRUE, all.y = TRUE)


### Fill constant variables
repeat_last = function(x, forward = TRUE, maxgap = Inf, na.rm = FALSE) {
  if (!forward) x = rev(x)           # reverse x twice if carrying backward
  ind = which(!is.na(x))             # get positions of nonmissing values
  if (is.na(x[1]) && !na.rm)         # if it begins with NA
    ind = c(1,ind)                 # add first pos
  rep_times = diff(                  # diffing the indices + length yields how often
    c(ind, length(x) + 1) )          # they need to be repeated
  if (maxgap < Inf) {
    exceed = rep_times - 1 > maxgap  # exceeding maxgap
    if (any(exceed)) {               # any exceed?
      ind = sort(c(ind[exceed] + 1, ind))      # add NA in gaps
      rep_times = diff(c(ind, length(x) + 1) ) # diff again
    }
  }
  x = rep(x[ind], times = rep_times) # repeat the values at these indices
  if (!forward) x = rev(x)           # second reversion
  x
}

vars <- c("countryname", "continent", "iso2c",
          "population"                               ,"population_density",                      
          "median_age"                               ,"aged_65_older"                            ,"aged_70_older"        ,                   
          "gdp_per_capita"                           ,"extreme_poverty"                          ,"cardiovasc_death_rate",                   
          "diabetes_prevalence"                      ,"female_smokers"                           ,"male_smokers"         ,                   
          "handwashing_facilities"                   ,"hospital_beds_per_thousand"               ,"life_expectancy"      ,                   
          "human_development_index"                  ,"stat_cap_all"                            ,
          "stat_cap_meth"                            ,"stat_cap_source"                         ,"stat_cap_time")

for(i in vars){
  corona_all.df[, i] <- ave(corona_all.df[, i], 
                            corona_all.df$iso_code,
                            FUN = function(x) repeat_last(x, forward = FALSE))
  corona_all.df[, i] <- ave(corona_all.df[, i], 
                            corona_all.df$iso_code,
                            FUN = function(x) repeat_last(x, forward = TRUE))
}


### Drop obs which are all NA in either of them
corona_all.df$tmp <- ave(corona_all.df$new_deaths,
                         corona_all.df$iso_code,
                         FUN = function(x) mean(x, na.rm = TRUE))
table(corona_all.df$countryname[which(is.nan(corona_all.df$tmp))])
corona_all.df <- corona_all.df[which(!is.nan(corona_all.df$tmp)),]

corona_all.df$tmp <- ave(corona_all.df$stringencyindex,
                         corona_all.df$iso_code,
                         FUN = function(x) mean(x, na.rm = TRUE))
table(corona_all.df$countryname[which(is.nan(corona_all.df$tmp))])
corona_all.df <- corona_all.df[which(!is.nan(corona_all.df$tmp)),]

corona_all.df$tmp <- NULL



####################################
#### Rename or create variables ####
####################################

# Rename
names(corona_all.df)[which(names(corona_all.df) == "new_cases")] <- "cases"
names(corona_all.df)[which(names(corona_all.df) == "new_deaths")] <- "deaths"

names(corona_all.df)[which(names(corona_all.df) == "new_cases_per_million")] <- "cases_c"
names(corona_all.df)[which(names(corona_all.df) == "new_deaths_per_million")] <- "deaths_c"

names(corona_all.df)[which(names(corona_all.df) == "new_tests")] <- "tests"

names(corona_all.df)[which(names(corona_all.df) == "stringencyindex")] <- "stringency"



# Create
corona_all.df$cases_t <- corona_all.df$cases/corona_all.df$tests
corona_all.df$deaths_t <- corona_all.df$deaths/corona_all.df$tests

corona_all.df$cases_ct <- corona_all.df$cases_c/corona_all.df$tests
corona_all.df$deaths_ct <- corona_all.df$deaths_c/corona_all.df$tests

# Generate log
corona_all.df$deaths_c_ln <- log(corona_all.df$deaths_c + 1)



##################################
#### Interpolate vaccinations ####
##################################

# Use cumulative count
corona_all.df$vaccinations_c <- corona_all.df$total_vaccinations_per_hundred
corona_all.df$fully_vaccinationed_c <- corona_all.df$people_fully_vaccinated_per_hundred


### Linear interpolation (only between non-NA values)
interpolfun <- function(x, varx, vary){
  v1 <- x[, varx]
  v2 <- x[, vary]
  nona <- which(!is.na(v2))
  if(length(nona) < 2){
    if(length(nona) == 0){
      res <- rep(0, length(v2))
    }else{
      res <- v2
      res[1:(nona-1)] <- 0
    }
    res2 <- res
  }else{
    # Interpolation
    res <- approx(x = v1, y = v2,
                  xout = v1, rule = 1)$y
    # Extrapolate
    nona <- which(!is.na(res))
    ext <- 1:(nona[1] - 1)
    lm <- lm(res ~ v1, data = data.frame(v1, res))
    pred <- predict(lm, newdata = data.frame(v1 = v1[ext]) )
    pred[pred < 0] <- 0
    
    res2 <- c(pred[ext], res[nona[1]:length(res)])
  }
  return(res2)
}

### Interpolate by id

# Order
corona_all.df$day <- as.numeric(as.Date(corona_all.df$date, format = "%Y%m%d"))
corona_all.df <- corona_all.df[order(corona_all.df$iso_code, corona_all.df$day), ]

# Interpolation, loop over vars
vars <- c("vaccinations_c", "fully_vaccinationed_c")
for(v in vars){
  imp <- by(corona_all.df, corona_all.df$iso_code,
            FUN = function(x) interpolfun(x = x, varx = "day", vary = v))
  corona_all.df[, paste0(v, "_imp")] <- unlist(imp)
}




### Save data
save(corona_all.df, file = "Corona_all_daily.RData")








##########################
#### Make weekly data ####
##########################

### Create day since day 2020
corona_all.df$date <- as.numeric(as.Date(corona_all.df$date, format = "%Y%m%d"))
corona_all.df$day <- corona_all.df$date - as.numeric(as.Date("2019-12-29", format = "%Y-%m-%d"))                              
                                 

### Week
d <- seq(1, max(corona_all.df$day, na.rm = TRUE)+7, 7)
corona_all.df$week <- cut(corona_all.df$day, breaks = d, include.lowest = TRUE, right = FALSE)
corona_all.df$week <- ave(corona_all.df$date, 
                          corona_all.df$week,
                          FUN = function(x) max(x))
corona_all.df$week <- as.Date(corona_all.df$week, origin = "1970-01-01")
corona_all.df$day <- NULL

### Average over variables
var_sum <- c("deaths",
             "deaths_c",
             "deaths_t",
             "deaths_ct",
             "cases",
             "cases_c",
             "cases_t",
             "cases_ct",
             "tests")

# Set negative values to NA
corona_all.df[, var_sum] <- apply(corona_all.df[, var_sum], 2,
                                  FUN = function(x){
                                    x <- ifelse(x < 0, NA, x)
                                    return(x)
                                  })

for(i in var_sum){
  v1 <- paste0(i, "_week")
  corona_all.df[, v1] <- ave(corona_all.df[, i],
                             corona_all.df$countryname, corona_all.df$week,
                             FUN = function(x) mean(x, na.rm = TRUE))
}


### Sum over variables
var_sum <- c("deaths",
             "deaths_c",
             "deaths_t",
             "deaths_ct",
             "cases",
             "cases_c",
             "cases_t",
             "cases_ct",
             "tests")

# Set negative values to NA
corona_all.df[, var_sum] <- apply(corona_all.df[, var_sum], 2,
                                  FUN = function(x){
                                    x <- ifelse(x < 0, NA, x)
                                    return(x)
                                  })

for(i in var_sum){
  v1 <- paste0(i, "_weeksum")
  corona_all.df[, v1] <- ave(corona_all.df[, i],
                             corona_all.df$countryname, corona_all.df$week,
                             FUN = function(x) sum(x, na.rm = TRUE))
}





### Mean over variables
var_mean <- c("stringency",
              "c1_schoolclosing",
              "c2_workplaceclosing",
              "c3_cancelpublicevents",
              "c4_restrictionsongatherings",
              "c5_closepublictransport",
              "c6_stayathomerequirements",
              "c7_restrictionsoninternalmovement",
              "c8_internationaltravelcontrols",
              "h1_publicinformationcampaigns",
              "h2_testingpolicy",
              "h3_contacttracing",
              "h1_publicinformationcampaigns",
              "h6_facialcoverings",
              "h7_vaccinationpolicy",
              "h8_protectionofelderlypeople",
              "new_vaccinations",
              "total_vaccinations_per_hundred",
              "people_vaccinated_per_hundred",
              "vaccinations_c_imp",
              "fully_vaccinationed_c_imp")

for(i in var_mean){
  v1 <- paste0(i, "_week")
  corona_all.df[, v1] <- ave(corona_all.df[, i],
                             corona_all.df$countryname, corona_all.df$week,
                             FUN = function(x) mean(x, na.rm = TRUE))
}
for(i in var_mean[-1]){
  v <- gsub("_.*", "_flag", i)
  v1 <- paste0(i, "_week_flag")
  corona_all.df[, v1] <- ave(corona_all.df[, i],
                             corona_all.df$countryname, corona_all.df$week,
                             FUN = function(x) mean(x, na.rm = TRUE))
}



### Reduce to weekly (use seq from day 5)
corona_all.df$date <- as.Date(corona_all.df$date, origin = "1970-01-01")
oo <- which(corona_all.df$date == corona_all.df$week)
corona_all.df <- corona_all.df[oo, ]



### Reduce to relevant variables
vars <- c("iso_code", "countryname", "continent", "date", "week",
         "population"                               ,"population_density",                      
         "median_age"                               ,"aged_65_older"                            ,"aged_70_older"        ,                   
         "gdp_per_capita"                           ,"extreme_poverty"                          ,"cardiovasc_death_rate",                   
         "diabetes_prevalence"                      ,"female_smokers"                           ,"male_smokers"         ,                   
         "handwashing_facilities"                   ,"hospital_beds_per_thousand"               ,"life_expectancy"      ,                   
         "human_development_index"                  ,"stat_cap_all"                            ,
         "stat_cap_meth"                            ,"stat_cap_source"                         ,"stat_cap_time")

vars2 <- names(corona_all.df)[which(grepl("week", names(corona_all.df)))]
corona_all.df <- corona_all.df[, c(vars, vars2)]

### Generate log
corona_all.df$deaths_c_week_ln <- log(corona_all.df$deaths_c_week + 1)

### Save data
save(corona_all.df, file = "Corona_all_weekly.RData")



