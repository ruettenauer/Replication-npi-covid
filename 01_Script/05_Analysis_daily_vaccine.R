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

library(plm)
library(splm)
library(feisr)
library(lfe)
library(texreg)
library(extrafont)
loadfonts()
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(colorspace)

library(gsynth)
library(panelView)
library(parallel)
library(doParallel)

# library(glmnet)
# library(caret)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}




setwd("Path/02_Data")


ncores <- parallel::detectCores(logical = FALSE)




###################
#### Load data ####
###################


load("Corona_all_daily_2.RData")




##############################################
#### Construct new variables for analysis ####
##############################################

### Numerical date
corona_all.df$date <- as.Date(corona_all.df$date, format = "%Y%m%d")
corona_all.df$date_num <- as.numeric(corona_all.df$date)

### date group (3 dates intervals)
corona_all.df$date_gr <- cut(corona_all.df$date_num, breaks = seq(min(corona_all.df$date_num), max(corona_all.df$date_num)+3, 3),
                             include.lowest = TRUE, right = FALSE)


### Christmas / new years dummy
ave.df <- aggregate(corona_all.df[, "deaths_c", drop = FALSE],
              by = list(date = corona_all.df$date),
              FUN = function(x) mean(x, na.rm = TRUE))
ggplot(ave.df[which(ave.df$date > "2020-12-01" & ave.df$date < "2021-01-15"), ],
       aes(x = date, y = deaths_c)) +
  geom_bar(stat = "identity") + coord_flip()

# 24/12 - 03.01
oo <- which(corona_all.df$date >= "2020-12-24" & corona_all.df$date <= "2021-01-03")
corona_all.df$christmas <- 0
corona_all.df$christmas[oo] <- 1

### Temperature squared
corona_all.df$t2m_sq <- corona_all.df$t2m^2


### NPI index
# This is similar to stringency index but uses only a subset of measures
# See codebook: https://github.com/OxCGRT/covid-policy-tracker/blob/master/documentation/index_methodology.md

# Subindizes
vars <- c("c1_schoolclosing"                 ,
          "c2_workplaceclosing"              ,
          "c3_cancelpublicevents"            ,
          "c4_restrictionsongatherings"      ,
          "c5_closepublictransport"          ,
          "c6_stayathomerequirements"        ,
          "c7_restrictionsoninternalmovement",
          "c8_internationaltravelcontrols"   ,
          "h8_protectionofelderlypeople"    ,
          "h2_testingpolicy"                 ,
          "h3_contacttracing"                ,
          "h6_facialcoverings")

for(v in vars){
  v1 <- paste0(substr(v, 1, 2), "_flag")
  x <- corona_all.df[, v]
  max <- max(corona_all.df[, v], na.rm = TRUE)
  if(v %in% c("c8_internationaltravelcontrols", "h2_testingpolicy", "h3_contacttracing")){ # No flag for international travel
    f <- rep(0, length(x))
  }else{
    f <- ifelse(corona_all.df[, v1] > 0.5, 1, 0)
  }
  
  F <- ifelse(max(f, na.rm = TRUE) > 0, 1, 0)
  F_f <- F - f
  F_f <- ifelse(x <= 0.5, x, F_f) # original F - f treated as zero if x is zero, <= take 0.5 here for rounding
  ind <- 100 * ((x - 0.5 * F_f) / max)
  
  # Subindex variable
  n1 <- paste0(substr(v, 1, 2), "_subindex")
  corona_all.df[, n1] <- ind
}

# Mean of subindizes
vars <- paste0(substr(vars, 1, 2), "_subindex")
corona_all.df$npi_index <- rowMeans(corona_all.df[, vars], na.rm = TRUE)

summary(corona_all.df$npi_index)
# hist(corona_all.df$npi_index)


### Dichotomize instances

vars <- c("c1_schoolclosing"                 ,
          "c2_workplaceclosing"              ,
          "c3_cancelpublicevents"            ,
          "c4_restrictionsongatherings"      ,
          "c5_closepublictransport"          ,
          "c6_stayathomerequirements"        ,
          "c7_restrictionsoninternalmovement",
          "c8_internationaltravelcontrols"   ,
          "h8_protectionofelderlypeople"    ,
          "h2_testingpolicy"                 ,
          "h3_contacttracing"                ,
          "h6_facialcoverings")

for(v in vars){
  n1 <- sub("_.*", "", v)
  n1 <- paste0(n1, "_impl")
  
  # use maximum value observe overall
  maxv <- max(corona_all.df[, v], na.rm = TRUE) 
  
  # # for masks, use category 3:
  # # 3 - Required in all shared/public spaces outside the home with other people present or all situations when social distancing not possible
  # if(v == "h6_facialcoverings_date"){
  #   maxv <- 3
  # }

  
  corona_all.df[, n1] <- NA
  corona_all.df[which(round(corona_all.df[, v], 0) < maxv), n1] <- 0
  corona_all.df[which(round(corona_all.df[, v], 0) >= maxv), n1] <- 1 # set one if rounded value == max
  
  print(n1)
  print(table(corona_all.df[, n1]))
}



#########################################################
#### Vaccination at above 50% of population received ####
#########################################################

oo <- which(corona_all.df$vaccinations_c_imp > 70)
corona_all.df$x1_impl <- ifelse(!is.na(corona_all.df$vaccinations_c_imp), 0, NA)
corona_all.df$x1_impl[oo] <- 1

# Use metric variable as subindex
corona_all.df$x1_subindex <- corona_all.df$fully_vaccinationed_c_imp


#######################
#### Reduce sample ####
#######################


### Exclude countries with zero cases
corona_all.df$meandeaths <- ave(corona_all.df$deaths_c,
                                corona_all.df$countryname,
                                FUN = function(x) mean(x, na.rm = TRUE))
summary(corona_all.df$meandeaths)
table(corona_all.df$iso_code[which(corona_all.df$meandeaths == 0 | is.na(corona_all.df$meandeaths))])

corona_sample.df <- corona_all.df[which(corona_all.df$meandeaths != 0), ]




### ----- ATTENTION: Set NA all dates with negative cases / deaths ----- ###

corona_sample.df$cases_c[which(corona_sample.df$cases_c < 0)] <- NA
corona_sample.df$deaths_c[which(corona_sample.df$deaths_c < 0)] <- NA






##########################################
#### Fill zero cases before reporting ####
##########################################




### Make assumption: If country starts with zero cases, then impute zero before that


### Cases
# First reported number
corona_sample.df$startcases <- ave(corona_sample.df$cases_c,
                                      corona_sample.df$iso_code,
                                      FUN = function(x) x[which(!is.na(x))][1])
# First reporting date
corona_sample.df$startdate <- corona_sample.df$date
corona_sample.df$startdate[which(is.na(corona_sample.df$cases_c))] <- NA
corona_sample.df$startdate <- ave(corona_sample.df$startdate,
                                     corona_sample.df$iso_code,
                                     FUN = function(x) min(x, na.rm = TRUE))

# Impute zeros
oo <- which(corona_sample.df$startdeths == 0 & corona_sample.df$date < corona_sample.df$startdate)

corona_sample.df$cases_s_sm <- corona_sample.df$cases_c
corona_sample.df$cases_s_sm[oo] <- 0


### Deaths
# First reported number
corona_sample.df$startdeaths <- ave(corona_sample.df$deaths_c,
                                       corona_sample.df$iso_code,
                                       FUN = function(x) x[which(!is.na(x))][1])
# First reporting date
corona_sample.df$startdate <- corona_sample.df$date
corona_sample.df$startdate[which(is.na(corona_sample.df$deaths_c))] <- NA
corona_sample.df$startdate <- ave(corona_sample.df$startdate,
                                     corona_sample.df$iso_code,
                                     FUN = function(x) min(x, na.rm = TRUE))

# Impute zeros
oo <- which(corona_sample.df$startdeths == 0 & corona_sample.df$date < corona_sample.df$startdate)

corona_sample.df$deaths_s_sm <- corona_sample.df$deaths_c
corona_sample.df$deaths_s_sm[oo] <- 0






##########################
#### Rolling averages ####
##########################

corona_sample.df <- corona_sample.df[order(corona_sample.df$iso_code, corona_sample.df$date),]


# 7 Rolling ave in deaths 
corona_sample.df$deaths_rolling <- ave(corona_sample.df$deaths_s_sm,
                                       corona_sample.df$iso_code,
                                       FUN = function(x) zoo::rollmean(x, 7, fill = NA, align = "right"))
# 7 Rolling ave in cases 
corona_sample.df$cases_rolling <- ave(corona_sample.df$cases_s_sm,
                                      corona_sample.df$iso_code,
                                      FUN = function(x) zoo::rollmean(x, 7, fill = NA, align = "right"))


### One week lag in 7d rolling average
corona_sample.df$l_deaths_rolling <- ave(corona_sample.df$deaths_rolling,
                                       corona_sample.df$iso_code,
                                       FUN = function(x) dplyr::lag(x, 7))
corona_sample.df$l_cases_rolling <- ave(corona_sample.df$cases_rolling,
                                        corona_sample.df$iso_code,
                                        FUN = function(x) dplyr::lag(x, 7))


### Rolling average of NPIs and lags

vars <- c("c1_schoolclosing"                 ,
          "c2_workplaceclosing"              ,
          "c3_cancelpublicevents"            ,
          "c4_restrictionsongatherings"      ,
          "c5_closepublictransport"          ,
          "c6_stayathomerequirements"        ,
          "c7_restrictionsoninternalmovement",
          "c8_internationaltravelcontrols"   ,
          "h8_protectionofelderlypeople"    ,
          "h2_testingpolicy"                 ,
          "h3_contacttracing"                ,
          "h6_facialcoverings",
          "x1_vaccination")

for(v in vars){
  n1 <- sub("_.*", "", v)
  i1 <- paste0(n1, "_subindex")
  
  r1 <- paste0(n1, "_subindex_rolling")
  
  l1 <- paste0(n1, "_subindex_rolling_lag1")
  l2 <- paste0(n1, "_subindex_rolling_lag2")
  l3 <- paste0(n1, "_subindex_rolling_lag3")
  
  # Rolling av
  corona_sample.df[, r1] <- ave(corona_sample.df[, i1],
                                corona_sample.df$iso_code,
                                FUN = function(x) zoo::rollmean(x, 7, fill = NA, align = "right"))
  # lags
  corona_sample.df[, l1] <- ave(corona_sample.df[, r1],
                                corona_sample.df$iso_code,
                                FUN = function(x) dplyr::lag(x, 7))
  corona_sample.df[, l2] <- ave(corona_sample.df[, r1],
                                corona_sample.df$iso_code,
                                FUN = function(x) dplyr::lag(x, 14))
  corona_sample.df[, l3] <- ave(corona_sample.df[, r1],
                                corona_sample.df$iso_code,
                                FUN = function(x) dplyr::lag(x, 21))
}






####---------------------------------------####
#### Impact functions of single indicators ####
####---------------------------------------####



###############################
#### Make event time clock ####
###############################

vars <- c("c1_schoolclosing"                 ,
          "c2_workplaceclosing"              ,
          "c3_cancelpublicevents"            ,
          "c4_restrictionsongatherings"      ,
          "c5_closepublictransport"          ,
          "c6_stayathomerequirements"        ,
          "c7_restrictionsoninternalmovement",
          "c8_internationaltravelcontrols"   ,
          "h8_protectionofelderlypeople"    ,
          "h2_testingpolicy"                 ,
          "h3_contacttracing"                ,
          "h6_facialcoverings",
          "x1_vaccination")



### Create indicator for running number of implementation and event time clock within each running number
corona_sample.df <- corona_sample.df[order(corona_sample.df$iso_code, corona_sample.df$date), ]
for(v in vars){
  n <- sub("_.*", "", v)
  v1 <- paste0(n, "_impl")
  n1 <- paste0(n, "_impl_rn")
  n2 <- paste0(n, "_impl_et")
  n3 <- paste0(n, "_impl_etneg")
  n4 <- paste0(n, "_impl_ettotal")
  n5 <- paste0(n, "_impl_nt")
  
  # Order
  corona_sample.df <- corona_sample.df[order(corona_sample.df$iso_code, corona_sample.df$date), ]
  
  # Replace single missing days where day before and after is known and equal
  lag <- ave(corona_sample.df[, v1],
             corona_sample.df$iso_code,
             FUN = function(x) dplyr::lag(x, 1))
  lead <- ave(corona_sample.df[, v1],
              corona_sample.df$iso_code,
              FUN = function(x) dplyr::lead(x, 1))
  corona_sample.df[, v1] <- ifelse(is.na(corona_sample.df[, v1]) & lag == lead, lag, corona_sample.df[, v1]) 
  
  # Count 0 and 1 until transition from 1 to 0, than count new running number
  corona_sample.df$tmp <- ave(corona_sample.df[, v1],
                           corona_sample.df$iso_code,
                           FUN = function(x) dplyr::lag(x, 1, default = 0) - x)
  corona_sample.df$tmp[which(corona_sample.df$tmp < 0)] <- 0
  corona_sample.df$tmp[is.na(corona_sample.df$tmp)] <- 0
  
  corona_sample.df[, n1] <- ave(corona_sample.df[, "tmp"],
                             corona_sample.df$iso_code,
                             FUN = function(x) cumsum(x))
  corona_sample.df$tmp <- NULL
  
  # Within each running number of implementation, use cumsum to count, start with date of implementation
  corona_sample.df[, n2] <- ave(corona_sample.df[, v1],
                             corona_sample.df$iso_code, corona_sample.df[, n1],
                             FUN = function(x) cumsum(x))
  
  ### Reverse for negative counts within each running number
  corona_sample.df <- corona_sample.df[order(corona_sample.df$iso_code, -as.numeric(corona_sample.df$date)),]
  
  # Count zeros within each running number (starting with one for period before treatment)
  corona_sample.df[, n3] <- -ave(1 - corona_sample.df[, v1],
                              corona_sample.df$iso_code, corona_sample.df[, n1],
                              FUN = function(x) cumsum(x))
  
  # Replace positives with zero
  corona_sample.df[which(corona_sample.df[, n3] > 0), n3] <- 0
  
  # If there is no treatment within the running number set negative counts to zero
  oo <- ave(corona_sample.df[, v1],
            corona_sample.df$iso_code, corona_sample.df[, n1],
            FUN = function(x) max(x, na.rm = TRUE))
  oo <- which(oo == 0)
  corona_sample.df[oo, n3] <- 0
  
  # Save as never treated indicator
  corona_sample.df[ , n5] <- 0
  corona_sample.df[oo, n5] <- 1
  
  # Combine to total event count
  corona_sample.df[ , n4] <- corona_sample.df[, n3]
  oo <- which(corona_sample.df[, n3] == 0 & corona_sample.df[, n2] > 0)
  corona_sample.df[oo, n4] <- corona_sample.df[oo, n2]
  
  # Sort again forwards
  corona_sample.df <- corona_sample.df[order(corona_sample.df$iso_code, as.numeric(corona_sample.df$date)),]
  
  
}


### Relevel with zero as reference category (include all before 5 dates)
for(v in vars){
  n <- sub("_.*", "", v)
  n1 <- paste0(n, "_impl_ettotal")
  n2 <- paste0(n, "_impl_ettotal2")
  n5 <- paste0(n, "_impl_nt")
  
  
  corona_sample.df[, n2] <- corona_sample.df[, n1]
  corona_sample.df[which(corona_sample.df[, n2] < -35), n2] <- 0
  
  # Never treated periods as own category
  corona_sample.df[which(corona_sample.df[, n5] == 1), n2] <- 999
  
  corona_sample.df[, n2] <- as.factor(corona_sample.df[, n2])
  corona_sample.df[, n2] <- relevel(corona_sample.df[, n2], ref = "0")
  
  print(table(corona_sample.df[, n2]))
}








####-----------------------####
#### Construct lag of NPIS ####
####-----------------------####


### Lag of NPIs
vars <- c("c1_schoolclosing"                 ,
          "c2_workplaceclosing"              ,
          # "c3_cancelpublicevents"            ,
          # "c4_restrictionsongatherings"      ,
          "c5_closepublictransport"          ,
          "c6_stayathomerequirements"        ,
          "c7_restrictionsoninternalmovement",
          "c8_internationaltravelcontrols"   ,
          "h8_protectionofelderlypeople"    ,
          "h2_testingpolicy"                 ,
          "h3_contacttracing"                ,
          "h6_facialcoverings",
          "x1_vaccination")

for(v in vars){
  n <- sub("_.*", "", v)
  v1 <- paste0(n, "_impl")
  n1 <- paste0(n, "_impl_rn")
  
  # Order
  corona_sample.df <- corona_sample.df[order(corona_sample.df$iso_code, corona_sample.df$date), ]
  
  # Contrstruct lags
  for(i in 7:35){
    nl <- paste0(n, "_npilag", i)
    corona_sample.df[, nl] <- ave(corona_sample.df[, v1],
                                  corona_sample.df$iso_code,
                                  FUN = function(x) dplyr::lag(x, i, default = 0))
  }

  
  
  
}






####--------------------------####
#### Select different subsets ####
####--------------------------####

# ### Full sample
corona_subsample.df <- corona_sample.df

### Fist wave
oo <- which(as.Date(corona_sample.df$date, format = "%Y%m%d") <= "2020-08-31")
corona_subsample.dfa <- corona_sample.df[oo,]

### Second wave
oo <- which(as.Date(corona_sample.df$date, format = "%Y%m%d") > "2020-06-30")
corona_subsample.dfb <- corona_sample.df[oo,]



# Save sample suffixes
# samples <- c("", "a", "b")
samples <- c("b")





####----------------------------------------------####
#### Make cases independent of past NPI treatment ####
####----------------------------------------------####


for(s in samples){
  
  corona_tmp.df <- NULL
  corona_tmp.df <- get(paste0("corona_subsample.df", s))
  
  vars <- c("c1_schoolclosing"                 ,
            "c2_workplaceclosing"              ,
            # "c3_cancelpublicevents"            ,
            # "c4_restrictionsongatherings"      ,
            "c5_closepublictransport"          ,
            "c6_stayathomerequirements"        ,
            "c7_restrictionsoninternalmovement",
            "c8_internationaltravelcontrols"   ,
            "h8_protectionofelderlypeople"    ,
            "h2_testingpolicy"                 ,
            "h3_contacttracing"                ,
            "h6_facialcoverings",
            "x1_vaccination")
  
  
  ### One overall model of lagged NPIs on cases
  
  vs <- names(corona_tmp.df)[grep("_npilag", names(corona_tmp.df))]
  fm <- paste0("cases_s_sm ~ ",
              paste(vs, collapse = " + "),
              " + as.factor(iso_code) + as.factor(date)")
  lm.mod <- lm(as.formula(fm), data = corona_tmp.df, na.action = na.exclude)
  corona_tmp.df$resid <- resid(lm.mod)
  
  
  
  ### Lags of resids
  
  # Order
  corona_tmp.df <- corona_tmp.df[order(corona_tmp.df$iso_code, corona_tmp.df$date), ]
  
  # 7d Rolling ave in residuals
  corona_tmp.df$resid_rolling <- ave(corona_tmp.df$resid,
                                        corona_tmp.df$iso_code,
                                        FUN = function(x) zoo::rollmean(x, 7, fill = NA, align = "right"))
  
  # One/two/three week lag in 7d rolling average
  corona_tmp.df$resid_rolling_lag1 <- ave(corona_tmp.df$resid_rolling,
                                             corona_tmp.df$iso_code,
                                             FUN = function(x) dplyr::lag(x, 7))
  corona_tmp.df$resid_rolling_lag2 <- ave(corona_tmp.df$resid_rolling,
                                             corona_tmp.df$iso_code,
                                             FUN = function(x) dplyr::lag(x, 14))
  corona_tmp.df$resid_rolling_lag3 <- ave(corona_tmp.df$resid_rolling,
                                             corona_tmp.df$iso_code,
                                             FUN = function(x) dplyr::lag(x, 21))
  corona_tmp.df$resid_rolling_lag4 <- ave(corona_tmp.df$resid_rolling,
                                             corona_tmp.df$iso_code,
                                             FUN = function(x) dplyr::lag(x, 28))
  corona_tmp.df$resid_rolling_lag5 <- ave(corona_tmp.df$resid_rolling,
                                             corona_tmp.df$iso_code,
                                             FUN = function(x) dplyr::lag(x, 35))
  corona_tmp.df$resid_rolling_lag6 <- ave(corona_tmp.df$resid_rolling,
                                             corona_tmp.df$iso_code,
                                             FUN = function(x) dplyr::lag(x, 42))
  
  
  ### Reassign name
  assign(paste0("corona_subsample.df", s), corona_tmp.df)
  
  
}









####---------------------------------####
#### Synthetic control methods daily ####
####---------------------------------####

# Define dependent variables
depvars <- c("deaths_s_sm")


### Loop over samples and dep vars
for(s in samples){
  
  corona_tmp.df <- NULL
  corona_tmp.df <- get(paste0("corona_subsample.df", s))
  
  
  # Over depvars
  for(dv in depvars){
    
    dvn <- sub("_.*", "", dv)
    
    # Axis ranges for plots
    if(dvn == "cases"){
      ymax <- 200
      ymin <- -300
    }
    if(dvn == "deaths"){
      ymax <- 2.5
      ymin <- -2.5
    }
    
    vars <- c("c1",
              "c2",
              # "c3",
              # "c4",
              "c5",
              "c6",
              "c7",
              "c8",
              "h8",
              "h2",
              "h3",
              "h6",
              "x1")
    names(vars) <- c("School closing",
                     "Workplace closing",
                     # "Cancel public events",
                     # "Restrict gatherings",
                     "Close public transport",
                     "Stay at home",
                     "Restrict internal movement",
                     "International travel",
                     "Protecting elderly",
                     "Testing policy",
                     "Contact tracing",
                     "Masks",
                     "Vaccinations")
    
    # Introduce backwards shift of treatment (in days)
    if(dvn == "cases"){
      shift <- NULL
    }
    if(dvn == "deaths"){
      shift <- NULL
    }
    
    
    
    ###############################################################
    #### Estimate models with cuts into id-intervention period ####
    ###############################################################
    
    
    ### Loop trough all the measures, while controlling for the others

    synthcont.list2 <- vector(mode = "list", length = length(vars))


    for(v in vars[11]){
      k <- which(vars == v)
      v1 <- paste0(v, "_impl")
      v2 <- paste0(v, "_impl_rn")
      vlist <- paste0(vars, "_impl")
      vlist2 <- paste0(vars, "_subindex")
      
      r1 <- paste0(v, "_resid")
      r2 <- paste0(v, "_resid_rolling")
      r3 <- paste0(v, "_resid_rolling_lag1")
      r4 <- paste0(v, "_resid_rolling_lag2")
      r5 <- paste0(v, "_resid_rolling_lag3")
      r6 <- paste0(v, "_resid_rolling_lag4")
      r7 <- paste0(v, "_resid_rolling_lag5")
      
      if(s == "a"){
        formula <- paste0(dv, " ~ ", v1, " + ", paste(vlist2[-k], collapse = " + "), " + t2m + t2m_sq + tcc + tp + q")
      }else{
        formula <- paste0(dv, " ~ ", v1, " + ", paste(vlist2[-k], collapse = " + "), " + t2m + t2m_sq + tcc + tp + q")
      }

      
      # Temporal ID for each iso + running number in intervention
      tmp_sample.df <- corona_tmp.df[which(!is.na(corona_tmp.df[, v2])),]
      tmp_sample.df$tmp_id <- paste0(tmp_sample.df$iso_code, "_", tmp_sample.df[, v2])
      
      # If deaths, add residualized cases as control
      if(dvn == "deaths"){
        
        for(i in 1:5){
          for(j in 2:3){
            tmp_sample.df[, paste0("resid_rolling_lag", i, "_p", j)] <- tmp_sample.df[, paste0("resid_rolling_lag", i)]^j
          }
          
        }
        
        formula <- paste0(formula, " + ",
                          "resid_rolling_lag1", " + ",
                          "resid_rolling_lag2", " + ",
                          "resid_rolling_lag3", " + ",
                          "resid_rolling_lag4", " + ",
                          "resid_rolling_lag5", " + ",
                          "resid_rolling_lag1_p2", " + ",
                          "resid_rolling_lag2_p2", " + ",
                          "resid_rolling_lag3_p2", " + ",
                          "resid_rolling_lag4_p2", " + ",
                          "resid_rolling_lag5_p2", " + ",
                          "resid_rolling_lag1_p3", " + ",
                          "resid_rolling_lag2_p3", " + ",
                          "resid_rolling_lag3_p3", " + ",
                          "resid_rolling_lag4_p3", " + ",
                          "resid_rolling_lag5_p3")


        
      }

      # Drop periods with less than 2 observations in complete cases
      tmp_sample.df <- tmp_sample.df[complete.cases(tmp_sample.df[, all.vars(terms(as.formula(formula)))]), ]
      obs <- ave(rep(1, nrow(tmp_sample.df)), tmp_sample.df$tmp_id, FUN = function(x) sum(x))
      tmp_sample.df <- tmp_sample.df[which(obs >= 2), ]
      
      # Shift treatment 
      if(!is.null(shift)){
        tmp_sample.df <- tmp_sample.df[order(tmp_sample.df$tmp_id, tmp_sample.df$date), ]
        tmp_sample.df[, v1] <- ave(tmp_sample.df[, v1], 
                                   tmp_sample.df$tmp_id, 
                                   FUN = function(x) dplyr::lag(x, shift, default = 0))
        shiftnumb <- shift
      }else{
        shiftnumb <- 0
      }
      

      # Drop time-points with only one observation
      obs <- ave(rep(1, nrow(tmp_sample.df)), tmp_sample.df$date, FUN = function(x) sum(x))
      tmp_sample.df <- tmp_sample.df[which(obs >= 2), ]
      
      # Add some random variation at last (anyway missing period) to christmas to circumvent gsynth variance test
      oo <- which(tmp_sample.df$date == max(tmp_sample.df$date))
      tmp_sample.df$christmas[oo] <- rnorm(length(oo), mean = 0, sd = 0.001) 

      # Drop time-points without control groups (to avoid error in gsynth)
      if(s == "a"){
        pt <- ave(tmp_sample.df[, v1], 
                  tmp_sample.df$tmp_id,
                  FUN = function(x) sum(ifelse(x == 0, 1, 0)))
        oo <- which(pt >= 35)
        dts <- unique(tmp_sample.df$date[oo])
        tmp_sample.df <- tmp_sample.df[tmp_sample.df$date %in% dts, ]
        
        cl = NULL # avoid problems with clustering of small n treated in wave a
      }else{
        cl = "iso_code"
      }

      ### Synth with cross-validation of factors and minimum 5 pre-treatment time periods
      sc_1.mod <- NULL
      sc_1.mod <- gsynth(as.formula(formula),
                         data = tmp_sample.df, na.rm = TRUE, estimator = "mc",
                         index = c("tmp_id", "date"), min.T0 = 35 + shiftnumb,
                         force = "two-way", CV = TRUE, r = c(0, 5), se = TRUE,
                         nlambda = 20, k = 10,
                         inference = "nonparametric", cl = cl, nboots = 1000, seed = 5451528,
                         parallel = TRUE, cores = ncores, tol = 1e-3, EM = FALSE)

      print(sc_1.mod)
      
      sc_1.mod <- sc_1.mod[which(names(sc_1.mod) %in% c("Y.dat", "Y",  "D",  "X", "W", "index",  "id", "time", 
                                                        "obs.missing", "id.tr", "id.co", "D.tr", "I.tr", "Y.tr", 
                                                        "Y.ct", "Y.co", "eff", "Y.bar", "Ntr", "Nco", "T0",  "tr", 
                                                        "pre", "post", "lambda.cv", "beta", "est.att", "att.boot",
                                                        "est.att.avg", "est.beta", "est.ind", "wgt.implied",
                                                        "factor", "lambda.co", "lambda.tr"))]
      

      synthcont.list2[[k]] <- sc_1.mod

    }
    save(synthcont.list2, file = paste0(s, "SynthControl_vacc_", dvn, "_cut.RData"))
    
    
    ### Create plot
    plot.list1 <- vector(mode = "list", length = length(vars))
    plot.list2 <- vector(mode = "list", length = length(vars))
    
    
    for(v in vars[11]){
      k <- which(vars == v)
      name <- names(vars)[which(vars == v)]
      sc.tmp <- synthcont.list2[[k]]
      
      # Reverse Shift of treatment again 
      if(!is.null(shift)){
        rownames(sc.tmp$est.att) <- as.numeric(rownames(sc.tmp$est.att)) + shift
      }
      
      # Dataframe
      coef.df <- data.frame(sc.tmp$est.att[which(rownames(sc.tmp$est.att) %in% c(-34:151)), ])
      coef.df$date <- as.numeric(rownames(coef.df)) -1
      oo <- which(is.na(coef.df$ATT))
      coef.df$CI.lower[oo] <- NA; coef.df$CI.upper[oo] <- NA
      
      # Pretreatment trend
      if(nrow(coef.df[coef.df$date < 0, ]) == 0 | all(is.na(coef.df$ATT))){ # NUll model if synth coefs are NA
        tmp <- lm(rep(0, 35) ~ 1 +  rnorm(35, 0, 0.0001))
      }else{
        tmp <- lm(ATT ~ date, coef.df[coef.df$date < 0, ])
      }
      
      
      # N treated
      coef.df$ntreated <- cut(coef.df$n.Treated, breaks = c(0, 5, 10, 20, 30,  max(max(coef.df$n.Treated), 31)), right = TRUE)
      coef.df$ribmin <- -Inf
      coef.df$ribmax <- ymin + abs(ymin*0.05)
      
      # Duplicate if only one date per value ad drop single day jumps
      rib.df <- coef.df[coef.df$date >= 0, ]
      rib.df$tmp <- dplyr::lag(as.numeric(rib.df$ntreated)) - as.numeric(rib.df$ntreated)
      rib.df$tmp[is.na(rib.df$tmp)] <- 0
      rib.df$tmp <- cumsum(abs(rib.df$tmp))
      rib.df$tmp[is.na(rib.df$ntreated)] <- 99
      
      levels(rib.df$ntreated)[5] <- "> 30"
      
      cols_rib <- c("#ffffcc",
                    "#a1dab4",
                    "#41b6c4",
                    "#2c7fb8",
                    "#253494")
      
      p1 <- ggplot(coef.df, aes(date, ATT))+
        geom_ribbon(data = rib.df, aes(x = date, ymin = ribmin, ymax = ribmax, fill = ntreated, group = tmp), 
                    alpha = 1, na.rm = TRUE, show.legend = TRUE) +
        geom_line(data = coef.df, lwd = 1.0)+
        geom_ribbon(data = coef.df,aes(ymin = CI.lower, ymax = CI.upper), alpha = 0.2)
      p1 <- p1 + geom_hline(yintercept = 0, colour = "black", lty = 3, lwd = 1, alpha = 0.7)+
        geom_abline(intercept = tmp$coefficients[1], slope = tmp$coefficients[2],
                    colour = "black", lty = 2, lwd = 1, alpha = 0.7)+
        geom_vline(xintercept = 0, colour = gray(1/2), lty = 1, lwd = 1)
      p1 <- p1 + scale_x_continuous(limits = c(-35, 150), breaks = seq(-30, 150, 30))
      p1 <- p1 + scale_fill_manual(values = cols_rib, drop = FALSE)
      p1 <- p1 + coord_cartesian(ylim = c(ymin, ymax))
      p1 <- p1 + labs(y = paste0("Change in COVID-19 ", dvn, " / mio capita"),
                      x = "Days since intervention") 
      p1 <- p1 + theme_bw()
      p1  <-  p1 + theme(text = element_text( size = 14),
                         legend.position = c(0.01,0.99), legend.justification = c(0.01,0.99),
                         legend.background = element_blank(),
                         legend.key = element_blank(),
                         axis.text.y = element_text(colour = "black"),
                         plot.title = element_text(hjust = 0.5),
                         axis.title.x = element_text(colour = "black", margin = margin(t = 10, r = 0 , b = 0, l = 0)),
                         axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 10 , b = 0, l = 0)),
                         axis.text.x = element_text(colour = "black"),
      )
      p1 <- p1 + ggtitle(names(vars)[which(vars == v)])
      p1 <- p1 + guides(fill = guide_legend(title = "N treated", reverse = T))
      p1
      
      png(file = paste0("../03_Output/", s, "Vaccination_", dvn, "_cut", ".png"), width = 7, height = 5,
          units = "in", res = 300)
      par(mar = c(0, 0, 0, 0))
      par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
      print(p1)
      dev.off()
      
      pdf(file = paste0("../03_Output/", s, "Vaccination_", dvn, "_cut", ".pdf"), width = 7, height = 5)
      par(mar = c(0, 0, 0, 0))
      par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
      print(p1)
      dev.off()
      
      
    }
    

    
    
  }
}


