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
library(texreg)
library(extrafont)
loadfonts()
library(ggplot2)
library(scales)
library(cowplot)
library(grid)
library(gridExtra)
library(ggridges)
library(ungeviz)
library(colorspace)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

setwd("Path/02_Data")



###################
#### Load data ####
###################


load("Corona_all_weekly_2.RData")




##############################################
#### Construct new variables for analysis ####
##############################################

### Numerical week
corona_all.df$week_num <- as.numeric(as.factor(corona_all.df$week))

### Temperature squared
corona_all.df$t2m_sq <- corona_all.df$t2m^2


### Week group (3 weeks intervals)
corona_all.df$week_gr <- cut(corona_all.df$week_num, breaks = seq(min(corona_all.df$week_num)-3, max(corona_all.df$week_num)+3, 3),
                             include.lowest = FALSE)


### NPI index
# This is similar to stringency index but uses only a subset of measures
# See codebook: https://github.com/OxCGRT/COVID-policy-tracker/blob/master/documentation/index_methodology.md

# Subindizes
vars <- c("c1_schoolclosing_week"                 ,
          "c2_workplaceclosing_week"              ,
          "c3_cancelpublicevents_week"            ,
          "c4_restrictionsongatherings_week"      ,
          "c5_closepublictransport_week"          ,
          "c6_stayathomerequirements_week"        ,
          "c7_restrictionsoninternalmovement_week",
          "c8_internationaltravelcontrols_week"   ,
          "h8_protectionofelderlypeople_week"    ,
          "h2_testingpolicy_week"                 ,
          "h3_contacttracing_week"                ,
          "h6_facialcoverings_week")

for(v in vars){
  v1 <- paste0(v, "_flag")
  x <- corona_all.df[, v]
  max <- max(corona_all.df[, v], na.rm = TRUE)
  f <- ifelse(corona_all.df[, v1] > 0.5, 1, 0)
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
hist(corona_all.df$npi_index)





####----------------------------------####
#### Descriptives of single countries ####
####----------------------------------####

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

### Line plot for single countries
ids <- unique(corona_all.df$countryname)
ids <- ids[order(ids)]

plot.list <- vector(mode = "list", length = length(ids))
names(plot.list) <- ids

c <- seq(1, length(plot.list), 16)
k <- 1
for(i in ids){
  tmp <- corona_all.df[corona_all.df$countryname == i , ]
  name <- unique(tmp$countryname)

  max1 <- max(tmp$deaths_c_week, na.rm = T)
  max2 <- max(tmp$npi_index, na.rm = T)
  r <- max1/100
  if(!is.infinite(max1) & !is.infinite(max2)){
    cols <- gg_color_hue(2)
    
    p1 <- ggplot(tmp, aes(x = week_num, y = deaths_c_week)) +
      geom_blank(aes(x = 0, y = 0)) +
      geom_line(aes(x = week_num, y = deaths_c_week, colour = "COVID Deaths"), lwd = 1.5) +
      geom_point(aes(x = week_num, y = deaths_c_week, colour = "COVID Deaths", shape = "COVID Deaths"), size  = 3) +
      geom_line(aes(x = week_num, y = npi_index*r, colour = "NPI Index"), lwd = 1.5) +
      geom_point(aes(x = week_num, y = npi_index*r, colour = "NPI Index", shape = "NPI Index"), size  = 3) +
      scale_color_manual(name = "", values = colorspace::darken(cols, 0.2), labels = c("COVID Deaths", "NPI Index")) + 
      scale_shape_manual(name = "", values = c(15, 17), labels = c("COVID Deaths", "NPI Index")) +
      scale_y_continuous(sec.axis = sec_axis(~./r, name = "Overall NPI Index")) +
      xlab("Week") + ylab("COVID deaths per mio capita") +
      ggtitle(name) +
      theme_bw() +
      theme(text = element_text(size = 20),
            axis.text.y = element_text(colour = "black"),
            axis.text.x = element_text(colour = "black"),
            axis.title.x = element_text(colour = "black",  margin = margin(t = 10, r = 0 , b = 0, l = 0)),
            axis.title.y = element_text(colour = "black",  margin = margin(t = 0, r = 15 , b = 0, l = 0)),
            axis.title.y.right = element_text(colour = "black",  margin = margin(t = 0, r = 0 , b = 0, l = 15)),
            plot.title = element_text(hjust = 0.5),
            legend.key = element_blank(), legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.position = c(0.01,0.93), legend.justification = 'left',
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black"),
            legend.spacing.y = unit(-0.1, "cm"), legend.key.size = unit(0.7,"line"))
    print(p1)
    plot.list[[i]] <- recordPlot() # use record, as saving in list messes with second axis scale
    dev.off()
    
    # # Export
    # pdf(file = paste0("./Output/", "Country_Tractory_", i, ".pdf"), width = 9, height = 7)
    # par(mar = c(0, 0, 0, 0))
    # par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
    # print(p1)
    # dev.off()
  }
  
}

# Export
n <- which(lapply(plot.list, FUN = function(x) !is.null(x)) == TRUE)
plot.list <- plot.list[n]
n <- names(plot.list)
pdf(file = paste0("../03_Output/", "Country_Tractories", ".pdf"), width = 9, height = 7)
par(mar = c(0, 0, 0, 0))
par(mfrow = c(4, 4), oma = c(0, 0, 0, 0))
for(i in 1:length(plot.list)){
  replayPlot(plot.list[[n[i]]])
}
dev.off()




### Line plot for single countries with temperature
ids <- unique(corona_all.df$countryname)
ids <- ids[order(ids)]

plot.list <- vector(mode = "list", length = length(ids))
names(plot.list) <- ids

c <- seq(1, length(plot.list), 16)
k <- 1
for(i in ids){
  tmp <- corona_all.df[corona_all.df$countryname == i , ]
  name <- unique(tmp$countryname)
  
  max1 <- max(tmp$deaths_c_week, na.rm = T)
  max2 <- max(tmp$t2m, na.rm = T)
  r <- max1/max2
  if(!is.infinite(max1) & !is.infinite(max2)){
    cols <- gg_color_hue(2)
    
    p1 <- ggplot(tmp, aes(x = week_num, y = deaths_c_week)) +
      geom_blank(aes(x = 0, y = 0)) +
      geom_line(aes(x = week_num, y = deaths_c_week, colour = "COVID Deaths"), lwd = 1.5) +
      geom_point(aes(x = week_num, y = deaths_c_week, colour = "COVID Deaths", shape = "COVID Deaths"), size  = 3) +
      geom_line(aes(x = week_num, y = t2m*r, colour = "Temperature"), lwd = 1.5) +
      geom_point(aes(x = week_num, y = t2m*r, colour = "Temperature", shape = ""), size  = 3) +
      scale_color_manual(name = "", values = colorspace::darken(cols, 0.2), labels = c("COVID Deaths", "Temperature")) + 
      scale_shape_manual(name = "", values = c(15, 17), labels = c("COVID Deaths", "Temperature")) +
      scale_y_continuous(sec.axis = sec_axis(~./r, name = "Temperature (degrees Celsius)")) +
      xlab("Week") + ylab("COVID deaths per mio capita") +
      ggtitle(name) +
      theme_bw() +
      theme(text = element_text(size = 20),
            axis.text.y = element_text(colour = "black"),
            axis.text.x = element_text(colour = "black"),
            axis.title.x = element_text(colour = "black",  margin = margin(t = 10, r = 0 , b = 0, l = 0)),
            axis.title.y = element_text(colour = "black",  margin = margin(t = 0, r = 15 , b = 0, l = 0)),
            axis.title.y.right = element_text(colour = "black",  margin = margin(t = 0, r = 0 , b = 0, l = 15)),
            plot.title = element_text(hjust = 0.5),
            legend.key = element_blank(), legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.position = c(0.01,0.93), legend.justification = 'left',
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black"),
            legend.spacing.y = unit(-0.1, "cm"), legend.key.size = unit(0.7,"line"))
    print(p1)
    plot.list[[i]] <- recordPlot() # use record, as saving in list messes with second axis scale
    dev.off()
    
    # # Export
    # pdf(file = paste0("./Output/", "Country_Tractory_", i, ".pdf"), width = 9, height = 7)
    # par(mar = c(0, 0, 0, 0))
    # par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
    # print(p1)
    # dev.off()
  }
  
}

# Export
n <- which(lapply(plot.list, FUN = function(x) !is.null(x)) == TRUE)
plot.list <- plot.list[n]
n <- names(plot.list)
pdf(file = paste0("../03_Output/", "Country_Tractories_temperature", ".pdf"), width = 9, height = 7)
par(mar = c(0, 0, 0, 0))
par(mfrow = c(4, 4), oma = c(0, 0, 0, 0))
for(i in 1:length(plot.list)){
  replayPlot(plot.list[[n[i]]])
}
dev.off()





#######################
#### Reduce sample ####
#######################


### Exclude countries with zero cases
corona_all.df$meandeaths <- ave(corona_all.df$deaths_week,
                                corona_all.df$countryname,
                                FUN = function(x) mean(x, na.rm = TRUE))
summary(corona_all.df$meandeaths)
table(corona_all.df$countryname[corona_all.df$meandeaths == 0])

corona_sample.df <- corona_all.df[which(corona_all.df$meandeaths != 0), ]








###################################
#### Single measures indicator ####
###################################

### Dichotomize instances

vars <- c("c1_schoolclosing_week"                 ,
          "c2_workplaceclosing_week"              ,
          "c3_cancelpublicevents_week"            ,
          "c4_restrictionsongatherings_week"      ,
          "c5_closepublictransport_week"          ,
          "c6_stayathomerequirements_week"        ,
          "c7_restrictionsoninternalmovement_week",
          "c8_internationaltravelcontrols_week"   ,
          "h8_protectionofelderlypeople_week"    ,
          "h2_testingpolicy_week"                 ,
          "h3_contacttracing_week"                ,
          "h6_facialcoverings_week",
          "h1_publicinformationcampaigns_week")

for(v in vars){
  n1 <- sub("_.*", "", v)
  n1 <- paste0(n1, "_impl")
  
  # use maximum value observe overall
  maxv <- max(corona_sample.df[, v], na.rm = TRUE) 
  
  # # for masks, use category 3:
  # # 3 - Required in all shared/public spaces outside the home with other people present or all situations when social distancing not possible
  # if(v == "h6_facialcoverings_week"){
  #   maxv <- 3
  # }
  
  corona_sample.df[, n1] <- NA
  corona_sample.df[which(round(corona_sample.df[, v], 0) < maxv), n1] <- 0
  corona_sample.df[which(round(corona_sample.df[, v], 0) >= maxv), n1] <- 1 # set one if rounded value == max
  
  print(n1)
  print(table(corona_sample.df[, n1]))
}




###############################
#### Make event time clock ####
###############################



### Create indicator for running number of implementation and event time clock within each running number
corona_sample.df <- corona_sample.df[order(corona_sample.df$iso_code, corona_sample.df$week_num), ]
for(v in vars){
  n <- sub("_.*", "", v)
  v1 <- paste0(n, "_impl")
  n1 <- paste0(n, "_impl_rn")
  n2 <- paste0(n, "_impl_et")
  n3 <- paste0(n, "_impl_etneg")
  n4 <- paste0(n, "_impl_ettotal")
  n5 <- paste0(n, "_impl_nt")
  
  # Order
  corona_sample.df <- corona_sample.df[order(corona_sample.df$iso_code, corona_sample.df$week_num), ]
  
  # Count 0 and 1 until transition from 1 to 0, than count new running number
  corona_sample.df$tmp <- ave(corona_sample.df[, v1],
                           corona_sample.df$iso_code,
                           FUN = function(x) dplyr::lag(x, 1, default = 0) - x)
  corona_sample.df$tmp[which(corona_sample.df$tmp < 0)] <- 0
  
  corona_sample.df[, n1] <- ave(corona_sample.df[, "tmp"],
                             corona_sample.df$iso_code,
                             FUN = function(x) cumsum(x))
  corona_sample.df$tmp <- NULL
  
  # Within each running number of implementation, use cumsum to count, start with week of implementation
  corona_sample.df[, n2] <- ave(corona_sample.df[, v1],
                             corona_sample.df$iso_code, corona_sample.df[, n1],
                             FUN = function(x) cumsum(x))
  
  ### Reverse for negative counts within each running number
  corona_sample.df <- corona_sample.df[order(corona_sample.df$iso_code, -as.numeric(corona_sample.df$week_num)),]
  
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
  corona_sample.df <- corona_sample.df[order(corona_sample.df$iso_code, as.numeric(corona_sample.df$week_num)),]
  
  
}


### Relevel with zero as reference category (include all before 5 weeks)
for(v in vars){
  n <- sub("_.*", "", v)
  n1 <- paste0(n, "_impl_ettotal")
  n2 <- paste0(n, "_impl_ettotal2")
  n5 <- paste0(n, "_impl_nt")
  
  
  corona_sample.df[, n2] <- corona_sample.df[, n1]
  corona_sample.df[which(corona_sample.df[, n2] < -5), n2] <- 0
  
  # Never treated periods as own category
  corona_sample.df[which(corona_sample.df[, n5] == 1), n2] <- 999
  
  corona_sample.df[, n2] <- as.factor(corona_sample.df[, n2])
  corona_sample.df[, n2] <- relevel(corona_sample.df[, n2], ref = "0")
  
  print(table(corona_sample.df[, n2]))
}









####-------------------------------------####
#### Deaths cases, and median stringency ####
####-------------------------------------####




country_av1.df <- aggregate(corona_sample.df[, c("cases_c_weeksum", "deaths_c_weeksum")],
                            by = list(iso_code = corona_sample.df$iso_code),
                            FUN = function(x) sum(x, na.rm = TRUE))
country_av2.df <- aggregate(corona_sample.df[, c("npi_index", "stringency_week")],
                            by = list(iso_code = corona_sample.df$iso_code),
                            FUN = function(x) median(x, na.rm = TRUE))

country_av.df <- merge(country_av1.df, country_av2.df, by = "iso_code")


### Order
country_av.df <- country_av.df[order(country_av.df$deaths_c_weeksum),]
country_av.df$iso_code <- factor(country_av.df$iso_code, levels = country_av.df$iso_code)


### Plot
cols <- rev(gg_color_hue(2))
max1 <- max(country_av.df$deaths_c_weeksum, na.rm = T)
max2 <- max(country_av.df$npi_index, na.rm = T)
r <- max1/max2

m1 <- mean(country_av.df$deaths_c_weeksum)
m2 <- mean(country_av.df$npi_index)

country_av.df$page <- 1
country_av.df$page[round((nrow(country_av.df)/2)):nrow(country_av.df)] <- 2


p1 <- ggplot(country_av.df[country_av.df$page == 1,], 
             aes(x = deaths_c_weeksum, y = iso_code)) +
  geom_point(aes(x = deaths_c_weeksum, colour = "Total COVID-19 deaths per mio capita", shape = "Total COVID-19 deaths per mio capita"), size = 3) +
  geom_point(aes(x = npi_index*r, colour = "Median COVID-19 NPI stringency index", shape = "Median COVID-19 NPI stringency index"), size  = 3) +
  scale_x_continuous(sec.axis = sec_axis(~./r, name = "Median COVID-19 NPI stringency index")) +
  xlab("COVID-19 deaths per mio capita") +
  scale_color_manual(name = "", values = cols, labels = rev(c("Total COVID-19 deaths per mio capita", 
                                                              "Median COVID-19 NPI stringency index")),
                     guide = guide_legend(reverse = TRUE)) + 
  scale_shape_manual(name = "", values = c(15, 17), labels = rev(c("Total COVID-19 deaths per mio capita", 
                                                                   "Median COVID-19 NPI stringency index")),
                     guide = guide_legend(reverse = TRUE)) +
  scale_y_discrete(name = NULL, expand = c(0, 0.5), guide = guide_axis(n.dodge = 2)) +
  geom_vline(xintercept = m1, color = cols[2], size = 1.5) +
  geom_vline(xintercept = m2*r, color = cols[1], size = 1.5) +
  geom_blank(aes(x = c(max1), y = c(1))) +
  geom_blank(aes(x = c(min(country_av.df$deaths_c_weeksum)), y = c(1))) +
  theme_bw() +
  theme(axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 14),
        axis.title.x = element_text(colour = "black",  size = 16),
        axis.title.y = element_text(colour = "black",  size = 16),
        legend.key = element_blank(), legend.title = element_blank(),
        legend.text = element_text(size = 14),
        # legend.position = c(0.05,0.95), legend.justification = 'left',
        legend.position = "none",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(-0.1, "cm"), legend.key.size = unit(0.7,"line")
  )
p1

p2 <- ggplot(country_av.df[country_av.df$page == 2,], 
             aes(x = deaths_c_weeksum, y = forcats::fct_reorder(iso_code, deaths_c_weeksum))) +
  geom_point(aes(x = deaths_c_weeksum, colour = "Total COVID-19 deaths per mio capita", shape = "Total COVID-19 deaths per mio capita"), size = 3) +
  geom_point(aes(x = npi_index*r, colour = "Median COVID-19 NPI stringency index", shape = "Median COVID-19 NPI stringency index"), size  = 3) +
  scale_x_continuous(sec.axis = sec_axis(~./r, name = "Median COVID-19 NPI stringency index")) +
  xlab("COVID-19 deaths per mio capita") +
  scale_color_manual(name = "", values = cols, labels = rev(c("Total COVID-19 deaths per mio capita", 
                                                          "Median COVID-19 NPI stringency index")),
                     guide = guide_legend(reverse = TRUE)) + 
  scale_shape_manual(name = "", values = c(15, 17), labels = rev(c("Total COVID-19 deaths per mio capita", 
                                                               "Median COVID-19 NPI stringency index")),
                     guide = guide_legend(reverse = TRUE)) +
  scale_y_discrete(name = NULL, expand = c(0, 0.5), guide = guide_axis(n.dodge = 2)) +
  geom_vline(xintercept = m1, color = cols[2], size = 1.5) +
  geom_vline(xintercept = m2*r, color = cols[1], size = 1.5) +
  geom_blank(aes(x = c(max1), y = c(1))) +
  geom_blank(aes(x = c(min(country_av.df$deaths_c_weeksum)), y = c(1))) +
  geom_blank(aes(x = m1, y = 1)) +
  theme_bw() +
  theme(axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 14),
        axis.title.x = element_text(colour = "black",  size = 16),
        axis.title.y = element_text(colour = "black",  size = 16),
        legend.key = element_blank(), legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = c(0.03,0.9), legend.justification = 'left',
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(-0.1, "cm"), legend.key.size = unit(0.7,"line")
  )
p2



### Export
pdf(file = paste0("../03_Output/", "Descriptives_Country_average", ".pdf"), width = 16, height = 9)
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
grid.arrange(p2, p1, ncol = 2)
dev.off()
### Export
png(file = paste0("../03_Output/", "Descriptives_Country_average", ".png"), width = 16, height = 9,
    units = "in", res = 300)
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
grid.arrange(p2, p1, ncol = 2)
dev.off()





####-----------------------------------------------####
#### Deaths cases, and median stringency over time ####
####-----------------------------------------------####


time_av1.df <- aggregate(corona_sample.df[, c("cases_c_week", "deaths_c_week")],
                            by = list(date = corona_sample.df$date),
                            FUN = function(x) sum(x, na.rm = TRUE))
time_av2.df <- aggregate(corona_sample.df[, c("npi_index", "stringency_week")],
                            by = list(date = corona_sample.df$date),
                            FUN = function(x) mean(x, na.rm = TRUE))

time_av.df <- merge(time_av1.df, time_av2.df, by = "date")


tmp <- data.frame(time = c(time_av.df$date, time_av.df$date),
                  x = c(time_av.df$deaths_c_week, time_av.df$npi_index),
                  value = rep(c(2, 1), each = nrow(time_av.df)))
tmp$time <- as.Date(tmp$time, format = "%Y%m%d")

tmp$value <- factor(tmp$value, levels = c(1, 2), labels = c("Weekly average COVID-19 NPI stringency index",
                                                            "Weekly average COVID-19 deaths per mio capita"))




### Plot 
cols <- gg_color_hue(2)
max1 <- max(time_av.df$deaths_c_week, na.rm = T)
max2 <- max(time_av.df$npi_index, na.rm = T)
r <- (max1/max2) / 1.318
tmp$x[which(tmp$value == "Weekly average COVID-19 NPI stringency index")] <- tmp$x[which(tmp$value == "Weekly average COVID-19 NPI stringency index")]*r

tmp$max <- ave(tmp$x, by = tmp$value, FUN = function(x) max(x))
max <- tmp[which(tmp$x == tmp$max),]


p1 <- ggplot(tmp, aes(x = time, y = x, color = value, fill = value)) + 
  geom_density_line(stat = "identity", lwd = 1.2) +
  scale_color_manual(
    values = darken(cols, 0),
    breaks = c("Weekly average COVID-19 deaths per mio capita", 
               "Weekly average COVID-19 NPI stringency index"),
    #guide = "none"
  ) +
  scale_fill_manual(
    values = scales::alpha(cols, 0.3),
    breaks = c("Weekly average COVID-19 deaths per mio capita", 
               "Weekly average COVID-19 NPI stringency index"),
    #guide = "none"
  ) +
  scale_y_continuous(sec.axis = sec_axis(~./r, name = "COVID-19 NPI stringency index", breaks = seq(0, 100, 20))) +
  scale_x_date(breaks = date_breaks("3 months"), date_labels = "%Y-%m") +
  geom_blank(mapping = aes(x = as.Date("2020-04-01", format = "Y%-%m-%d"), y = 500)) +
  coord_cartesian(clip = "off") +
  xlab("Time") + ylab("COVID-19 deaths per mio capita") +
  theme_bw() +
  theme(axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 14),
        axis.title.x = element_text(colour = "black",  size = 16),
        axis.title.y = element_text(colour = "black",  size = 16),
        legend.key = element_blank(), legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = c(0.01,0.95), legend.justification = 'left',
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(-0.1, "cm"), legend.key.size = unit(0.7,"line")
  )
p1



### Export
pdf(file = paste0("../03_Output/", "Descriptives_Time_average", ".pdf"), width = 12, height = 7)
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
p1
dev.off()
### Export
png(file = paste0("../03_Output/", "Descriptives_Time_average", ".png"), width = 12, height = 7,
    units = "in", res = 300)
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
p1
dev.off()







####---------------------------####
#### Sequence plot of measures ####
####---------------------------####

vars <- c("c1_schoolclosing_week"                 ,
          "c2_workplaceclosing_week"              ,
          "c3_cancelpublicevents_week"            ,
          "c4_restrictionsongatherings_week"      ,
          "c5_closepublictransport_week"          ,
          "c6_stayathomerequirements_week"        ,
          "c7_restrictionsoninternalmovement_week",
          "c8_internationaltravelcontrols_week"   ,
          "h8_protectionofelderlypeople_week"    ,
          "h2_testingpolicy_week"                 ,
          "h3_contacttracing_week"                ,
          "h6_facialcoverings_week",
          "h1_publicinformationcampaigns_week")
names(vars) <- c("School closing",
                 "Workplace closing",
                 "Cancel public events",
                 "Restrict gatherings",
                 "Close public transport",
                 "Stay at home",
                 "Restrict internal movement",
                 "International travel",
                 "Protecting elderly",
                 "Testing policy",
                 "Contact tracing",
                 "Masks",
                 "Public information")

corona_sample.df$time <- as.Date(corona_sample.df$date, format = "%Y%m%d")
corona_sample.df$iso <- factor(corona_sample.df$iso_code, levels = rev(sort(unique(corona_sample.df$iso_code))))

n <- rev(sort(unique(corona_sample.df$iso_code)))
n <- n[1:(length(n)/2)]
corona_sample.df$sep <- ifelse(corona_sample.df$iso_code %in% n, 1, 2)

corona_sample.df$nr <- as.numeric(corona_sample.df$iso) - ave(as.numeric(corona_sample.df$iso),
                                                              corona_sample.df$sep, 
                                                     FUN = function(x) (min(x) - 1))

### Order times
times <- unique(corona_sample.df$time)
times <- format(times,"%Y-%m")
# oo <- c(1:length(times))[-seq(1, length(times), 12)]
qrt <- c("2020-03", "2020-09", "2021-03", "2021-08")
oo <- which(!times %in% qrt)
times[oo] <- ""
times[duplicated(times)] <- ""
ticks <- which(times != "")

corona_sample.df$time2 <- as.numeric(as.factor(corona_sample.df$time))

plot.list <- vector(mode = "list", length = length(vars))

for(v in vars){
  k <- which(vars == v)
  n <- sub("_.*", "", v)
  v1 <- paste0(n, "_impl")
  n1 <- paste0(n, "_impl_rn")
  
  ### Df for single cuts
  cuts.df <- corona_sample.df[, c("iso", "time2", "sep", "nr", v, n1)]
  cuts.df$pnr <- ave(rep(1, nrow(cuts.df)),
                     by = paste0(cuts.df$iso, "_", cuts.df[, n1]),
                     FUN = function(x) cumsum(x))
  cuts.df <- cuts.df[cuts.df$pnr == 1, ]
  
  ### Df for treatment within cuts
  treatment.df <- corona_sample.df[, c("iso", "time2", "sep", "nr", v, v1, n1)]
  treatment.df$pnr <- ave(rep(1, nrow(treatment.df)),
                     by = paste0(treatment.df$iso, "_", treatment.df[, n1], treatment.df[, v1]),
                     FUN = function(x) cumsum(x))
  treatment.df <- treatment.df[treatment.df$pnr == 1 & treatment.df[, v1] == 1, ]
  treatment.df$nr <- as.numeric(treatment.df$iso) - ave(as.numeric(treatment.df$iso),
                                              treatment.df$sep, FUN = function(x) (min(x) - 1))

  
  ### Plot
  
  p1 <- ggplot(data = corona_sample.df[corona_sample.df$sep == 1, ],
               mapping = aes_string(x = "time2", y = "iso", fill = v)) +
    geom_tile() +
    scale_fill_viridis_c(
      option = "A", begin = 0.02, end = 0.98,
      name = names(vars)[k],
      guide = guide_colorbar(
        direction = "horizontal",
        label.position = "bottom",
        title.position = "top",
        ticks = FALSE,
        barwidth = grid::unit(3.5, "in"),
        barheight = grid::unit(0.2, "in")
      )
    ) +
    # geom_segment(data = cuts.df, mapping = aes(x = time, y = iso -0.5, xend = time, yend = iso + 0.5), col = "red") + 
    scale_y_discrete(name = NULL, position = "left", guide = guide_axis(n.dodge = 2),
                     limits = levels(droplevels(corona_sample.df[corona_sample.df$sep == 1, "iso"]))) +
    scale_x_continuous(breaks = ticks,
                       labels = times[which(times != "")],
                       name = "Time", 
                       expand = c(0,0)) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(colour = "black", size = 12),
      axis.text.x = element_text(colour = "black", size = 14),
      axis.title.x = element_text(colour = "black",  size = 16),
      axis.title.y = element_text(colour = "black",  size = 16),
      legend.key = element_blank(), legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      legend.background = element_blank(),
      legend.box.background = element_blank()
    ) +
  geom_vpline(data = cuts.df[cuts.df$sep == 1, ],
             mapping = aes(x = time2 -0.5, y = iso, col = "col1"),
              height = 1.5, lwd = 1.1, show.legend = TRUE)  +
  geom_point(data = treatment.df[treatment.df$sep == 1, ],
                 mapping = aes(x = time2 , y = iso, col = "col2"),
                 pch = 4, size = 1.3, show.legend = TRUE) +
    scale_color_manual(name = "cols", values = c("#4d9221", "#2b8cbe"),
                       labels = c("Period split", "Treatment")) + 
    guides(colour = guide_legend(override.aes = list(size = 5, shape = c("|", "x"), linetype = 0)))
    
    
  p1
  
  p2 <- ggplot(corona_sample.df[corona_sample.df$sep == 2, ],
               aes_string(x = "time2", y = "iso", fill = v)) +
    geom_tile() +
    scale_fill_viridis_c(
      option = "A", begin = 0.05, end = 0.98,
      name = names(vars)[k],
      guide = guide_colorbar(
        direction = "horizontal",
        label.position = "bottom",
        title.position = "top",
        ticks = FALSE,
        barwidth = grid::unit(3.5, "in"),
        barheight = grid::unit(0.2, "in")
      )
    ) +
    scale_y_discrete(name = NULL, position = "left", guide = guide_axis(n.dodge = 2),
                     limits = levels(droplevels(corona_sample.df[corona_sample.df$sep == 2, "iso"]))) +
    scale_x_continuous(breaks = ticks,
                       labels = times[which(times != "")],
                       name = "Time", 
                       expand = c(0,0)) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(colour = "black", size = 12),
      axis.text.x = element_text(colour = "black", size = 14),
      axis.title.x = element_text(colour = "black",  size = 16),
      axis.title.y = element_text(colour = "black",  size = 16),
      legend.key = element_blank(), legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      legend.background = element_blank(),
      legend.box.background = element_blank()
    ) +
    geom_vpline(data = cuts.df[cuts.df$sep == 2, ],
                mapping = aes(x = time2  -0.5, y = iso, col = "col1"),
                height = 1.5, lwd = 1.1, show.legend = TRUE)  +
    geom_point(data = treatment.df[treatment.df$sep == 2, ],
               mapping = aes(x = time2 , y = iso, col = "col2"),
               pch = 4, size = 1.3, show.legend = TRUE) +
    scale_color_manual(name = "cols", values = c("#4d9221", "#2b8cbe"),
                       labels = c("Period split", "Treatment")) + 
    guides(colour = guide_legend(override.aes = list(size = 5, shape = c("|", "x"), linetype = 0)))
  p2
  
  get_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  p2_legend <- get_legend(p1)
  
  ### Export
  pdf(file = paste0("../03_Output/", "Descriptives_Seq_", k, ".pdf"), width = 12, height = 7)
  par(mar = c(0, 0, 0, 0))
  par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
  grid.arrange(arrangeGrob(p2 + theme(legend.position="none"), 
                           p1 + theme(legend.position="none"), ncol = 2), 
               p2_legend, 
               nrow = 2, heights = c(10, 1),
               top = textGrob(names(vars)[which(vars == v)], gp = gpar(fontsize = 20)))
  dev.off()
  ### Export
  png(file = paste0("../03_Output/", "Descriptives_Seq_", k, ".png"), width = 12, height = 7,
      units = "in", res = 300)
  par(mar = c(0, 0, 0, 0))
  par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
  grid.arrange(arrangeGrob(p2 + theme(legend.position="none"), 
                           p1 + theme(legend.position="none"), ncol = 2), 
               p2_legend, 
               nrow = 2, heights = c(10, 1),
               top = textGrob(names(vars)[which(vars == v)], gp = gpar(fontsize = 20)))
  dev.off()
  
  
}  

















