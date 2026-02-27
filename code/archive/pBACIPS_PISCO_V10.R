#This code was written by Emily Donham to analyze kelp forest survey data using pBACI analysis

source("code/0_libraries.r")
rm(list = ls())

####################################################################################################
##Importing size frequency data to use for biomass estimates########################################
##Note that his size frequency dataset is used for resampling for all biomass estimates across programs
SizeFreq = read.csv(paste('data/ALL_sizefreq_2024.csv')) #Import PISCO/VRG/KFM/LTER size freq data

#Subset only taxa of interest, note that the codes differ across programs
SizeFreq.south <- subset(SizeFreq, campus == "VRG" & classcode == "MESFRA" | campus == "UCSB" & classcode == "MESFRA" |campus == "KFM" & classcode == "Strongylocentrotus franciscanus"
                         | campus == "VRG" & classcode == "STRPUR"| campus == "UCSB" & classcode == "STRPUR" | campus == "KFM" & classcode == "Strongylocentrotus purpuratus"
                         | campus == "UCSB" & classcode == "PANINT" | campus == "VRG" & classcode == "PANINT" 
                         | campus == "LTER" & classcode == "SFL" | campus == "LTER" & classcode == "SFS"
                         | campus == "LTER" & classcode == "SPS" | campus == "LTER" & classcode == "SPL")

#Remove KFM sites that we won't be using because 1) they are not a control or impact site OR 2) because they're only short time series and long time series exists
SizeFreq.south <- subset(SizeFreq.south, site != "a-k-21" & site != "a-k-05"  & site != "a-k-07" & site != "a-k-26"
                         & site != "a-k-27" & site != "a-k-28" & site != "a-k-29" & site != "a-k-30" 
                         & site != "a-k-36" & site != "a-k-37" & site != "a-k-21" & site != "a-k-12" & site != "a-k-13"
                         & site != "a-k-31" & site != "a-k-35" & site != "a-k-36" & site != "a-k-37")

#Remove sites that are not good references where ideal reference comparisons exist (based on conversation with Scott Hamilton)
SizeFreq.south <- subset(SizeFreq.south, site != "ANACAPA_BLACK_SEA_BASS" & site != "SMI_BAY_POINT" & site != "SCAI_SHIP_ROCK" 
                         & site != "SMI_CROOK_POINT_E" & site != "SMI_CROOK_POINT_W" & site != "SBI_SUTIL"& site != "SCI_YELLOWBANKS_CEN" & site != "SCI_YELLOWBANKS_W"
                         & site != "SRI_FORD_POINT" & site != "SMI_PRINCE_ISLAND_CEN" & site != "SMI_PRINCE_ISLAND_N" & site != "	SMI_HARE_ROCK"
                         & site != "SMI_TYLER_BIGHT_E" & site != "SMI_TYLER_BIGHT_W" & site != "SBI_GRAVEYARD_CANYON" & site != "SBI_GRAVEYARD_CANYON_N")

#Change name of KFM and LTER codes to be consistent with PISCO
SizeFreq <- mutate(SizeFreq.south, classcode = case_when(
  classcode == "Strongylocentrotus purpuratus"  ~ "STRPUR",
  classcode == "Strongylocentrotus franciscanus"    ~ "MESFRA",
  classcode == "SFL" ~ "MESFRA",
  classcode == "SFS" ~ "MESFRA",
  classcode == "SPS" ~ "STRPUR",
  classcode == "SPL" ~ "STRPUR",
  TRUE ~ classcode))



####################################################################################################

####################################################################################################
#######Import and wrangle PISCO swath data###########################################################
Swath = read.csv(paste('data/PISCO/MLPA_kelpforest_swath_2024.csv')) #Import PISCO/VRG swath
Swath.edit = subset(Swath, campus == "VRG" | campus == "UCSB" & zone != "DEEP" & zone != "SPECIALO") #Remove deep site, only surveyed by VRG. Also not always kelp forest and so aren't comparable/of interest

sites = read.csv(paste('data/PISCO/master_site_table_Emilyedit.csv')) #Edited table so that the PISCO references that had been changed during their decadal review run with this script (i.e. a couple of sites needed to be duplicated as they are controls for mult. MPAs)
sites.short <- sites %>%
  filter(!duplicated(cbind(site,CA_MPA_Name_Short))) #Remove duplicated Sites x MPA combos

##Transform count into stipes for macrocystis (size = stipes, count = # ind, so mult together to get total stipes)
Swath.edit$count <- ifelse(Swath.edit$classcode == "MACPYRAD",Swath.edit$count*Swath.edit$size, Swath.edit$count)

##Create list of all MPA*Status*Year*Month combos so we know if there was actually a survey to accurately say if there were zeros in the data
Site.yr <- unique(Swath.edit[,c("site","year","month","zone","transect","campus")])
Classcodes <- unique(Swath.edit$classcode) #unique species
Site.All.Classes <- merge(x = Site.yr, y = Classcodes, by = NULL) #merge to create full list of species for every survey 

##Merge our big "empty" dataset with the data to create dataframe with counts for ALL possible species
Swath.site.PISCO <- merge(Site.All.Classes, Swath.edit, by.x = c("site","year","month", "zone","transect","y","campus"), 
                              by.y = c("site","year","month","zone","transect","classcode","campus"),
                              all = TRUE)
Swath.site.PISCO[is.na(Swath.site.PISCO)] <- 0 #add zeros if no data

#Subset only summer since this is consistent across survey methods
Swath.site.PISCO.sub <- subset(Swath.site.PISCO, month <= 10 & month >= 5) 

####################################################################################################

#####################################################################################################
##Summarize Swath data to annual mean################################################################
##First, Sum the observations for species at each transect (this really just combines lobster and macro to 
#get a true sum since there are multiple observations due to sizes, all other species are a single row)
Swath.site.PISCO.sum.int <- Swath.site.PISCO.sub %>%
  group_by(year, month, site, zone, transect, y) %>% 
  summarise_at(c("count"), sum) %>%
  ungroup()

##Now take the average for each year since there are some years that multiple surveys were conducted
Swath.site.PISCO.sum <-  Swath.site.PISCO.sum.int %>%
  group_by(year, site, zone, transect, y) %>%
  summarize_at(c("count"), mean) 

##Ave transects across zones to have ave counts per zone
Swath.ave.zone <- Swath.site.PISCO.sum %>%
  group_by(year, site, y, zone) %>% 
  summarise_at(c("count"), mean)

##Now average the zones across site to have ave counts per site
Swath.ave <- Swath.ave.zone %>%
  group_by(year, site, y) %>% 
  summarise_at(c("count"), mean)

##Join average swath data to site table so we have mpas and other info
Swath.ave.site <- merge(Swath.ave, sites.short, by.x = c("site"), by.y = c("site"),
                        all.x = TRUE)

##Remove any sites that are not associated with an MPA according to PISCO table (some sites are duplicated since there are duplicate
#controls. This was just how I decided to work through the data as MPAs)
Swath.ave.site <- Swath.ave.site[complete.cases(Swath.ave.site$CA_MPA_Name_Short), ]

##Subset columns since we don't need all of them and this is easier when I want to look at things
Swath.ave.site.sub <- Swath.ave.site[, colnames(Swath.ave.site)[c(1:4,6,11:13,16)]]

Swath.ave.site.PISCO <- subset(Swath.ave.site.sub, BaselineRegion == "SOUTH") #Only sites from SoCal, not interested in Cen/North
Swath.ave.site.PISCO.subset <- subset(Swath.ave.site.PISCO, y == "STRPURAD" | y == "MESFRAAD" | y == "MACPYRAD" | y == "PANINT") #Only taxa of interest

Swath.ave.site.PISCO.subset$Density <- Swath.ave.site.PISCO.subset$count/60 #Convert to ind/m2 since transects are 60m2 (30*2)
Swath.ave.site.PISCO.subset$biomass <- NA
#####################################################################################################

#####################################################################################################
#################Wrangle size frequency data for lobsters############################################

#Subset only lobster data that has NOT been averaged by transect/zone so that we can also calculate biomass from size, since it's per individual in the data
Swath.site.PISCO.PANINT <- subset(Swath.site.PISCO.sub,   y == "PANINT")

PANINT <- merge(Swath.site.PISCO.PANINT, sites.short, by.x = c("site"), by.y = c("site"),
                     all.x = TRUE)

#Subset only PISCO data where we have sizes (2010 onward); VRG sizes are in a separate table and need to be analyzed differently
PANINT.UCSB <- subset(PANINT, year > 2009 & campus.x == "UCSB") #2010 is when PISCO started measuring carapace length

#Convert carapace length to biomass for lobsters
PANINT.UCSB$biomass <- 0.001352821 * ((PANINT.UCSB$size*10) ^ 2.913963113) * PANINT.UCSB$count #from LTER biomass relationships, input should be in mm, units are now grams

PANINT.UCSB.sub <- PANINT.UCSB[, colnames(PANINT.UCSB)[c(1:6,8,11:12,24:26,29,33)]] #Get rid of columns we don't need

PANINT.UCSB.null <- subset(PANINT.UCSB.sub, count > 0 & biomass == 0) #Find sites where lobster was found on transect but no size was recorded (only 6 observations)

PANINT.UCSB.sub.sub <- PANINT.UCSB.sub

#Remove sites where there is null data since we can't use the data from this year for biomass, wouldn't be accurate
for(i in 1:nrow(PANINT.UCSB.null)) {
  PANINT.UCSB.sub.sub <- subset(PANINT.UCSB.sub.sub, site != PANINT.UCSB.null$site[i] | year != PANINT.UCSB.null$year[i])
}

#Now sum across year, month, site, zone, transect so that we have a transect level total
PANINT.All.sub.zone.sum <- PANINT.UCSB.sub.sub %>%
  group_by(year, month, site, CA_MPA_Name_Short, site_designation, site_status, BaselineRegion, zone, transect, y) %>% 
  summarise_at(c("count", "biomass"), sum)

##Now average the transects across zones
PANINT.All.sub.zone.mean <- PANINT.All.sub.zone.sum %>%
  group_by(year, site, CA_MPA_Name_Short, site_designation, site_status, BaselineRegion, zone, y) %>% 
  summarise_at(c("count", "biomass"), mean)

#Now ave zones to have ave counts/biomass per site
PANINT.All.sub.mean <- PANINT.All.sub.zone.mean %>%
  group_by(year, site, CA_MPA_Name_Short, site_designation, site_status, BaselineRegion, y) %>% 
  summarise_at(c("count", "biomass"), mean)

PANINT.All.sub.mean$Density <- PANINT.All.sub.mean$count/60
PANINT.All.sub.mean$y <- "PANINT"
PANINT.All.sub.mean$campus <- "UCSB"

#####################################################################################################
#Now we need to use size frequency data to get sizes of lobsters from VRG surveys####################
PANINT.VRG <- subset(PANINT, year > 2010 & campus.x == "VRG") #2011 is when Vantuna Research Group started measuring carapace length on transects and roving divers, but size was not attached to the transect counts

#Merge size frequency table with site table and subset only lobsters
SizeFreq.PANINT.VRG <- merge(SizeFreq, sites.short, by.x = c("site", "site_new", "campus"), by.y = c("site", "site_new", "campus"),
                       all.x = TRUE)
SizeFreq.PANINT.VRG <- subset(SizeFreq.PANINT.VRG, classcode == "PANINT")

#Create dataframe with all unique zones across years
u <- unique(PANINT.VRG[,c('year','site','CA_MPA_Name_Short','site_designation','site_status','BaselineRegion','zone')])

#Create empty dataframes
VRG.PANINT.site <- data.frame(site = NA, year = NA, zone = NA, transect = NA, biomass = NA, count = NA)
VRG.PANINT.SF <- data.frame(site = NA, year = NA, zone = NA, transect = NA, biomass = NA, count = NA)

#Create function to calculate biomass of lobsters equation from LTER
bio_lobster <- function(size) {
  0.001352821 * (size) ^ 2.913963113
}

#Select site/year/zone if there are lobsters then look up same site/year/zone from size freq table since this is where all size data are.
#Create lobster "population" from size freq table and draw with replacement the number of lobsters from swath data table if lobsters were counted on transect
#However, if no lobsters then equals zero, if there were lobsters but no size freq data exists, these will be marked as NA and removed later since they can't be included
for(i in 1:nrow(u)) {
t <- which(PANINT.VRG$site == u$site[i] & PANINT.VRG$year == u$year[i] & PANINT.VRG$zone == u$zone[i])
n <- sum(PANINT.VRG$count[t])
t2 <- which(SizeFreq.PANINT.VRG$CA_MPA_Name_Short == u$CA_MPA_Name_Short[i] & SizeFreq.PANINT.VRG$year == u$year[i] &
          SizeFreq.PANINT.VRG$site_status == u$site_status[i])

if (n != 0 & length(t2) != 0) { #This means a lobster was observed on the transect and in the size freq table
  a = NA
  s = data.frame(matrix(NA, nrow = 1000, ncol = n)) #resampling dataframe
  l <- length(t2) #number of lobsters to resample
  for(j in 1:l){
      list <- as.data.frame(rep(SizeFreq.PANINT.VRG$size[t2[j]], SizeFreq.PANINT.VRG$count[t2[j]] ))
      a <- rbind(a,list)
      a <- na.omit(a) }
  a <- as.numeric(as.character(unlist(a))) #change list to numeric so that it can be sampled from, this is each lobster size that we saw in the transects and roving diver surveys
  #randomly sample the lobster "population" with replacement; do this 1000 times
  for(k in 1:1000){
  s[k,] <- sample(a,n, replace = TRUE) 
  }
  #Convert CL to biomass for all individuals
  s_bio <- s %>% 
    mutate(across(starts_with("X"), ~bio_lobster(.)))
  s_bio <- s_bio %>% 
    mutate(sum = rowSums(.))
  ave <- mean(s_bio$sum)
  VRG.PANINT.SF$site <- u$site[i]
  VRG.PANINT.SF$year <- u$year[i]
  VRG.PANINT.SF$zone <- u$zone[i]
  VRG.PANINT.SF$transect <- length(t)
  VRG.PANINT.SF$biomass <- ave
  VRG.PANINT.SF$count <- n
  VRG.PANINT.site <- rbind(VRG.PANINT.site,VRG.PANINT.SF)
  
} else if (n == 0 & length(t2) == 0){#This means there were no lobsters on transect and none in the size freq data either
  VRG.PANINT.SF$site <- u$site[i]
  VRG.PANINT.SF$year <- u$year[i]
  VRG.PANINT.SF$zone <- u$zone[i]
  VRG.PANINT.SF$transect <- length(t)
  VRG.PANINT.SF$biomass <- 0
  VRG.PANINT.SF$count <- 0
  VRG.PANINT.site <- rbind(VRG.PANINT.site,VRG.PANINT.SF)
  
} else if (n != 0 & length(t2) == 0){ #If there are lobsters on transects, but they weren't sized we will mark size as NA
  VRG.PANINT.SF$site <- u$site[i]
  VRG.PANINT.SF$year <- u$year[i]
  VRG.PANINT.SF$zone <- u$zone[i]
  VRG.PANINT.SF$transect <- length(t)
  VRG.PANINT.SF$biomass <- NA
  VRG.PANINT.SF$count <- n
  VRG.PANINT.site <- rbind(VRG.PANINT.site,VRG.PANINT.SF)
  
} else {
  VRG.PANINT.SF$site <- u$site[i]
  VRG.PANINT.SF$year <- u$year[i]
  VRG.PANINT.SF$zone <- u$zone[i]
  VRG.PANINT.SF$transect <- length(t)
  VRG.PANINT.SF$biomass <- 0
  VRG.PANINT.SF$count <- 0
  VRG.PANINT.site <- rbind(VRG.PANINT.site,VRG.PANINT.SF)
}}

#Merge sampling data with site table 
VRG.PANINT.site.merge <- merge(VRG.PANINT.site, sites.short, by.x = c("site"), by.y = c("site"),
                all.x = TRUE)

VRG.PANINT.site.merge <- VRG.PANINT.site.merge[, colnames(VRG.PANINT.site.merge)[c(2,1,13:15,18,3:6)]]
VRG.PANINT.site.merge[is.na(VRG.PANINT.site.merge)] <- 0

#Find sites where lobster was found on transect but no size was recorded
PANINT.VRG.null <- subset(VRG.PANINT.site.merge, count > 0 & biomass == 0) 

VRG.PANINT.site.merge.sub <- VRG.PANINT.site.merge

#Remove sites where lobsters were seen but not measured since we can't get an accurate biomass estimate
for(i in 1:nrow(PANINT.VRG.null)) {
  VRG.PANINT.site.merge.sub <- subset(VRG.PANINT.site.merge.sub, site != PANINT.VRG.null$site[i] | year != PANINT.VRG.null$year[i])
}

VRG.PANINT.site.all.TransectSum <- VRG.PANINT.site.merge.sub %>% #summing across transects
  group_by(year, site, CA_MPA_Name_Short, site_designation, site_status, BaselineRegion, zone, transect) %>% 
  summarise_at(c("count", "biomass"), sum, na.rm = TRUE)

#Need to account for the different number of transects at the sites, since they sometimes vary
VRG.PANINT.site.all.TransectSum$count <- VRG.PANINT.site.all.TransectSum$count/VRG.PANINT.site.all.TransectSum$transect
VRG.PANINT.site.all.TransectSum$biomass <- VRG.PANINT.site.all.TransectSum$biomass/VRG.PANINT.site.all.TransectSum$transect

VRG.PANINT.site.all <- VRG.PANINT.site.all.TransectSum %>% #averaging across zones
  group_by(year, site, CA_MPA_Name_Short, site_designation, site_status, BaselineRegion) %>% 
  summarise_at(c("count", "biomass"), mean)

VRG.PANINT.site.all$Density <- VRG.PANINT.site.all$count/60
VRG.PANINT.site.all$y <- "PANINT"

VRG.PANINT <- VRG.PANINT.site.all

VRG.PANINT <- VRG.PANINT[, colnames(VRG.PANINT)[c(1:6,10,7:9)]]
VRG.PANINT$campus <- "VRG"
#####################################################################################################

#####Merge lobster datasets##########################################################################
PANINT.all <- rbind(VRG.PANINT, PANINT.All.sub.mean)
PANINT.all <- PANINT.all[, colnames(PANINT.all)[c(2,1,7,8,11,3:6,10,9)]]

Swath.PISCO.subset <- rbind(PANINT.all, Swath.ave.site.PISCO.subset) 

Swath.PISCO <- Swath.PISCO.subset[!is.na(Swath.PISCO.subset$count),] #Delete any site years where there were transects with lobsters that weren't measured, so we just can't get sizes
###################################################################################################################
#Now we need to use size frequency data to get sizes/biomass of urchins###################################
Swath.Urchin.Bio <- subset(Swath.site.PISCO.sum.int, y == "STRPURAD" | y == "MESFRAAD") 
URCHINS <- merge(Swath.Urchin.Bio, sites.short, by.x = c("site"), by.y = c("site"),
                all.x = TRUE)
SizeFreq.Urch <- subset(SizeFreq, classcode == "MESFRA" | classcode == "STRPUR") #Subset only urchins
SizeFreq.Urch <- mutate(SizeFreq.Urch, classcode = case_when(
  classcode == "STRPUR"  ~ "STRPURAD",
  classcode == "MESFRA"    ~ "MESFRAAD",
  TRUE ~ classcode))

SizeFreq.Urch <- merge(SizeFreq.Urch, sites.short, by.x = c("site", "site_new", "campus"), by.y = c("site", "site_new", "campus"),
                 all.x = TRUE)

SizeFreq.Urch.OG <- SizeFreq.Urch

#Only subset urchins 25mm or larger since PISCO only counts urchins this size or larger
SizeFreq.Urch <- subset(SizeFreq.Urch, size >= 25)

#Create function to calculate biomass of urchins for red
bio_redurch <- function(size) {
  0.00059 * (size) ^ 2.917
}
#and purps
bio_purpurch <- function(size) {
  0.00059 * (size) ^ 2.870
}

#Find unique zones since we'll average the transects to zone
u <- unique(URCHINS[,c('year','site','site_new','CA_MPA_Name_Short','site_designation','site_status','BaselineRegion','zone','y')])

#Create empty dataframes
Urchin.site <- data.frame(site = NA, CA_MPA_Name_Short = NA, site_status = NA, year = NA, 
                          zone = NA, transect = NA, count = NA, biomass = NA, y = NA, site_new = NA)
Urchin.SF <- data.frame(site = NA, CA_MPA_Name_Short = NA, site_status = NA, year = NA, 
                          zone = NA, transect = NA, count = NA, biomass = NA, y = NA, site_new = NA)

#Select site/year/zone/species if there are urchins then look up same site/year/zone from size freq table since this is where 
#all size data are create urchin "population" from size freq table and draw with replacement the number of urchins 
#from swath data table if urchins were measured on transect; however, if no urchins then equals zero, if there were 
#urchins but no size freq these will be removed, so marked as NA

for(i in 1:nrow(u)) {
  t <- which(URCHINS$site == u$site[i] & URCHINS$CA_MPA_Name_Short == u$CA_MPA_Name_Short[i] & URCHINS$year == u$year[i] 
             & URCHINS$zone == u$zone[i] & URCHINS$y == u$y[i])
  n <- sum(URCHINS$count[t])
  t2 <- which(SizeFreq.Urch$CA_MPA_Name_Short == u$CA_MPA_Name_Short[i] & SizeFreq.Urch$site_status == u$site_status[i] &
                SizeFreq.Urch$year == u$year[i] & SizeFreq.Urch$classcode == u$y[i])

  if (n != 0 & length(t2) != 0) { #This means a lobster was observed on the transect and in the size freq table
    a = NA
    s = data.frame(matrix(NA, nrow = 1000, ncol = n))
    l <- length(t2)
    for(j in 1:l){
      list <- as.data.frame(rep(SizeFreq.Urch$size[t2[j]], SizeFreq.Urch$count[t2[j]] ))
      a <- rbind(a,list)
      a <- na.omit(a) }
      a <- as.numeric(as.character(unlist(a))) #change list to numeric so that it can be sampled from, this is each urchin size 
      #randomly sample the urchin "population" with replacement; do this 1000 times
    for(k in 1:1000){
      s[k,] <- sample(a,n, replace = TRUE) 
    }
    #Convert CL to biomass for all individuals
    if(u$y[i] == "MESFRAAD"){
       s_bio <- s %>% 
         mutate(across(starts_with("X"), ~bio_redurch(.)))
       s_bio <- s_bio %>% 
         mutate(sum = rowSums(.))
      ave <- mean(s_bio$sum)
      }else {
       s_bio <- s %>% 
         mutate(across(starts_with("X"), ~bio_purpurch(.)))
       s_bio <- s_bio %>% 
        mutate(sum = rowSums(.))
       ave <- mean(s_bio$sum)
       }
    
    Urchin.SF$site <- u$site[i]
    Urchin.SF$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
    Urchin.SF$site_status <- u$site_status[i]
    Urchin.SF$year <- u$year[i]
    Urchin.SF$zone <- u$zone[i]
    Urchin.SF$transect <- length(t)
    Urchin.SF$count <- n
    Urchin.SF$biomass <- ave
    Urchin.SF$y <- u$y[i]
    Urchin.SF$site_new <- u$site_new[i]
    Urchin.site <- rbind(Urchin.site, Urchin.SF)
    
  } else if (n == 0 & length(t2) == 0){#This means there were no lobsters on transect and none in the size freq data either
    Urchin.SF$site <- u$site[i]
    Urchin.SF$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
    Urchin.SF$site_status <- u$site_status[i]
    Urchin.SF$year <- u$year[i]
    Urchin.SF$zone <- u$zone[i]
    Urchin.SF$transect <- length(t)
    Urchin.SF$count <- 0
    Urchin.SF$biomass <- 0
    Urchin.SF$y <- u$y[i]
    Urchin.SF$site_new <- u$site_new[i]
    Urchin.site <- rbind(Urchin.site,Urchin.SF)
    
  } else if (n != 0 & length(t2) == 0){ #If there are lobsters on transects, but they weren't sized we will mark size as NA
    Urchin.SF$site <- u$site[i]
    Urchin.SF$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
    Urchin.SF$site_status <- u$site_status[i]
    Urchin.SF$year <- u$year[i]
    Urchin.SF$zone <- u$zone[i]
    Urchin.SF$transect <- length(t)
    Urchin.SF$count <- n
    Urchin.SF$biomass <- NA
    Urchin.SF$y <- u$y[i]
    Urchin.SF$site_new <- u$site_new[i]
    Urchin.site <- rbind(Urchin.site,Urchin.SF)
    
  } else {
    Urchin.SF$site <- u$site[i]
    Urchin.SF$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
    Urchin.SF$site_status <- u$site_status[i]
    Urchin.SF$year <- u$year[i]
    Urchin.SF$zone <- u$zone[i]
    Urchin.SF$transect <- length(t)
    Urchin.SF$count <- 0
    Urchin.SF$biomass <- 0
    Urchin.SF$y <- u$y[i]
    Urchin.SF$site_new <- u$site_new[i]
    Urchin.site <- rbind(Urchin.site,Urchin.SF)
  }}

PISCO.Urchin.site.merge <- merge(Urchin.site, sites.short, by.x = c("site", "CA_MPA_Name_Short", "site_new"), by.y = c("site", "CA_MPA_Name_Short", "site_new"),
                               all.x = TRUE)

PISCO.Urchin.site.merge <- PISCO.Urchin.site.merge[, colnames(PISCO.Urchin.site.merge)[c(2,1,3,17,5:10,16,20)]]

PISCO.Urchin.site.merge[is.na(PISCO.Urchin.site.merge)] <- 0

#Find sites where urchin was found on transect but no size was recorded
PISCO.Urchin.null <- subset(PISCO.Urchin.site.merge, count > 0 & biomass == 0) 

PISCO.Urchin.site.merge.sub <- PISCO.Urchin.site.merge

#Remove sites where urchins were seen but not measured since we can't get an accurate biomass estimate
for(i in 1:nrow(PISCO.Urchin.null)) {
  PISCO.Urchin.site.merge.sub <- subset(PISCO.Urchin.site.merge.sub, site != PISCO.Urchin.null$site[i] | year != PISCO.Urchin.null$year[i])
}

names(PISCO.Urchin.site.merge.sub)[names(PISCO.Urchin.site.merge.sub) == "site_status.y"] <- "site_status"

PISCO.Urchin.site.transectSum <- PISCO.Urchin.site.merge.sub %>% #summing across transects
  group_by(year, site, CA_MPA_Name_Short, site_designation, site_status, BaselineRegion, zone, transect, y) %>% 
  summarise_at(c("count", "biomass"), sum, na.rm = TRUE)

#Need to account for the different number of transects at the sites
PISCO.Urchin.site.transectSum$count <- PISCO.Urchin.site.transectSum$count/PISCO.Urchin.site.transectSum$transect
PISCO.Urchin.site.transectSum$biomass <- PISCO.Urchin.site.transectSum$biomass/PISCO.Urchin.site.transectSum$transect

PISCO.Urchin.site.all <- PISCO.Urchin.site.transectSum %>% #averaging across zones
  group_by(year, site, CA_MPA_Name_Short, site_designation, site_status, BaselineRegion, y) %>% 
  summarise_at(c("count", "biomass"), mean)

PISCO.Urchin.site.all$Density <- PISCO.Urchin.site.all$count/60

Swath.PISCO.allbio <- merge(Swath.PISCO, PISCO.Urchin.site.all, by.x = c("site", "year", "y", "CA_MPA_Name_Short"), 
              by.y = c("site", "year", "y", "CA_MPA_Name_Short"),
              all.x = TRUE )

Swath.PISCO.allbio$biomass.x <- ifelse(is.na(Swath.PISCO.allbio$biomass.x) == TRUE,  
                                          Swath.PISCO.allbio$biomass.y, Swath.PISCO.allbio$biomass.x)


Swath.PISCO <- Swath.PISCO.allbio
Swath.PISCO <- Swath.PISCO[, colnames(Swath.PISCO)[c(1:11,15)]]
names(Swath.PISCO)[names(Swath.PISCO) == "biomass.x"] <- "biomass"
names(Swath.PISCO)[names(Swath.PISCO) == "count.x"] <- "count"
names(Swath.PISCO)[names(Swath.PISCO) == "Density.x"] <- "Density"
names(Swath.PISCO)[names(Swath.PISCO) == "site_designation.x"] <- "site_designation"
names(Swath.PISCO)[names(Swath.PISCO) == "site_status.x"] <- "site_status"
names(Swath.PISCO)[names(Swath.PISCO) == "BaselineRegion.x"] <- "BaselineRegion"

Purp.sub <- subset(Swath.PISCO, y == "STRPURAD")
Purp.sub$diff <- Purp.sub$count - Purp.sub$count.y

#####################################################################################################
###############################Now Import and wrangle PISCO fish data################################
FishAtt = read.csv(paste('data/PISCO/spp_attribute_table_downloaded_9-13-22_SHSPUL.csv')) # Attribute table that has parameters to convert size to biomass, using Scott Hamiltons parameters for sheephead
Fish = read.csv(paste('data/PISCO/MLPA_kelpforest_fish_2024.csv')) # PISCO fish database
Fish = subset(Fish, campus == "UCSB" | campus == "VRG")

#Merge the attributes to the fish
Fish.merged <- merge(Fish, FishAtt, by.x = "classcode", by.y = "pisco_classcode",
                        all.x = TRUE)

Fish.merged <- Fish.merged[complete.cases(Fish.merged$WL_input_length), ]

#Calculate biomass from length. Some are in TL and some are in SL, so need to flip equation depending on SL vs TL
l = nrow(Fish.merged)
for(i in 1:l) {
if (Fish.merged$WL_input_length[i] == "TL") { 
  Fish.merged$biomass[i] = Fish.merged$count[i]*(Fish.merged$WL_a[i]*(Fish.merged$fish_tl[i]^Fish.merged$WL_b[i])) #bio = N * (a * TL^b)
} else if (Fish.merged$WL_input_length[i] == "SL" & Fish.merged$LC_type_for_WL[i] == "TYPICAL") {
  Fish.merged$biomass[i] = Fish.merged$count[i]*(Fish.merged$WL_a[i]*
                           ((Fish.merged$fish_tl[i]*Fish.merged$LC.a._for_WL[i]+Fish.merged$LC.b._for_WL[i])^Fish.merged$WL_b[i])) #bio = N * a * (TL * LCa + (LCb))^b
} else if (Fish.merged$WL_input_length[i] == "SL" & Fish.merged$LC_type_for_WL[i] == "REVERSE") {
  Fish.merged$biomass[i] = Fish.merged$count[i]*(Fish.merged$WL_a[i]*
                            (((Fish.merged$fish_tl[i]-Fish.merged$LC.b._for_WL[i])/Fish.merged$LC.a._for_WL[i])^Fish.merged$WL_b[i]))#bio = N * a * ((TL - LCb)/LCa)^b
} else if (Fish.merged$WL_input_length[i] == "DW") {
    Fish.merged$biomass[i] = NA
} else  {
    Fish.merged$biomass[i] = NA
}}

#Remove excess columns we don't need
Fish.sub <- select(Fish.merged, classcode, year, month, site, zone, level, transect, count, fish_tl, TrophicGroup, BroadTrophic, Targeted, biomass)

##Join fish data to site table 
Fish.sub.site <- merge(Fish.sub, sites.short, by.x = c("site"), by.y = c("site"),
                        all.x = TRUE)

##Subset columns since we don't need all of them and this is easier to look at
Fish.sub.site <- Fish.sub.site[, colnames(Fish.sub.site)[c(1:13,20:22,25)]]

#Remove canopy since this is irrelevant to sheephead
Fish.sub.site.NOcanopy <-  subset(Fish.sub.site, level == "BOT")
Fish.sub.site.NOcanopy <-  subset(Fish.sub.site, BaselineRegion == "SOUTH") #Pull out just study region

#Sum biomass/counts across levels for each mpa/reference site
Fish.sum <- Fish.sub.site.NOcanopy %>%
  group_by(CA_MPA_Name_Short, site, year, zone, transect, classcode, site_status, site_designation) %>% 
  summarise_at(c("biomass","count"), sum)


##Create list of all MPA*Status*Year*Month combos so we know if there was actually a survey to accurately say if there were zeros in the data
Site.yr <- unique(Fish.sum[,c("site","year","zone","transect","CA_MPA_Name_Short","site_status", "site_designation")])
Classcodes <- unique(Fish.sum[,c("classcode")])
Site.All.Classes <- merge(x = Site.yr, y = Classcodes, by = NULL)

##Merge our big "empty" dataset with the data to create dataframe with counts for ALL possible species
PISCO.fish <- merge(Site.All.Classes, Fish.sum, by.x = c("site","year","zone","transect","CA_MPA_Name_Short","site_status","classcode","site_designation"), 
                    by.y = c("site","year","zone","transect","CA_MPA_Name_Short","site_status","classcode","site_designation"),
                    all = TRUE)

##Replace NA with zero
PISCO.fish[is.na(PISCO.fish)] <- 0

#Calculate mean biomass and density over time at site level for each mpa/reference site
Fish.mean <- PISCO.fish %>%
  group_by(CA_MPA_Name_Short, site, year, zone, classcode, site_status, site_designation) %>% 
  summarise_at(c("biomass","count"), mean) %>%
  ungroup()

Fish.mean.final <- Fish.mean %>%
  group_by(CA_MPA_Name_Short, site, year, classcode, site_status, site_designation) %>% 
  summarise_at(c("biomass","count"), mean) %>%
  ungroup()

##########################################################################################################################
#######Subset dataset for sheephead#######################################################################################

Fish.SPUL.All <- subset(Fish.mean.final, classcode == "SPUL") #Subset averaged data where zeroes have been added
Fish.SPUL.All$density <- Fish.SPUL.All$count/60
colnames(Fish.SPUL.All)[4] <- "y"

##################Combine all data from PISCO surveys######################################################################
Swath.PISCO.sub <- Swath.PISCO[, colnames(Swath.PISCO)[c(4,1,2,3,8,7,11,5)]]
Swath.PISCO.sub$density <- Swath.PISCO.sub$count/60

All_PISCO <- rbind(Swath.PISCO.sub, Fish.SPUL.All)

x <- unique(All_PISCO$site)

l <- length(x)
All_PISCO.short <- All_PISCO
for(i in 1:l) {
  idx <- which(All_PISCO$site == x[i] )
  x2 <- unique(All_PISCO[idx,c("site","year")])
  l2 <- nrow(x2)
  if (l2 < 5) { #Must have 5 or more years worth of data
    All_PISCO.short <- subset(All_PISCO.short, All_PISCO.short$site != x[i])
  } else  {
  }}

###################Remove sites based on being good references in channel islands##########################################
All_PISCO.short.sub <- All_PISCO.short
#Remove sites that are not good references where ideal reference comparisons exist (based on conversation with Scott Hamilton)
All_PISCO.short.sub <- subset(All_PISCO.short.sub, site != "ANACAPA_BLACK_SEA_BASS" & site != "SMI_BAY_POINT" & site != "SCAI_SHIP_ROCK" 
                              & site != "SMI_CROOK_POINT_E" & site != "SMI_CROOK_POINT_W" & site != "SBI_SUTIL"& site != "SCI_YELLOWBANKS_CEN" & site != "SCI_YELLOWBANKS_W"
                              & site != "SRI_FORD_POINT" & site != "SMI_PRINCE_ISLAND_CEN" & site != "SMI_PRINCE_ISLAND_N" & site != "	SMI_HARE_ROCK"
                  & site != "SMI_TYLER_BIGHT_E" & site != "SMI_TYLER_BIGHT_W" & site != "SBI_GRAVEYARD_CANYON" & site != "SBI_GRAVEYARD_CANYON_N")

#Remove mpas that don't meet criteria for inclusion
All_PISCO.short.sub <- subset(All_PISCO.short.sub,  CA_MPA_Name_Short != "Carrington Point SMR" & CA_MPA_Name_Short != "N/A"
                              & CA_MPA_Name_Short != "Arrow Point to Lion Head Point SMCA" & CA_MPA_Name_Short != "Crystal Cove SMCA" &
                                CA_MPA_Name_Short != "Laguna Beach SMR" & CA_MPA_Name_Short != "South La Jolla SMR" &
                                CA_MPA_Name_Short != "Vandenberg SMR" & CA_MPA_Name_Short != "Point Conception SMR" &
                                CA_MPA_Name_Short != "Anacapa Island SMR 1978" & CA_MPA_Name_Short != "Painted Cave SMCA" &
                                CA_MPA_Name_Short != "Anacapa Island SMCA")

Site = read.csv(paste('data/Site_List_All.csv'))
All_PISCO.short.sub$time <- 0
x <- unique(All_PISCO.short.sub$CA_MPA_Name_Short) # to redo since sites have been removed
l = length(x)
for(i in 1:l) {
  idx <-which(x[i] == All_PISCO.short.sub$CA_MPA_Name_Short)
  j <-which(x[i] == Site$Site)
  m <- length(idx)
  if (Site$MPA_Start[j] == 2003) {
    for (n in 1:m) {
      if (All_PISCO.short.sub$year[idx[n]] <=2003) {
        All_PISCO.short.sub$time[idx[n]] = 0
      } else if (All_PISCO.short.sub$year[idx[n]] == 2004) {
        All_PISCO.short.sub$time[idx[n]] = 1
      } else if (All_PISCO.short.sub$year[idx[n]] == 2005) {
        All_PISCO.short.sub$time[idx[n]] = 2
      } else if (All_PISCO.short.sub$year[idx[n]] == 2006) {
        All_PISCO.short.sub$time[idx[n]] = 3
      } else if (All_PISCO.short.sub$year[idx[n]] == 2007) {
        All_PISCO.short.sub$time[idx[n]] = 4
      } else if (All_PISCO.short.sub$year[idx[n]] == 2008) {
        All_PISCO.short.sub$time[idx[n]] = 5
      } else if (All_PISCO.short.sub$year[idx[n]] == 2009) {
        All_PISCO.short.sub$time[idx[n]] = 6
      } else if (All_PISCO.short.sub$year[idx[n]] == 2010) {
        All_PISCO.short.sub$time[idx[n]] = 7
      } else if (All_PISCO.short.sub$year[idx[n]] == 2011) {
        All_PISCO.short.sub$time[idx[n]] = 8
      } else if (All_PISCO.short.sub$year[idx[n]] == 2012) {
        All_PISCO.short.sub$time[idx[n]] = 9
      } else if (All_PISCO.short.sub$year[idx[n]] == 2013) {
        All_PISCO.short.sub$time[idx[n]] = 10
      } else if (All_PISCO.short.sub$year[idx[n]] == 2014) {
        All_PISCO.short.sub$time[idx[n]] = 11
      } else if (All_PISCO.short.sub$year[idx[n]] == 2015) {
        All_PISCO.short.sub$time[idx[n]] = 12
      } else if (All_PISCO.short.sub$year[idx[n]] == 2016) {
        All_PISCO.short.sub$time[idx[n]] = 13
      } else if (All_PISCO.short.sub$year[idx[n]] == 2017) {
        All_PISCO.short.sub$time[idx[n]] = 14
      } else if (All_PISCO.short.sub$year[idx[n]] == 2018) {
        All_PISCO.short.sub$time[idx[n]] = 15
      } else if (All_PISCO.short.sub$year[idx[n]] == 2019) {
        All_PISCO.short.sub$time[idx[n]] = 16
      } else if (All_PISCO.short.sub$year[idx[n]] == 2020) {
        All_PISCO.short.sub$time[idx[n]] = 17
      } else if (All_PISCO.short.sub$year[idx[n]] == 2021) {
        All_PISCO.short.sub$time[idx[n]] = 18
      } else if (All_PISCO.short.sub$year[idx[n]] == 2022) {
        All_PISCO.short.sub$time[idx[n]] = 19
      } else if (All_PISCO.short.sub$year[idx[n]] == 2023) {
        All_PISCO.short.sub$time[idx[n]] = 20
      }}
    
  } else if (Site$MPA_Start[j] == 2012) {
    for (n in 1:m) {
      if (All_PISCO.short.sub$year[idx[n]] <= 2012) {
        All_PISCO.short.sub$time[idx[n]] = 0
      } else if (All_PISCO.short.sub$year[idx[n]] == 2013) {
        All_PISCO.short.sub$time[idx[n]] = 1
      } else if (All_PISCO.short.sub$year[idx[n]] == 2014) {
        All_PISCO.short.sub$time[idx[n]] = 2
      } else if (All_PISCO.short.sub$year[idx[n]] == 2015) {
        All_PISCO.short.sub$time[idx[n]] = 3
      } else if (All_PISCO.short.sub$year[idx[n]] == 2016) {
        All_PISCO.short.sub$time[idx[n]] = 4
      } else if (All_PISCO.short.sub$year[idx[n]] == 2017) {
        All_PISCO.short.sub$time[idx[n]] = 5
      } else if (All_PISCO.short.sub$year[idx[n]] == 2018) {
        All_PISCO.short.sub$time[idx[n]] = 6
      } else if (All_PISCO.short.sub$year[idx[n]] == 2019) {
        All_PISCO.short.sub$time[idx[n]] = 7
      } else if (All_PISCO.short.sub$year[idx[n]] == 2020) {
        All_PISCO.short.sub$time[idx[n]] = 8
      } else if (All_PISCO.short.sub$year[idx[n]] == 2021) {
        All_PISCO.short.sub$time[idx[n]] = 9
      } else if (All_PISCO.short.sub$year[idx[n]] == 2022) {
        All_PISCO.short.sub$time[idx[n]] = 10
      } else if (All_PISCO.short.sub$year[idx[n]] == 2023) {
        All_PISCO.short.sub$time[idx[n]] = 11
      }}}
  else {}
}

################################################################################################################
#Calculating biomass of macrocystis#############################################################################
################################################################################################################
#We're going to calculate biomass of macro based on equations from LTER
aveslope <- mean(c(0.10386, 0.10103, 0.09267,	0.09204,0.08054, 0.08505)) #Slopes from May-Oct from LTER
All_PISCO.short.sub$biomass <- ifelse(All_PISCO.short.sub$y == "MACPYRAD",All_PISCO.short.sub$density*aveslope*1000, All_PISCO.short.sub$biomass)

#################################################################################################################
#Average across mpa and reference sites so we have one value per MPA#############################################
#################################################################################################################
All_PISCO.mpa <- All_PISCO.short.sub %>%
  group_by(CA_MPA_Name_Short, year, y, site_status, site_designation) %>% 
  summarise_at(c("biomass", "count", "density"), mean, na.rm = TRUE) %>%
  ungroup()
#################################################################################################################

#################################################################################################################
##################Create short format for log response ratio calculations########################################
#################################################################################################################
All.bio <- select(All_PISCO.mpa, CA_MPA_Name_Short, year, y, site_status, site_designation, biomass)
All.bio <- subset(All.bio, site_status != "")
All.bio <- subset(All.bio, site_designation != "N/A")
All.bio <- subset(All.bio, biomass != "N/A") #removing NA since there is no data from that year
All.bio <- subset(All.bio, biomass != "NaN") #removing NA since there is no data from that year

#####Start with converting data to proportions
All.bio$Prop <- NA
All.bio$PropCorr <- NA
x = unique(All.bio$CA_MPA_Name_Short)
t = unique(All.bio$y)
l = length(x)
p = length(t)

#Turn maximum annual biomass into proportion of time series maximum biomass
for(i in 1:l) {
  for(j in 1:p) {
    idx <-which(x[i] == All.bio$CA_MPA_Name_Short & t[j] == All.bio$y)
    m <- max(All.bio$biomass[idx]) #Find the max of the mpa specific time series to use to create a proportion of the max
    All.bio$Prop[idx] <- All.bio$biomass[idx]/m #Calculate proportion of max for every data point in the time series
    All.bio$PropCorr[idx] <- All.bio$Prop[idx] + 0.01 #Add 1% so that there are no zeroes in the dataset
}}


# Create short form to easily line up and know whether we have data in both mpa and reference for a given year and calc diff
All.bio.sub <- All.bio[, colnames(All.bio)[c(1:5,8)]] #First subset only columns we need
Short.bio <- All.bio.sub %>% #Create short format
  spread(site_status, PropCorr)
Short.bio.diff <- Short.bio[complete.cases(Short.bio), ] #Remove data if there are not both mpa and reference
Short.bio.diff$Diff <- Short.bio.diff$mpa/Short.bio.diff$reference #Calculate the response ratio
Short.bio.diff$lnDiff <- log(Short.bio.diff$Diff) #Log response ratio

All.bio.diff.final <- Short.bio.diff
All.bio.diff.final$resp <- "Bio"


#Now do the same for density
All.den <- select(All_PISCO.mpa, CA_MPA_Name_Short, year, y, site_status, site_designation, density)
All.den <- subset(All.den, site_status != "")
All.den$Prop <- 0
All.den$PropCorr <- 0
x = unique(All.den$CA_MPA_Name_Short)
t = unique(All.den$y)
l = length(x)
p = length(t)

#Turn maximum annual density into proportion of time series maximum density
for(i in 1:l) {
  for(j in 1:p) {
  idx <-which(x[i] == All.den$CA_MPA_Name_Short & t[j] == All.den$y)
  m <- max(All.den$density[idx]) 
  All.den$Prop[idx] <- All.den$density[idx]/m
  All.den$PropCorr[idx] <- All.den$Prop[idx] + 0.01 
}}

All.den.sub <- All.den[, colnames(All.den.SPUL)[c(1:5,8)]] #subset only columns needed
Short.den <- All.den.sub %>% #short format to calculate diff
  spread(site_status, PropCorr)
Short.den.diff <- Short.den[complete.cases(Short.den), ] #remove any site/years where there is not data in both mpa and reference
Short.den.diff$Diff <- Short.den.diff$mpa/Short.den.diff$reference #Calculate the response ratio
Short.den.diff$lnDiff <- log(Short.den.diff$Diff) #Log response ratio

All.den.diff.final <- Short.den.diff
All.den.diff.final$resp <- "Den"



#################################################################################################################
#Create dataframe with all raw average density/biomass inside/outside mpas#######################################

All_PISCO.edit <- All_PISCO.mpa[, colnames(All_PISCO.mpa)[c(4,1,2,8,3)]]
All_PISCO.edit <- All_PISCO.edit[complete.cases(All_PISCO.edit), ] #Delete any rows where there isn't an mpa and reference
All_PISCO.edit <- subset(All_PISCO.edit, site_status != "")
PISCO.den <- All_PISCO.edit %>% 
  spread(site_status, density)
PISCO.den$source <- "PISCO"
PISCO.den <- PISCO.den[complete.cases(PISCO.den), ] #Delete any rows where there isn't an mpa and reference
PISCO.den.long <- gather(PISCO.den, status, value, 4:5)
PISCO.den.long$resp <- 'Den'

All_PISCO.mpa$Bio <- ifelse(All_PISCO.mpa$y == "MACPYRAD", All_PISCO.mpa$biomass, All_PISCO.mpa$biomass/60)
All_PISCO.edit <- All_PISCO.mpa[, colnames(All_PISCO.mpa)[c(4,1,2,9,3)]]
All_PISCO.edit <- All_PISCO.edit[complete.cases(All_PISCO.edit), ] #Delete any rows where there isn't an mpa and reference
All_PISCO.edit <- subset(All_PISCO.edit, site_status != "")
PISCO.bio <- All_PISCO.edit %>% 
  spread(site_status, Bio)
PISCO.bio$source <- "PISCO"
PISCO.bio <- PISCO.bio[complete.cases(PISCO.bio), ] #Delete any rows where there isn't an mpa and reference
PISCO.bio.long <- gather(PISCO.bio, status, value, 4:5)
PISCO.bio.long$resp <- 'Bio'

PISCO.resp <- rbind(PISCO.den.long, PISCO.bio.long)
#################################################################################################################

#################################################################################################################
######Combine PISCO response ratio data##########################################################################
Swath.join <- rbind(All.bio.diff.final, All.den.diff.final)
Swath.join.sub <- Swath.join %>% arrange(CA_MPA_Name_Short, year)
x <- unique(Swath.join.sub$CA_MPA_Name_Short) # redo since sites have been removed
l = length(x)
Swath.join.sub$time <- 0
Swath.join.sub <- data.frame(Swath.join.sub)

for(i in 1:l) {
  idx <-which(x[i] == Swath.join.sub$CA_MPA_Name_Short)
  j <-which(x[i] == Site$Site)
  m <- length(idx)
  if (Site$MPA_Start[j] == 2003) {
    for (n in 1:m) {
      if (Swath.join.sub$year[idx[n]] <=2003) {
        Swath.join.sub$time[idx[n]] = 0
      } else if (Swath.join.sub$year[idx[n]] == 2004) {
        Swath.join.sub$time[idx[n]] = 1
      } else if (Swath.join.sub$year[idx[n]] == 2005) {
        Swath.join.sub$time[idx[n]] = 2
      } else if (Swath.join.sub$year[idx[n]] == 2006) {
        Swath.join.sub$time[idx[n]] = 3
      } else if (Swath.join.sub$year[idx[n]] == 2007) {
        Swath.join.sub$time[idx[n]] = 4
      } else if (Swath.join.sub$year[idx[n]] == 2008) {
        Swath.join.sub$time[idx[n]] = 5
      } else if (Swath.join.sub$year[idx[n]] == 2009) {
        Swath.join.sub$time[idx[n]] = 6
      } else if (Swath.join.sub$year[idx[n]] == 2010) {
        Swath.join.sub$time[idx[n]] = 7
      } else if (Swath.join.sub$year[idx[n]] == 2011) {
        Swath.join.sub$time[idx[n]] = 8
      } else if (Swath.join.sub$year[idx[n]] == 2012) {
        Swath.join.sub$time[idx[n]] = 9
      } else if (Swath.join.sub$year[idx[n]] == 2013) {
        Swath.join.sub$time[idx[n]] = 10
      } else if (Swath.join.sub$year[idx[n]] == 2014) {
        Swath.join.sub$time[idx[n]] = 11
      } else if (Swath.join.sub$year[idx[n]] == 2015) {
        Swath.join.sub$time[idx[n]] = 12
      } else if (Swath.join.sub$year[idx[n]] == 2016) {
        Swath.join.sub$time[idx[n]] = 13
      } else if (Swath.join.sub$year[idx[n]] == 2017) {
        Swath.join.sub$time[idx[n]] = 14
      } else if (Swath.join.sub$year[idx[n]] == 2018) {
        Swath.join.sub$time[idx[n]] = 15
      } else if (Swath.join.sub$year[idx[n]] == 2019) {
        Swath.join.sub$time[idx[n]] = 16
      } else if (Swath.join.sub$year[idx[n]] == 2020) {
        Swath.join.sub$time[idx[n]] = 17
      } else if (Swath.join.sub$year[idx[n]] == 2021) {
        Swath.join.sub$time[idx[n]] = 18
      } else if (Swath.join.sub$year[idx[n]] == 2022) {
        Swath.join.sub$time[idx[n]] = 19
      } else if (Swath.join.sub$year[idx[n]] == 2023) {
        Swath.join.sub$time[idx[n]] = 20
      }}
        
    } else if (Site$MPA_Start[j] == 2012) {
      for (n in 1:m) {
        if (Swath.join.sub$year[idx[n]] <= 2012) {
            Swath.join.sub$time[idx[n]] = 0
          } else if (Swath.join.sub$year[idx[n]] == 2013) {
            Swath.join.sub$time[idx[n]] = 1
          } else if (Swath.join.sub$year[idx[n]] == 2014) {
            Swath.join.sub$time[idx[n]] = 2
          } else if (Swath.join.sub$year[idx[n]] == 2015) {
            Swath.join.sub$time[idx[n]] = 3
          } else if (Swath.join.sub$year[idx[n]] == 2016) {
            Swath.join.sub$time[idx[n]] = 4
          } else if (Swath.join.sub$year[idx[n]] == 2017) {
            Swath.join.sub$time[idx[n]] = 5
          } else if (Swath.join.sub$year[idx[n]] == 2018) {
            Swath.join.sub$time[idx[n]] = 6
          } else if (Swath.join.sub$year[idx[n]] == 2019) {
            Swath.join.sub$time[idx[n]] = 7
          } else if (Swath.join.sub$year[idx[n]] == 2020) {
            Swath.join.sub$time[idx[n]] = 8
          } else if (Swath.join.sub$year[idx[n]] == 2021) {
            Swath.join.sub$time[idx[n]] = 9
          } else if (Swath.join.sub$year[idx[n]] == 2022) {
            Swath.join.sub$time[idx[n]] = 10
          } else if (Swath.join.sub$year[idx[n]] == 2023) {
            Swath.join.sub$time[idx[n]] = 11
          }}
} else if (Site$MPA_Start[j] == 1978) {
      for (n in 1:m) {
        if (Swath.join.sub$year[idx[n]] <= 1978) {
            Swath.join.sub$time[idx[n]] = 0
          } else if (Swath.join.sub$year[idx[n]] == 2007) {
            Swath.join.sub$time[idx[n]] = 29
          } else if (Swath.join.sub$year[idx[n]] == 2008) {
            Swath.join.sub$time[idx[n]] = 30
          } else if (Swath.join.sub$year[idx[n]] == 2009) {
            Swath.join.sub$time[idx[n]] = 31
          } else if (Swath.join.sub$year[idx[n]] == 2010) {
            Swath.join.sub$time[idx[n]] = 32
          } else if (Swath.join.sub$year[idx[n]] == 2011) {
            Swath.join.sub$time[idx[n]] = 33
          } else if (Swath.join.sub$year[idx[n]] == 2012) {
            Swath.join.sub$time[idx[n]] = 34
          } else if (Swath.join.sub$year[idx[n]] == 2013) {
            Swath.join.sub$time[idx[n]] = 35
          } else if (Swath.join.sub$year[idx[n]] == 2014) {
            Swath.join.sub$time[idx[n]] = 36
          } else if (Swath.join.sub$year[idx[n]] == 2015) {
            Swath.join.sub$time[idx[n]] = 37
          } else if (Swath.join.sub$year[idx[n]] == 2016) {
            Swath.join.sub$time[idx[n]] = 38
          } else if (Swath.join.sub$year[idx[n]] == 2017) {
            Swath.join.sub$time[idx[n]] = 39
          } else if (Swath.join.sub$year[idx[n]] == 2018) {
            Swath.join.sub$time[idx[n]] = 40
          } else if (Swath.join.sub$year[idx[n]] == 2019) {
            Swath.join.sub$time[idx[n]] = 41
          } else if (Swath.join.sub$year[idx[n]] == 2020) {
            Swath.join.sub$time[idx[n]] = 42
          } else if (Swath.join.sub$year[idx[n]] == 2021) {
            Swath.join.sub$time[idx[n]] = 43
          } else if (Swath.join.sub$year[idx[n]] == 2022) {
            Swath.join.sub$time[idx[n]] = 44
          } else if (Swath.join.sub$year[idx[n]] == 2023) {
            Swath.join.sub$time[idx[n]] = 45
          }}}
    else {}
  }
   
for(i in 1:l) {
  idx <-which(x[i] == Swath.join.sub$CA_MPA_Name_Short)
  j <-which(x[i] == Site$Site)
  m <- length(idx)
  if (Site$MPA_Start[j] == 2003) {
    for (n in 1:m) {
      if (Swath.join.sub$year[idx[n]] <=2003) {
        Swath.join.sub$time[idx[n]] = 0
      } else if (Swath.join.sub$year[idx[n]] == 2004) {
        Swath.join.sub$time[idx[n]] = 1
      } else if (Swath.join.sub$year[idx[n]] == 2005) {
        Swath.join.sub$time[idx[n]] = 2
      } else if (Swath.join.sub$year[idx[n]] == 2006) {
        Swath.join.sub$time[idx[n]] = 3
      } else if (Swath.join.sub$year[idx[n]] == 2007) {
        Swath.join.sub$time[idx[n]] = 4
      } else if (Swath.join.sub$year[idx[n]] == 2008) {
        Swath.join.sub$time[idx[n]] = 5
      } else if (Swath.join.sub$year[idx[n]] == 2009) {
        Swath.join.sub$time[idx[n]] = 6
      } else if (Swath.join.sub$year[idx[n]] == 2010) {
        Swath.join.sub$time[idx[n]] = 7
      } else if (Swath.join.sub$year[idx[n]] == 2011) {
        Swath.join.sub$time[idx[n]] = 8
      } else if (Swath.join.sub$year[idx[n]] == 2012) {
        Swath.join.sub$time[idx[n]] = 9
      } else if (Swath.join.sub$year[idx[n]] == 2013) {
        Swath.join.sub$time[idx[n]] = 10
      } else if (Swath.join.sub$year[idx[n]] == 2014) {
        Swath.join.sub$time[idx[n]] = 11
      } else if (Swath.join.sub$year[idx[n]] == 2015) {
        Swath.join.sub$time[idx[n]] = 12
      } else if (Swath.join.sub$year[idx[n]] == 2016) {
        Swath.join.sub$time[idx[n]] = 13
      } else if (Swath.join.sub$year[idx[n]] == 2017) {
        Swath.join.sub$time[idx[n]] = 14
      } else if (Swath.join.sub$year[idx[n]] == 2018) {
        Swath.join.sub$time[idx[n]] = 15
      } else if (Swath.join.sub$year[idx[n]] == 2019) {
        Swath.join.sub$time[idx[n]] = 16
      } else if (Swath.join.sub$year[idx[n]] == 2020) {
        Swath.join.sub$time[idx[n]] = 17
      } else if (Swath.join.sub$year[idx[n]] == 2021) {
        Swath.join.sub$time[idx[n]] = 18
      } else if (Swath.join.sub$year[idx[n]] == 2022) {
        Swath.join.sub$time[idx[n]] = 19
      } else if (Swath.join.sub$year[idx[n]] == 2023) {
        Swath.join.sub$time[idx[n]] = 20
      }}
        
    } else if (Site$MPA_Start[j] == 2012) {
      for (n in 1:m) {
        if (Swath.join.sub$year[idx[n]] <= 2012) {
            Swath.join.sub$time[idx[n]] = 0
          } else if (Swath.join.sub$year[idx[n]] == 2013) {
            Swath.join.sub$time[idx[n]] = 1
          } else if (Swath.join.sub$year[idx[n]] == 2014) {
            Swath.join.sub$time[idx[n]] = 2
          } else if (Swath.join.sub$year[idx[n]] == 2015) {
            Swath.join.sub$time[idx[n]] = 3
          } else if (Swath.join.sub$year[idx[n]] == 2016) {
            Swath.join.sub$time[idx[n]] = 4
          } else if (Swath.join.sub$year[idx[n]] == 2017) {
            Swath.join.sub$time[idx[n]] = 5
          } else if (Swath.join.sub$year[idx[n]] == 2018) {
            Swath.join.sub$time[idx[n]] = 6
          } else if (Swath.join.sub$year[idx[n]] == 2019) {
            Swath.join.sub$time[idx[n]] = 7
          } else if (Swath.join.sub$year[idx[n]] == 2020) {
            Swath.join.sub$time[idx[n]] = 8
          } else if (Swath.join.sub$year[idx[n]] == 2021) {
            Swath.join.sub$time[idx[n]] = 9
          } else if (Swath.join.sub$year[idx[n]] == 2022) {
            Swath.join.sub$time[idx[n]] = 10
          } else if (Swath.join.sub$year[idx[n]] == 2023) {
            Swath.join.sub$time[idx[n]] = 11
          }}
} else if (Site$MPA_Start[j] == 1978) {
      for (n in 1:m) {
        if (Swath.join.sub$year[idx[n]] <= 1978) {
            Swath.join.sub$time[idx[n]] = 0
          } else if (Swath.join.sub$year[idx[n]] == 2007) {
            Swath.join.sub$time[idx[n]] = 29
          } else if (Swath.join.sub$year[idx[n]] == 2008) {
            Swath.join.sub$time[idx[n]] = 30
          } else if (Swath.join.sub$year[idx[n]] == 2009) {
            Swath.join.sub$time[idx[n]] = 31
          } else if (Swath.join.sub$year[idx[n]] == 2010) {
            Swath.join.sub$time[idx[n]] = 32
          } else if (Swath.join.sub$year[idx[n]] == 2011) {
            Swath.join.sub$time[idx[n]] = 33
          } else if (Swath.join.sub$year[idx[n]] == 2012) {
            Swath.join.sub$time[idx[n]] = 34
          } else if (Swath.join.sub$year[idx[n]] == 2013) {
            Swath.join.sub$time[idx[n]] = 35
          } else if (Swath.join.sub$year[idx[n]] == 2014) {
            Swath.join.sub$time[idx[n]] = 36
          } else if (Swath.join.sub$year[idx[n]] == 2015) {
            Swath.join.sub$time[idx[n]] = 37
          } else if (Swath.join.sub$year[idx[n]] == 2016) {
            Swath.join.sub$time[idx[n]] = 38
          } else if (Swath.join.sub$year[idx[n]] == 2017) {
            Swath.join.sub$time[idx[n]] = 39
          } else if (Swath.join.sub$year[idx[n]] == 2018) {
            Swath.join.sub$time[idx[n]] = 40
          } else if (Swath.join.sub$year[idx[n]] == 2019) {
            Swath.join.sub$time[idx[n]] = 41
          } else if (Swath.join.sub$year[idx[n]] == 2020) {
            Swath.join.sub$time[idx[n]] = 42
          } else if (Swath.join.sub$year[idx[n]] == 2021) {
            Swath.join.sub$time[idx[n]] = 43
          } else if (Swath.join.sub$year[idx[n]] == 2022) {
            Swath.join.sub$time[idx[n]] = 44
          } else if (Swath.join.sub$year[idx[n]] == 2023) {
            Swath.join.sub$time[idx[n]] = 45
          }}}
    else {}
  }
#Import site info that has MPA established date
Site.size = read.csv(paste('data/MPAfeatures_subset.csv'))
Site = read.csv(paste('data/Site_List_All.csv'))

#Merge site data tables so that we have size info about the respective MPA
Site <- merge(Site, Site.size, by.x = c("Site"), 
              by.y = c("NAME"),
              all = TRUE)

#Remove rows that are missing data
Site <- Site[!is.na(Site$Lat),] 
sites.short.edit = subset(Site, select = c(CA_MPA_Name_Short, type, Location, Hectares) )
Swath.join.sub <- merge(Swath.join.sub, sites.short.edit, by.x = c("CA_MPA_Name_Short"), 
                        by.y = c("CA_MPA_Name_Short"),
                        all = TRUE)

Swath.join.sub <- Swath.join.sub[complete.cases(Swath.join.sub$year), ]
Swath.join.sub <- subset(Swath.join.sub, CA_MPA_Name_Short != "Anacapa Island SMR 1978" | 
                           CA_MPA_Name_Short != "Anacapa Island SMCA" |
                           CA_MPA_Name_Short != "Point Dume SMR") # Anacapa only have data 20 years + after implementation, PT. Dume SMR only a since data point
Swath.join.sub$source <- "PISCO"
Swath.join.sub$BA <- NA

l = nrow(Swath.join.sub)

for(i in 1:l) {
  if (Swath.join.sub$time[i] == 0) {
        Swath.join.sub$BA[i] = "Before"}
  else {
    Swath.join.sub$BA[i] = "After"
  }}


####################################################################################################################

####################################################################################################################
####################MBON_Synthesis##################################################################################
MBON = read.csv(paste('data/MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv'))
Sites2 = read.csv(paste('data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv'))
Sites2$site_id <- as.factor(Sites2$site_id)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings
if (class(MBON$data_source)!="factor") MBON$data_source<- as.factor(MBON$data_source)
if (class(MBON$sample_method)!="factor") MBON$sample_method<- as.factor(MBON$sample_method)                                   
# convert MBON$date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1date<-as.Date(MBON$date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1date) == length(tmp1date[!is.na(tmp1date)])){MBON$date <- tmp1date } else {print("Date conversion failed for MBON$date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1date) 
if (class(MBON$site_id)!="factor") MBON$site_id<- as.factor(MBON$site_id)
if (class(MBON$subsite_id)!="factor") MBON$subsite_id<- as.factor(MBON$subsite_id)
if (class(MBON$transect_id)!="factor") MBON$transect_id<- as.factor(MBON$transect_id)
if (class(MBON$replicate_id)!="factor") MBON$replicate_id<- as.factor(MBON$replicate_id)
if (class(MBON$proj_taxon_id)!="factor") MBON$proj_taxon_id<- as.factor(MBON$proj_taxon_id)
if (class(MBON$area)=="factor") MBON$area <-as.numeric(levels(MBON$area))[as.integer(MBON$area) ]               
if (class(MBON$area)=="character") MBON$area <-as.numeric(MBON$area)
if (class(MBON$count)=="factor") MBON$count <-as.numeric(levels(MBON$count))[as.integer(MBON$count) ]               
if (class(MBON$count)=="character") MBON$count <-as.numeric(MBON$count)
if (class(MBON$auth_taxon_id)!="factor") MBON$auth_taxon_id<- as.factor(MBON$auth_taxon_id)
if (class(MBON$auth_name)!="factor") MBON$auth_name<- as.factor(MBON$auth_name)
if (class(MBON$taxon_name)!="factor") MBON$taxon_name<- as.factor(MBON$taxon_name)
if (class(MBON$site_name)!="factor") MBON$site_name<- as.factor(MBON$site_name)
if (class(MBON$subsite_name)!="factor") MBON$subsite_name<- as.factor(MBON$subsite_name)
if (class(MBON$latitude)=="factor") MBON$latitude <-as.numeric(levels(MBON$latitude))[as.integer(MBON$latitude) ]               
if (class(MBON$latitude)=="character") MBON$latitude <-as.numeric(MBON$latitude)
if (class(MBON$longitude)=="factor") MBON$longitude <-as.numeric(levels(MBON$longitude))[as.integer(MBON$longitude) ]               
if (class(MBON$longitude)=="character") MBON$longitude <-as.numeric(MBON$longitude)

# Convert Missing Values to NA for non-dates

MBON$count <- ifelse((trimws(as.character(MBON$count))==trimws(".")),NA,MBON$count)               
suppressWarnings(MBON$count <- ifelse(!is.na(as.numeric(".")) & (trimws(as.character(MBON$count))==as.character(as.numeric("."))),NA,MBON$count))
MBON$auth_taxon_id <- as.factor(ifelse((trimws(as.character(MBON$auth_taxon_id))==trimws(".")),NA,as.character(MBON$auth_taxon_id)))
MBON$auth_name <- as.factor(ifelse((trimws(as.character(MBON$auth_name))==trimws(".")),NA,as.character(MBON$auth_name)))
MBON$subsite_name <- as.factor(ifelse((trimws(as.character(MBON$subsite_name))==trimws(".")),NA,as.character(MBON$subsite_name)))
MBON$latitude <- ifelse((trimws(as.character(MBON$latitude))==trimws(".")),NA,MBON$latitude)               
suppressWarnings(MBON$latitude <- ifelse(!is.na(as.numeric(".")) & (trimws(as.character(MBON$latitude))==as.character(as.numeric("."))),NA,MBON$latitude))
MBON$longitude <- ifelse((trimws(as.character(MBON$longitude))==trimws(".")),NA,MBON$longitude)               
suppressWarnings(MBON$longitude <- ifelse(!is.na(as.numeric(".")) & (trimws(as.character(MBON$longitude))==as.character(as.numeric("."))),NA,MBON$longitude))

kfm <- subset(MBON, data_source == "kfm")
lter <- subset(MBON, data_source == "lter")
###########################################################################################

###########################################################################################
######Processing KFM data#################################################################

###Start with SWATH data##################################################################
##Join average swath/quad data to site table###########################################
#######################################################################################
kfm.site <- merge(kfm, Sites2, by = c("site_id"))
kfm.site$year <- year(kfm.site$date)

kfm.site <- subset(kfm.site, CA_MPA_Name_Short != "")
#Remove sites where there are long time series so don't need short time series that wasn't a part of the before comparison; 2) wrong side of island and better reference exists
kfm.site <- subset(kfm.site, site_id != "a-k-35" & site_id != "a-k-36" & site_id != "a-k-37" & site_id != "a-k-21" & site_id != "a-k-05" 
                         & site_id != "a-k-30" & site_id != "a-k-29" & site_id != "a-k-28" & site_id != "a-k-27" & site_id != "a-k-26" 
                         & site_id != "a-k-21" & site_id != "a-k-07" & site_id != "a-k-12" & site_id != "a-k-13" & site_id != "a-k-31"
                         & site_id != "a-k-35" & site_id != "a-k-36" & site_id != "a-k-37")

kfm.site$time <- 0
kfm.site$year <- year(kfm.site$date)
x <- unique(kfm.site$CA_MPA_Name_Short) # to redo since sites have been removed
l = length(x)
for(i in 1:l) {
  idx <-which(x[i] == kfm.site$CA_MPA_Name_Short)
  j <-which(x[i] == Sites2$CA_MPA_Name_Short)
  m <- length(idx)
  if (Sites2$MPA_Start[j] == 2003) {
    for (n in 1:m) {
      if (kfm.site$year[idx[n]] <=2003) {
        kfm.site$time[idx[n]] = 0
      } else if (kfm.site$year[idx[n]] == 2004) {
        kfm.site$time[idx[n]] = 1
      } else if (kfm.site$year[idx[n]] == 2005) {
        kfm.site$time[idx[n]] = 2
      } else if (kfm.site$year[idx[n]] == 2006) {
        kfm.site$time[idx[n]] = 3
      } else if (kfm.site$year[idx[n]] == 2007) {
        kfm.site$time[idx[n]] = 4
      } else if (kfm.site$year[idx[n]] == 2008) {
        kfm.site$time[idx[n]] = 5
      } else if (kfm.site$year[idx[n]] == 2009) {
        kfm.site$time[idx[n]] = 6
      } else if (kfm.site$year[idx[n]] == 2010) {
        kfm.site$time[idx[n]] = 7
      } else if (kfm.site$year[idx[n]] == 2011) {
        kfm.site$time[idx[n]] = 8
      } else if (kfm.site$year[idx[n]] == 2012) {
        kfm.site$time[idx[n]] = 9
      } else if (kfm.site$year[idx[n]] == 2013) {
        kfm.site$time[idx[n]] = 10
      } else if (kfm.site$year[idx[n]] == 2014) {
        kfm.site$time[idx[n]] = 11
      } else if (kfm.site$year[idx[n]] == 2015) {
        kfm.site$time[idx[n]] = 12
      } else if (kfm.site$year[idx[n]] == 2016) {
        kfm.site$time[idx[n]] = 13
      } else if (kfm.site$year[idx[n]] == 2017) {
        kfm.site$time[idx[n]] = 14
      } else if (kfm.site$year[idx[n]] == 2018) {
        kfm.site$time[idx[n]] = 15
      } else if (kfm.site$year[idx[n]] == 2019) {
        kfm.site$time[idx[n]] = 16
      } else if (kfm.site$year[idx[n]] == 2020) {
        kfm.site$time[idx[n]] = 17
      } else if (kfm.site$year[idx[n]] == 2021) {
        kfm.site$time[idx[n]] = 18
      } else if (kfm.site$year[idx[n]] == 2022) {
        kfm.site$time[idx[n]] = 19
      } else if (kfm.site$year[idx[n]] == 2023) {
        kfm.site$time[idx[n]] = 20
      }}
    
  } else if (Sites2$MPA_Start[j] == 1978) {
    for (n in 1:m) {
      if (kfm.site$year[idx[n]] <= 1978) {
        kfm.site$time[idx[n]] = 0
      } else if (kfm.site$year[idx[n]] == 1982) {
        kfm.site$time[idx[n]] = 4
      } else if (kfm.site$year[idx[n]] == 1983) {
        kfm.site$time[idx[n]] = 5
      } else if (kfm.site$year[idx[n]] == 1984) {
        kfm.site$time[idx[n]] = 6
      } else if (kfm.site$year[idx[n]] == 1985) {
        kfm.site$time[idx[n]] = 7
      } else if (kfm.site$year[idx[n]] == 1986) {
        kfm.site$time[idx[n]] = 8
      } else if (kfm.site$year[idx[n]] == 1987) {
        kfm.site$time[idx[n]] = 9
      } else if (kfm.site$year[idx[n]] == 1988) {
        kfm.site$time[idx[n]] = 10
      } else if (kfm.site$year[idx[n]] == 1989) {
        kfm.site$time[idx[n]] = 11
      } else if (kfm.site$year[idx[n]] == 1990) {
        kfm.site$time[idx[n]] = 12
      } else if (kfm.site$year[idx[n]] == 1991) {
        kfm.site$time[idx[n]] = 13
      } else if (kfm.site$year[idx[n]] == 1992) {
        kfm.site$time[idx[n]] = 14
      } else if (kfm.site$year[idx[n]] == 1993) {
        kfm.site$time[idx[n]] = 15
      } else if (kfm.site$year[idx[n]] == 1994) {
        kfm.site$time[idx[n]] = 16
      } else if (kfm.site$year[idx[n]] == 1995) {
        kfm.site$time[idx[n]] = 17
      } else if (kfm.site$year[idx[n]] == 1996) {
        kfm.site$time[idx[n]] = 18
      } else if (kfm.site$year[idx[n]] == 1997) {
        kfm.site$time[idx[n]] = 19
      } else if (kfm.site$year[idx[n]] == 1998) {
        kfm.site$time[idx[n]] = 20
      } else if (kfm.site$year[idx[n]] == 1999) {
        kfm.site$time[idx[n]] = 21
      } else if (kfm.site$year[idx[n]] == 2000) {
        kfm.site$time[idx[n]] = 22
      } else if (kfm.site$year[idx[n]] == 2001) {
        kfm.site$time[idx[n]] = 23
      } else if (kfm.site$year[idx[n]] == 2002) {
        kfm.site$time[idx[n]] = 24
      } else if (kfm.site$year[idx[n]] == 2003) {
        kfm.site$time[idx[n]] = 25
      } else if (kfm.site$year[idx[n]] == 2004) {
        kfm.site$time[idx[n]] = 26
      } else if (kfm.site$year[idx[n]] == 2005) {
        kfm.site$time[idx[n]] = 27
      } else if (kfm.site$year[idx[n]] == 2006) {
        kfm.site$time[idx[n]] = 28
      } else if (kfm.site$year[idx[n]] == 2007) {
        kfm.site$time[idx[n]] = 29
      } else if (kfm.site$year[idx[n]] == 2008) {
        kfm.site$time[idx[n]] = 30
      } else if (kfm.site$year[idx[n]] == 2009) {
        kfm.site$time[idx[n]] = 31
      } else if (kfm.site$year[idx[n]] == 2010) {
        kfm.site$time[idx[n]] = 32
      } else if (kfm.site$year[idx[n]] == 2011) {
        kfm.site$time[idx[n]] = 33
      } else if (kfm.site$year[idx[n]] == 2012) {
        kfm.site$time[idx[n]] = 34
      } else if (kfm.site$year[idx[n]] == 2013) {
        kfm.site$time[idx[n]] = 35
      } else if (kfm.site$year[idx[n]] == 2014) {
        kfm.site$time[idx[n]] = 36
      } else if (kfm.site$year[idx[n]] == 2015) {
        kfm.site$time[idx[n]] = 37
      } else if (kfm.site$year[idx[n]] == 2016) {
        kfm.site$time[idx[n]] = 38
      } else if (kfm.site$year[idx[n]] == 2017) {
        kfm.site$time[idx[n]] = 39
      } else if (kfm.site$year[idx[n]] == 2018) {
        kfm.site$time[idx[n]] = 40
      } else if (kfm.site$year[idx[n]] == 2019) {
        kfm.site$time[idx[n]] = 41
      } else if (kfm.site$year[idx[n]] == 2020) {
        kfm.site$time[idx[n]] = 42
      } else if (kfm.site$year[idx[n]] == 2021) {
        kfm.site$time[idx[n]] = 43
      } else if (kfm.site$year[idx[n]] == 2022) {
        kfm.site$time[idx[n]] = 44
      } else if (kfm.site$year[idx[n]] == 2023) {
        kfm.site$time[idx[n]] = 45
      }}}
  else {}
}

kfm.site <- subset(kfm.site, taxon_name == "Strongylocentrotus purpuratus" | taxon_name == "Mesocentrotus franciscanus"
                   | taxon_name == "Macrocystis pyrifera" | taxon_name == "Panulirus interruptus")
x <- unique(kfm.site$CA_MPA_Name_Short) # to redo since sites have been removed
MPA_implement <- filter(Sites2, CA_MPA_Name_Short %in% x)
MPA_implement.sub <- filter(MPA_implement, CA_MPA_Name_Short %in% x)

#Remove juvenile giant kelp since this is not done by any other programs and therefore is not comparable
kfm.site <- subset(kfm.site, taxon_name != "Macrocystis pyrifera" | proj_taxon_id != "t-k-002")

#Need to sum across different Macro names since we're considering all 
kfm.sum <- kfm.site %>%
  group_by(site_id, taxon_name, status, CA_MPA_Name_Short, area, ChannelIsland, replicate_id, MPA_Start, time, year) %>% 
  summarise_at(c("count"), sum) %>%
  ungroup()

#Transform to DENSITY since there are different survey methods
kfm.sum$den <- kfm.sum$count/kfm.sum$area

#Ave across sites and replicates
kfm.ave <- kfm.sum %>%
  group_by(site_id, taxon_name, status, CA_MPA_Name_Short, area, ChannelIsland, MPA_Start, time, year) %>% 
  summarise_at(c("den"), mean) %>%
  ungroup()
##########################################################################################

##########################################################################################
########Calculating biomass of urchins from size freq distributions#######################
##########################################################################################
URCHINS <- subset(kfm.site, taxon_name == "Strongylocentrotus purpuratus" | taxon_name == "Mesocentrotus franciscanus")
URCHINS <- mutate(URCHINS, taxon_name = case_when(
  taxon_name == "Strongylocentrotus purpuratus"  ~ "STRPURAD",
  taxon_name == "Mesocentrotus franciscanus"    ~ "MESFRAAD",
  TRUE ~ taxon_name))

u <- unique(URCHINS[,c('year','site_id','CA_MPA_Name_Short','status','area','taxon_name')])

#Create empty dataframes
Urchin.site <- data.frame(site = NA, CA_MPA_Name_Short = NA,  site_status = NA, year = NA, transect = NA, count = NA, biomass = NA, y = NA, area = NA)
Urchin.SF <- data.frame(site = NA, CA_MPA_Name_Short = NA,  site_status = NA, year = NA, transect = NA, count = NA, biomass = NA, y = NA, area = NA)

SizeFreq.Urch <- SizeFreq.Urch.OG
#Select site/year/zone if there are lobsters then look up same site/year/zone from size freq table since this is where all size data are
#create lobster "population" from size freq table and draw with replacement the number of lobsters from swath data table 
#if lobsters were measured on transect; however, if no lobsters then equals zero, if there were lobsters but no size freq
#these will be removed and marked as NA

for(i in 1:nrow(u)) {
  t <- which(URCHINS$site_id == u$site_id[i] & URCHINS$CA_MPA_Name_Short == u$CA_MPA_Name_Short[i] & URCHINS$year == u$year[i] 
             & URCHINS$taxon_name == u$taxon_name[i])
  n <- sum(URCHINS$count[t])
  t2 <- which(SizeFreq.Urch$CA_MPA_Name_Short == u$CA_MPA_Name_Short[i] & SizeFreq.Urch$site_status == u$status[i] &
                SizeFreq.Urch$year == u$year[i] & SizeFreq.Urch$classcode == u$taxon_name[i])
  
  if (n != 0 & length(t2) != 0) { #This means an urchin was observed on the transect and in the size freq table
    a = NA
    s = data.frame(matrix(NA, nrow = 1000, ncol = n))
    l <- length(t2)
    for(j in 1:l){
      list <- as.data.frame(rep(SizeFreq.Urch$size[t2[j]], SizeFreq.Urch$count[t2[j]] ))
      a <- rbind(a,list)
      a <- na.omit(a) }
      a <- as.numeric(as.character(unlist(a))) #change list to numeric so that it can be sampled from, this is each urchin size 
      #randomly sample the urchin "population" with replacement; do this 1000 times
      for(k in 1:1000){
        s[k,] <- sample(a,n, replace = TRUE) 
      }
      #Convert CL to biomass for all individuals
      if(u$taxon_name[i] == "MESFRAAD"){
        s_bio <- s %>% 
          mutate(across(starts_with("X"), ~bio_redurch(.)))
        s_bio <- s_bio %>% 
          mutate(sum = rowSums(.))
        ave <- mean(s_bio$sum)
      }else {
        s_bio <- s %>% 
          mutate(across(starts_with("X"), ~bio_purpurch(.)))
        s_bio <- s_bio %>% 
          mutate(sum = rowSums(.))
        ave <- mean(s_bio$sum)
      }
      
      Urchin.SF$site <- u$site_id[i]
      Urchin.SF$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
      Urchin.SF$site_status <- u$status[i]
      Urchin.SF$year <- u$year[i]
      Urchin.SF$area <- u$area[i]
      Urchin.SF$transect <- length(t)
      Urchin.SF$count <- n
      Urchin.SF$biomass <- ave
      Urchin.SF$y <- u$taxon_name[i]
      Urchin.site <- rbind(Urchin.site, Urchin.SF)
      
    } else if (n == 0 & length(t2) == 0){#This means there were no lobsters on transect and none in the size freq data either
      Urchin.SF$site <- u$site_id[i]
      Urchin.SF$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
      Urchin.SF$site_status <- u$status[i]
      Urchin.SF$year <- u$year[i]
      Urchin.SF$area <- u$area[i]
      Urchin.SF$transect <- length(t)
      Urchin.SF$count <- 0
      Urchin.SF$biomass <- 0
      Urchin.SF$y <- u$taxon_name[i]
      Urchin.site <- rbind(Urchin.site,Urchin.SF)
      
    } else if (n != 0 & length(t2) == 0){ #If there are lobsters on transects, but they weren't sized we will mark size as NA
      Urchin.SF$site <- u$site_id[i]
      Urchin.SF$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
      Urchin.SF$site_status <- u$status[i]
      Urchin.SF$year <- u$year[i]
      Urchin.SF$area <- u$area[i]
      Urchin.SF$transect <- length(t)
      Urchin.SF$count <- n
      Urchin.SF$biomass <- NA
      Urchin.SF$y <- u$taxon_name[i]
      Urchin.site <- rbind(Urchin.site,Urchin.SF)
      
    } else {
      Urchin.SF$site <- u$site_id[i]
      Urchin.SF$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
      Urchin.SF$site_status <- u$status[i]
      Urchin.SF$year <- u$year[i]
      Urchin.SF$area <- u$area[i]
      Urchin.SF$transect <- length(t)
      Urchin.SF$count <- 0
      Urchin.SF$biomass <- 0
      Urchin.SF$y <- u$taxon_name[i]
      Urchin.site <- rbind(Urchin.site,Urchin.SF)
    }}

KFM.Urchin.site.merge <- merge(Urchin.site, sites.short, by.x = c("site"), by.y = c("site"),
                                 all.x = TRUE)

KFM.Urchin.site.merge <- KFM.Urchin.site.merge[, colnames(KFM.Urchin.site.merge)[c(2,1,17,18,4:9)]]

#Remove Old Anacapa MPA and the SMCA since they don't meet criteria for inclusion
KFM.Urchin.site.merge <- subset(KFM.Urchin.site.merge, CA_MPA_Name_Short.x != "Anacapa Island SMCA" & 
                                  CA_MPA_Name_Short.x != "Anacapa Island SMR 1978")

KFM.Urchin.site.merge[is.na(KFM.Urchin.site.merge)] <- 0

#Find sites where urchin was found on transect but no size was recorded
KFM.Urchin.null <- subset(KFM.Urchin.site.merge, count > 0 & biomass == 0) 

KFM.Urchin.site.merge.sub <- KFM.Urchin.site.merge

#Remove sites where urchins were seen but not measured since we can't get an accurate biomass estimate
for(i in 1:nrow(KFM.Urchin.null)) {
  KFM.Urchin.site.merge.sub <- subset(KFM.Urchin.site.merge.sub, site != KFM.Urchin.null$site[i] | year != KFM.Urchin.null$year[i])
}

names(KFM.Urchin.site.merge.sub)[names(KFM.Urchin.site.merge.sub) == "CA_MPA_Name_Short.x"] <- "CA_MPA_Name_Short"
names(KFM.Urchin.site.merge.sub)[names(KFM.Urchin.site.merge.sub) == "site_status.y"] <- "site_status"

#Need to account for the different number of transects at the sites
KFM.Urchin.site.merge.sub$count.ave <- KFM.Urchin.site.merge.sub$count/KFM.Urchin.site.merge.sub$transect
KFM.Urchin.site.merge.sub$biomass.ave <- KFM.Urchin.site.merge.sub$biomass/KFM.Urchin.site.merge.sub$transect

KFM.Urchin.site.all <- KFM.Urchin.site.merge.sub 

#Quads are 2m2, so need to divide by 2 to get per m2
KFM.Urchin.site.all$Density <- KFM.Urchin.site.all$count.ave/2
KFM.Urchin.site.all$bio.m2 <- KFM.Urchin.site.all$biomass.ave/2

KFM.Urchin.site.all <- mutate(KFM.Urchin.site.all, y = case_when(
  y == "STRPURAD"  ~ "Strongylocentrotus purpuratus",
  y == "MESFRAAD"    ~ "Mesocentrotus franciscanus",
  TRUE ~ y))


######Start with stipes####################################################################

#Import raw stipe counts
kfm.stipe = read.csv(paste('data/MBON/KFM_Macrocystis_RawData_1984-2023.csv'))
kfm.stipe <- kfm.stipe[!duplicated(kfm.stipe),] #Data issues caused duplicates

#Create site_id that is the same as in site table
kfm.stipe$site_id <- ifelse(kfm.stipe$SiteNumber > 9,  paste0("a-k-", kfm.stipe$SiteNumber), paste0("a-k-0", kfm.stipe$SiteNumber))
kfm.stipe.site <- merge(kfm.stipe, Sites2, by = c("site_id"))
kfm.stipe.site <- kfm.stipe.site %>% 
  mutate(date = as.Date(SurveyDate, format = "%m/%d/%Y"))

kfm.stipe.site$year <- year(kfm.stipe.site$date)

#Remove sites where 1) there are long time series so don't need short time series that wasn't a part of the before comparison; 2) wrong side of island and better reference exists
kfm.stipe.site <- subset(kfm.stipe.site, site_id != "a-k-35" & site_id != "a-k-36" & site_id != "a-k-37" & site_id != "a-k-21" & site_id != "a-k-05" 
                         & site_id != "a-k-30" & site_id != "a-k-29" & site_id != "a-k-28" & site_id != "a-k-27" & site_id != "a-k-26" 
                         & site_id != "a-k-21" & site_id != "a-k-07" & site_id != "a-k-12" & site_id != "a-k-13" & site_id != "a-k-31"
                         & site_id != "a-k-35" & site_id != "a-k-36" & site_id != "a-k-37")
###########################################################################################
#Calc ave number of stipes per site########################################################
kfm.stipe.mean <- kfm.stipe.site %>%
  group_by(site_id, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, year) %>% 
  summarise_at(c("Stipe_Count"), mean)


#Calc ave number of stipes per MPA/Reference
kfm.stipe.max <- kfm.stipe.mean %>%
  group_by(status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, year) %>% 
  summarise_at(c("Stipe_Count"), mean)


kfm.stipe.max$time <- 0
kfm.stipe.ave <- subset(kfm.stipe.max, CA_MPA_Name_Short != "")
x <- unique(kfm.stipe.ave$CA_MPA_Name_Short) # to redo since sites have been removed
l = length(x)
for(i in 1:l) {
  idx <-which(x[i] == kfm.stipe.ave$CA_MPA_Name_Short)
  j <-which(x[i] == Sites2$CA_MPA_Name_Short)
  m <- length(idx)
  if (Sites2$MPA_Start[j] == 2003) {
    for (n in 1:m) {
      if (kfm.stipe.ave$year[idx[n]] <=2003) {
        kfm.stipe.ave$time[idx[n]] = 0
      } else if (kfm.stipe.ave$year[idx[n]] == 2004) {
        kfm.stipe.ave$time[idx[n]] = 1
      } else if (kfm.stipe.ave$year[idx[n]] == 2005) {
        kfm.stipe.ave$time[idx[n]] = 2
      } else if (kfm.stipe.ave$year[idx[n]] == 2006) {
        kfm.stipe.ave$time[idx[n]] = 3
      } else if (kfm.stipe.ave$year[idx[n]] == 2007) {
        kfm.stipe.ave$time[idx[n]] = 4
      } else if (kfm.stipe.ave$year[idx[n]] == 2008) {
        kfm.stipe.ave$time[idx[n]] = 5
      } else if (kfm.stipe.ave$year[idx[n]] == 2009) {
        kfm.stipe.ave$time[idx[n]] = 6
      } else if (kfm.stipe.ave$year[idx[n]] == 2010) {
        kfm.stipe.ave$time[idx[n]] = 7
      } else if (kfm.stipe.ave$year[idx[n]] == 2011) {
        kfm.stipe.ave$time[idx[n]] = 8
      } else if (kfm.stipe.ave$year[idx[n]] == 2012) {
        kfm.stipe.ave$time[idx[n]] = 9
      } else if (kfm.stipe.ave$year[idx[n]] == 2013) {
        kfm.stipe.ave$time[idx[n]] = 10
      } else if (kfm.stipe.ave$year[idx[n]] == 2014) {
        kfm.stipe.ave$time[idx[n]] = 11
      } else if (kfm.stipe.ave$year[idx[n]] == 2015) {
        kfm.stipe.ave$time[idx[n]] = 12
      } else if (kfm.stipe.ave$year[idx[n]] == 2016) {
        kfm.stipe.ave$time[idx[n]] = 13
      } else if (kfm.stipe.ave$year[idx[n]] == 2017) {
        kfm.stipe.ave$time[idx[n]] = 14
      } else if (kfm.stipe.ave$year[idx[n]] == 2018) {
        kfm.stipe.ave$time[idx[n]] = 15
      } else if (kfm.stipe.ave$year[idx[n]] == 2019) {
        kfm.stipe.ave$time[idx[n]] = 16
      } else if (kfm.stipe.ave$year[idx[n]] == 2020) {
        kfm.stipe.ave$time[idx[n]] = 17
      } else if (kfm.stipe.ave$year[idx[n]] == 2021) {
        kfm.stipe.ave$time[idx[n]] = 18
      } else if (kfm.stipe.ave$year[idx[n]] == 2022) {
        kfm.stipe.ave$time[idx[n]] = 19
      } else if (kfm.stipe.ave$year[idx[n]] == 2023) {
        kfm.stipe.ave$time[idx[n]] = 20
      }}
    
  } else if (Sites2$MPA_Start[j] == 1978) {
    for (n in 1:m) {
      if (kfm.stipe.ave$year[idx[n]] <= 1978) {
        kfm.stipe.ave$time[idx[n]] = 0
      } else if (kfm.stipe.ave$year[idx[n]] == 1982) {
        kfm.stipe.ave$time[idx[n]] = 4
      } else if (kfm.stipe.ave$year[idx[n]] == 1983) {
        kfm.stipe.ave$time[idx[n]] = 5
      } else if (kfm.stipe.ave$year[idx[n]] == 1984) {
        kfm.stipe.ave$time[idx[n]] = 6
      } else if (kfm.stipe.ave$year[idx[n]] == 1985) {
        kfm.stipe.ave$time[idx[n]] = 7
      } else if (kfm.stipe.ave$year[idx[n]] == 1986) {
        kfm.stipe.ave$time[idx[n]] = 8
      } else if (kfm.stipe.ave$year[idx[n]] == 1987) {
        kfm.stipe.ave$time[idx[n]] = 9
      } else if (kfm.stipe.ave$year[idx[n]] == 1988) {
        kfm.stipe.ave$time[idx[n]] = 10
      } else if (kfm.stipe.ave$year[idx[n]] == 1989) {
        kfm.stipe.ave$time[idx[n]] = 11
      } else if (kfm.stipe.ave$year[idx[n]] == 1990) {
        kfm.stipe.ave$time[idx[n]] = 12
      } else if (kfm.stipe.ave$year[idx[n]] == 1991) {
        kfm.stipe.ave$time[idx[n]] = 13
      } else if (kfm.stipe.ave$year[idx[n]] == 1992) {
        kfm.stipe.ave$time[idx[n]] = 14
      } else if (kfm.stipe.ave$year[idx[n]] == 1993) {
        kfm.stipe.ave$time[idx[n]] = 15
      } else if (kfm.stipe.ave$year[idx[n]] == 1994) {
        kfm.stipe.ave$time[idx[n]] = 16
      } else if (kfm.stipe.ave$year[idx[n]] == 1995) {
        kfm.stipe.ave$time[idx[n]] = 17
      } else if (kfm.stipe.ave$year[idx[n]] == 1996) {
        kfm.stipe.ave$time[idx[n]] = 18
      } else if (kfm.stipe.ave$year[idx[n]] == 1997) {
        kfm.stipe.ave$time[idx[n]] = 19
      } else if (kfm.stipe.ave$year[idx[n]] == 1998) {
        kfm.stipe.ave$time[idx[n]] = 20
      } else if (kfm.stipe.ave$year[idx[n]] == 1999) {
        kfm.stipe.ave$time[idx[n]] = 21
      } else if (kfm.stipe.ave$year[idx[n]] == 2000) {
        kfm.stipe.ave$time[idx[n]] = 22
      } else if (kfm.stipe.ave$year[idx[n]] == 2001) {
        kfm.stipe.ave$time[idx[n]] = 23
      } else if (kfm.stipe.ave$year[idx[n]] == 2002) {
        kfm.stipe.ave$time[idx[n]] = 24
      } else if (kfm.stipe.ave$year[idx[n]] == 2003) {
        kfm.stipe.ave$time[idx[n]] = 25
      } else if (kfm.stipe.ave$year[idx[n]] == 2004) {
        kfm.stipe.ave$time[idx[n]] = 26
      } else if (kfm.stipe.ave$year[idx[n]] == 2005) {
        kfm.stipe.ave$time[idx[n]] = 27
      } else if (kfm.stipe.ave$year[idx[n]] == 2006) {
        kfm.stipe.ave$time[idx[n]] = 28
      } else if (kfm.stipe.ave$year[idx[n]] == 2007) {
        kfm.stipe.ave$time[idx[n]] = 29
      } else if (kfm.stipe.ave$year[idx[n]] == 2008) {
        kfm.stipe.ave$time[idx[n]] = 30
      } else if (kfm.stipe.ave$year[idx[n]] == 2009) {
        kfm.stipe.ave$time[idx[n]] = 31
      } else if (kfm.stipe.ave$year[idx[n]] == 2010) {
        kfm.stipe.ave$time[idx[n]] = 32
      } else if (kfm.stipe.ave$year[idx[n]] == 2011) {
        kfm.stipe.ave$time[idx[n]] = 33
      } else if (kfm.stipe.ave$year[idx[n]] == 2012) {
        kfm.stipe.ave$time[idx[n]] = 34
      } else if (kfm.stipe.ave$year[idx[n]] == 2013) {
        kfm.stipe.ave$time[idx[n]] = 35
      } else if (kfm.stipe.ave$year[idx[n]] == 2014) {
        kfm.stipe.ave$time[idx[n]] = 36
      } else if (kfm.stipe.ave$year[idx[n]] == 2015) {
        kfm.stipe.ave$time[idx[n]] = 37
      } else if (kfm.stipe.ave$year[idx[n]] == 2016) {
        kfm.stipe.ave$time[idx[n]] = 38
      } else if (kfm.stipe.ave$year[idx[n]] == 2017) {
        kfm.stipe.ave$time[idx[n]] = 39
      } else if (kfm.stipe.ave$year[idx[n]] == 2018) {
        kfm.stipe.ave$time[idx[n]] = 40
      } else if (kfm.stipe.ave$year[idx[n]] == 2019) {
        kfm.stipe.ave$time[idx[n]] = 41
      } else if (kfm.stipe.ave$year[idx[n]] == 2020) {
        kfm.stipe.ave$time[idx[n]] = 42
      } else if (kfm.stipe.ave$year[idx[n]] == 2021) {
        kfm.stipe.ave$time[idx[n]] = 43
      } else if (kfm.stipe.ave$year[idx[n]] == 2022) {
        kfm.stipe.ave$time[idx[n]] = 44
      } else if (kfm.stipe.ave$year[idx[n]] == 2023) {
        kfm.stipe.ave$time[idx[n]] = 45
      }}}
  else {}
}

########Calculating biomass of macro stipes from size freq distributions#######################

MACRO <- subset(kfm.sum, taxon_name == "Macrocystis pyrifera")
#Remove Old Anacapa MPA and the SMCA since they don't meet criteria for inclusion
MACRO <- subset(MACRO, CA_MPA_Name_Short != "Anacapa Island SMCA" & 
                  CA_MPA_Name_Short != "Anacapa Island SMR 1978")

u <- unique(MACRO[,c('year','site_id','CA_MPA_Name_Short','status','area','taxon_name')])

#Create empty dataframes
Macro.site <- data.frame(site = NA, CA_MPA_Name_Short = NA,  year = NA,  transect = NA, biomass = NA, count = NA, y = NA, area = NA, ind = NA)
Meta <- data.frame(site = NA, CA_MPA_Name_Short = NA,  year = NA,  transect = NA, biomass = NA, count = NA, y = NA, area = NA, ind = NA)

kfm.stipe.OG <- kfm.stipe.site

bio_macro <- function(stipe) {
  stipe * aveslope * 1000
}

for(i in 1:nrow(u)) {
  t <- which(MACRO$site_id == u$site_id[i] &  MACRO$year == u$year[i] & MACRO$area == u$area[i] & MACRO$taxon_name == u$taxon_name[i])
  n <- sum(MACRO$count[t])
  t2 <- which(kfm.stipe.site$site_id == u$site_id[i] &
                kfm.stipe.site$SurveyYear == u$year[i] & kfm.stipe.site$ScientificName == u$taxon_name[i])
  
  
  if (n != 0 & length(t2) != 0) { #This means an urchin was observed on the transect and in the size freq table
    a = NA
    s = data.frame(matrix(NA, nrow = 1000, ncol = n))
    a = kfm.stipe.site[t2,]
    for(k in 1:1000){
      s[k,] <- sample(a$Stipe_Count,n, replace = TRUE) #randomly sample the kelp "population"
    }
    s_bio <- s %>%
      mutate(across(starts_with("X"), ~bio_macro(.)))
    s_bio <- s_bio %>% 
      mutate(sum = rowSums(.))
    s <- s %>% 
      mutate(sum = rowSums(.))
    aveB <- mean(s_bio$sum)
    aveD <- mean(s$sum)
    
    Meta$site <- u$site_id[i]
    Meta$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
    Meta$year <- u$year[i]
    Meta$transect <- length(t)
    Meta$area <- u$area[i]
    Meta$y <- u$taxon_name[i]
    Meta$count <- aveD
    Meta$biomass <- aveB
    Meta$ind <- n
    Macro.site <- rbind(Macro.site,Meta)
    
  } else if (n == 0 & length(t2) == 0){#This means there were no urchins on transect and none in the size freq data either
    Meta$site <- u$site_id[i]
    Meta$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
    Meta$year <- u$year[i]
    Meta$transect <- length(t)
    Meta$area <- u$area[i]
    Meta$y <- u$taxon_name[i]
    Meta$count <- 0
    Meta$biomass <- 0
    Meta$ind <- 0
    Macro.site <- rbind(Macro.site,Meta)
    
  } else if (n != 0 & length(t2) == 0){ #If there are urchins on transects, but they weren't sized we will mark size as NA
    Meta$site <- u$site_id[i]
    Meta$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
    Meta$year <- u$year[i]
    Meta$transect <- length(t)
    Meta$area <- u$area[i]
    Meta$y <- u$taxon_name[i]
    Meta$count <- NA
    Meta$biomass <- NA
    Meta$ind <- n
    Macro.site <- rbind(Macro.site,Meta)
    
  } else {
    Meta$site <- u$site_id[i]
    Meta$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
    Meta$year <- u$year[i]
    Meta$transect <- length(t)
    Meta$area <- u$area[i]
    Meta$y <- u$taxon_name[i]
    Meta$count <- 0
    Meta$biomass <- 0
    Meta$ind <- n
    Macro.site <- rbind(Macro.site,Meta)
  }
}

##########################################################################################
KFM.Macro.site.merge <- merge(Macro.site, sites.short, by.x = c("site"), by.y = c("site"),
                              all.x = TRUE)
KFM.Macro.site.merge <- KFM.Macro.site.merge[, colnames(KFM.Macro.site.merge)[c(2,1,3,12,17,18,4:9)]]

#Convert to m2 by dividing by # transects and area of quads
KFM.Macro.site.merge$biomassCorr <- KFM.Macro.site.merge$biomass/KFM.Macro.site.merge$area/KFM.Macro.site.merge$transect
KFM.Macro.site.merge$densityCorr <- KFM.Macro.site.merge$count/KFM.Macro.site.merge$area/KFM.Macro.site.merge$transect
KFM.Macro.site.merge$densityIndCorr <- KFM.Macro.site.merge$ind/KFM.Macro.site.merge$area/KFM.Macro.site.merge$transect
#######################################################################################

#############################################################################################
####Wrangling urchin, macro KFM data###############################
KFM.Macro.site.all.sub <- KFM.Macro.site.merge[, colnames(KFM.Macro.site.merge)[c(1,2,5,6,3,7,11,9,8,10,14,13,15)]]
names(KFM.Macro.site.all.sub)[names(KFM.Macro.site.all.sub) == "biomass"] <- "biomassRaw"
names(KFM.Macro.site.all.sub)[names(KFM.Macro.site.all.sub) == "biomassCorr"] <- "biomass"
names(KFM.Macro.site.all.sub)[names(KFM.Macro.site.all.sub) == "densityCorr"] <- "density"
names(KFM.Macro.site.all.sub)[names(KFM.Macro.site.all.sub) == "CA_MPA_Name_Short.x"] <- "CA_MPA_Name_Short"

KFM.Urchin.site.all.sub <- KFM.Urchin.site.all[, colnames(KFM.Urchin.site.all)[c(1:6,10,7:9,13,14)]]
names(KFM.Urchin.site.all.sub)[names(KFM.Urchin.site.all.sub) == "Density"] <- "density"
names(KFM.Urchin.site.all.sub)[names(KFM.Urchin.site.all.sub) == "biomass"] <- "biomassRaw"
names(KFM.Urchin.site.all.sub)[names(KFM.Urchin.site.all.sub) == "bio.m2"] <- "biomass"
KFM.Urchin.site.all.sub$densityIndCorr <- NA

KFM.bio.site.all <- rbind(KFM.Macro.site.all.sub, KFM.Urchin.site.all.sub)

#After merging den = density of individuals, whereas density means density of ind for everything except kelp, 
#which is actually the density of stipes, biomass estimates have been divided by number of transects and area
KFM.allbio <- merge(kfm.ave, KFM.bio.site.all, by.x = c("site_id",
                      "status", "year", "taxon_name", "CA_MPA_Name_Short", "area"), 
                    by.y = c("site","site_status", "year", "y", "CA_MPA_Name_Short", "area"),
                    all.x = TRUE )

Swath.KFM <- KFM.allbio
names(Swath.KFM)[names(Swath.KFM) == "taxon_name"] <- "y"
names(Swath.KFM)[names(Swath.KFM) == "count"] <- "countRaw"
names(Swath.KFM)[names(Swath.KFM) == "density"] <- "den.y"
names(Swath.KFM)[names(Swath.KFM) == "site_id"] <- "site"

Swath.KFM.Corr <- Swath.KFM
#Make density the stipe density, not density of ind for macro, this is the comparison across programs
Swath.KFM.Corr$den <-  ifelse(Swath.KFM.Corr$y == "Macrocystis pyrifera", Swath.KFM.Corr$den.y, Swath.KFM.Corr$den)
  
kfm.ave.ave <- Swath.KFM.Corr %>% #Calculate the mean across sites
group_by(y, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, time, year, area) %>% 
summarise_at(c("den","biomass"), mean)

#Remove MPAs that don't meet our criteria for inclusion
kfm.ave.ave <- subset(kfm.ave.ave, CA_MPA_Name_Short != "Anacapa Island SMCA" & CA_MPA_Name_Short != "Anacapa Island SMR 1978"
                      & CA_MPA_Name_Short != "Carrington Point SMR"
                      & CA_MPA_Name_Short != "Painted Cave SMCA")
###############################################################################################

###############################################################################################
##Transform all data into proportions########################################################## 
######KFM density##############################################################################
All.den <- kfm.ave.ave
All.den <- subset(All.den, den != "N/A") #removing NA since there is no data from that year
All.den$Prop <- 0
All.den$PropCorr <- 0
x = unique(All.den$CA_MPA_Name_Short)
t = unique(All.den$y)
l = length(x)
p = length(t)

#Turn maximum annual density into proportion of time series maximum den
for(i in 1:l) {
  for(j in 1:p) {
    idx <-which(x[i] == All.den$CA_MPA_Name_Short & t[j] == All.den$y)
    m <- max(All.den$den[idx]) #Find the max of the mpa specific time series to use to create a proportion of the max
    All.den$Prop[idx] <- All.den$den[idx]/m #Calculate proportion of max for every data point in the time series
    All.den$PropCorr[idx] <- All.den$Prop[idx] + 0.01 #Add 1% so that there are no zeroes in the dataset
  }}


# Create short form to easily line up and know whether we have data in both mpa and reference for a given year and calc diff
All.den.sub <- All.den[, colnames(All.den)[c(3,7,1,2,8,12)]] #First subset only columns we need
Short.den <- All.den.sub %>% #Create short format
  spread(status, PropCorr)
Short.den.diff <- Short.den[complete.cases(Short.den), ] #Remove data if there are not both mpa and reference
Short.den.diff$Diff <- Short.den.diff$mpa/Short.den.diff$reference #Calculate the response ratio
Short.den.diff$lnDiff <- log(Short.den.diff$Diff) #Log response ratio
Short.den.diff$resp <- "Den"

######KFM biomass##############################################################################
All.bio <- kfm.ave.ave
All.bio <- subset(All.bio, biomass != "N/A") #removing NA since there is no data from that year
All.bio$Prop <- 0
All.bio$PropCorr <- 0
x = unique(All.bio$CA_MPA_Name_Short)
t = unique(All.bio$y)
l = length(x)
p = length(t)

#Turn maximum annual biomass into proportion of time series maximum biomass
for(i in 1:l) {
  for(j in 1:p) {
    idx <-which(x[i] == All.bio$CA_MPA_Name_Short & t[j] == All.bio$y)
    m <- max(All.bio$biomass[idx]) #Find the max of the mpa specific time series to use to create a proportion of the max
    All.bio$Prop[idx] <- All.bio$biomass[idx]/m #Calculate proportion of max for every data point in the time series
    All.bio$PropCorr[idx] <- All.bio$Prop[idx] + 0.01 #Add 1% so that there are no zeroes in the dataset
  }}


# Create short form to easily line up and know whether we have data in both mpa and reference for a given year and calc diff
All.bio.sub <- All.bio[, colnames(All.bio)[c(3,7,1,2,8,12)]] #First subset only columns we need
Short.bio <- All.bio.sub %>% #Create short format
  spread(status, PropCorr)
Short.bio.diff <- Short.bio[complete.cases(Short.bio), ] #Remove data if there are not both mpa and reference
Short.bio.diff$Diff <- Short.bio.diff$mpa/Short.bio.diff$reference #Calculate the response ratio
Short.bio.diff$lnDiff <- log(Short.bio.diff$Diff) #Log response ratio
Short.bio.diff$resp <- "Bio"

KFM.join.ave <- rbind(Short.den.diff, Short.bio.diff)
#####################################################################################################

#####################################################################################################
#Create dataframe with all average density/biomass inside/outside mpas KFM macro, urchins, lob#######

#Only want 10m area for Macro, though 2m appears similar
kfm.edit.den <- subset(kfm.ave.ave, y != "Macrocystis pyrifera" | area != 2)
kfm.edit.den <- subset(kfm.edit.den, y != "Macrocystis pyrifera" | area != 1)

kfm.edit.den <- kfm.edit.den[, colnames(kfm.edit.den)[c(2,3,5,7,9,1)]]
kfm.edit.den <- kfm.edit.den[complete.cases(kfm.edit.den), ] #Delete any rows where there isn't an mpa and reference
KFM.den <- kfm.edit.den %>% 
  spread(status, den)
KFM.den$source <- "KFM"
KFM.den <- KFM.den[complete.cases(KFM.den), ] #Delete any rows where there isn't an mpa and reference
KFM.den.long <- gather(KFM.den, status, value, 5:6)
KFM.den.long$resp <- 'Den'

kfm.edit.bio <- subset(kfm.ave.ave, y != "Macrocystis pyrifera" | area != 2)
kfm.edit.bio <- subset(kfm.edit.bio, y != "Macrocystis pyrifera" | area != 1)

kfm.edit.bio <- kfm.edit.bio[, colnames(kfm.edit.bio)[c(2,3,5,7,10,1)]]
kfm.edit.bio <- kfm.edit.bio[complete.cases(kfm.edit.bio), ] #Delete any rows where there isn't an mpa and reference
KFM.bio <- kfm.edit.bio %>% 
  spread(status, biomass)
KFM.bio$source <- "KFM"
KFM.bio <- KFM.bio[complete.cases(KFM.bio), ] #Delete any rows where there isn't an mpa and reference
KFM.bio.long <- gather(KFM.bio, status, value, 5:6)
KFM.bio.long$resp <- 'Bio'

KFM.resp <- rbind(KFM.den.long, KFM.bio.long)
#################################################################################

#################################################################################
##############KFM data wrangling#################################################
KFM.join.ave$time <- 0
x <- unique(KFM.join.ave$CA_MPA_Name_Short) # to redo since sites have been removed
l = length(x)
for(i in 1:l) {
  idx <-which(x[i] == KFM.join.ave$CA_MPA_Name_Short)
  j <-which(x[i] == Sites2$CA_MPA_Name_Short)
  m <- length(idx)
  if (Sites2$MPA_Start[j] == 2003) {
    for (n in 1:m) {
      if (KFM.join.ave$year[idx[n]] <=2003) {
        KFM.join.ave$time[idx[n]] = 0
      } else if (KFM.join.ave$year[idx[n]] == 2004) {
        KFM.join.ave$time[idx[n]] = 1
      } else if (KFM.join.ave$year[idx[n]] == 2005) {
        KFM.join.ave$time[idx[n]] = 2
      } else if (KFM.join.ave$year[idx[n]] == 2006) {
        KFM.join.ave$time[idx[n]] = 3
      } else if (KFM.join.ave$year[idx[n]] == 2007) {
        KFM.join.ave$time[idx[n]] = 4
      } else if (KFM.join.ave$year[idx[n]] == 2008) {
        KFM.join.ave$time[idx[n]] = 5
      } else if (KFM.join.ave$year[idx[n]] == 2009) {
        KFM.join.ave$time[idx[n]] = 6
      } else if (KFM.join.ave$year[idx[n]] == 2010) {
        KFM.join.ave$time[idx[n]] = 7
      } else if (KFM.join.ave$year[idx[n]] == 2011) {
        KFM.join.ave$time[idx[n]] = 8
      } else if (KFM.join.ave$year[idx[n]] == 2012) {
        KFM.join.ave$time[idx[n]] = 9
      } else if (KFM.join.ave$year[idx[n]] == 2013) {
        KFM.join.ave$time[idx[n]] = 10
      } else if (KFM.join.ave$year[idx[n]] == 2014) {
        KFM.join.ave$time[idx[n]] = 11
      } else if (KFM.join.ave$year[idx[n]] == 2015) {
        KFM.join.ave$time[idx[n]] = 12
      } else if (KFM.join.ave$year[idx[n]] == 2016) {
        KFM.join.ave$time[idx[n]] = 13
      } else if (KFM.join.ave$year[idx[n]] == 2017) {
        KFM.join.ave$time[idx[n]] = 14
      } else if (KFM.join.ave$year[idx[n]] == 2018) {
        KFM.join.ave$time[idx[n]] = 15
      } else if (KFM.join.ave$year[idx[n]] == 2019) {
        KFM.join.ave$time[idx[n]] = 16
      } else if (KFM.join.ave$year[idx[n]] == 2020) {
        KFM.join.ave$time[idx[n]] = 17
      } else if (KFM.join.ave$year[idx[n]] == 2021) {
        KFM.join.ave$time[idx[n]] = 18
      } else if (KFM.join.ave$year[idx[n]] == 2022) {
        KFM.join.ave$time[idx[n]] = 19
      } else if (KFM.join.ave$year[idx[n]] == 2023) {
        KFM.join.ave$time[idx[n]] = 20
      }}
    
  } else if (Sites2$MPA_Start[j] == 1978) {
    for (n in 1:m) {
      if (KFM.join.ave$year[idx[n]] <= 1978) {
        KFM.join.ave$time[idx[n]] = 0
      } else if (KFM.join.ave$year[idx[n]] == 1982) {
        KFM.join.ave$time[idx[n]] = 4
      } else if (KFM.join.ave$year[idx[n]] == 1983) {
        KFM.join.ave$time[idx[n]] = 5
      } else if (KFM.join.ave$year[idx[n]] == 1984) {
        KFM.join.ave$time[idx[n]] = 6
      } else if (KFM.join.ave$year[idx[n]] == 1985) {
        KFM.join.ave$time[idx[n]] = 7
      } else if (KFM.join.ave$year[idx[n]] == 1986) {
        KFM.join.ave$time[idx[n]] = 8
      } else if (KFM.join.ave$year[idx[n]] == 1987) {
        KFM.join.ave$time[idx[n]] = 9
      } else if (KFM.join.ave$year[idx[n]] == 1988) {
        KFM.join.ave$time[idx[n]] = 10
      } else if (KFM.join.ave$year[idx[n]] == 1989) {
        KFM.join.ave$time[idx[n]] = 11
      } else if (KFM.join.ave$year[idx[n]] == 1990) {
        KFM.join.ave$time[idx[n]] = 12
      } else if (KFM.join.ave$year[idx[n]] == 1991) {
        KFM.join.ave$time[idx[n]] = 13
      } else if (KFM.join.ave$year[idx[n]] == 1992) {
        KFM.join.ave$time[idx[n]] = 14
      } else if (KFM.join.ave$year[idx[n]] == 1993) {
        KFM.join.ave$time[idx[n]] = 15
      } else if (KFM.join.ave$year[idx[n]] == 1994) {
        KFM.join.ave$time[idx[n]] = 16
      } else if (KFM.join.ave$year[idx[n]] == 1995) {
        KFM.join.ave$time[idx[n]] = 17
      } else if (KFM.join.ave$year[idx[n]] == 1996) {
        KFM.join.ave$time[idx[n]] = 18
      } else if (KFM.join.ave$year[idx[n]] == 1997) {
        KFM.join.ave$time[idx[n]] = 19
      } else if (KFM.join.ave$year[idx[n]] == 1998) {
        KFM.join.ave$time[idx[n]] = 20
      } else if (KFM.join.ave$year[idx[n]] == 1999) {
        KFM.join.ave$time[idx[n]] = 21
      } else if (KFM.join.ave$year[idx[n]] == 2000) {
        KFM.join.ave$time[idx[n]] = 22
      } else if (KFM.join.ave$year[idx[n]] == 2001) {
        KFM.join.ave$time[idx[n]] = 23
      } else if (KFM.join.ave$year[idx[n]] == 2002) {
        KFM.join.ave$time[idx[n]] = 24
      } else if (KFM.join.ave$year[idx[n]] == 2003) {
        KFM.join.ave$time[idx[n]] = 25
      } else if (KFM.join.ave$year[idx[n]] == 2004) {
        KFM.join.ave$time[idx[n]] = 26
      } else if (KFM.join.ave$year[idx[n]] == 2005) {
        KFM.join.ave$time[idx[n]] = 27
      } else if (KFM.join.ave$year[idx[n]] == 2006) {
        KFM.join.ave$time[idx[n]] = 28
      } else if (KFM.join.ave$year[idx[n]] == 2007) {
        KFM.join.ave$time[idx[n]] = 29
      } else if (KFM.join.ave$year[idx[n]] == 2008) {
        KFM.join.ave$time[idx[n]] = 30
      } else if (KFM.join.ave$year[idx[n]] == 2009) {
        KFM.join.ave$time[idx[n]] = 31
      } else if (KFM.join.ave$year[idx[n]] == 2010) {
        KFM.join.ave$time[idx[n]] = 32
      } else if (KFM.join.ave$year[idx[n]] == 2011) {
        KFM.join.ave$time[idx[n]] = 33
      } else if (KFM.join.ave$year[idx[n]] == 2012) {
        KFM.join.ave$time[idx[n]] = 34
      } else if (KFM.join.ave$year[idx[n]] == 2013) {
        KFM.join.ave$time[idx[n]] = 35
      } else if (KFM.join.ave$year[idx[n]] == 2014) {
        KFM.join.ave$time[idx[n]] = 36
      } else if (KFM.join.ave$year[idx[n]] == 2015) {
        KFM.join.ave$time[idx[n]] = 37
      } else if (KFM.join.ave$year[idx[n]] == 2016) {
        KFM.join.ave$time[idx[n]] = 38
      } else if (KFM.join.ave$year[idx[n]] == 2017) {
        KFM.join.ave$time[idx[n]] = 39
      } else if (KFM.join.ave$year[idx[n]] == 2018) {
        KFM.join.ave$time[idx[n]] = 40
      } else if (KFM.join.ave$year[idx[n]] == 2019) {
        KFM.join.ave$time[idx[n]] = 41
      } else if (KFM.join.ave$year[idx[n]] == 2020) {
        KFM.join.ave$time[idx[n]] = 42
      } else if (KFM.join.ave$year[idx[n]] == 2021) {
        KFM.join.ave$time[idx[n]] = 43
      } else if (KFM.join.ave$year[idx[n]] == 2022) {
        KFM.join.ave$time[idx[n]] = 44
      } else if (KFM.join.ave$year[idx[n]] == 2023) {
        KFM.join.ave$time[idx[n]] = 45
      }}}
  else {}
}

KFM.join.ave <- merge(KFM.join.ave, sites.short.edit, by.x = c("CA_MPA_Name_Short"), 
                        by.y = c("CA_MPA_Name_Short"),
                        all = TRUE)
KFM.join.ave <- KFM.join.ave[complete.cases(KFM.join.ave$year), ]
KFM.join.ave$source <- "KFM"
for(i in 1:l) {
  idx <-which(x[i] == KFM.join.ave$CA_MPA_Name_Short)
  j <-which(x[i] == Sites2$CA_MPA_Name_Short)
  m <- length(idx)
if (Sites2$MPA_Start[j] == 2003) {
  for (n in 1:m) {
    if (KFM.join.ave$year[idx[n]] <=2003) {
      KFM.join.ave$BA[idx[n]] = "Before"
    } else if (KFM.join.ave$year[idx[n]] == 2004) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2005) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2006) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2007) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2008) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2009) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2010) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2011) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2012) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2013) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2014) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2015) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2016) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2017) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2018) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2019) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2020) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2021) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2022) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2023) {
      KFM.join.ave$BA[idx[n]] = "After"
    }}
  
} else if (Sites2$MPA_Start[j] == 1978) {
  for (n in 1:m) {
    if (KFM.join.ave$year[idx[n]] <= 1978) {
      KFM.join.ave$BA[idx[n]] = "Before"
    } else if (KFM.join.ave$year[idx[n]] == 1982) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1983) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1984) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1985) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1986) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1987) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1988) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1989) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1990) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1991) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1992) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1993) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1994) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1995) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1996) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1997) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1998) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 1999) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2000) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2001) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2002) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2003) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2004) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2005) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2006) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2007) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2008) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2009) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2010) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2011) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2012) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2013) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2014) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2015) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2016) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2017) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2018) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2019) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2020) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2021) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2022) {
      KFM.join.ave$BA[idx[n]] = "After"
    } else if (KFM.join.ave$year[idx[n]] == 2023) {
      KFM.join.ave$BA[idx[n]] = "After"
    }}}
else {}
}

#################################################################################
##Now lets analyze the LTER data##################################################################################
##Join average swath data to site table
lter.site <- merge(lter, Sites2, by = c("site_id"))
lter.site <- subset(lter.site, CA_MPA_Name_Short != "" & site_id != "a-l-02") #This is arroyo honda where kelp is essentially absent, Bob Miller said to remove not a good reference
lter.site$time <- 0
lter.site$year <- year(lter.site$date)
x <- unique(lter.site$CA_MPA_Name_Short) # to redo since sites have been removed
l = length(x)
for(i in 1:l) {
  idx <-which(x[i] == lter.site$CA_MPA_Name_Short)
  j <-which(x[i] == Sites2$CA_MPA_Name_Short)
  m <- length(idx)
  if (Sites2$MPA_Start[j] == 2003) {
    for (n in 1:m) {
      if (lter.site$year[idx[n]] <=2003) {
        lter.site$time[idx[n]] = 0
      } else if (lter.site$year[idx[n]] == 2004) {
        lter.site$time[idx[n]] = 1
      } else if (lter.site$year[idx[n]] == 2005) {
        lter.site$time[idx[n]] = 2
      } else if (lter.site$year[idx[n]] == 2006) {
        lter.site$time[idx[n]] = 3
      } else if (lter.site$year[idx[n]] == 2007) {
        lter.site$time[idx[n]] = 4
      } else if (lter.site$year[idx[n]] == 2008) {
        lter.site$time[idx[n]] = 5
      } else if (lter.site$year[idx[n]] == 2009) {
        lter.site$time[idx[n]] = 6
      } else if (lter.site$year[idx[n]] == 2010) {
        lter.site$time[idx[n]] = 7
      } else if (lter.site$year[idx[n]] == 2011) {
        lter.site$time[idx[n]] = 8
      } else if (lter.site$year[idx[n]] == 2012) {
        lter.site$time[idx[n]] = 9
      } else if (lter.site$year[idx[n]] == 2013) {
        lter.site$time[idx[n]] = 10
      } else if (lter.site$year[idx[n]] == 2014) {
        lter.site$time[idx[n]] = 11
      } else if (lter.site$year[idx[n]] == 2015) {
        lter.site$time[idx[n]] = 12
      } else if (lter.site$year[idx[n]] == 2016) {
        lter.site$time[idx[n]] = 13
      } else if (lter.site$year[idx[n]] == 2017) {
        lter.site$time[idx[n]] = 14
      } else if (lter.site$year[idx[n]] == 2018) {
        lter.site$time[idx[n]] = 15
      } else if (lter.site$year[idx[n]] == 2019) {
        lter.site$time[idx[n]] = 16
      } else if (lter.site$year[idx[n]] == 2020) {
        lter.site$time[idx[n]] = 17
      } else if (lter.site$year[idx[n]] == 2021) {
        lter.site$time[idx[n]] = 18
      } else if (lter.site$year[idx[n]] == 2022) {
        lter.site$time[idx[n]] = 19
      } else if (lter.site$year[idx[n]] == 2023) {
        lter.site$time[idx[n]] = 20
      }}
    
  } else if (Sites2$MPA_Start[j] == 2012) {
    for (n in 1:m) {
      if (lter.site$year[idx[n]] <= 2012) {
        lter.site$time[idx[n]] = 0
      } else if (lter.site$year[idx[n]] == 2013) {
        lter.site$time[idx[n]] = 1
      } else if (lter.site$year[idx[n]] == 2014) {
        lter.site$time[idx[n]] = 2
      } else if (lter.site$year[idx[n]] == 2015) {
        lter.site$time[idx[n]] = 3
      } else if (lter.site$year[idx[n]] == 2016) {
        lter.site$time[idx[n]] = 4
      } else if (lter.site$year[idx[n]] == 2017) {
        lter.site$time[idx[n]] = 5
      } else if (lter.site$year[idx[n]] == 2018) {
        lter.site$time[idx[n]] = 6
      } else if (lter.site$year[idx[n]] == 2019) {
        lter.site$time[idx[n]] = 7
      } else if (lter.site$year[idx[n]] == 2020) {
        lter.site$time[idx[n]] = 8
      } else if (lter.site$year[idx[n]] == 2021) {
        lter.site$time[idx[n]] = 9
      } else if (lter.site$year[idx[n]] == 2022) {
        lter.site$time[idx[n]] = 10
      } else if (lter.site$year[idx[n]] == 2023) {
        lter.site$time[idx[n]] = 11
      }}}
  else {}
}

#Subset only Campus Point and Naples
lter.site <- subset(lter.site, CA_MPA_Name_Short == "Campus Point SMCA" | CA_MPA_Name_Short == "Naples SMCA")
x <- unique(lter.site$CA_MPA_Name_Short) # to redo since sites have been removed
MPA_implement <- filter(Sites2, CA_MPA_Name_Short %in% x)
MPA_implement.sub <- filter(MPA_implement, CA_MPA_Name_Short %in% x)

lter.sub <- subset(lter.site, taxon_name == "Strongylocentrotus purpuratus" & sample_method == "quad" | taxon_name == "Mesocentrotus franciscanus" & sample_method == "quad" |
                       taxon_name == "Macrocystis pyrifera" & sample_method == "swath")
lter.sub$den <- lter.sub$count/lter.sub$area

#Calculate site means
lter.ave <- lter.sub %>%
  group_by(site_id,  transect_id, taxon_name, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, time, year) %>% 
  summarise_at(c("den"), mean)

lter.ave <- lter.ave %>%
  group_by(site_id,  taxon_name, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, time, year) %>% 
  summarise_at(c("den"), mean)

###########################################################################################
##########################################################################################
########Calculating biomass of urchins from size freq distributions#######################
URCHINS <- subset(lter.site, taxon_name == "Strongylocentrotus purpuratus" | taxon_name == "Mesocentrotus franciscanus")
URCHINS <- mutate(URCHINS, taxon_name = case_when(
  taxon_name == "Strongylocentrotus purpuratus"  ~ "STRPURAD",
  taxon_name == "Mesocentrotus franciscanus"    ~ "MESFRAAD",
  TRUE ~ taxon_name))

URCHINS$replicate_id <- as.numeric(URCHINS$replicate_id)

u <- unique(URCHINS[,c('year','site_id','CA_MPA_Name_Short','status','area','taxon_name','transect_id')])

#Create empty dataframes
Urchin.site <- data.frame(site = NA, CA_MPA_Name_Short = NA, site_status = NA, year = NA,  transect = NA, quad = NA, biomass = NA, count = NA, y = NA, area = NA)
Urchin.SF <- data.frame(site = NA, CA_MPA_Name_Short = NA, site_status = NA, year = NA,  transect = NA, quad = NA, biomass = NA, count = NA, y = NA, area = NA)

SizeFreq.Urch <- SizeFreq.Urch.OG
SizeFreq.Urch <- SizeFreq.Urch[complete.cases(SizeFreq.Urch$size), ]

#Select site/year/zone if there are lobsters then look up same site/year/zone from size freq table since this is where all size data are
#create lobster "population" from size freq table and draw with replacement the number of lobsters from swath data table 
#if lobsters were measured on transect; however, if no lobsters then equals zero, if there were lobsters but no size freq
#these will be removed and marked as NA


for(i in 1:nrow(u)) {
  t <- which(URCHINS$site_id == u$site_id[i]  & URCHINS$year == u$year[i]  & URCHINS$transect_id == u$transect_id[i] &
               URCHINS$taxon_name == u$taxon_name[i])
  n <- sum(URCHINS$count[t])
  t2 <- which(SizeFreq.Urch$CA_MPA_Name_Short == u$CA_MPA_Name_Short[i] & SizeFreq.Urch$site_status == u$status[i] &
                SizeFreq.Urch$year == u$year[i] & SizeFreq.Urch$classcode == u$taxon_name[i])
  
  if (n != 0 & length(t2) != 0) { #This means an urchin was observed on the transect and in the size freq table
    a = NA
    s = data.frame(matrix(NA, nrow = 1000, ncol = n))
    l <- length(t2)
    for(j in 1:l){
      list <- as.data.frame(rep(SizeFreq.Urch$size[t2[j]], SizeFreq.Urch$count[t2[j]] ))
      a <- rbind(a,list)
      a <- na.omit(a) }
    a <- as.numeric(as.character(unlist(a))) #change list to numeric so that it can be sampled from, this is each urchin size 
    #randomly sample the urchin "population" with replacement; do this 1000 times
    for(k in 1:1000){
      s[k,] <- sample(a,n, replace = TRUE) 
    }
    #Convert CL to biomass for all individuals
    if(u$taxon_name[i] == "MESFRAAD"){
      s_bio <- s %>% 
        mutate(across(starts_with("X"), ~bio_redurch(.)))
      s_bio <- s_bio %>% 
        mutate(sum = rowSums(.))
      ave <- mean(s_bio$sum)
    }else {
      s_bio <- s %>% 
        mutate(across(starts_with("X"), ~bio_purpurch(.)))
      s_bio <- s_bio %>% 
        mutate(sum = rowSums(.))
      ave <- mean(s_bio$sum)
    }
    Urchin.SF$site <- u$site[i]
    Urchin.SF$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
    Urchin.SF$year <- u$year[i]
    Urchin.SF$transect <- u$transect_id[i]
    Urchin.SF$quad <- length(t)
    Urchin.SF$site_status <- u$status[i]
    Urchin.SF$area <- u$area[i]
    Urchin.SF$y <- u$taxon_name[i]
    Urchin.SF$biomass <- ave
    Urchin.SF$count <- n
    Urchin.site <- rbind(Urchin.site,Urchin.SF)
    
  } else if (n == 0 & length(t2) == 0){#This means there were no urchins on transect and none in the size freq data either
    Urchin.SF$site <- u$site[i]
    Urchin.SF$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
    Urchin.SF$year <- u$year[i]
    Urchin.SF$transect <- u$transect_id[i]
    Urchin.SF$quad <- length(t)
    Urchin.SF$site_status <- u$status[i]
    Urchin.SF$area <- u$area[i]
    Urchin.SF$y <- u$taxon_name[i]
    Urchin.SF$biomass <- 0
    Urchin.SF$count <- 0
    Urchin.site <- rbind(Urchin.site,Urchin.SF)
    
  } else if (n != 0 & length(t2) == 0){ #If there are urchins on transects, but they weren't sized we will mark size as NA
    Urchin.SF$site <- u$site[i]
    Urchin.SF$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
    Urchin.SF$year <- u$year[i]
    Urchin.SF$transect <- u$transect_id[i]
    Urchin.SF$quad <- length(t)
    Urchin.SF$site_status <- u$status[i]
    Urchin.SF$area <- u$area[i]
    Urchin.SF$y <- u$taxon_name[i]
    Urchin.SF$biomass <- NA
    Urchin.SF$count <- n
    Urchin.site <- rbind(Urchin.site,Urchin.SF)
    
  } else {
    Urchin.SF$site <- u$site[i]
    Urchin.SF$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
    Urchin.SF$year <- u$year[i]
    Urchin.SF$transect <- u$transect_id[i]
    Urchin.SF$quad <- length(t)
    Urchin.SF$site_status <- u$status[i]
    Urchin.SF$area <- u$area[i]
    Urchin.SF$y <- u$taxon_name[i]
    Urchin.SF$biomass <- 0
    Urchin.SF$count <- 0
    Urchin.site <- rbind(Urchin.site,Urchin.SF)
  }
}

####################################################################################
LTER.Urchin.site.merge <- merge(Urchin.site, sites.short, by.x = c("site", "site_status", "CA_MPA_Name_Short"), by.y = c("site", "site_status", "CA_MPA_Name_Short"),
                               all.x = TRUE)

LTER.Urchin.site.merge <- LTER.Urchin.site.merge[, colnames(LTER.Urchin.site.merge)[c(3,1,2,4,17,20,5:10)]]
LTER.Urchin.site.merge[is.na(LTER.Urchin.site.merge)] <- 0

#Find sites where urchin was found on transect but no size was recorded
LTER.Urchin.null <- subset(LTER.Urchin.site.merge, count > 0 & biomass == 0) 

LTER.Urchin.site.merge.sub <- LTER.Urchin.site.merge

#Remove sites where lobsters were seen but not measured since we can't get an accurate biomass estimate
for(i in 1:nrow(LTER.Urchin.null)) {
  LTER.Urchin.site.merge.sub <- subset(LTER.Urchin.site.merge.sub, site != LTER.Urchin.null$site[i] | year != LTER.Urchin.null$year[i])
}

LTER.Urchin.site.transectSum <- LTER.Urchin.site.merge.sub 

#Need to account for the different number of quads on transects at the sites (only occured once)
LTER.Urchin.site.transectSum$count <- LTER.Urchin.site.transectSum$count/LTER.Urchin.site.transectSum$quad
LTER.Urchin.site.transectSum$biomass <- LTER.Urchin.site.transectSum$biomass/LTER.Urchin.site.transectSum$quad

LTER.Urchin.site.all <- LTER.Urchin.site.transectSum %>% #averaging across transects
  group_by(year, site, CA_MPA_Name_Short, site_designation, site_status, BaselineRegion, y, area) %>% 
  summarise_at(c("count", "biomass"), mean)

LTER.Urchin.site.all$Density <- LTER.Urchin.site.all$count/LTER.Urchin.site.all$area

lter.ave <- mutate(lter.ave, taxon_name = case_when(
  taxon_name == "Strongylocentrotus purpuratus"  ~ "STRPURAD",
  taxon_name == "Mesocentrotus franciscanus"    ~ "MESFRAAD",
  TRUE ~ taxon_name))

LTER.allbio <- merge(lter.ave, LTER.Urchin.site.all, by.x = c("site_id", "status", "year", "taxon_name", "CA_MPA_Name_Short"), 
                    by.y = c("site", "site_status", "year", "y", "CA_MPA_Name_Short"), all = T)

Swath.LTER <- LTER.allbio
Swath.LTER <- Swath.LTER[, colnames(Swath.LTER)[c(1:5,12,6:10,13:15)]]
names(Swath.LTER)[names(Swath.LTER) == "Density"] <- "count.y"
names(Swath.LTER)[names(Swath.LTER) == "status"] <- "site_status"
names(Swath.LTER)[names(Swath.LTER) == "den"] <- "Density"
names(Swath.LTER)[names(Swath.LTER) == "site_id"] <- "site"
names(Swath.LTER)[names(Swath.LTER) == "taxon_name"] <- "y"

###########################################################################################
###########################################################################################
###########################################################################################
lter.ave.ave <- Swath.LTER %>% #Calculate the mean across sites
  group_by(y, site_status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, time, year) %>% 
  summarise_at(vars("Density", "biomass"), mean)

lter.ave.ave <- mutate(lter.ave.ave, y = case_when(
  y == "STRPURAD"  ~ "Strongylocentrotus purpuratus",
  y == "MESFRAAD"    ~ "Mesocentrotus franciscanus",
  TRUE ~ y))
names(lter.ave.ave)[names(lter.ave.ave) == "Density"] <- "den"

########LTER density of kelp individuals###################################################################
All.den <- lter.ave.ave
All.den <- All.den[complete.cases(All.den$den), ]
All.den$Prop <- 0
All.den$PropCorr <- 0
x = unique(All.den$CA_MPA_Name_Short)
t = unique(All.den$y)
l = length(x)
p = length(t)

#Turn maximum annual density into proportion of time series maximum den
for(i in 1:l) {
  for(j in 1:p) {
    idx <-which(x[i] == All.den$CA_MPA_Name_Short & t[j] == All.den$y)
    m <- max(All.den$den[idx]) #Find the max of the mpa specific time series to use to create a proportion of the max
    All.den$Prop[idx] <- All.den$den[idx]/m #Calculate proportion of max for every data point in the time series
    All.den$PropCorr[idx] <- All.den$Prop[idx] + 0.01 #Add 1% so that there are no zeroes in the dataset
  }}


# Create short form to easily line up and know whether we have data in both mpa and reference for a given year and calc diff
All.den.sub <- All.den[, colnames(All.den)[c(3,7,1,2,11)]] #First subset only columns we need
Short.den <- All.den.sub %>% #Create short format
  spread(site_status, PropCorr)
Short.den.diff <- Short.den[complete.cases(Short.den), ] #Remove data if there are not both mpa and reference
Short.den.diff$Diff <- Short.den.diff$mpa/Short.den.diff$reference #Calculate the response ratio
Short.den.diff$lnDiff <- log(Short.den.diff$Diff) #Log response ratio
Short.den.diff$resp <- "Den"

########LTER biomass of individuals###################################################################
All.bio <- lter.ave.ave
All.bio <- All.bio[complete.cases(All.bio$biomass), ]
All.bio$Prop <- 0
All.bio$PropCorr <- 0
x = unique(All.bio$CA_MPA_Name_Short)
t = unique(All.bio$y)
l = length(x)
p = length(t)

#Turn maximum annual density into proportion of time series maximum den
for(i in 1:l) {
  for(j in 1:p) {
    idx <-which(x[i] == All.bio$CA_MPA_Name_Short & t[j] == All.bio$y)
    m <- max(All.bio$biomass[idx]) #Find the max of the mpa specific time series to use to create a proportion of the max
    All.bio$Prop[idx] <- All.bio$biomass[idx]/m #Calculate proportion of max for every data point in the time series
    All.bio$PropCorr[idx] <- All.bio$Prop[idx] + 0.01 #Add 1% so that there are no zeroes in the dataset
  }}


# Create short form to easily line up and know whether we have data in both mpa and reference for a given year and calc diff
All.bio.sub <- All.bio[, colnames(All.bio)[c(3,7,1,2,11)]] #First subset only columns we need
Short.bio <- All.bio.sub %>% #Create short format
  spread(site_status, PropCorr)
Short.bio.diff <- Short.bio[complete.cases(Short.bio), ] #Remove data if there are not both mpa and reference
Short.bio.diff$Diff <- Short.bio.diff$mpa/Short.bio.diff$reference #Calculate the response ratio
Short.bio.diff$lnDiff <- log(Short.bio.diff$Diff) #Log response ratio
Short.bio.diff$resp <- "Bio"

LTER.join.ave <- rbind(Short.den.diff,Short.bio.diff)
#####################################################################################################


#####################################################################################################
#Create dataframe with all average density/biomass inside/outside mpas LTER kelp fronds##############
#####################################################################################################
lter.edit.den <- lter.ave.ave[, colnames(lter.ave.ave)[c(2,3,5,7,8,1)]]
lter.edit.den <- lter.edit.den[complete.cases(lter.edit.den), ] #Delete any rows where there isn't an mpa and reference
LTER.den <- lter.edit.den %>% 
  spread(site_status, den)
LTER.den <- subset(LTER.den, y != "Macrocystis pyrifera")
LTER.den$source <- "LTER"
LTER.den <- LTER.den[complete.cases(LTER.den), ] #Delete any rows where there isn't an mpa and reference
LTER.den.long <- gather(LTER.den, site_status, value, 5:6)
LTER.den.long$resp <- 'Den'

lter.edit.bio <- lter.ave.ave[, colnames(lter.ave.ave)[c(2,3,5,7,9,1)]]
lter.edit.bio <- lter.edit.bio[complete.cases(lter.edit.bio), ] #Delete any rows where there isn't an mpa and reference
LTER.bio <- lter.edit.bio %>% 
  spread(site_status, biomass)
LTER.bio$source <- "LTER"
LTER.bio <- LTER.bio[complete.cases(LTER.bio), ] #Delete any rows where there isn't an mpa and reference
LTER.bio.long <- gather(LTER.bio, site_status, value, 5:6)
LTER.bio.long$resp <- 'Bio'

LTER.resp <- rbind(LTER.bio.long, LTER.den.long)

#################################################################################
##############Wrangling of LTER data for analyses################################
LTER.join.ave$time <- 0
x <- unique(LTER.join.ave$CA_MPA_Name_Short) # to redo since sites have been removed
l = length(x)
for(i in 1:l) {
  idx <-which(x[i] == LTER.join.ave$CA_MPA_Name_Short)
  j <-which(x[i] == Sites2$CA_MPA_Name_Short)
  m <- length(idx)
  if (Sites2$MPA_Start[j] == 2012) {
    for (n in 1:m) {
      if (LTER.join.ave$year[idx[n]] <=2012) {
        LTER.join.ave$time[idx[n]] = 0
      } else if (LTER.join.ave$year[idx[n]] == 2013) {
        LTER.join.ave$time[idx[n]] = 1
      } else if (LTER.join.ave$year[idx[n]] == 2014) {
        LTER.join.ave$time[idx[n]] = 2
      } else if (LTER.join.ave$year[idx[n]] == 2015) {
        LTER.join.ave$time[idx[n]] = 3
      } else if (LTER.join.ave$year[idx[n]] == 2016) {
        LTER.join.ave$time[idx[n]] = 4
      } else if (LTER.join.ave$year[idx[n]] == 2017) {
        LTER.join.ave$time[idx[n]] = 5
      } else if (LTER.join.ave$year[idx[n]] == 2018) {
        LTER.join.ave$time[idx[n]] = 6
      } else if (LTER.join.ave$year[idx[n]] == 2019) {
        LTER.join.ave$time[idx[n]] = 7
      } else if (LTER.join.ave$year[idx[n]] == 2020) {
        LTER.join.ave$time[idx[n]] = 8
      } else if (LTER.join.ave$year[idx[n]] == 2021) {
        LTER.join.ave$time[idx[n]] = 9
      } else if (LTER.join.ave$year[idx[n]] == 2022) {
        LTER.join.ave$time[idx[n]] = 10
      } else if (LTER.join.ave$year[idx[n]] == 2023) {
        LTER.join.ave$time[idx[n]] = 11
      } else if (LTER.join.ave$year[idx[n]] == 2024) {
        LTER.join.ave$time[idx[n]] = 12
      }}}
  else {}
}


LTER.join.ave <- merge(LTER.join.ave, sites.short.edit, by.x = c("CA_MPA_Name_Short"), 
                      by.y = c("CA_MPA_Name_Short"),
                      all = TRUE)
LTER.join.ave <- LTER.join.ave[complete.cases(LTER.join.ave$year), ]
LTER.join.ave$source <- "LTER"
for(i in 1:l) {
  idx <-which(x[i] == LTER.join.ave$CA_MPA_Name_Short)
  j <-which(x[i] == Sites2$CA_MPA_Name_Short)
  m <- length(idx)
  if (Sites2$MPA_Start[j] == 2012) {
    for (n in 1:m) {
      if (LTER.join.ave$year[idx[n]] <=2012) {
        LTER.join.ave$BA[idx[n]] = "Before"
      } else if (LTER.join.ave$year[idx[n]] == 2013) {
        LTER.join.ave$BA[idx[n]] = "After"
      } else if (LTER.join.ave$year[idx[n]] == 2014) {
        LTER.join.ave$BA[idx[n]] = "After"
      } else if (LTER.join.ave$year[idx[n]] == 2015) {
        LTER.join.ave$BA[idx[n]] = "After"
      } else if (LTER.join.ave$year[idx[n]] == 2016) {
        LTER.join.ave$BA[idx[n]] = "After"
      } else if (LTER.join.ave$year[idx[n]] == 2017) {
        LTER.join.ave$BA[idx[n]] = "After"
      } else if (LTER.join.ave$year[idx[n]] == 2018) {
        LTER.join.ave$BA[idx[n]] = "After"
      } else if (LTER.join.ave$year[idx[n]] == 2019) {
        LTER.join.ave$BA[idx[n]] = "After"
      } else if (LTER.join.ave$year[idx[n]] == 2020) {
        LTER.join.ave$BA[idx[n]] = "After"
      } else if (LTER.join.ave$year[idx[n]] == 2021) {
        LTER.join.ave$BA[idx[n]] = "After"
      } else if (LTER.join.ave$year[idx[n]] == 2022) {
        LTER.join.ave$BA[idx[n]] = "After"
      } else if (LTER.join.ave$year[idx[n]] == 2023) {
        LTER.join.ave$BA[idx[n]] = "After"
      } else if (LTER.join.ave$year[idx[n]] == 2024) {
        LTER.join.ave$BA[idx[n]] = "After"
      }}}
  else {}
}

##############################################################################
########Lobster specific LTER surveys#########################################
#Now we need to bring in size data for lobsters and also LTER fish data
lter.lob = read.csv(paste('data/LTER/Lobster_Abundance_All_Years_20230831.csv'))
lter.lob.site <- merge(lter.lob, Sites2, by.x = c("SITE"), by.y = c("site_id")) #Combine data with site table
lter.lob.site <- subset(lter.lob.site, CA_MPA_Name_Short != "") #Remove sites that aren't associated with MPA

lter.lob.site$biomass <- lter.lob.site$COUNT*(0.001352821 * (lter.lob.site$SIZE ^ 2.913963113)) #Using LTER equation to convert to biomass
lter.lob.site$biomass[is.na(lter.lob.site$biomass)] <- 0

#Now all biomass/count
lter.lob.site.sum <- lter.lob.site %>%
  group_by(YEAR, SITE, TRANSECT, REPLICATE, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>% #Panulirus only 
  summarise_at(c("COUNT", "biomass"), sum)

lter.lob.site.ave <- lter.lob.site.sum %>%
  group_by(YEAR, SITE, TRANSECT, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>% #Panulirus only
  summarise_at(c("COUNT", "biomass"), mean)

lter.lob.site.ave <- lter.lob.site.ave %>%
  group_by(YEAR, SITE, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>% 
  summarise_at(c("COUNT", "biomass"), mean)

#####################################################################################################
lter.lob.site.max <- lter.lob.site.ave %>%
  group_by(status, CA_MPA_Name_Short, MPA_Start, YEAR) %>% 
  summarise_at(c("COUNT", "biomass"), mean)
#####################################################################################################
#Now turn data into proportions######################################################################
lter.lob.site.max$Prop <- 0
lter.lob.site.max$PropCorr <- 0
x = unique(lter.lob.site.max$CA_MPA_Name_Short)
l = length(x)

##########Turn maximum annual biomass into proportion of time series maximum count###################
for(i in 1:l) {
  idx <-which(x[i] == lter.lob.site.max$CA_MPA_Name_Short)
  m <- max(lter.lob.site.max$biomass[idx]) #Find the max of the mpa specific time series to use to create a proportion of the max
  lter.lob.site.max$Prop[idx] <- lter.lob.site.max$biomass[idx]/m
  lter.lob.site.max$PropCorr[idx] <- lter.lob.site.max$Prop[idx] + 0.01 #Add the min value so that there are no zeroes in the dataset
}
lter.lob.site.max$taxon_name <- "Panulirus interruptus"

All.bio.panLTERall.sub <- lter.lob.site.max[, colnames(lter.lob.site.max)[c(2,4,9,1,8)]]
Short.bio.panLTERall <- All.bio.panLTERall.sub %>% 
  spread(status, PropCorr)
Short.bio.panLTERall.diff <- Short.bio.panLTERall[complete.cases(Short.bio.panLTERall), ]
Short.bio.panLTERall.diff$Diff <- Short.bio.panLTERall.diff$mpa/Short.bio.panLTERall.diff$reference #Calculate the response ratio
Short.bio.panLTERall.diff$lnDiff <- log(Short.bio.panLTERall.diff$Diff) #Log response ratio
Short.bio.panLTERall.diff$resp <- "Bio"

#########Now all lter lobster max density##########################################################
for(i in 1:l) {
  idx <-which(x[i] == lter.lob.site.max$CA_MPA_Name_Short)
  m <- max(lter.lob.site.max$COUNT[idx]) #Find the max of the mpa specific time series to use to create a proportion of the max
  lter.lob.site.max$Prop[idx] <- lter.lob.site.max$COUNT[idx]/m
  lter.lob.site.max$PropCorr[idx] <- lter.lob.site.max$Prop[idx] + 0.01 #Add the min value so that there are no zeroes in the dataset
}

All.den.panLTERall.sub <- lter.lob.site.max[, colnames(lter.lob.site.max)[c(2,4,9,1,8)]]
Short.den.panLTERall <- All.den.panLTERall.sub %>% 
  spread(status, PropCorr)
Short.den.panLTERall.diff <- Short.den.panLTERall[complete.cases(Short.den.panLTERall), ]
Short.den.panLTERall.diff$Diff <- Short.den.panLTERall.diff$mpa/Short.den.panLTERall.diff$reference #Calculate the response ratio
Short.den.panLTERall.diff$lnDiff <- log(Short.den.panLTERall.diff$Diff) #Log response ratio
Short.den.panLTERall.diff$resp <- "Den"

LTER.lob <- rbind(Short.bio.panLTERall.diff, Short.den.panLTERall.diff)
colnames(LTER.lob)[2] <- "year"
colnames(LTER.lob)[3] <- "y"

LTER.lob$time <- 0
x <- unique(LTER.lob$CA_MPA_Name_Short) # to redo since sites have been removed
l = length(x)
for(i in 1:l) {
  idx <-which(x[i] == LTER.lob$CA_MPA_Name_Short)
  j <-which(x[i] == Sites2$CA_MPA_Name_Short)
  m <- length(idx)
  if (Sites2$MPA_Start[j] == 2012) {
    for (n in 1:m) {
      if (LTER.lob$year[idx[n]] <=2012) {
        LTER.lob$time[idx[n]] = 0
      } else if (LTER.lob$year[idx[n]] == 2013) {
        LTER.lob$time[idx[n]] = 1
      } else if (LTER.lob$year[idx[n]] == 2014) {
        LTER.lob$time[idx[n]] = 2
      } else if (LTER.lob$year[idx[n]] == 2015) {
        LTER.lob$time[idx[n]] = 3
      } else if (LTER.lob$year[idx[n]] == 2016) {
        LTER.lob$time[idx[n]] = 4
      } else if (LTER.lob$year[idx[n]] == 2017) {
        LTER.lob$time[idx[n]] = 5
      } else if (LTER.lob$year[idx[n]] == 2018) {
        LTER.lob$time[idx[n]] = 6
      } else if (LTER.lob$year[idx[n]] == 2019) {
        LTER.lob$time[idx[n]] = 7
      } else if (LTER.lob$year[idx[n]] == 2020) {
        LTER.lob$time[idx[n]] = 8
      } else if (LTER.lob$year[idx[n]] == 2021) {
        LTER.lob$time[idx[n]] = 9
      } else if (LTER.lob$year[idx[n]] == 2022) {
        LTER.lob$time[idx[n]] = 10
      } else if (LTER.lob$year[idx[n]] == 2023) {
        LTER.lob$time[idx[n]] = 11
      } else if (LTER.lob$year[idx[n]] == 2024) {
        LTER.lob$time[idx[n]] = 12
      }}}
  else {}
}

LTER.lob <- merge(LTER.lob, sites.short.edit, by.x = c("CA_MPA_Name_Short"), 
                       by.y = c("CA_MPA_Name_Short"),
                       all = TRUE)
LTER.lob <- LTER.lob[complete.cases(LTER.lob$year), ]
LTER.lob$source <- "LTER lob surveys"
LTER.lob$BA <- "N/A"

for(i in 1:l) {
  idx <-which(x[i] == LTER.lob$CA_MPA_Name_Short)
  j <-which(x[i] == Sites2$CA_MPA_Name_Short)
  m <- length(idx)
  if (Sites2$MPA_Start[j] == 2012) {
    for (n in 1:m) {
      if (LTER.lob$year[idx[n]] <=2012) {
        LTER.lob$BA[idx[n]] = "Before"
      } else if (LTER.lob$year[idx[n]] == 2013) {
        LTER.lob$BA[idx[n]] = "After"
      } else if (LTER.lob$year[idx[n]] == 2014) {
        LTER.lob$BA[idx[n]] = "After"
      } else if (LTER.lob$year[idx[n]] == 2015) {
        LTER.lob$BA[idx[n]] = "After"
      } else if (LTER.lob$year[idx[n]] == 2016) {
        LTER.lob$BA[idx[n]] = "After"
      } else if (LTER.lob$year[idx[n]] == 2017) {
        LTER.lob$BA[idx[n]] = "After"
      } else if (LTER.lob$year[idx[n]] == 2018) {
        LTER.lob$BA[idx[n]] = "After"
      } else if (LTER.lob$year[idx[n]] == 2019) {
        LTER.lob$BA[idx[n]] = "After"
      } else if (LTER.lob$year[idx[n]] == 2020) {
        LTER.lob$BA[idx[n]] = "After"
      } else if (LTER.lob$year[idx[n]] == 2021) {
        LTER.lob$BA[idx[n]] = "After"
      } else if (LTER.lob$year[idx[n]] == 2022) {
        LTER.lob$BA[idx[n]] = "After"
      } else if (LTER.lob$year[idx[n]] == 2023) {
        LTER.lob$BA[idx[n]] = "After"
      } else if (LTER.lob$year[idx[n]] == 2024) {
        LTER.lob$BA[idx[n]] = "After"
      }}}
  else {}
}

lter.lob.site.max$Den <- lter.lob.site.max$COUNT/300 #plots are 300m2
lter.lob.site.max$Bio <- lter.lob.site.max$biomass/300 

#####################################################################################################
#Create dataframe with all average density/biomass inside/outside mpas LTER kelp fronds##############

lter.lob.edit.den <- lter.lob.site.max[, colnames(lter.lob.site.max)[c(1:4,10,9)]]
LTER.lob.den <- lter.lob.edit.den %>% 
  spread(status, Den)
LTER.lob.den$source <- "LTER"
LTER.lob.den.sub <- LTER.lob.den[complete.cases(LTER.lob.den), ] #Delete any rows where there isn't an mpa and reference
LTER.lob.den.long <- gather(LTER.lob.den.sub, status, value, 5:6)
LTER.lob.den.long$resp <- 'Den'

lter.lob.edit.bio <- lter.lob.site.max[, colnames(lter.lob.site.max)[c(1:4,11,9)]]
LTER.lob.bio <- lter.lob.edit.bio %>% 
  spread(status, Bio)
LTER.lob.bio$source <- "LTER"
LTER.lob.bio.sub <- LTER.lob.bio[complete.cases(LTER.lob.bio), ] #Delete any rows where there isn't an mpa and reference
LTER.lob.bio.long <- gather(LTER.lob.bio.sub, status, value, 5:6)
LTER.lob.bio.long$resp <- 'Bio'

LTER.lob.resp <- rbind(LTER.lob.den.long, LTER.lob.bio.long)
#################################################################################

#################################################################################
#Now we'll bring in LTER kelp stipes#############################################
lter.macro = read.csv(paste('data/LTER/Annual_Kelp_All_Years_20240305.csv'))
lter.macro.site <- merge(lter.macro, Sites2, by.x = c("SITE"), by.y = c("site_id"))

lter.macro.site <- subset(lter.macro.site, CA_MPA_Name_Short != "" & SITE != "AHND") #This is arroyo honda where kelp is essentially absent, Bob Miller said to remove not a good reference
lter.macro.site <- subset(lter.macro.site, FRONDS != -99999)
lter.macro.site <- subset(lter.macro.site, CA_MPA_Name_Short == "Campus Point SMCA" | 
                            CA_MPA_Name_Short == "Naples SMCA")
lter.macro.site$frondDen <- lter.macro.site$FRONDS/lter.macro.site$AREA
lter.macro.site$biomass <- lter.macro.site$frondDen*aveslope*1000 #aveslope is from line 674

#Summarize macro fronds
lter.macro.site.sum <- lter.macro.site %>%
  group_by(YEAR, SITE, TRANSECT, QUAD, SIDE, AREA, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>% #Panulirus only 
  summarise_at(c("FRONDS","biomass"), sum)

lter.macro.site.sum$den <- lter.macro.site.sum$FRONDS/lter.macro.site.sum$AREA

lter.macro.site.ave <- lter.macro.site.sum %>%
  group_by(YEAR, TRANSECT, SITE, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>% 
  summarise_at(c("den","biomass"), mean)

lter.macro.site.ave <- lter.macro.site.ave %>%
  group_by(YEAR, SITE, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>% 
  summarise_at(c("den","biomass"), mean)

lter.macro.site.max <- lter.macro.site.ave %>%
  group_by(status, CA_MPA_Name_Short, MPA_Start, YEAR) %>% 
  summarise_at(c("den","biomass"), mean)

#################################################################################
#Now turn macrocystis frond data into proportions################################
lter.macro.site.max$Prop <- 0
lter.macro.site.max$PropCorr <- 0
x = unique(lter.macro.site.max$CA_MPA_Name_Short)
l = length(x)

#Turn mean annual fronds into proportion of time series maximum count
for(i in 1:l) {
  idx <-which(x[i] == lter.macro.site.max$CA_MPA_Name_Short)
  m <- max(lter.macro.site.max$den[idx]) #Find the max of the mpa specific time series to use to create a proportion of the max
  lter.macro.site.max$Prop[idx] <- lter.macro.site.max$den[idx]/m
  lter.macro.site.max$PropCorr[idx] <- lter.macro.site.max$Prop[idx] + 0.01 #Add the min value so that there are no zeroes in the dataset
}
lter.macro.site.max$taxon_name <- "Macrocystis pyrifera"

All.lter.macro.frond.sub <- lter.macro.site.max[, colnames(lter.macro.site.max)[c(2,4,9,1,8)]]
Short.lter.macro.frond <- All.lter.macro.frond.sub %>% 
  spread(status, PropCorr)
Short.lter.macro.frond.diff <- Short.lter.macro.frond[complete.cases(Short.lter.macro.frond), ]
Short.lter.macro.frond.diff$Diff <- Short.lter.macro.frond.diff$mpa/Short.lter.macro.frond.diff$reference #Calculate the response ratio
Short.lter.macro.frond.diff$lnDiff <- log(Short.lter.macro.frond.diff$Diff) #Log response ratio
Short.lter.macro.frond.diff$resp <- "Den"

colnames(Short.lter.macro.frond.diff)[2] <- "year"
colnames(Short.lter.macro.frond.diff)[3] <- "y"

for(i in 1:l) {
  idx <-which(x[i] == lter.macro.site.max$CA_MPA_Name_Short)
  m <- max(lter.macro.site.max$biomass[idx]) #Find the max of the mpa specific time series to use to create a proportion of the max
  lter.macro.site.max$Prop[idx] <- lter.macro.site.max$biomass[idx]/m
  lter.macro.site.max$PropCorr[idx] <- lter.macro.site.max$Prop[idx] + 0.01 #Add the min value so that there are no zeroes in the dataset
}
lter.macro.site.max$taxon_name <- "Macrocystis pyrifera"

All.lter.macro.bio.sub <- lter.macro.site.max[, colnames(lter.macro.site.max)[c(2,4,9,1,8)]]
Short.lter.macro.bio <- All.lter.macro.bio.sub %>% 
  spread(status, PropCorr)
Short.lter.macro.bio.diff <- Short.lter.macro.bio[complete.cases(Short.lter.macro.bio), ]
Short.lter.macro.bio.diff$Diff <- Short.lter.macro.bio.diff$mpa/Short.lter.macro.bio.diff$reference #Calculate the response ratio
Short.lter.macro.bio.diff$lnDiff <- log(Short.lter.macro.bio.diff$Diff) #Log response ratio
Short.lter.macro.bio.diff$resp <- "Bio"

colnames(Short.lter.macro.bio.diff)[2] <- "year"
colnames(Short.lter.macro.bio.diff)[3] <- "y"

Short.lter.macro <- rbind(Short.lter.macro.frond.diff,Short.lter.macro.bio.diff)

Short.lter.macro$time <- 0
x <- unique(Short.lter.macro$CA_MPA_Name_Short) # to redo since sites have been removed
l = length(x)
for(i in 1:l) {
  idx <-which(x[i] == Short.lter.macro$CA_MPA_Name_Short)
  j <-which(x[i] == Sites2$CA_MPA_Name_Short)
  m <- length(idx)
  if (Sites2$MPA_Start[j] == 2012) {
    for (n in 1:m) {
      if (Short.lter.macro$year[idx[n]] <=2012) {
        Short.lter.macro$time[idx[n]] = 0
      } else if (Short.lter.macro$year[idx[n]] == 2013) {
        Short.lter.macro$time[idx[n]] = 1
      } else if (Short.lter.macro$year[idx[n]] == 2014) {
        Short.lter.macro$time[idx[n]] = 2
      } else if (Short.lter.macro$year[idx[n]] == 2015) {
        Short.lter.macro$time[idx[n]] = 3
      } else if (Short.lter.macro$year[idx[n]] == 2016) {
        Short.lter.macro$time[idx[n]] = 4
      } else if (Short.lter.macro$year[idx[n]] == 2017) {
        Short.lter.macro$time[idx[n]] = 5
      } else if (Short.lter.macro$year[idx[n]] == 2018) {
        Short.lter.macro$time[idx[n]] = 6
      } else if (Short.lter.macro$year[idx[n]] == 2019) {
        Short.lter.macro$time[idx[n]] = 7
      } else if (Short.lter.macro$year[idx[n]] == 2020) {
        Short.lter.macro$time[idx[n]] = 8
      } else if (Short.lter.macro$year[idx[n]] == 2021) {
        Short.lter.macro$time[idx[n]] = 9
      } else if (Short.lter.macro$year[idx[n]] == 2022) {
        Short.lter.macro$time[idx[n]] = 10
      } else if (Short.lter.macro$year[idx[n]] == 2023) {
        Short.lter.macro$time[idx[n]] = 11
      } else if (Short.lter.macro$year[idx[n]] == 2024) {
        Short.lter.macro$time[idx[n]] = 12
      }}}
  else {}
}


Short.lter.macro <- merge(Short.lter.macro, sites.short.edit, by.x = c("CA_MPA_Name_Short"), 
                  by.y = c("CA_MPA_Name_Short"),
                  all = TRUE)
Short.lter.macro <- Short.lter.macro[complete.cases(Short.lter.macro$year), ]
Short.lter.macro$source <- "LTER macro surveys"

for(i in 1:l) {
  idx <-which(x[i] == Short.lter.macro$CA_MPA_Name_Short)
  j <-which(x[i] == Sites2$CA_MPA_Name_Short)
  m <- length(idx)
  if (Sites2$MPA_Start[j] == 2012) {
    for (n in 1:m) {
      if (Short.lter.macro$year[idx[n]] <=2012) {
        Short.lter.macro$BA[idx[n]] = "Before"
      } else if (Short.lter.macro$year[idx[n]] == 2013) {
        Short.lter.macro$BA[idx[n]] = "After"
      } else if (Short.lter.macro$year[idx[n]] == 2014) {
        Short.lter.macro$BA[idx[n]] = "After"
      } else if (Short.lter.macro$year[idx[n]] == 2015) {
        Short.lter.macro$BA[idx[n]] = "After"
      } else if (Short.lter.macro$year[idx[n]] == 2016) {
        Short.lter.macro$BA[idx[n]] = "After"
      } else if (Short.lter.macro$year[idx[n]] == 2017) {
        Short.lter.macro$BA[idx[n]] = "After"
      } else if (Short.lter.macro$year[idx[n]] == 2018) {
        Short.lter.macro$BA[idx[n]] = "After"
      } else if (Short.lter.macro$year[idx[n]] == 2019) {
        Short.lter.macro$BA[idx[n]] = "After"
      } else if (Short.lter.macro$year[idx[n]] == 2020) {
        Short.lter.macro$BA[idx[n]] = "After"
      } else if (Short.lter.macro$year[idx[n]] == 2021) {
        Short.lter.macro$BA[idx[n]] = "After"
      } else if (Short.lter.macro$year[idx[n]] == 2022) {
        Short.lter.macro$BA[idx[n]] = "After"
      } else if (Short.lter.macro$year[idx[n]] == 2023) {
        Short.lter.macro$BA[idx[n]] = "After"
      } else if (Short.lter.macro$year[idx[n]] == 2024) {
        Short.lter.macro$BA[idx[n]] = "After"
      }}}
  else {}
}

#####################################################################################################
#Create dataframe with all average density/biomass inside/outside mpas LTER kelp fronds##############

lter.macro.edit.den <- lter.macro.site.max[, colnames(lter.macro.site.max)[c(1:4,5,9)]]
LTER.macro.den <- lter.macro.edit.den %>% 
  spread(status, den)
LTER.macro.den$source <- "LTER"
LTER.macro.den.sub <- LTER.macro.den[complete.cases(LTER.macro.den), ] #Delete any rows where there isn't an mpa and reference
LTER.macro.den.long <- gather(LTER.macro.den.sub, status, value, 5:6)
LTER.macro.den.long$resp <- 'Den'

lter.macro.edit.bio <- lter.macro.site.max[, colnames(lter.macro.site.max)[c(1:4,6,9)]]
LTER.macro.bio <- lter.macro.edit.bio %>% 
  spread(status, biomass)
LTER.macro.bio$source <- "LTER"
LTER.macro.bio.sub <- LTER.macro.bio[complete.cases(LTER.macro.bio), ] #Delete any rows where there isn't an mpa and reference
LTER.macro.bio.long <- gather(LTER.macro.bio.sub, status, value, 5:6)
LTER.macro.bio.long$resp <- 'Bio'

LTER.macro.resp <- rbind(LTER.macro.den.long, LTER.macro.bio.long)


#################################################################################
#Now LTER fish###################################################################
lter.fish = read.csv(paste('data/LTER/Annual_fish_comb_20240307.csv'))
lter.fish.sub <- subset(lter.fish, SCIENTIFIC_NAME == "Semicossyphus pulcher")

lter.fish.sub.site <- merge(lter.fish.sub, Sites2, by.x = c("SITE"), by.y = c("site_id"))
lter.fish.sub.site <- subset(lter.fish.sub.site, CA_MPA_Name_Short != "" & SITE != "AHND")
lter.fish.sub.site <- subset(lter.fish.sub.site, CA_MPA_Name_Short == "Campus Point SMCA" | 
                            CA_MPA_Name_Short == "Naples SMCA")
lter.fish.sub.site <- subset(lter.fish.sub.site, COUNT != -99999)
lter.fish.sub.site$biomass <- 0.0144*(lter.fish.sub.site$SIZE^3.04)
lter.fish.sub.site$biomass[is.na(lter.fish.sub.site$biomass)] <- 0

#Summarize counts of SPUL
lter.fish.sum <- lter.fish.sub.site %>%
  group_by(YEAR, SITE, TRANSECT, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>% 
  summarise_at(c("COUNT", "biomass"), sum)

lter.fish.ave <- lter.fish.sum %>%
  group_by(YEAR, SITE, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>% 
  summarise_at(c("COUNT", "biomass"), mean)

lter.fish.max <- lter.fish.ave %>%
  group_by(status, CA_MPA_Name_Short, MPA_Start, YEAR) %>% 
  summarise_at(c("COUNT", "biomass"), mean)

#Now turn fish data into proportions#################################################################

lter.fish.max$Prop <- 0
lter.fish.max$PropCorr <- 0
x = unique(lter.fish.max$CA_MPA_Name_Short)
l = length(x)

#Turn maximum annual fronds into proportion of time series maximum count
for(i in 1:l) {
  idx <-which(x[i] == lter.fish.max$CA_MPA_Name_Short)
  m <- max(lter.fish.max$COUNT[idx]) #Find the max of the mpa specific time series to use to create a proportion of the max
  lter.fish.max$Prop[idx] <- lter.fish.max$COUNT[idx]/m
  lter.fish.max$PropCorr[idx] <- lter.fish.max$Prop[idx] + 0.01 #Add the min value so that there are no zeroes in the dataset
}
lter.fish.max$taxon_name <- "SPUL"

All.lter.fish.count.sub <- lter.fish.max[, colnames(lter.fish.max)[c(2,4,9,1,8)]]
Short.lter.fish.count <- All.lter.fish.count.sub %>% 
  spread(status, PropCorr)
Short.lter.fish.count.diff <- Short.lter.fish.count[complete.cases(Short.lter.fish.count), ]
Short.lter.fish.count.diff$Diff <- Short.lter.fish.count.diff$mpa/Short.lter.fish.count.diff$reference #Calculate the response ratio
Short.lter.fish.count.diff$lnDiff <- log(Short.lter.fish.count.diff$Diff) #Log response ratio
Short.lter.fish.count.diff$resp <- "Den"

for(i in 1:l) {
  idx <-which(x[i] == lter.fish.max$CA_MPA_Name_Short)
  m <- max(lter.fish.max$biomass[idx]) #Find the max of the mpa specific time series to use to create a proportion of the max
  lter.fish.max$Prop[idx] <- lter.fish.max$biomass[idx]/m
  lter.fish.max$PropCorr[idx] <- lter.fish.max$Prop[idx] + 0.01 #Add the min value so that there are no zeroes in the dataset
}
lter.fish.max$taxon_name <- "SPUL"

All.lter.spul.bio.sub <- lter.fish.max[, colnames(lter.fish.max)[c(2,4,9,1,8)]]
Short.lter.spul.bio <- All.lter.spul.bio.sub %>% 
  spread(status, PropCorr)
Short.lter.spul.bio.diff <- Short.lter.spul.bio[complete.cases(Short.lter.spul.bio), ]
Short.lter.spul.bio.diff$Diff <- Short.lter.spul.bio.diff$mpa/Short.lter.spul.bio.diff$reference #Calculate the response ratio
Short.lter.spul.bio.diff$lnDiff <- log(Short.lter.spul.bio.diff$Diff) #Log response ratio
Short.lter.spul.bio.diff$resp <- "Bio"

LTER.fish <- rbind(Short.lter.fish.count.diff, Short.lter.spul.bio.diff)

colnames(LTER.fish)[2] <- "year"
colnames(LTER.fish)[3] <- "y"

LTER.fish$time <- 0
x <- unique(LTER.fish$CA_MPA_Name_Short) # to redo since sites have been removed
l = length(x)
for(i in 1:l) {
  idx <-which(x[i] == LTER.fish$CA_MPA_Name_Short)
  j <-which(x[i] == Sites2$CA_MPA_Name_Short)
  m <- length(idx)
  if (Sites2$MPA_Start[j] == 2012) {
    for (n in 1:m) {
      if (LTER.fish$year[idx[n]] <=2012) {
        LTER.fish$time[idx[n]] = 0
      } else if (LTER.fish$year[idx[n]] == 2013) {
        LTER.fish$time[idx[n]] = 1
      } else if (LTER.fish$year[idx[n]] == 2014) {
        LTER.fish$time[idx[n]] = 2
      } else if (LTER.fish$year[idx[n]] == 2015) {
        LTER.fish$time[idx[n]] = 3
      } else if (LTER.fish$year[idx[n]] == 2016) {
        LTER.fish$time[idx[n]] = 4
      } else if (LTER.fish$year[idx[n]] == 2017) {
        LTER.fish$time[idx[n]] = 5
      } else if (LTER.fish$year[idx[n]] == 2018) {
        LTER.fish$time[idx[n]] = 6
      } else if (LTER.fish$year[idx[n]] == 2019) {
        LTER.fish$time[idx[n]] = 7
      } else if (LTER.fish$year[idx[n]] == 2020) {
        LTER.fish$time[idx[n]] = 8
      } else if (LTER.fish$year[idx[n]] == 2021) {
        LTER.fish$time[idx[n]] = 9
      } else if (LTER.fish$year[idx[n]] == 2022) {
        LTER.fish$time[idx[n]] = 10
      } else if (LTER.fish$year[idx[n]] == 2023) {
        LTER.fish$time[idx[n]] = 11
      } else if (LTER.fish$year[idx[n]] == 2024) {
        LTER.fish$time[idx[n]] = 12
      }}}
  else {}
}


LTER.fish <- merge(LTER.fish, sites.short.edit, by.x = c("CA_MPA_Name_Short"), 
                                     by.y = c("CA_MPA_Name_Short"),
                                     all = TRUE)
LTER.fish <- LTER.fish[complete.cases(LTER.fish$year), ]
LTER.fish$source <- "LTER"

for(i in 1:l) {
  idx <-which(x[i] == LTER.fish$CA_MPA_Name_Short)
  j <-which(x[i] == Sites2$CA_MPA_Name_Short)
  m <- length(idx)
  if (Sites2$MPA_Start[j] == 2012) {
    for (n in 1:m) {
      if (LTER.fish$year[idx[n]] <=2012) {
        LTER.fish$BA[idx[n]] = "Before"
      } else if (LTER.fish$year[idx[n]] == 2013) {
        LTER.fish$BA[idx[n]] = "After"
      } else if (LTER.fish$year[idx[n]] == 2014) {
        LTER.fish$BA[idx[n]] = "After"
      } else if (LTER.fish$year[idx[n]] == 2015) {
        LTER.fish$BA[idx[n]] = "After"
      } else if (LTER.fish$year[idx[n]] == 2016) {
        LTER.fish$BA[idx[n]] = "After"
      } else if (LTER.fish$year[idx[n]] == 2017) {
        LTER.fish$BA[idx[n]] = "After"
      } else if (LTER.fish$year[idx[n]] == 2018) {
        LTER.fish$BA[idx[n]] = "After"
      } else if (LTER.fish$year[idx[n]] == 2019) {
        LTER.fish$BA[idx[n]] = "After"
      } else if (LTER.fish$year[idx[n]] == 2020) {
        LTER.fish$BA[idx[n]] = "After"
      } else if (LTER.fish$year[idx[n]] == 2021) {
        LTER.fish$BA[idx[n]] = "After"
      } else if (LTER.fish$year[idx[n]] == 2022) {
        LTER.fish$BA[idx[n]] = "After"
      } else if (LTER.fish$year[idx[n]] == 2023) {
        LTER.fish$BA[idx[n]] = "After"
      } else if (LTER.fish$year[idx[n]] == 2024) {
        LTER.fish$BA[idx[n]] = "After"
      }}}
  else {}
}

lter.fish.max$Bio <- lter.fish.max$biomass/80
lter.fish.max$Den <- lter.fish.max$COUNT/80
#####################################################################################################
#Create dataframe with all average density/biomass inside/outside mpas LTER FISH#####################
lter.fish.edit.den <- lter.fish.max[, colnames(lter.fish.max)[c(1:4,11,9)]]
LTER.fish.den <- lter.fish.edit.den %>% 
  spread(status, Den)
LTER.fish.den$source <- "LTER"
LTER.fish.den.sub <- LTER.fish.den[complete.cases(LTER.fish.den), ] #Delete any rows where there isn't an mpa and reference
LTER.fish.den.long <- gather(LTER.fish.den.sub, status, value, 5:6)
LTER.fish.den.long$resp <- 'Den'

lter.fish.edit.bio <- lter.fish.max[, colnames(lter.fish.max)[c(1:4,10,9)]]
LTER.fish.bio <- lter.fish.edit.bio %>% 
  spread(status, Bio)
LTER.fish.bio$source <- "LTER"
LTER.fish.bio.sub <- LTER.fish.bio[complete.cases(LTER.fish.den), ] #Delete any rows where there isn't an mpa and reference
LTER.fish.bio.long <- gather(LTER.fish.bio.sub, status, value, 5:6)
LTER.fish.bio.long$resp <- 'Bio'

LTER.fish.resp <- rbind(LTER.fish.den.long, LTER.fish.bio.long)

#####################################################################################################
#####KFM fish processing#####################################################################
mbon.fish = read.csv(paste('data/MBON/SBCMBON_kelp_forest_integrated_fish_20231022.csv'))
mbon.fish = subset(mbon.fish, data_source == "kfm")
unique(mbon.fish$sample_method)
mbon.fish.spul = subset(mbon.fish, proj_taxon_id == "t-k-127" | proj_taxon_id == "t-k-233")
tmpDateFormat<-"%Y-%m-%d"
tmp1date<-as.Date(mbon.fish.spul$date,format=tmpDateFormat)
if(length(tmp1date) == length(tmp1date[!is.na(tmp1date)])){mbon.fish.spul$date <- tmp1date } else {print("Date conversion failed for MBON$date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1date) 
mbon.fish.spul$year <- year(mbon.fish.spul$date)
mbon.fish.spul.sub = subset(mbon.fish.spul, sample_method == "rdfc" & year >= 2003 | sample_method == "visualfish")
mbon.spul.site <- merge(mbon.fish.spul.sub, Sites2, by = c("site_id"))
mbon.spul.site$count <- as.numeric(mbon.spul.site$count)
mbon.spul.site <- subset(mbon.spul.site, CA_MPA_Name_Short != "")
mbon.spul.site$den <- mbon.spul.site$count/mbon.spul.site$area #transform to density

#Calc sheephead per MPA/Reference
mbon.spul.sum <- mbon.spul.site %>%
  group_by(site_id, sample_method, transect_id, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, year) %>% 
  summarise_at(c("den"), mean)

mbon.spul.ave <- mbon.spul.sum %>%
  group_by(site_id, sample_method, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, year) %>% 
  summarise_at(c("den"), mean)

mbon.spul.max <- mbon.spul.ave %>%
  group_by(status, sample_method, CA_MPA_Name_Short, ChannelIsland, MPA_Start, year) %>% 
  summarise_at(c("den"), mean)

mbon.spul.max <- subset(mbon.spul.max, CA_MPA_Name_Short != "Anacapa Island SMR 1978")
colnames(mbon.spul.max)[colnames(mbon.spul.max) == 'den'] <- 'count'

mbon.spul.max$Prop <- 0
mbon.spul.max$PropCorr <- 0
x = unique(mbon.spul.max$CA_MPA_Name_Short)
l = length(x)
mbon.spul.max$taxon_name <- "SPUL"

mbon.spul.max.rd <- subset(mbon.spul.max, sample_method == "rdfc")
mbon.spul.max.vf <- subset(mbon.spul.max, sample_method == "visualfish")

#####################################################################################################
#Create dataframe with all average density/biomass inside/outside mpas KFM FISH######################
mbon.fish.edit <- mbon.spul.max.vf[, colnames(mbon.spul.max.vf)[c(1,3,5,6,7,10)]]
KFM.fish.den <- mbon.fish.edit %>% 
  spread(status, count)
KFM.fish.den$source <- "KFM"
KFM.fish.den.sub <- KFM.fish.den[complete.cases(KFM.fish.den), ] #Delete any rows where there isn't an mpa and reference
KFM.fish.den.long <- gather(KFM.fish.den.sub, status, value, 5:6)
KFM.fish.den.long$resp <- 'Den'

##########################################################################################
#Turn maximum annual fish data from KFM into proportion of time series maximum density####
for(i in 1:l) {
  idx <-which(x[i] == mbon.spul.max.rd$CA_MPA_Name_Short)
  m <- max(mbon.spul.max.rd$count[idx]) #Find the max of the mpa specific time series to use to create a proportion of the max
  mbon.spul.max.rd$Prop[idx] <- mbon.spul.max.rd$count[idx]/m
  mbon.spul.max.rd$PropCorr[idx] <- mbon.spul.max.rd$Prop[idx] + 0.01 #Add the min value so that there are no zeroes in the dataset
}

All.den.spul.kfm.sub.rd <- mbon.spul.max.rd[, colnames(mbon.spul.max.rd)[c(3,2,6,10,1,9)]]
Short.den.spul.kfm.sub.rd <- All.den.spul.kfm.sub.rd %>% 
  spread(status, PropCorr)
Short.den.spul.kfm.sub.rd.diff <- Short.den.spul.kfm.sub.rd[complete.cases(Short.den.spul.kfm.sub.rd), ]
Short.den.spul.kfm.sub.rd.diff$Diff <- Short.den.spul.kfm.sub.rd.diff$mpa/Short.den.spul.kfm.sub.rd.diff$reference #Calculate the response ratio
Short.den.spul.kfm.sub.rd.diff$lnDiff <- log(Short.den.spul.kfm.sub.rd.diff$Diff) #Log response ratio
colnames(Short.den.spul.kfm.sub.rd.diff)[4] <- "y"
Short.den.spul.kfm.sub.rd.diff$resp <- "RD"

#Turn maximum annual vf into proportion of time series maximum biomass
for(i in 1:l) {
  idx <-which(x[i] == mbon.spul.max.vf$CA_MPA_Name_Short)
  m <- max(mbon.spul.max.vf$count[idx]) #Find the max of the mpa specific time series to use to create a proportion of the max
  mbon.spul.max.vf$Prop[idx] <- mbon.spul.max.vf$count[idx]/m
  mbon.spul.max.vf$PropCorr[idx] <- mbon.spul.max.vf$Prop[idx] + 0.01 #Add the min value so that there are no zeroes in the dataset
}

All.den.spul.kfm.sub.vf <- mbon.spul.max.vf[, colnames(mbon.spul.max.vf)[c(3,2,6,10,1,9)]]
Short.den.spul.kfm.sub.vf <- All.den.spul.kfm.sub.vf %>% 
  spread(status, PropCorr)
Short.den.spul.kfm.sub.vf.diff <- Short.den.spul.kfm.sub.vf[complete.cases(Short.den.spul.kfm.sub.vf), ]
Short.den.spul.kfm.sub.vf.diff$Diff <- Short.den.spul.kfm.sub.vf.diff$mpa/Short.den.spul.kfm.sub.vf.diff$reference #Calculate the response ratio
Short.den.spul.kfm.sub.vf.diff$lnDiff <- log(Short.den.spul.kfm.sub.vf.diff$Diff) #Log response ratio
colnames(Short.den.spul.kfm.sub.vf.diff)[4] <- "y"
Short.den.spul.kfm.sub.vf.diff$resp <- "Den"

Short.den.spul.kfm.sub.diff <- rbind(Short.den.spul.kfm.sub.vf.diff, Short.den.spul.kfm.sub.rd.diff)



Short.den.spul.kfm.sub.diff$time <- 0
x <- unique(Short.den.spul.kfm.sub.diff$CA_MPA_Name_Short) # to redo since sites have been removed
l = length(x)
for(i in 1:l) {
  idx <-which(x[i] == Short.den.spul.kfm.sub.diff$CA_MPA_Name_Short)
  j <-which(x[i] == Sites2$CA_MPA_Name_Short)
  m <- length(idx)
  if (Sites2$MPA_Start[j] == 2003) {
    for (n in 1:m) {
      if (Short.den.spul.kfm.sub.diff$year[idx[n]] <=2003) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 0
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2004) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 1
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2005) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 2
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2006) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 3
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2007) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 4
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2008) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 5
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2009) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 6
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2010) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 7
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2011) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 8
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2012) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 9
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2013) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 10
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2014) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 11
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2015) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 12
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2016) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 13
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2017) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 14
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2018) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 15
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2019) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 16
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2020) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 17
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2021) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 18
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2022) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 19
      } else if (Short.den.spul.kfm.sub.diff$year[idx[n]] == 2023) {
        Short.den.spul.kfm.sub.diff$time[idx[n]] = 20
      }}}
  else {}
}



kfm.fish <- merge(Short.den.spul.kfm.sub.diff, sites.short.edit, by.x = c("CA_MPA_Name_Short"), 
                   by.y = c("CA_MPA_Name_Short"),
                   all = TRUE)
kfm.fish <- kfm.fish[complete.cases(kfm.fish$year), ]
kfm.fish$source <- "KFM"
kfm.fish$BA <- NA
x <- unique(kfm.fish$CA_MPA_Name_Short) # to redo since sites have been removed
l = length(x)
for(i in 1:l) {
  idx <-which(x[i] == kfm.fish$CA_MPA_Name_Short)
  j <-which(x[i] == Sites2$CA_MPA_Name_Short)
  m <- length(idx)
  if (Sites2$MPA_Start[j] == 2003) {
    for (n in 1:m) {
      if (kfm.fish$year[idx[n]] <=2003) {
        kfm.fish$BA[idx[n]] = "Before"
      } else if (kfm.fish$year[idx[n]] == 2004) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2005) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2006) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2007) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2008) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2009) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2010) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2011) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2012) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2013) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2014) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2015) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2016) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2017) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2018) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2019) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2020) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2021) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2022) {
        kfm.fish$BA[idx[n]] = "After"
      } else if (kfm.fish$year[idx[n]] == 2023) {
        kfm.fish$BA[idx[n]] = "After"
      }}}
  else {}
}

###############################################################################################

###############################################################################################
##########Final wrangling before join of all data streams######################################
Swath.join.sub <- Swath.join.sub[, colnames(Swath.join.sub)[c(1:3,5:15)]]
kfm.fish <- kfm.fish[, colnames(kfm.fish)[c(1,3:15)]]
KFM.join.ave <- KFM.join.ave[, colnames(KFM.join.ave)[c(1:3,5:15)]]
#Create dataframe with all Response Ratio data
All.RR <- rbind(LTER.join.ave, Swath.join.sub, KFM.join.ave, LTER.lob, Short.lter.macro, LTER.fish, kfm.fish)

All.RR.sub <- subset(All.RR, CA_MPA_Name_Short != 'Blue Cavern Onshore SMCA' & CA_MPA_Name_Short != 'Painted Cave SMCA' &
                       CA_MPA_Name_Short != 'Dana Point SMCA' & CA_MPA_Name_Short != 'Farnsworth Onshore SMCA' & 
                       CA_MPA_Name_Short != 'Point Dume SMR' & CA_MPA_Name_Short != 'Cat Harbor SMCA' &
                       CA_MPA_Name_Short != 'Swamis SMCA' & CA_MPA_Name_Short != "Anacapa Island SMCA" &
                       CA_MPA_Name_Short != 'Long Point SMR' & CA_MPA_Name_Short != 'Point Dume SMCA' &
                       CA_MPA_Name_Short != 'Santa Barbara Island SMR')
All.spul <- subset(All.RR, CA_MPA_Name_Short == 'Blue Cavern Onshore SMCA' & y == 'SPUL' | 
                     CA_MPA_Name_Short == 'Dana Point SMCA' & y == 'SPUL' |
                     CA_MPA_Name_Short == 'Farnsworth Onshore SMCA' & y == 'SPUL' |
                     CA_MPA_Name_Short == 'Swamis SMCA' & y == 'SPUL' |
                     CA_MPA_Name_Short == 'Cat Harbor SMCA' & y == 'SPUL' |
                     CA_MPA_Name_Short == 'Long Point SMR' & y == 'SPUL' |
                     CA_MPA_Name_Short == 'Santa Barbara Island SMR' & y == 'SPUL' & source == 'PISCO' |
                     CA_MPA_Name_Short == 'Santa Barbara Island SMR' & source == 'KFM')               
  
All.RR.sub <- rbind(All.RR.sub, All.spul)      
All.RR.sub.trans <- All.RR.sub
All.RR.sub.trans$y <- as.character(All.RR.sub.trans$y)
All.RR.sub.trans$y[All.RR.sub.trans$y == "PANINT"] <- "Panulirus interruptus"

All.RR.sub.trans$y[All.RR.sub.trans$y == "SPUL"] <- "Semicossyphus pulcher"

All.RR.sub.trans$y[All.RR.sub.trans$y == "MACPYRAD"] <- "Macrocystis pyrifera"
All.RR.sub.trans$y[All.RR.sub.trans$y == "MESFRAAD"] <- "Mesocentrotus franciscanus"
All.RR.sub.trans$y[All.RR.sub.trans$y == "STRPURAD"] <- "Strongylocentrotus purpuratus"

All.RR.sub.trans$source[All.RR.sub.trans$source == "LTER macro surveys"] <- "LTER"
All.RR.sub.trans$source[All.RR.sub.trans$source == "LTER lob surveys"] <- "LTER"

#########################################################################################################
#Create dataframe with all processed mean (not RR) density/biomass data#################
PISCO.resp$y[PISCO.resp$y == "PANINT"] <- "Panulirus interruptus"
PISCO.resp$y[PISCO.resp$y == "SPUL"] <- "Semicossyphus pulcher"
PISCO.resp$y[PISCO.resp$y == "MACPYRAD"] <- "Macrocystis pyrifera"
PISCO.resp$y[PISCO.resp$y == "MESFRAAD"] <- "Mesocentrotus franciscanus"
PISCO.resp$y[PISCO.resp$y == "STRPURAD"] <- "Strongylocentrotus purpuratus"
colnames(PISCO.resp)[3] <- 'taxon_name'

KFM.fish.den.long <- KFM.fish.den.long[, colnames(KFM.fish.den.long)[c(1,3:8)]]
KFM.fish.den.long$taxon_name[KFM.fish.den.long$taxon_name == "SPUL"] <- "Semicossyphus pulcher"

KFM.resp <- KFM.resp[, colnames(KFM.resp)[c(1,3:8)]]
colnames(KFM.resp)[3] <- 'taxon_name'
KFM.resp$taxon_name[KFM.resp$taxon_name == "PANINT"] <- "Panulirus interruptus"
KFM.resp$taxon_name[KFM.resp$taxon_name == "MACPYRAD"] <- "Macrocystis pyrifera"
KFM.resp$taxon_name[KFM.resp$taxon_name == "MESFRAAD"] <- "Mesocentrotus franciscanus"
KFM.resp$taxon_name[KFM.resp$taxon_name == "STRPURAD"] <- "Strongylocentrotus purpuratus"

LTER.fish.resp <- LTER.fish.resp[, colnames(LTER.fish.resp)[c(1,3:8)]]
LTER.fish.resp$taxon_name[LTER.fish.resp$taxon_name == "SPUL"] <- "Semicossyphus pulcher"

colnames(LTER.fish.resp)[2] <- 'year'
LTER.macro.resp <- LTER.macro.resp[, colnames(LTER.macro.resp)[c(1,3:8)]]
colnames(LTER.macro.resp)[2] <- 'year'
LTER.lob.resp <- LTER.lob.resp[, colnames(LTER.lob.resp)[c(1,3:8)]]
colnames(LTER.lob.resp)[2] <- 'year'
LTER.resp <- LTER.resp[, colnames(LTER.resp)[c(1,3:8)]]

All.Resp <- rbind(KFM.fish.den.long, LTER.fish.resp, LTER.macro.resp, LTER.lob.resp, LTER.resp,
                  KFM.resp, PISCO.resp)

All.Resp.sub <- subset(All.Resp, CA_MPA_Name_Short != 'Blue Cavern Onshore SMCA' & CA_MPA_Name_Short != 'Painted Cave SMCA' &
                       CA_MPA_Name_Short != 'Dana Point SMCA' & CA_MPA_Name_Short != 'Farnsworth Onshore SMCA' & 
                       CA_MPA_Name_Short != 'Point Dume SMR' & CA_MPA_Name_Short != 'Cat Harbor SMCA' &
                       CA_MPA_Name_Short != 'Swamis SMCA' & CA_MPA_Name_Short != "Anacapa Island SMCA" &
                       CA_MPA_Name_Short != 'Long Point SMR' & CA_MPA_Name_Short != 'Point Dume SMCA' &
                       CA_MPA_Name_Short != 'Santa Barbara Island SMR')
All.Resp.spul <- subset(All.Resp, CA_MPA_Name_Short == 'Blue Cavern Onshore SMCA' & taxon_name == 'Semicossyphus pulcher' | 
                     CA_MPA_Name_Short == 'Dana Point SMCA' & taxon_name == 'Semicossyphus pulcher' |
                     CA_MPA_Name_Short == 'Farnsworth Onshore SMCA' & taxon_name == 'Semicossyphus pulcher' |
                     CA_MPA_Name_Short == 'Swamis SMCA' & taxon_name == 'Semicossyphus pulcher' |
                     CA_MPA_Name_Short == 'Cat Harbor SMCA' & taxon_name == 'Semicossyphus pulcher' |
                     CA_MPA_Name_Short == 'Long Point SMR' & taxon_name == 'Semicossyphus pulcher' |
                     CA_MPA_Name_Short == 'Santa Barbara Island SMR' & taxon_name == 'SPUL' & source == 'PISCO' |
                     CA_MPA_Name_Short == 'Santa Barbara Island SMR' & source == 'KFM')               
All.Resp.sub <- rbind(All.Resp.sub, All.Resp.spul)

x <- unique(All.Resp.sub$CA_MPA_Name_Short) # to redo since sites have been removed
l = length(x)
for(i in 1:l) {
  idx <-which(x[i] == All.Resp.sub$CA_MPA_Name_Short)
  j <-which(x[i] == Site$CA_MPA_Name_Short)
  m <- length(idx)
  if (Site$MPA_Start[j] == 2003) {
    for (n in 1:m) {
      if (All.Resp.sub$year[idx[n]] <=2003) {
        All.Resp.sub$BA[idx[n]] = "Before"
      } else if (All.Resp.sub$year[idx[n]] == 2004) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2005) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2006) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2007) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2008) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2009) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2010) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2011) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2012) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2013) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2014) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2015) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2016) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2017) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2018) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2019) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2020) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2021) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2022) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2023) {
        All.Resp.sub$BA[idx[n]] = "After"
      }}
    
  } else if (Site$MPA_Start[j] == 2012) {
    for (n in 1:m) {
      if (All.Resp.sub$year[idx[n]] <= 2012) {
        All.Resp.sub$BA[idx[n]] = "Before"
      } else if (All.Resp.sub$year[idx[n]] == 2013) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2014) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2015) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2016) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2017) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2018) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2019) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2020) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2021) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2022) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2023) {
        All.Resp.sub$BA[idx[n]] = "After"
      } else if (All.Resp.sub$year[idx[n]] == 2024) {
        All.Resp.sub$BA[idx[n]] = "After"
      }}}}

x <- unique(All.Resp.sub$CA_MPA_Name_Short) # to redo since sites have been removed
l = length(x)
for(i in 1:l) {
  idx <-which(x[i] == All.Resp.sub$CA_MPA_Name_Short)
  j <-which(x[i] == Site$CA_MPA_Name_Short)
  m <- length(idx)
  if (Site$MPA_Start[j] == 2003) {
    for (n in 1:m) {
      if (All.Resp.sub$year[idx[n]] <=2003) {
        All.Resp.sub$time[idx[n]] = 0
      } else if (All.Resp.sub$year[idx[n]] == 2004) {
        All.Resp.sub$time[idx[n]] = 1
      } else if (All.Resp.sub$year[idx[n]] == 2005) {
        All.Resp.sub$time[idx[n]] = 2
      } else if (All.Resp.sub$year[idx[n]] == 2006) {
        All.Resp.sub$time[idx[n]] = 3
      } else if (All.Resp.sub$year[idx[n]] == 2007) {
        All.Resp.sub$time[idx[n]] = 4
      } else if (All.Resp.sub$year[idx[n]] == 2008) {
        All.Resp.sub$time[idx[n]] = 5
      } else if (All.Resp.sub$year[idx[n]] == 2009) {
        All.Resp.sub$time[idx[n]] = 6
      } else if (All.Resp.sub$year[idx[n]] == 2010) {
        All.Resp.sub$time[idx[n]] = 7
      } else if (All.Resp.sub$year[idx[n]] == 2011) {
        All.Resp.sub$time[idx[n]] = 8
      } else if (All.Resp.sub$year[idx[n]] == 2012) {
        All.Resp.sub$time[idx[n]] = 9
      } else if (All.Resp.sub$year[idx[n]] == 2013) {
        All.Resp.sub$time[idx[n]] = 10
      } else if (All.Resp.sub$year[idx[n]] == 2014) {
        All.Resp.sub$time[idx[n]] = 11
      } else if (All.Resp.sub$year[idx[n]] == 2015) {
        All.Resp.sub$time[idx[n]] = 12
      } else if (All.Resp.sub$year[idx[n]] == 2016) {
        All.Resp.sub$time[idx[n]] = 13
      } else if (All.Resp.sub$year[idx[n]] == 2017) {
        All.Resp.sub$time[idx[n]] = 14
      } else if (All.Resp.sub$year[idx[n]] == 2018) {
        All.Resp.sub$time[idx[n]] = 15
      } else if (All.Resp.sub$year[idx[n]] == 2019) {
        All.Resp.sub$time[idx[n]] = 16
      } else if (All.Resp.sub$year[idx[n]] == 2020) {
        All.Resp.sub$time[idx[n]] = 17
      } else if (All.Resp.sub$year[idx[n]] == 2021) {
        All.Resp.sub$time[idx[n]] = 18
      } else if (All.Resp.sub$year[idx[n]] == 2022) {
        All.Resp.sub$time[idx[n]] = 19
      } else if (All.Resp.sub$year[idx[n]] == 2023) {
        All.Resp.sub$time[idx[n]] = 20
      }}
    
  } else if (Site$MPA_Start[j] == 2012) {
    for (n in 1:m) {
      if (All.Resp.sub$year[idx[n]] <= 2012) {
        All.Resp.sub$time[idx[n]] = 0
      } else if (All.Resp.sub$year[idx[n]] == 2013) {
        All.Resp.sub$time[idx[n]] = 1
      } else if (All.Resp.sub$year[idx[n]] == 2014) {
        All.Resp.sub$time[idx[n]] = 2
      } else if (All.Resp.sub$year[idx[n]] == 2015) {
        All.Resp.sub$time[idx[n]] = 3
      } else if (All.Resp.sub$year[idx[n]] == 2016) {
        All.Resp.sub$time[idx[n]] = 4
      } else if (All.Resp.sub$year[idx[n]] == 2017) {
        All.Resp.sub$time[idx[n]] = 5
      } else if (All.Resp.sub$year[idx[n]] == 2018) {
        All.Resp.sub$time[idx[n]] = 6
      } else if (All.Resp.sub$year[idx[n]] == 2019) {
        All.Resp.sub$time[idx[n]] = 7
      } else if (All.Resp.sub$year[idx[n]] == 2020) {
        All.Resp.sub$time[idx[n]] = 8
      } else if (All.Resp.sub$year[idx[n]] == 2021) {
        All.Resp.sub$time[idx[n]] = 9
      } else if (All.Resp.sub$year[idx[n]] == 2022) {
        All.Resp.sub$time[idx[n]] = 10
      } else if (All.Resp.sub$year[idx[n]] == 2023) {
        All.Resp.sub$time[idx[n]] = 11
      } else if (All.Resp.sub$year[idx[n]] == 2024) {
        All.Resp.sub$time[idx[n]] = 12
      }}}}

All.Resp.sub <- All.Resp.sub[, colnames(All.Resp.sub)[c(1:7,10:11)]]



AveResponse <- summarySE(data = All.Resp.sub, measurevar = "value", groupvars = c("taxon_name", "source", "resp"))
write_csv(AveResponse, "AveResponses.csv")


##########################################################################################################

##########################################################################################################
#Plot combined time series of RR across sources###############################################
All.purps <- subset(All.RR.sub.trans, y == "Strongylocentrotus purpuratus" & resp == "Den")
x <- unique(All.purps$CA_MPA_Name_Short) # to redo since sites have been removed
MPA_implement <- filter(Site, Site %in% x)
MPA_implement.sub <- filter(MPA_implement, Site %in% x)
PurpAll <- ggplot(data=All.purps, aes(x=year, y=lnDiff, color = source, group_by = source, linetype = source)) +
  #geom_line(size = 1) + 
  geom_point(size = 3) + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  facet_wrap(~CA_MPA_Name_Short, scales = "free") +
  ylab(expression("lnRR("*italic("S. purpuratus")~" Density)")) + 
  xlab("Year") + 
  theme_classic() +
  scale_linetype_manual(values=c(1,4,3)) +  
  #scale_color_manual(name = "Status", values = c("#d8b365", "#5ab4ac"),labels = c("MPA", "Reference")) +
  theme(plot.title = element_text(hjust = 0.5), strip.text = element_text(size = 14)) +
  theme(legend.position = "top", legend.key.width = unit(1,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank(), 
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold")) 
PurpAll

ggsave(plot = PurpAll, file = "plots/All_TimeSeries_PurpDen.png", 
       device = "png",  bg = "white",
       width = 35, height = 30, units = "cm", dpi = 300)

All.purps <- subset(All.RR.sub.trans, y == "Strongylocentrotus purpuratus" & resp == "Bio")

PurpAll <- ggplot(data=All.purps, aes(x=year, y=lnDiff, color = source, group_by = source, linetype = source)) +
  #geom_line(size = 1) + 
  geom_point(size = 3) + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  facet_wrap(~CA_MPA_Name_Short, scales = "free") +
  ylab(expression("lnRR("*italic("S. purpuratus")~" Biomass)")) + 
  xlab("Year") + 
  theme_classic() +
  scale_linetype_manual(values=c(1,4,3)) +  
  #scale_color_manual(name = "Status", values = c("#d8b365", "#5ab4ac"),labels = c("MPA", "Reference")) +
  theme(plot.title = element_text(hjust = 0.5), strip.text = element_text(size = 14)) +
  theme(legend.position = "top", legend.key.width = unit(1,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank(), 
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold")) 
PurpAll

ggsave(plot = PurpAll, file = "plots/All_TimeSeries_PurpBio.png", 
       device = "png",  bg = "white",
       width = 35, height = 30, units = "cm", dpi = 300)

All.red <- subset(All.RR.sub.trans, y == "Mesocentrotus franciscanus" & resp == "Den")
x <- unique(All.red$CA_MPA_Name_Short) # to redo since sites have been removed
MPA_implement <- filter(Site, Site %in% x)
MPA_implement.sub <- filter(MPA_implement, Site %in% x)
RedAll <- ggplot(data=All.red, aes(x=year, y=lnDiff, color = source, group_by = source, linetype = source)) +
  #geom_line(size = 1) + 
  geom_point(size = 3) + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  facet_wrap(~CA_MPA_Name_Short, scales = "free") +
  ylab(expression("lnRR("*italic("M. franciscanus")~" Density)")) + 
  xlab("Year") + 
  theme_classic() +
  scale_linetype_manual(values=c(1,4,3)) +  
  theme(plot.title = element_text(hjust = 0.5), strip.text = element_text(size = 14)) +
  theme(legend.position = "top", legend.key.width = unit(1,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank(), 
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold")) 
RedAll

ggsave(plot = RedAll, file = "plots/All_TimeSeries_RedDen.png", 
       device = "png",  bg = "white",
       width = 35, height = 30, units = "cm", dpi = 300)

All.red <- subset(All.RR.sub.trans, y == "Mesocentrotus franciscanus" & resp == "Bio")

RedAll <- ggplot(data=All.red, aes(x=year, y=lnDiff, color = source, group_by = source, linetype = source)) +
  #geom_line(size = 1) + 
  geom_point(size = 3) + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  facet_wrap(~CA_MPA_Name_Short, scales = "free") +
  ylab(expression("lnRR("*italic("M. franciscanus")~" Biomass)")) + 
  xlab("Year") + 
  theme_classic() +
  scale_linetype_manual(values=c(1,4,3)) +  
  theme(plot.title = element_text(hjust = 0.5), strip.text = element_text(size = 14)) +
  theme(legend.position = "top", legend.key.width = unit(1,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank(), 
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold")) 
RedAll

ggsave(plot = RedAll, file = "plots/All_TimeSeries_RedBio.png", 
       device = "png",  bg = "white",
       width = 35, height = 30, units = "cm", dpi = 300)

All.macro <- subset(All.RR.sub.trans, y == "Macrocystis pyrifera" & resp == "Bio")
x <- unique(All.macro$CA_MPA_Name_Short) # to redo since sites have been removed
MPA_implement <- filter(Site, Site %in% x)
MPA_implement.sub <- filter(MPA_implement, Site %in% x)
MacroAll <- ggplot(data=All.macro, aes(x=year, y=lnDiff, color = source, group_by = source, linetype = source)) +
  geom_point(size = 3) + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  facet_wrap(~CA_MPA_Name_Short, scales = "free") +
  ylab(expression("lnRR("*italic("M. pyrifera")~" Density)")) + 
  xlab("Year") + 
  theme_classic() +
  scale_linetype_manual(values=c(1,4,3)) +  
  theme(plot.title = element_text(hjust = 0.5), strip.text = element_text(size = 14)) +
  theme(legend.position = "top", legend.key.width = unit(1,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank(), 
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold")) 
MacroAll

ggsave(plot = MacroAll, file = "plots/All_TimeSeries_MacroStipe.png", 
       device = "png",  bg = "white",
       width = 35, height = 30, units = "cm", dpi = 300)


All.spul <- subset(All.RR.sub.trans, y == "Semicossyphus pulcher" & resp == "Den")
x <- unique(All.spul$CA_MPA_Name_Short) # to redo since sites have been removed
MPA_implement <- filter(Site, Site %in% x)
SpulAll <- ggplot(data=All.spul, aes(x=year, y=lnDiff, color = source, group_by = source, linetype = source)) +
  geom_point(size = 3) + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  facet_wrap(~CA_MPA_Name_Short, scales = "free") +
  ylab(expression("lnRR("*italic("S. pulcher")~" Density)")) + 
  xlab("Year") + 
  theme_classic() +
  scale_linetype_manual(values=c(1,4,3)) +  
  theme(plot.title = element_text(hjust = 0.5), strip.text = element_text(size = 14)) +
  theme(legend.position = "top", legend.key.width = unit(1,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank(), 
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold")) 
SpulAll

ggsave(plot = SpulAll, file = "plots/All_TimeSeries_SheepheadDen.png", 
       device = "png",  bg = "white",
       width = 35, height = 30, units = "cm", dpi = 300)

All.spul <- subset(All.RR.sub.trans, y == "Semicossyphus pulcher" & resp == "Bio")

SpulAll <- ggplot(data=All.spul, aes(x=year, y=lnDiff, color = source, group_by = source, linetype = source)) +
  geom_point(size = 3) + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  facet_wrap(~CA_MPA_Name_Short, scales = "free") +
  ylab(expression("lnRR("*italic("S. pulcher")~" Biomass)")) + 
  xlab("Year") + 
  theme_classic() +
  scale_linetype_manual(values=c(1,4,3)) +  
  theme(plot.title = element_text(hjust = 0.5), strip.text = element_text(size = 14)) +
  theme(legend.position = "top", legend.key.width = unit(1,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank(), 
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold")) 
SpulAll

ggsave(plot = SpulAll, file = "plots/All_TimeSeries_SheepheadBio.png", 
       device = "png",  bg = "white",
       width = 35, height = 30, units = "cm", dpi = 300)

All.lob <- subset(All.RR.sub.trans, y == "Panulirus interruptus" & resp == "Den")
x <- unique(All.lob$CA_MPA_Name_Short) # to redo since sites have been removed
MPA_implement <- filter(Site, Site %in% x)
LobAll <- ggplot(data=All.lob, aes(x=year, y=lnDiff, color = source, group_by = source, linetype = source)) +
  geom_point(size = 3) + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  facet_wrap(~CA_MPA_Name_Short, scales = "free") +
  ylab(expression("lnRR("*italic("P. interruptus")~" Density)")) + 
  xlab("Year") + 
  theme_classic() +
  scale_linetype_manual(values=c(1,4,3)) +  
  theme(plot.title = element_text(hjust = 0.5), strip.text = element_text(size = 14)) +
  theme(legend.position = "top", legend.key.width = unit(1,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank(), 
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold")) 
LobAll

ggsave(plot = LobAll, file = "plots/All_TimeSeries_LobsterDen.png", 
       device = "png",  bg = "white",
       width = 35, height = 30, units = "cm", dpi = 300)

All.lob <- subset(All.RR.sub.trans, y == "Panulirus interruptus" & resp == "Bio")

LobAll <- ggplot(data=All.lob, aes(x=year, y=lnDiff, color = source, group_by = source, linetype = source)) +
  geom_point(size = 3) + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  facet_wrap(~CA_MPA_Name_Short, scales = "free") +
  ylab(expression("lnRR("*italic("P. interruptus")~" Biomass)")) + 
  xlab("Year") + 
  theme_classic() +
  scale_linetype_manual(values=c(1,4,3)) +  
  theme(plot.title = element_text(hjust = 0.5), strip.text = element_text(size = 14)) +
  theme(legend.position = "top", legend.key.width = unit(1,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank(), 
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold")) 
LobAll

ggsave(plot = LobAll, file = "plots/All_TimeSeries_LobsterBio.png", 
       device = "png",  bg = "white",
       width = 35, height = 30, units = "cm", dpi = 300)

##########################################################################################################

##########################################################################################################
#Now start analyses!######################################################################################
##########################################################################################################
#LTER first###############################################################################################
cols <- c("Sheephead" = "#377eb8", "Kelp" = "#4daf4a", "Lobs" = "#ff7f00", "Purps" = "#984ea3", "Reds" = "#e41a1c")
MPA_implement <- Site
names(MPA_implement)[names(MPA_implement) == "Site"] <- "Site_ID"
names(MPA_implement)[names(MPA_implement) == "CA_MPA_Name_Short"] <- "Site"
##Create a dataframe for all effect sizes#################################################################
SumStats = data.frame(Taxa = NA, MPA = NA, Mean = NA, SE=NA, SD = NA, CI = NA, Model = NA, Source = NA, Resp = NA, BA = NA, Primary = NA, Type = NA, LinearBefore = NA)  #Initialize empty dataframe

###Functions##############################################################################################
source("code/MEE_Supplement_example.r") #Import pBACIPS code
## Create an asymptotic function
myASYfun<-function(delta, time.model)
{
  funAsy<-function(parS, time.model)	(parS$M * time.model) / (parS$L + time.model) + parS$B
  residFun<-function(p, observed, time.model) observed + funAsy(p,time.model)
  parStart <- list(M=mean(delta[time.model.of.impact:length(time.true)]), B=mean(delta[1:time.model.of.impact]), L=1)
  nls_ASY_out <- nls.lm(par=parStart, fn= residFun, observed=delta, time.model=time.model, control = nls.lm.control(maxfev = integer(), maxiter = 1000))
  foAsy<-delta~(M * time.model) / (L + time.model) + B
  startPar<-c(-coef(nls_ASY_out)[1], coef(nls_ASY_out)[2], coef(nls_ASY_out)[3])
  asyFit<-nls2(foAsy, start=startPar, algorithm="brute-force") # nls2 enables to calculate AICc
  asyFit
}
## Create a sigmoid function
mySIGfun<-function(delta, time.model)
{
  funSIG<-function(parS, time.model)	(parS$M * (time.model/parS$L)^parS$K) / (1 + (time.model/parS$L) ^ parS$K) + parS$B
  residFun<-function(p, observed, time.model) observed + funSIG(p,time.model)
  parStart <- list(M=mean(delta[time.model.of.impact:length(time.true)]), B=mean(delta[1:time.model.of.impact]), L=mean(time.model), K=5)
  nls_SIG_out <- nls.lm(par=parStart, fn= residFun, observed=delta, time.model=time.model, control = nls.lm.control(maxfev = integer(), maxiter = 1000))
  foSIG<-delta~(M * (time.model/L) ^ K) / (1 + (time.model/L) ^ K) + B
  startPar<-c(-coef(nls_SIG_out)[1],-coef(nls_SIG_out)[2],coef(nls_SIG_out)[3],coef(nls_SIG_out)[4])
  sigFit<-nls2(foSIG, start=startPar, algorithm="brute-force") # nls2 enables to calculate AICc
  sigFit
}
#Start with LTER ###############################################################################
##Purps Density#################################################################################
LTER.purps <- subset(All.RR.sub.trans, source == "LTER" & y == "Strongylocentrotus purpuratus" & resp == "Den")
LTER.purps$time.model <- NA #time.model sets all before to 0 and then all after sequential starting with 1
LTER.purps$time.true <- NA #time.true all time points are sequential
x <- unique(LTER.purps$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)
LTER.purps <- LTER.purps[order(as.numeric(LTER.purps$year),decreasing = FALSE), ]

#Create time.model and time.true columns
for(i in 1:l) {
  idx <- which(LTER.purps$CA_MPA_Name_Short == x[i])
  min.id <- min(LTER.purps$year[idx])
  max.id <- max(LTER.purps$year[idx])
  LTER.purps$time.model[idx] <-c(rep(0,length(which(LTER.purps$BA[idx]=="Before"))),seq(0:length(which(LTER.purps$BA[idx]=="After"))))
  LTER.purps$time.true[idx] <-seq(1:length(idx))
}

LinearBefore = data.frame(pInt = NA, pSlope = NA) #initialize dataframe to hold pvalues for linear regressions

for (j in x) {
  temp <- subset(LTER.purps, LTER.purps$CA_MPA_Name_Short == j, select = 
                   c(CA_MPA_Name_Short, mpa, reference, lnDiff, time.true, time.model, BA))
  temp = subset(temp, temp$BA == "Before")
  lmBefore <- lm(lnDiff ~ time.true, data = temp)
  pl <- summary(lmBefore)$coefficients[,4]
  LinearBefore <- rbind(LinearBefore, pl)
}

LinearBefore <- na.omit(LinearBefore) #Remove first line with no data
LinearBefore$site <- x #Join mpa info so we know which one we're looking at

###Linear in the before, so exclude both density and biomass 
###this is a conservative approach, but since there isn't enough data in before for biomass to evaluate
###it seems highly likely that if den was linear, so would biomass

#Now Reds#########################################################################################
LTER.reds <- subset(All.RR.sub.trans, source == "LTER" & y == "Mesocentrotus franciscanus" & resp == 'Den')
x <- unique(LTER.reds$CA_MPA_Name_Short) #Need to redo since sites have been removed


#Create time.model and time.true columns
for(i in 1:l) {
  idx <- which(LTER.reds$CA_MPA_Name_Short == x[i])
  min.id <- min(LTER.reds$year[idx])
  max.id <- max(LTER.reds$year[idx])
  LTER.reds$time.model[idx] <-c(rep(0,length(which(LTER.reds$BA[idx]=="Before"))),seq(0:length(which(LTER.reds$BA[idx]=="After"))))
  LTER.reds$time.true[idx] <-seq(1:length(idx))
}

LinearBefore = data.frame(pInt = NA, pSlope = NA) #initialize dataframe to hold pvalues for linear regressions

for (j in x) {
  temp <- subset(LTER.reds, LTER.reds$CA_MPA_Name_Short == j, select = 
                   c(CA_MPA_Name_Short, mpa, reference, lnDiff, time.true, time.model, BA))
  temp = subset(temp, temp$BA == "Before")
  lmBefore <- lm(lnDiff ~ time.true, data = temp)
  pl <- summary(lmBefore)$coefficients[,4]
  LinearBefore <- rbind(LinearBefore, pl)
}

LinearBefore <- na.omit(LinearBefore) #Remove first line with no data
LinearBefore$site <- x #Join mpa info so we know which one we're looking at

likelihood.poly <- data.frame(matrix(ncol = 5, nrow = 0)) #Initialize dataframe to put likelihood values from loop
pvals.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues
Norm.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues

LTER.reds$lnmpa <- log(LTER.reds$mpa)
LTER.reds$lnreference <- log(LTER.reds$reference)

for (j in x) {
  temp <- subset(LTER.reds, LTER.reds$CA_MPA_Name_Short == j, select = c(CA_MPA_Name_Short, lnmpa, 
                                                                           lnreference, lnDiff, time.true, time.model))
  ProgressiveChangeBACIPS(control=temp$lnreference,
                          impact=temp$lnmpa,
                          time.true=temp$time.true,
                          time.model=temp$time.model)
  ps <- summary(step.Model)[[1]][["Pr(>F)"]][1]
  pl <- summary(linear.Model)$coefficients[,4]
  psi <- summary(sigmoid.Model)$coefficients[,4]
  pa <- summary(asymptotic.Model)$coefficients[,4]
  pval <- c(ps, pl, psi, pa)
  likelihood.poly <- rbind(likelihood.poly, w)
  pvals.poly <- rbind(pvals.poly,pval)
  t <- shapiro.test(temp$lnDiff)
  Norm.poly <- rbind(Norm.poly, t$p.value)
}

#Rename column headers
colnames(likelihood.poly) <- c("Step","Linear","Asymptotic","Sigmoid")
colnames(pvals.poly) <- c("Step","Linear","Asymptotic M","Asymptotic B","Asymptotic L","Sigmoid M","Sigmoid B","Sigmoid L","Sigmoid K")
likelihood.poly <- cbind(likelihood.poly, x)

#Only doing Naples, since Campus Point is linear in before; step was best fit model

LTER.redsNap <- subset(LTER.reds, CA_MPA_Name_Short == "Naples SMCA")

mod1 <-lm(data = LTER.redsNap, lnDiff ~ BA)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
ols_plot_resid_qq(mod1)
acf(mod1$residuals, type = "correlation")
leveneTest(lnDiff ~ as.factor(BA), data = LTER.redsNap)

##Calculate effect size
mfran.lter.nap <- tidy(emmeans(mod1, ~ BA))
mfran.lter.nap$stdev <- mfran.lter.nap$std.error*sqrt(mfran.lter.nap$df+2)
pSD <- sqrt((((mfran.lter.nap$df[1]+1)*(mfran.lter.nap$stdev[1]^2))+
               ((mfran.lter.nap$df[2]+1)*(mfran.lter.nap$stdev[2]^2)))
            /(mfran.lter.nap$df[1]+2))
pSE <- pSD*sqrt((1/(mfran.lter.nap$df[1]+2)) + (1/(mfran.lter.nap$df[2]+2)))
pCI <- pSE*1.96
mean <- mfran.lter.nap$estimate[1]-mfran.lter.nap$estimate[2] #We need to calculate the effect size at each year

SumStats[1,] <- rbind('M. franciscanus', 'Naples SMCA', mean, pSE, pSD, pCI, 'Step', 'LTER', 'Den', 'Y', 'Y','BACI', 'N')


#Now Reds biomass#########################################################################################
LTER.reds <- subset(All.RR.sub.trans, source == "LTER" & y == "Mesocentrotus franciscanus" & resp == 'Bio')
LTER.reds$time.model <- NA #time.model sets all before to 0 and then all after sequential starting with 1
LTER.reds$time.true <- NA #time.true all time points are sequential
x <- unique(LTER.reds$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)
LTER.reds <- LTER.reds[order(as.numeric(LTER.reds$year),decreasing = FALSE), ]

#Create time.model and time.true columns
for(i in 1:l) {
  idx <- which(LTER.reds$CA_MPA_Name_Short == x[i])
  min.id <- min(LTER.reds$year[idx])
  max.id <- max(LTER.reds$year[idx])
  LTER.reds$time.model[idx] <-c(rep(0,length(which(LTER.reds$BA[idx]=="Before"))),seq(0:length(which(LTER.reds$BA[idx]=="After"))))
  LTER.reds$time.true[idx] <-seq(1:length(idx))
}

#Naples only since Campus Point was linear in density in the before
LTER.nap.red <- filter(LTER.reds, CA_MPA_Name_Short == 'Naples SMCA') 

## Fit a linear model
linear.Model<-lm(data = LTER.nap.red, lnDiff ~ time.model)
anova(linear.Model)
summary(linear.Model)

LTER.nap.red <- filter(LTER.reds, CA_MPA_Name_Short == 'Naples SMCA' & time > 0) 

CP.mean <- summarySE(LTER.nap.red, measurevar = "lnDiff" )
SumStats[2,] <- rbind('M. franciscanus', 'Naples SMCA', CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'LTER', 'Bio', 'N', 'Y','CI', 'NA')

#Now biomass of macro stipes#############################################################
##Note, we only do biomass since density is the exact same proportional data since bio is a mulitplier of density
LTER.macro <- subset(All.RR.sub.trans, source == "LTER" & y == "Macrocystis pyrifera" & resp == "Bio")

LTER.macro$time.model <- NA #time.model sets all before to 0 and then all after sequential starting with 1
LTER.macro$time.true <- NA #time.true all time points are sequential
x <- unique(LTER.macro$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)
LTER.macro <- LTER.macro[order(as.numeric(LTER.macro$year),decreasing = FALSE), ]

#Create time.model and time.true columns
for(i in 1:l) {
  idx <- which(LTER.macro$CA_MPA_Name_Short == x[i])
  min.id <- min(LTER.macro$year[idx])
  max.id <- max(LTER.macro$year[idx])
  LTER.macro$time.model[idx] <-c(rep(0,length(which(LTER.macro$BA[idx]=="Before"))),seq(0:length(which(LTER.macro$BA[idx]=="After"))))
  LTER.macro$time.true[idx] <-seq(1:length(idx))
}

LinearBefore = data.frame(pInt = NA, pSlope = NA) #initialize dataframe to hold pvalues for linear regressions

for (j in x) {
  temp <- subset(LTER.macro, LTER.macro$CA_MPA_Name_Short == j, select = 
                   c(CA_MPA_Name_Short, mpa, reference, lnDiff, time.true, time.model, BA))
  temp = subset(temp, temp$BA == "Before")
  lmBefore <- lm(lnDiff ~ time.true, data = temp)
  pl <- summary(lmBefore)$coefficients[,4]
  LinearBefore <- rbind(LinearBefore, pl)
}

LinearBefore <- na.omit(LinearBefore) #Remove first line with no data
LinearBefore$site <- x #Join mpa info so we know which one we're looking at

likelihood.poly <- data.frame(matrix(ncol = 5, nrow = 0)) #Initialize dataframe to put likelihood values from loop
pvals.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues
Norm.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues

LTER.macro$lnmpa <- log(LTER.macro$mpa)
LTER.macro$lnreference <- log(LTER.macro$reference)

for (j in x) {
  temp <- subset(LTER.macro, LTER.macro$CA_MPA_Name_Short == j, select = c(CA_MPA_Name_Short, lnmpa, 
                                                                         lnreference, lnDiff, time.true, time.model))
  ProgressiveChangeBACIPS(control=temp$lnreference,
                          impact=temp$lnmpa,
                          time.true=temp$time.true,
                          time.model=temp$time.model)
  ps <- summary(step.Model)[[1]][["Pr(>F)"]][1]
  pl <- summary(linear.Model)$coefficients[,4]
  psi <- summary(sigmoid.Model)$coefficients[,4]
  pa <- summary(asymptotic.Model)$coefficients[,4]
  pval <- c(ps, pl, psi, pa)
  likelihood.poly <- rbind(likelihood.poly, w)
  pvals.poly <- rbind(pvals.poly,pval)
  t <- shapiro.test(temp$lnDiff)
  Norm.poly <- rbind(Norm.poly, t$p.value)
}

#Rename column headers
colnames(likelihood.poly) <- c("Step","Linear","Asymptotic","Sigmoid")
colnames(pvals.poly) <- c("Step","Linear","Asymptotic M","Asymptotic B","Asymptotic L","Sigmoid M","Sigmoid B","Sigmoid L","Sigmoid K")
likelihood.poly <- cbind(likelihood.poly, x)

#Naples 
LTER.nap.macro <- filter(LTER.macro, CA_MPA_Name_Short == 'Naples SMCA') 

mod1 <-lm(data = LTER.nap.macro, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
ols_plot_resid_qq(mod1)
acf(mod1$residuals, type = "correlation")
qqPlot(mod1) #Plot residuals
leveneTest(lnDiff ~ as.factor(BA), data = LTER.nap.macro)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year

SumStats[3,] <- rbind('M. pyrifera', 'Naples SMCA', mean, pSE, pSD, pCI, 'Linear', 'LTER', 'Bio', 'Y', 'Y','pBACIPS', 'N')

#Campus Point
LTER.cp.macro <- filter(LTER.macro, CA_MPA_Name_Short == 'Campus Point SMCA') 

mod1 <-lm(data = LTER.cp.macro, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
ols_plot_resid_qq(mod1)
acf(mod1$residuals, type = "correlation")
qqPlot(mod1) #Plot residuals
leveneTest(lnDiff ~ as.factor(BA), data = LTER.cp.macro)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year

SumStats[4,] <- rbind('M. pyrifera', 'Campus Point SMCA', mean, pSE, pSD, pCI, 'Linear', 'LTER', 'Bio', 'Y', 'Y','pBACIPS', 'N')



#Now density of SPUL#############################################################
LTER.SPUL.den <- subset(All.RR.sub.trans, source == "LTER" & y == "Semicossyphus pulcher" & resp == "Den")
LTER.SPUL.den$time.model <- NA #time.model sets all before to 0 and then all after sequential starting with 1
LTER.SPUL.den$time.true <- NA #time.true all time points are sequential
x <- unique(LTER.SPUL.den$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)
LTER.SPUL.den <- LTER.SPUL.den[order(as.numeric(LTER.SPUL.den$year),decreasing = FALSE), ]

#Create time.model and time.true columns
for(i in 1:l) {
  idx <- which(LTER.SPUL.den$CA_MPA_Name_Short == x[i])
  min.id <- min(LTER.SPUL.den$year[idx])
  max.id <- max(LTER.SPUL.den$year[idx])
  LTER.SPUL.den$time.model[idx] <-c(rep(0,length(which(LTER.SPUL.den$BA[idx]=="Before"))),seq(0:length(which(LTER.SPUL.den$BA[idx]=="After"))))
  LTER.SPUL.den$time.true[idx] <-seq(1:length(idx))
}

LinearBefore = data.frame(pInt = NA, pSlope = NA) #initialize dataframe to hold pvalues for linear regressions

for (j in x) {
  temp <- subset(LTER.SPUL.den, LTER.SPUL.den$CA_MPA_Name_Short == j, select = 
                   c(CA_MPA_Name_Short, mpa, reference, lnDiff, time.true, time.model, BA))
  temp = subset(temp, temp$BA == "Before")
  lmBefore <- lm(lnDiff ~ time.true, data = temp)
  pl <- summary(lmBefore)$coefficients[,4]
  LinearBefore <- rbind(LinearBefore, pl)
}

LinearBefore <- na.omit(LinearBefore) #Remove first line with no data
LinearBefore$site <- x #Join mpa info so we know which one we're looking at

likelihood.poly <- data.frame(matrix(ncol = 5, nrow = 0)) #Initialize dataframe to put likelihood values from loop
pvals.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues
Norm.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues

LTER.SPUL.den$lnmpa <- log(LTER.SPUL.den$mpa)
LTER.SPUL.den$lnreference <- log(LTER.SPUL.den$reference)

for (j in x) {
  temp <- subset(LTER.SPUL.den, LTER.SPUL.den$CA_MPA_Name_Short == j, select = c(CA_MPA_Name_Short, lnmpa, 
                                                                                         lnreference, lnDiff, time.true, time.model))
  ProgressiveChangeBACIPS(control=temp$lnreference,
                          impact=temp$lnmpa,
                          time.true=temp$time.true,
                          time.model=temp$time.model)
  ps <- summary(step.Model)[[1]][["Pr(>F)"]][1]
  pl <- summary(linear.Model)$coefficients[,4]
  psi <- summary(sigmoid.Model)$coefficients[,4]
  pa <- summary(asymptotic.Model)$coefficients[,4]
  pval <- c(ps, pl, psi, pa)
  likelihood.poly <- rbind(likelihood.poly, w)
  pvals.poly <- rbind(pvals.poly,pval)
  t <- shapiro.test(temp$lnDiff)
  Norm.poly <- rbind(Norm.poly, t$p.value)
}

#Rename column headers
colnames(likelihood.poly) <- c("Step","Linear","Asymptotic","Sigmoid")
colnames(pvals.poly) <- c("Step","Linear","Asymptotic M","Asymptotic B","Asymptotic L","Sigmoid M","Sigmoid B","Sigmoid L","Sigmoid K")
likelihood.poly <- cbind(likelihood.poly, x)

#Naples biomass was linear in the before

#Campus Point
LTER.cp.SPUL.den <- filter(LTER.SPUL.den, CA_MPA_Name_Short == 'Campus Point SMCA') 

mod1 <-lm(data = LTER.cp.SPUL.den, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
ols_plot_resid_qq(mod1)
acf(mod1$residuals, type = "correlation")
qqPlot(mod1) #Plot residuals
leveneTest(lnDiff ~ as.factor(BA), data = LTER.cp.SPUL.den)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year

SumStats[5,] <- rbind('S. pulcher', 'Campus Point SMCA', mean, pSE, pSD, pCI, 'Linear', 'LTER', 'Den', 'Y', 'Y','pBACIPS', 'N')

#Now biomass of  SPUL#####################################################################
LTER.SPUL.bio <- subset(All.RR.sub.trans, source == "LTER" & y == "Semicossyphus pulcher" & resp == "Bio")
LTER.SPUL.bio$time.model <- NA #time.model sets all before to 0 and then all after sequential starting with 1
LTER.SPUL.bio$time.true <- NA #time.true all time points are sequential
x <- unique(LTER.SPUL.bio$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)
LTER.SPUL.bio <- LTER.SPUL.bio[order(as.numeric(LTER.SPUL.bio$year),decreasing = FALSE), ]

#Create time.model and time.true columns
for(i in 1:l) {
  idx <- which(LTER.SPUL.bio$CA_MPA_Name_Short == x[i])
  min.id <- min(LTER.SPUL.bio$year[idx])
  max.id <- max(LTER.SPUL.bio$year[idx])
  LTER.SPUL.bio$time.model[idx] <-c(rep(0,length(which(LTER.SPUL.bio$BA[idx]=="Before"))),seq(0:length(which(LTER.SPUL.bio$BA[idx]=="After"))))
  LTER.SPUL.bio$time.true[idx] <-seq(1:length(idx))
}

LinearBefore = data.frame(pInt = NA, pSlope = NA) #initialize dataframe to hold pvalues for linear regressions

for (j in x) {
  temp <- subset(LTER.SPUL.bio, LTER.SPUL.bio$CA_MPA_Name_Short == j, select = 
                   c(CA_MPA_Name_Short, mpa, reference, lnDiff, time.true, time.model, BA))
  temp = subset(temp, temp$BA == "Before")
  lmBefore <- lm(lnDiff ~ time.true, data = temp)
  pl <- summary(lmBefore)$coefficients[,4]
  LinearBefore <- rbind(LinearBefore, pl)
}

LinearBefore <- na.omit(LinearBefore) #Remove first line with no data
LinearBefore$site <- x #Join mpa info so we know which one we're looking at

likelihood.poly <- data.frame(matrix(ncol = 5, nrow = 0)) #Initialize dataframe to put likelihood values from loop
pvals.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues
Norm.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues

x <- unique(LTER.SPUL.bio$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)

LTER.SPUL.bio$lnmpa <- log(LTER.SPUL.bio$mpa)
LTER.SPUL.bio$lnreference <- log(LTER.SPUL.bio$reference)

for (j in x) {
  temp <- subset(LTER.SPUL.bio, LTER.SPUL.bio$CA_MPA_Name_Short == j, select = c(CA_MPA_Name_Short, lnmpa, 
                                                                                             lnreference, lnDiff, time.true, time.model))
  ProgressiveChangeBACIPS(control=temp$lnreference,
                          impact=temp$lnmpa,
                          time.true=temp$time.true,
                          time.model=temp$time.model)
  ps <- summary(step.Model)[[1]][["Pr(>F)"]][1]
  pl <- summary(linear.Model)$coefficients[,4]
  psi <- summary(sigmoid.Model)$coefficients[,4]
  pa <- summary(asymptotic.Model)$coefficients[,4]
  pval <- c(ps, pl, psi, pa)
  likelihood.poly <- rbind(likelihood.poly, w)
  pvals.poly <- rbind(pvals.poly,pval)
  t <- shapiro.test(temp$lnDiff)
  Norm.poly <- rbind(Norm.poly, t$p.value)
}

#Rename column headers
colnames(likelihood.poly) <- c("Step","Linear","Asymptotic","Sigmoid")
colnames(pvals.poly) <- c("Step","Linear","Asymptotic M","Asymptotic B","Asymptotic L","Sigmoid M","Sigmoid B","Sigmoid L","Sigmoid K")
likelihood.poly <- cbind(likelihood.poly, x)

#Campus Point
LTER.cp.SPUL.bio <- filter(LTER.SPUL.bio, CA_MPA_Name_Short == 'Campus Point SMCA') 

mod1 <-lm(data = LTER.cp.SPUL.bio, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
ols_plot_resid_qq(mod1)
acf(mod1$residuals, type = "correlation")
qqPlot(mod1) #Plot residuals
leveneTest(lnDiff ~ as.factor(BA), data = LTER.cp.SPUL.bio)

##Calculate effect size
spul.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
spul.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
spul.lter.before$stdev <- spul.lter.before$std.error*sqrt(spul.lter.before$df+2)
spul.lter.after$stdev <- spul.lter.after$std.error*sqrt(spul.lter.after$df+2)

pSD <- sqrt((((spul.lter.before$df[1]+1)*(spul.lter.before$stdev[1]^2))+
               ((spul.lter.after$df[1]+1)*(spul.lter.after$stdev[1]^2)))
            /(spul.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(spul.lter.before$df[1]+2)) + (1/(spul.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- spul.lter.after$estimate[1]-spul.lter.before$estimate[1] #We need to calculate the effect size at each year

SumStats[6,] <- rbind('S. pulcher', 'Campus Point SMCA', mean, pSE, pSD, pCI, 'Linear', 'LTER', 'Bio', 'Y', 'Y','pBACIPS', 'N')


##########Naples biomass of S. pulcher was linear in the before


#Now density of PANINT only CI data available unfortunately######################################################
LTER.PANINT.den <- subset(All.RR.sub.trans, source == "LTER" & y == "Panulirus interruptus" & resp == "Den")
x <- unique(LTER.PANINT.den$CA_MPA_Name_Short) #Need to redo since sites have been removed

#Campus Point
LTER.cp.LOB <- filter(LTER.PANINT.den, CA_MPA_Name_Short == 'Campus Point SMCA') 

mod1 <-lm(data = LTER.cp.LOB, lnDiff ~ time)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
ols_plot_resid_qq(mod1)
acf(mod1$residuals, type = "correlation")
qqPlot(mod1) #Plot residuals
leveneTest(lnDiff ~ as.factor(BA), data = LTER.cp.SPUL.bio)

##Calculate effect size
CP.mean <- summarySE(LTER.cp.LOB, measurevar = "lnDiff" )
SumStats[7,] <- rbind('P. interruptus', 'Campus Point SMCA', CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'LTER', 'Den', 'N', 'Y','CI', 'NA')

#Naples
LTER.nap.LOB <- filter(LTER.PANINT.den, CA_MPA_Name_Short == 'Naples SMCA') 

mod1 <-lm(data = LTER.nap.LOB, lnDiff ~ time)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
ols_plot_resid_qq(mod1)
acf(mod1$residuals, type = "correlation")
qqPlot(mod1) #Plot residuals
leveneTest(lnDiff ~ as.factor(BA), data = LTER.cp.SPUL.bio)

##Calculate effect size
CP.mean <- summarySE(LTER.nap.LOB, measurevar = "lnDiff" )
SumStats[8,] <- rbind('P. interruptus', 'Naples SMCA', CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'LTER', 'Den', 'N', 'Y','CI', 'NA')



#Now biomass of PANINT#####################################################################
LTER.PANINT.bio <- subset(All.RR.sub.trans, source == "LTER" & y == "Panulirus interruptus" & resp == "Bio")

#Campus Point
LTER.cp.LOB <- filter(LTER.PANINT.bio, CA_MPA_Name_Short == 'Campus Point SMCA') 

mod1 <-lm(data = LTER.cp.LOB, lnDiff ~ time)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
ols_plot_resid_qq(mod1)
acf(mod1$residuals, type = "correlation")
qqPlot(mod1) #Plot residuals

##Calculate effect size
lter.before <- tidy(emmeans(mod1, ~ time, at = list(time = 0) )) 
lter.after <- tidy(emmeans(mod1, ~ time, at = list(time = 11) ))
lter.before$stdev <- lter.before$std.error*sqrt(lter.before$df+2)
lter.after$stdev <- lter.after$std.error*sqrt(lter.after$df+2)

pSD <- sqrt((((lter.before$df[1]+1)*(lter.before$stdev[1]^2))+
               ((lter.after$df[1]+1)*(lter.after$stdev[1]^2)))
            /(lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(lter.before$df[1]+2)) + (1/(lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- lter.after$estimate[1]-lter.before$estimate[1] #We need to calculate the effect size at each year

SumStats[9,] <- rbind('P. interruptus', 'Campus Point SMCA', mean, pSE, pSD, pCI, 'Linear', 'LTER', 'Bio', 'N', 'Y','CI', 'NA')


#Naples
LTER.nap.LOB <- filter(LTER.PANINT.bio, CA_MPA_Name_Short == 'Naples SMCA') 

mod1 <-lm(data = LTER.nap.LOB, lnDiff ~ time)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
ols_plot_resid_qq(mod1)
acf(mod1$residuals, type = "correlation")
qqPlot(mod1) #Plot residuals

#Calculating just the mean effec
CP.mean <- summarySE(LTER.nap.LOB, measurevar = "lnDiff" )
SumStats[10,] <- rbind('P. interruptus', 'Naples SMCA', CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'LTER', 'Bio', 'N', 'Y','CI', 'NA')



###Now on to PISCO########################################################################################
####Reds##################################################################################################
######Start with Density##################################################################################
Mes <- subset(All.RR.sub.trans, y == "Mesocentrotus franciscanus" & source == "PISCO" & resp == "Den")

x <- unique(Mes$CA_MPA_Name_Short)
Coef.table.Mes.PISCO <- data.frame(Site = character(), Mes.Slope = double(), Mes.p = double(), Mes.r = double(), 
                         em0 = double(), l.CL0 = double(), u.CL0 = double(), em11 = double(), l.CL11 = double(), u.CL11 = double())
l <- length(x)
for(i in 1:l) {
  idx <- which(Mes$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = Mes[idx,])
  p <- summary(Lm.Ab)$coefficients[2,4]
  slope <- summary(Lm.Ab)$coefficients[2,1]
  r <- summary(Lm.Ab)$adj.r.squared
  slope = ifelse(p <= 0.05, slope, NA)
  Coef.table.Mes.PISCO[i,1] <-  x[i]
  Coef.table.Mes.PISCO[i,2] <- slope
  Coef.table.Mes.PISCO[i,3] <- p
  Coef.table.Mes.PISCO[i,4] <- r
  mod.emm <- emmeans(Lm.Ab, ~ time, at = list(time = c(0, 11)))
  Coef.table.Mes.PISCO[i,5] <- summary(mod.emm)$emmean[1]
  Coef.table.Mes.PISCO[i,6] <- summary(mod.emm)$lower.CL[1]
  Coef.table.Mes.PISCO[i,7] <- summary(mod.emm)$upper.CL[1]
  Coef.table.Mes.PISCO[i,8] <- summary(mod.emm)$emmean[2]
  Coef.table.Mes.PISCO[i,9] <- summary(mod.emm)$lower.CL[2]
  Coef.table.Mes.PISCO[i,10] <- summary(mod.emm)$upper.CL[2]
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(Mes[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
  SumStats[t+1,] <- rbind('M. franciscanus', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Den', 'N', 'Y', 'CI', 'NA')
  } else {
  SumStats[t+1,] <- rbind('M. franciscanus', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Den', 'N', 'N', 'CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('M. franciscanus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Den', 'N', 'N', 'CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('M. franciscanus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Den', 'N', 'Y', 'CI', 'NA')
  }

}

#####Now Biomass#######################################################################################

Mes <- subset(All.RR.sub.trans, y == "Mesocentrotus franciscanus" & source == "PISCO" & resp == "Bio")
Mes <- subset(Mes, CA_MPA_Name_Short != 'Matlahuayl SMR')
x <- unique(Mes$CA_MPA_Name_Short)
Coef.table.Mes.PISCO <- data.frame(Site = character(), Mes.Slope = double(), Mes.p = double(), Mes.r = double(), 
                                   em0 = double(), l.CL0 = double(), u.CL0 = double(), em11 = double(), l.CL11 = double(), u.CL11 = double())
l <- length(x)
for(i in 1:l) {
  idx <- which(Mes$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = Mes[idx,])
  p <- summary(Lm.Ab)$coefficients[2,4]
  slope <- summary(Lm.Ab)$coefficients[2,1]
  r <- summary(Lm.Ab)$adj.r.squared
  slope = ifelse(p <= 0.05, slope, NA)
  Coef.table.Mes.PISCO[i,1] <-  x[i]
  Coef.table.Mes.PISCO[i,2] <- slope
  Coef.table.Mes.PISCO[i,3] <- p
  Coef.table.Mes.PISCO[i,4] <- r
  mod.emm <- emmeans(Lm.Ab, ~ time, at = list(time = c(0, 11)))
  Coef.table.Mes.PISCO[i,5] <- summary(mod.emm)$emmean[1]
  Coef.table.Mes.PISCO[i,6] <- summary(mod.emm)$lower.CL[1]
  Coef.table.Mes.PISCO[i,7] <- summary(mod.emm)$upper.CL[1]
  Coef.table.Mes.PISCO[i,8] <- summary(mod.emm)$emmean[2]
  Coef.table.Mes.PISCO[i,9] <- summary(mod.emm)$lower.CL[2]
  Coef.table.Mes.PISCO[i,10] <- summary(mod.emm)$upper.CL[2]
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(Mes[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('M. franciscanus', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Bio', 'N', 'Y','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('M. franciscanus', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Bio', 'N', 'N','CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('M. franciscanus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Bio', 'N', 'N','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('M. franciscanus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Bio', 'N', 'Y','CI', 'NA')
  }  
}


########################################################################################################
#########Purps##########################################################################################
#Density##
Str <- subset(All.RR.sub.trans, y == "Strongylocentrotus purpuratus" & source == "PISCO" & resp == "Den")
x <- unique(Str$CA_MPA_Name_Short)
Coef.table.Str.PISCO <- data.frame(Site = character(), Mes.Slope = double(), Mes.p = double(), Mes.r = double(), 
                         em0 = double(), l.CL0 = double(), u.CL0 = double(), em11 = double(), l.CL11 = double(), u.CL11 = double())
l <- length(x)
for(i in 1:l) {
  idx <- which(Str$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = Str[idx,])
  p <- summary(Lm.Ab)$coefficients[2,4]
  slope <- summary(Lm.Ab)$coefficients[2,1]
  r <- summary(Lm.Ab)$adj.r.squared
  slope = ifelse(p <= 0.05, slope, NA)
  Coef.table.Str.PISCO[i,1] <-  x[i]
  Coef.table.Str.PISCO[i,2] <- slope
  Coef.table.Str.PISCO[i,3] <- p
  Coef.table.Str.PISCO[i,4] <- r
  mod.emm <- emmeans(Lm.Ab, ~ time, at = list(time = c(0, 11)))
  Coef.table.Str.PISCO[i,5] <- summary(mod.emm)$emmean[1]
  Coef.table.Str.PISCO[i,6] <- summary(mod.emm)$lower.CL[1]
  Coef.table.Str.PISCO[i,7] <- summary(mod.emm)$upper.CL[1]
  Coef.table.Str.PISCO[i,8] <- summary(mod.emm)$emmean[2]
  Coef.table.Str.PISCO[i,9] <- summary(mod.emm)$lower.CL[2]
  Coef.table.Str.PISCO[i,10] <- summary(mod.emm)$upper.CL[2]
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(Str[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Den', 'N', 'Y','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Den', 'N', 'N','CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Den', 'N', 'N','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Den', 'N', 'Y','CI', 'NA')
  }    }


#Biomass##
Str <- subset(All.RR.sub.trans, y == "Strongylocentrotus purpuratus" & source == "PISCO" & resp == "Bio")
Str <- subset(Str, CA_MPA_Name_Short != 'Matlahuayl SMR')
x <- unique(Str$CA_MPA_Name_Short)
Coef.table.Str.PISCO <- data.frame(Site = character(), Mes.Slope = double(), Mes.p = double(), Mes.r = double(), 
                                   em0 = double(), l.CL0 = double(), u.CL0 = double(), em11 = double(), l.CL11 = double(), u.CL11 = double())
l <- length(x)
for(i in 1:l) {
  idx <- which(Str$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = Str[idx,])
  p <- summary(Lm.Ab)$coefficients[2,4]
  slope <- summary(Lm.Ab)$coefficients[2,1]
  r <- summary(Lm.Ab)$adj.r.squared
  slope = ifelse(p <= 0.05, slope, NA)
  Coef.table.Str.PISCO[i,1] <-  x[i]
  Coef.table.Str.PISCO[i,2] <- slope
  Coef.table.Str.PISCO[i,3] <- p
  Coef.table.Str.PISCO[i,4] <- r
  mod.emm <- emmeans(Lm.Ab, ~ time, at = list(time = c(0, 11)))
  Coef.table.Str.PISCO[i,5] <- summary(mod.emm)$emmean[1]
  Coef.table.Str.PISCO[i,6] <- summary(mod.emm)$lower.CL[1]
  Coef.table.Str.PISCO[i,7] <- summary(mod.emm)$upper.CL[1]
  Coef.table.Str.PISCO[i,8] <- summary(mod.emm)$emmean[2]
  Coef.table.Str.PISCO[i,9] <- summary(mod.emm)$lower.CL[2]
  Coef.table.Str.PISCO[i,10] <- summary(mod.emm)$upper.CL[2]
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(Str[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Bio', 'N', 'Y','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Bio', 'N', 'N','CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Bio', 'N', 'N','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Bio', 'N', 'Y','CI', 'NA')
  }}


##################################################################################################################
###############Macrocystis##################################################################################
Mac <- subset(All.RR.sub.trans, y == "Macrocystis pyrifera" & source == "PISCO" & resp == "Bio")
x <- unique(Mac$CA_MPA_Name_Short)
Coef.table.Mac.PISCO <- data.frame(Site = character(), Mes.Slope = double(), Mes.p = double(), Mes.r = double(), 
                                   em0 = double(), l.CL0 = double(), u.CL0 = double(), em11 = double(), l.CL11 = double(), u.CL11 = double())
l <- length(x)
for(i in 1:l) {
  idx <- which(Mac$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = Mac[idx,])
  p <- summary(Lm.Ab)$coefficients[2,4]
  slope <- summary(Lm.Ab)$coefficients[2,1]
  r <- summary(Lm.Ab)$adj.r.squared
  slope = ifelse(p <= 0.05, slope, NA)
  Coef.table.Mac.PISCO[i,1] <-  x[i]
  Coef.table.Mac.PISCO[i,2] <- slope
  Coef.table.Mac.PISCO[i,3] <- p
  Coef.table.Mac.PISCO[i,4] <- r
  mod.emm <- emmeans(Lm.Ab, ~ time, at = list(time = c(0, 11)))
  Coef.table.Mac.PISCO[i,5] <- summary(mod.emm)$emmean[1]
  Coef.table.Mac.PISCO[i,6] <- summary(mod.emm)$lower.CL[1]
  Coef.table.Mac.PISCO[i,7] <- summary(mod.emm)$upper.CL[1]
  Coef.table.Mac.PISCO[i,8] <- summary(mod.emm)$emmean[2]
  Coef.table.Mac.PISCO[i,9] <- summary(mod.emm)$lower.CL[2]
  Coef.table.Mac.PISCO[i,10] <- summary(mod.emm)$upper.CL[2]
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(Mac[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('M. pyrifera', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Bio', 'N', 'Y','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('M. pyrifera', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Bio', 'N', 'N','CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('M. pyrifera', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Bio', 'N', 'N','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('M. pyrifera', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Bio', 'N', 'Y','CI', 'NA')
  }  }

##################################################################################################################
########Lobster##############################################################################################
Pan.den <- subset(All.RR.sub.trans, y == "Panulirus interruptus" & source == "PISCO" & resp == "Den")
x <- unique(Pan.den$CA_MPA_Name_Short)
Coef.table.Pan.den.PISCO <- data.frame(Site = character(), Mes.Slope = double(), Mes.p = double(), Mes.r = double(), 
                                   em0 = double(), l.CL0 = double(), u.CL0 = double(), em11 = double(), l.CL11 = double(), u.CL11 = double())
l <- length(x)
for(i in 1:l) {
  idx <- which(Pan.den$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = Pan.den[idx,])
  p <- summary(Lm.Ab)$coefficients[2,4]
  slope <- summary(Lm.Ab)$coefficients[2,1]
  r <- summary(Lm.Ab)$adj.r.squared
  slope = ifelse(p <= 0.05, slope, NA)
  Coef.table.Pan.den.PISCO[i,1] <-  x[i]
  Coef.table.Pan.den.PISCO[i,2] <- slope
  Coef.table.Pan.den.PISCO[i,3] <- p
  Coef.table.Pan.den.PISCO[i,4] <- r
  mod.emm <- emmeans(Lm.Ab, ~ time, at = list(time = c(0, 11)))
  Coef.table.Pan.den.PISCO[i,5] <- summary(mod.emm)$emmean[1]
  Coef.table.Pan.den.PISCO[i,6] <- summary(mod.emm)$lower.CL[1]
  Coef.table.Pan.den.PISCO[i,7] <- summary(mod.emm)$upper.CL[1]
  Coef.table.Pan.den.PISCO[i,8] <- summary(mod.emm)$emmean[2]
  Coef.table.Pan.den.PISCO[i,9] <- summary(mod.emm)$lower.CL[2]
  Coef.table.Pan.den.PISCO[i,10] <- summary(mod.emm)$upper.CL[2]
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(Pan.den[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('P. interruptus', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Den', 'N', 'Y','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('P. interruptus', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Den', 'N', 'N','CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('P. interruptus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Den', 'N', 'N','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('P. interruptus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Den', 'N', 'Y','CI', 'NA')
  }    }


Pan.bio <- subset(All.RR.sub.trans, y == "Panulirus interruptus" & resp == "Bio" & source == "PISCO")
Pan.bio <- subset(Pan.bio, CA_MPA_Name_Short != 'Gull Island SMR' & CA_MPA_Name_Short != 'Harris Point SMR' &
                    CA_MPA_Name_Short != 'South Point SMR' & CA_MPA_Name_Short != 'Scorpion SMR' &
                    CA_MPA_Name_Short != 'Anacapa Island SMR 2003')
Pan.bio <- subset(Pan.bio, CA_MPA_Name_Short != "Santa Barbara Island SMR" & CA_MPA_Name_Short != "Point Dume SMCA")
x <- unique(Pan.bio$CA_MPA_Name_Short)
Coef.table.Pan.bio.all.PISCO <- data.frame(Site = character(), Mes.Slope = double(), Mes.p = double(), Mes.r = double(), 
                                       em0 = double(), l.CL0 = double(), u.CL0 = double(), em11 = double(), l.CL11 = double(), u.CL11 = double())
l <- length(x)
for(i in 1:l) {
  idx <- which(Pan.bio$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = Pan.bio[idx,])
  p <- summary(Lm.Ab)$coefficients[2,4]
  slope <- summary(Lm.Ab)$coefficients[2,1]
  r <- summary(Lm.Ab)$adj.r.squared
  slope = ifelse(p <= 0.05, slope, NA)
  Coef.table.Pan.bio.all.PISCO[i,1] <-  x[i]
  Coef.table.Pan.bio.all.PISCO[i,2] <- slope
  Coef.table.Pan.bio.all.PISCO[i,3] <- p
  Coef.table.Pan.bio.all.PISCO[i,4] <- r
  mod.emm <- emmeans(Lm.Ab, ~ time, at = list(time = c(0, 11)))
  Coef.table.Pan.bio.all.PISCO[i,5] <- summary(mod.emm)$emmean[1]
  Coef.table.Pan.bio.all.PISCO[i,6] <- summary(mod.emm)$lower.CL[1]
  Coef.table.Pan.bio.all.PISCO[i,7] <- summary(mod.emm)$upper.CL[1]
  Coef.table.Pan.bio.all.PISCO[i,8] <- summary(mod.emm)$emmean[2]
  Coef.table.Pan.bio.all.PISCO[i,9] <- summary(mod.emm)$lower.CL[2]
  Coef.table.Pan.bio.all.PISCO[i,10] <- summary(mod.emm)$upper.CL[2]
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(Pan.bio[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('P. interruptus', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Bio', 'N', 'Y','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('P. interruptus', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Bio', 'N', 'N','CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('P. interruptus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Bio', 'N', 'N','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('P. interruptus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Bio', 'N', 'Y','CI', 'NA')
  }   }

##################################################################################################################
########Sheephead##############################################################################################

Sheep.den <- subset(All.RR.sub.trans, y == "Semicossyphus pulcher" & resp == "Den" & source == "PISCO")
Sheep.den$CA_MPA_Name_Short <- as.character(Sheep.den$CA_MPA_Name_Short)
x <- unique(Sheep.den$CA_MPA_Name_Short)
l <- length(x)
for(i in 1:l) {
  idx <- which(Sheep.den$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = Sheep.den[idx,])
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(Sheep.den[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('S. pulcher', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Den', 'N', 'Y','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('S. pulcher', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Den', 'N', 'N','CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('S. pulcher', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Den', 'N', 'N','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('S. pulcher', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Den', 'N', 'Y','CI', 'NA')
  } }


Sheep.bio <- subset(All.RR.sub.trans, y == "Semicossyphus pulcher" & resp == "Bio" & source == "PISCO")
Sheep.bio$CA_MPA_Name_Short <- as.character(Sheep.bio$CA_MPA_Name_Short)
x <- unique(Sheep.bio$CA_MPA_Name_Short)
l <- length(x)
for(i in 1:l) {
  idx <- which(Sheep.bio$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = Sheep.bio[idx,])
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(Sheep.bio[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('S. pulcher', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Bio', 'N', 'Y','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('S. pulcher', x[i], mean, pSE, pSD, pCI, 'Linear', 'PISCO', 'Bio', 'N', 'N','CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('S. pulcher', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Bio', 'N', 'N','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('S. pulcher', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'PISCO', 'Bio', 'N', 'Y','CI', 'NA')
  } }


##Now on to KFM###########################################################################################
#Start with Purps#########################################################################################
####Density#######
KFM.purps <- subset(All.RR.sub, source == "KFM" & y == "Strongylocentrotus purpuratus" & resp == "Den")
KFM.purps$time.model <- NA #time.model sets all before to 0 and then all after sequential starting with 1
KFM.purps$time.true <- NA #time.true all time points are sequential
x <- unique(KFM.purps$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)

#Create time.model and time.true columns
for(i in 1:l) {
  idx <- which(KFM.purps$CA_MPA_Name_Short == x[i])
  min.id <- min(KFM.purps$year[idx])
  max.id <- max(KFM.purps$year[idx])
  KFM.purps$time.model[idx] <-c(rep(0,length(which(KFM.purps$BA[idx]=="Before"))),seq(0:length(which(KFM.purps$BA[idx]=="After"))))
  KFM.purps$time.true[idx] <-seq(1:length(idx))
}

LinearBefore = data.frame(pInt = NA, pSlope = NA) #initialize dataframe to hold pvalues for linear regressions

#These sites have BA data, so will run using pBACIPS
KFM.purps.BA <- subset(KFM.purps, CA_MPA_Name_Short == "Scorpion SMR" | CA_MPA_Name_Short == "Santa Barbara Island SMR" | 
                         CA_MPA_Name_Short == "Harris Point SMR" | CA_MPA_Name_Short == "Gull Island SMR")
#These sites have CI data, so will do linear regression type analyses
KFM.purps.CI <- subset(KFM.purps, CA_MPA_Name_Short == "Anacapa Island SMR 2003" | CA_MPA_Name_Short == "South Point SMR")

x <- unique(KFM.purps.BA$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)
for (j in x) {
  temp <- subset(KFM.purps.BA, KFM.purps.BA$CA_MPA_Name_Short == j, select = 
                   c(CA_MPA_Name_Short, mpa, reference, lnDiff, time.true, time.model, BA))
  temp = subset(temp, temp$BA == "Before")
  lmBefore <- lm(lnDiff ~ time.true, data = temp)
  pl <- summary(lmBefore)$coefficients[,4]
  LinearBefore <- rbind(LinearBefore, pl)
}

LinearBefore <- na.omit(LinearBefore) #Remove first line with no data
LinearBefore$site <- x #Join mpa info so we know which one we're looking at

#Santa Barbara Island SMR is only site not linear in before period########################
KFM.SB.purp <- filter(KFM.purps, CA_MPA_Name_Short == 'Santa Barbara Island SMR') 

## Fit a linear model
linear.Model<-lm(data = KFM.SB.purp, lnDiff ~ time.model)
anova(linear.Model)
summary(linear.Model)

step.Model<-lm(data = KFM.SB.purp, lnDiff ~ BA)
anova(linear.Model)
summary(linear.Model)

AIC.test<-AIC(linear.Model, step.Model) # could not fit Asymptotic or Sigmoid functions
AIC.test

mod1 <-lm(data = KFM.SB.purp, lnDiff ~ BA)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
ols_plot_resid_qq(mod1)
acf(mod1$residuals, type = "correlation")
qqPlot(mod1) #Plot residuals
leveneTest(lnDiff ~ as.factor(BA), data = KFM.SB.purp)

##Calculate effect size
spurp.lter.nap <- tidy(emmeans(mod1, ~ BA))
spurp.lter.nap$stdev <- spurp.lter.nap$std.error*sqrt(spurp.lter.nap$df+2)
pSD <- sqrt((((spurp.lter.nap$df[1]+1)*(spurp.lter.nap$stdev[1]^2))+
               ((spurp.lter.nap$df[2]+1)*(spurp.lter.nap$stdev[2]^2)))
            /(spurp.lter.nap$df[1]+2))
pSE <- pSD*sqrt((1/(spurp.lter.nap$df[1]+2)) + (1/(spurp.lter.nap$df[2]+2)))
pCI <- pSE*1.96
mean <- spurp.lter.nap$estimate[1]-spurp.lter.nap$estimate[2] #We need to calculate the effect size at each year
SumStats[201,] <- rbind('S. purpuratus', 'Santa Barbara Island SMR', mean, pSE, pSD, pCI, 'Step', 'KFM', 'Den', 'Y', 'Y','BACI', 'N')


#Now just CI for Purp den###################
x <- unique(KFM.purps.CI$CA_MPA_Name_Short)
l <- length(x)
for(i in 1:l) {
  idx <- which(KFM.purps.CI$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = KFM.purps.CI[idx,])
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(KFM.purps.CI[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Den', 'N', 'Y','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Den', 'N', 'N','CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'KFM', 'Den', 'N', 'N','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'KFM', 'Den', 'N', 'Y','CI', 'NA')
  }}


#Biomass#################
KFM.purps <- subset(All.RR.sub, source == "KFM" & y == "Strongylocentrotus purpuratus" & resp == "Bio")
KFM.purps$time.model <- NA #time.model sets all before to 0 and then all after sequential starting with 1
KFM.purps$time.true <- NA #time.true all time points are sequential
x <- unique(KFM.purps$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)
KFM.purps <- KFM.purps[order(as.numeric(KFM.purps$year),decreasing = FALSE), ]

#Create time.model and time.true columns
for(i in 1:l) {
  idx <- which(KFM.purps$CA_MPA_Name_Short == x[i])
  min.id <- min(KFM.purps$year[idx])
  max.id <- max(KFM.purps$year[idx])
  KFM.purps$time.model[idx] <-c(rep(0,length(which(KFM.purps$BA[idx]=="Before"))),seq(0:length(which(KFM.purps$BA[idx]=="After"))))
  KFM.purps$time.true[idx] <-seq(1:length(idx))
}

#These sites have CI data, so will do linear regression type analyses, note that all other sites with before data were linear in before for density so removed
KFM.purps.CI <- subset(KFM.purps, CA_MPA_Name_Short == "Anacapa Island SMR 2003" | CA_MPA_Name_Short == "South Point SMR")

#KFM SPurp CI###############
#Now just CI for KFM spurp

x <- unique(KFM.purps.CI$CA_MPA_Name_Short)
l <- length(x)
for(i in 1:l) {
  idx <- which(KFM.purps.CI$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = KFM.purps.CI[idx,])
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(KFM.purps.CI[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Bio', 'N', 'Y','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Bio', 'N', 'N','CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'KFM', 'Bio', 'N', 'N','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('S. purpuratus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'KFM', 'Bio', 'N', 'Y','CI', 'NA')
  }}




#Now reds#########################################################################################
#########Density##########
KFM.reds <- subset(All.RR.sub, source == "KFM" & y == "Mesocentrotus franciscanus" & resp == "Den")
KFM.reds$time.model <- NA #time.model sets all before to 0 and then all after sequential starting with 1
KFM.reds$time.true <- NA #time.true all time points are sequential
x <- unique(KFM.reds$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)
KFM.reds <- KFM.reds[order(as.numeric(KFM.reds$year),decreasing = FALSE), ]

#Create time.model and time.true columns
for(i in 1:l) {
  idx <- which(KFM.reds$CA_MPA_Name_Short == x[i])
  min.id <- min(KFM.reds$year[idx])
  max.id <- max(KFM.reds$year[idx])
  KFM.reds$time.model[idx] <-c(rep(0,length(which(KFM.reds$BA[idx]=="Before"))),seq(0:length(which(KFM.reds$BA[idx]=="After"))))
  KFM.reds$time.true[idx] <-seq(1:length(idx))
}

LinearBefore = data.frame(pInt = NA, pSlope = NA) #initialize dataframe to hold pvalues for linear regressions
#These sites have BA data, so will run using pBACIPS
KFM.reds.BA <- subset(KFM.reds, CA_MPA_Name_Short == "Scorpion SMR" | CA_MPA_Name_Short == "Santa Barbara Island SMR" | 
                         CA_MPA_Name_Short == "Harris Point SMR" | CA_MPA_Name_Short == "Gull Island SMR")
#These sites have CI data, so will do linear regression type analyses
KFM.reds.CI <- subset(KFM.reds, CA_MPA_Name_Short == "Anacapa Island SMR 2003" | CA_MPA_Name_Short == "South Point SMR")

x <- unique(KFM.reds.BA$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)


for (j in x) {
  temp <- subset(KFM.reds.BA, KFM.reds.BA$CA_MPA_Name_Short == j, select = 
                   c(CA_MPA_Name_Short, mpa, reference, lnDiff, time.true, time.model, BA))
  temp = subset(temp, temp$BA == "Before")
  lmBefore <- lm(lnDiff ~ time.true, data = temp)
  pl <- summary(lmBefore)$coefficients[,4]
  LinearBefore <- rbind(LinearBefore, pl)
}

LinearBefore <- na.omit(LinearBefore) #Remove first line with no data
LinearBefore$site <- x #Join mpa info so we know which one we're looking at


#KFM MFran Scorpion SMR###############################
KFM.mfran <- filter(KFM.reds, CA_MPA_Name_Short == 'Scorpion SMR') 

mod1 <-lm(data = KFM.mfran, lnDiff ~ BA)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = KFM.mfran)

##Calculate effect size
spurp.lter.nap <- tidy(emmeans(mod1, ~ BA))
spurp.lter.nap$stdev <- spurp.lter.nap$std.error*sqrt(spurp.lter.nap$df+2)
pSD <- sqrt((((spurp.lter.nap$df[1]+1)*(spurp.lter.nap$stdev[1]^2))+
               ((spurp.lter.nap$df[2]+1)*(spurp.lter.nap$stdev[2]^2)))
            /(spurp.lter.nap$df[1]+2))
pSE <- pSD*sqrt((1/(spurp.lter.nap$df[1]+2)) + (1/(spurp.lter.nap$df[2]+2)))
pCI <- pSE*1.96
mean <- spurp.lter.nap$estimate[1]-spurp.lter.nap$estimate[2] #We need to calculate the effect size at each year
SumStats[210,] <- rbind('M. franciscanus', 'Scorpion SMR', mean, pSE, pSD, pCI, 'Step', 'KFM', 'Den', 'Y', 'Y','BACI', 'N')



#KFM MFran Gull Island SMR######################
KFM.mfran <- filter(KFM.reds, CA_MPA_Name_Short == 'Gull Island SMR') 

step.Model<-lm(data = KFM.mfran, lnDiff ~ BA)
anova(linear.Model)
summary(linear.Model)

##Calculate effect size
spurp.lter.nap <- tidy(emmeans(step.Model, ~ BA))
spurp.lter.nap$stdev <- spurp.lter.nap$std.error*sqrt(spurp.lter.nap$df+2)
pSD <- sqrt((((spurp.lter.nap$df[1]+1)*(spurp.lter.nap$stdev[1]^2))+
               ((spurp.lter.nap$df[2]+1)*(spurp.lter.nap$stdev[2]^2)))
            /(spurp.lter.nap$df[1]+2))
pSE <- pSD*sqrt((1/(spurp.lter.nap$df[1]+2)) + (1/(spurp.lter.nap$df[2]+2)))
pCI <- pSE*1.96
mean <- spurp.lter.nap$estimate[1]-spurp.lter.nap$estimate[2] #We need to calculate the effect size at each year

SumStats[211,] <- rbind('M. franciscanus', 'Gull Island SMR', mean, pSE, pSD, pCI, 'Step', 'KFM', 'Den', 'Y', 'Y','BACI', 'N')


#KFM MFran Harris Point SMR#################
KFM.mfran <- subset(KFM.reds, CA_MPA_Name_Short == "Harris Point SMR")

## Fit a linear model
linear.Model<-lm(data = KFM.mfran, lnDiff ~ time.model)
anova(linear.Model)
summary(linear.Model)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(linear.Model, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(linear.Model, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year

SumStats[212,] <- rbind('M. franciscanus', 'Harris Point SMR', mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Den', 'Y', 'Y','pBACIPS', 'N')

#Now just CI for KFM spurp

x <- unique(KFM.reds.CI$CA_MPA_Name_Short)
l <- length(x)
for(i in 1:l) {
  idx <- which(KFM.reds.CI$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = KFM.reds.CI[idx,])
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(KFM.reds.CI[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('M. franciscanus', x[i], mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Den', 'N', 'Y','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('M. franciscanus', x[i], mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Den', 'N', 'N','CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('M. franciscanus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'KFM', 'Den', 'N', 'N','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('M. franciscanus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'KFM', 'Den', 'N', 'Y','CI', 'NA')
  }}


#####Biomass MFran######################
KFM.reds <- subset(All.RR.sub, source == "KFM" & y == "Mesocentrotus franciscanus" & resp == "Bio")
KFM.reds$time.model <- NA #time.model sets all before to 0 and then all after sequential starting with 1
KFM.reds$time.true <- NA #time.true all time points are sequential
x <- unique(KFM.reds$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)
KFM.reds <- KFM.reds[order(as.numeric(KFM.reds$year),decreasing = FALSE), ]

#Create time.model and time.true columns
for(i in 1:l) {
  idx <- which(KFM.reds$CA_MPA_Name_Short == x[i])
  min.id <- min(KFM.reds$year[idx])
  max.id <- max(KFM.reds$year[idx])
  KFM.reds$time.model[idx] <-c(rep(0,length(which(KFM.reds$BA[idx]=="Before"))),seq(0:length(which(KFM.reds$BA[idx]=="After"))))
  KFM.reds$time.true[idx] <-seq(1:length(idx))
}

LinearBefore = data.frame(pInt = NA, pSlope = NA) #initialize dataframe to hold pvalues for linear regressions
#These sites have BA data, so will run using pBACIPS
KFM.reds.BA <- subset(KFM.reds, CA_MPA_Name_Short == "Scorpion SMR" | CA_MPA_Name_Short == "Santa Barbara Island SMR" | 
                        CA_MPA_Name_Short == "Harris Point SMR" | CA_MPA_Name_Short == "Gull Island SMR")
#These sites have CI data, so will do linear regression type analyses
KFM.reds.CI <- subset(KFM.reds, CA_MPA_Name_Short == "Anacapa Island SMR 2003" | CA_MPA_Name_Short == "South Point SMR")

x <- unique(KFM.reds.BA$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)


for (j in x) {
  temp <- subset(KFM.reds.BA, KFM.reds.BA$CA_MPA_Name_Short == j, select = 
                   c(CA_MPA_Name_Short, mpa, reference, lnDiff, time.true, time.model, BA))
  temp = subset(temp, temp$BA == "Before")
  lmBefore <- lm(lnDiff ~ time.true, data = temp)
  pl <- summary(lmBefore)$coefficients[,4]
  LinearBefore <- rbind(LinearBefore, pl)
}

LinearBefore <- na.omit(LinearBefore) #Remove first line with no data
LinearBefore$site <- x #Join mpa info so we know which one we're looking at

KFM.reds.BA$lnmpa <- log(KFM.reds.BA$mpa)
KFM.reds.BA$lnreference <- log(KFM.reds.BA$reference)

x <- unique(KFM.reds.BA$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)
likelihood.poly <- data.frame(matrix(ncol = 5, nrow = 0)) #Initialize dataframe to put likelihood values from loop
pvals.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues
Norm.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues

for (j in x) {
  temp <- subset(KFM.reds.BA, KFM.reds.BA$CA_MPA_Name_Short == j, select = c(CA_MPA_Name_Short, lnmpa, 
                                                                             lnreference, lnDiff, time.true, time.model))
  ProgressiveChangeBACIPS(control=temp$lnreference,
                          impact=temp$lnmpa,
                          time.true=temp$time.true,
                          time.model=temp$time.model)
  ps <- summary(step.Model)[[1]][["Pr(>F)"]][1]
  pl <- summary(linear.Model)$coefficients[,4]
  psi <- summary(sigmoid.Model)$coefficients[,4]
  pa <- summary(asymptotic.Model)$coefficients[,4]
  pval <- c(ps, pl, psi, pa)
  likelihood.poly <- rbind(likelihood.poly, w)
  pvals.poly <- rbind(pvals.poly,pval)
  t <- shapiro.test(temp$lnDiff)
  Norm.poly <- rbind(Norm.poly, t$p.value)
}

#Rename column headers
colnames(likelihood.poly) <- c("Step","Linear","Asymptotic","Sigmoid")
colnames(pvals.poly) <- c("Step","Linear","Asymptotic M","Asymptotic B","Asymptotic L","Sigmoid M","Sigmoid B","Sigmoid L","Sigmoid K")
likelihood.poly <- cbind(likelihood.poly, x)

#KFM MFran Scorpion SMR#####################
KFM.mfran <- filter(KFM.reds.BA, CA_MPA_Name_Short == 'Scorpion SMR') 
mod1 <-lm(data = KFM.mfran, lnDiff ~ BA)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
ols_plot_resid_qq(mod1)
acf(mod1$residuals, type = "correlation")
qqPlot(mod1) #Plot residuals
leveneTest(lnDiff ~ as.factor(BA), data = KFM.mfran)

##Calculate effect size
spurp.lter.nap <- tidy(emmeans(mod1, ~ BA))
spurp.lter.nap$stdev <- spurp.lter.nap$std.error*sqrt(spurp.lter.nap$df+2)
pSD <- sqrt((((spurp.lter.nap$df[1]+1)*(spurp.lter.nap$stdev[1]^2))+
               ((spurp.lter.nap$df[2]+1)*(spurp.lter.nap$stdev[2]^2)))
            /(spurp.lter.nap$df[1]+2))
pSE <- pSD*sqrt((1/(spurp.lter.nap$df[1]+2)) + (1/(spurp.lter.nap$df[2]+2)))
pCI <- pSE*1.96
mean <- spurp.lter.nap$estimate[1]-spurp.lter.nap$estimate[2] #We need to calculate the effect size at each year

SumStats[217,] <- rbind('M. franciscanus', 'Scorpion SMR', mean, pSE, pSD, pCI, 'Step', 'KFM', 'Bio', 'Y', 'Y','BACI', 'N')


#KFM MFran Gull Island SMR####################
KFM.mfran <- filter(KFM.reds.BA, CA_MPA_Name_Short == 'Gull Island SMR') 
mod1 <-lm(data = KFM.mfran, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
ols_plot_resid_qq(mod1)
acf(mod1$residuals, type = "correlation")
qqPlot(mod1) #Plot residuals
leveneTest(lnDiff ~ as.factor(BA), data = KFM.mfran)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year

SumStats[218,] <- rbind('M. franciscanus', 'Gull Island SMR', mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Bio', 'Y', 'Y','pBACIPS', 'N')


#KFM MFran Harris Point SMR##############
KFM.mfran <- subset(KFM.reds.BA, CA_MPA_Name_Short == "Harris Point SMR")

#Need these to fit sigmoid function above
time.model.of.impact=max(which(KFM.mfran$time.model==0))
time.true = KFM.mfran$time.true

sigmoid.Model<-mySIGfun(delta=KFM.mfran$lnDiff,time.model=KFM.mfran$time.model)
Scoef <- coef(sigmoid.Model) #These are the coefficients for the model
summary(sigmoid.Model)
res <- nlsResiduals(sigmoid.Model) # Extract residuals to test for normality and autocorrelation
test.nlsResiduals(res) # Test for autocorrelation

#Gives fitted value and standard error around estimates
interval <- data.frame(predFit(sigmoid.Model, newdata = KFM.mfran, interval = "prediction"))
interval$se <- abs((interval$lwr - interval$fit)/1.96)
interval$time <- KFM.mfran$time.model

KFM.mfran.plot <- ggplot(interval, aes(x = time, y = fit)) + 
  geom_line(size=1) +
  geom_point(data = KFM.mfran, aes(x = time.model, y = lnDiff), size = 1.5, color = "purple") +
  labs(y="Log(RR)", x="Time elaspsed (year)") +
  theme_classic() 
KFM.mfran.plot 

##Calculate effect size
interval$stdev <- interval$se*sqrt(32)

pSD <- sqrt((((31)*(interval$stdev[1]^2))+
               ((31)*(interval$stdev[29]^2)))
            /(32))
pSE <- pSD*sqrt((1/(32)) + (1/(32)))
pCI <- pSE*1.96
mean <- interval$fit[29]-interval$fit[1] #We need to calculate the effect size at each year

SumStats[219,] <- rbind('M. franciscanus', 'Harris Point SMR', mean, pSE, pSD, pCI, 'Sigmoid', 'KFM', 'Bio', 'Y', 'Y','pBACIPS', 'N')



#KFM MFran CI###################
x <- unique(KFM.reds.CI$CA_MPA_Name_Short)
l <- length(x)
for(i in 1:l) {
  idx <- which(KFM.reds.CI$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = KFM.reds.CI[idx,])
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(KFM.reds.CI[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('M. franciscanus', x[i], mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Bio', 'N', 'Y','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('M. franciscanus', x[i], mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Bio', 'N', 'N','CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('M. franciscanus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'KFM', 'Bio', 'N', 'N','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('M. franciscanus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'KFM', 'Bio', 'N', 'Y','CI', 'NA')
  }}


##Now on to Macrocystis pyrifera Stipe################################################################
KFM.macro <- subset(All.RR.sub.trans, source == "KFM" & y == "Macrocystis pyrifera" & resp == "Bio")
KFM.macro$time.model <- NA #time.model sets all before to 0 and then all after sequential starting with 1
KFM.macro$time.true <- NA #time.true all time points are sequential
x <- unique(KFM.macro$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)
KFM.macro <- KFM.macro[order(as.numeric(KFM.macro$year),decreasing = FALSE), ]

#Create time.model and time.true columns
for(i in 1:l) {
  idx <- which(KFM.macro$CA_MPA_Name_Short == x[i])
  min.id <- min(KFM.macro$year[idx])
  max.id <- max(KFM.macro$year[idx])
  KFM.macro$time.model[idx] <-c(rep(0,length(which(KFM.macro$BA[idx]=="Before"))),seq(0:length(which(KFM.macro$BA[idx]=="After"))))
  KFM.macro$time.true[idx] <-seq(1:length(idx))
}

LinearBefore = data.frame(pInt = NA, pSlope = NA) #initialize dataframe to hold pvalues for linear regressions

#These sites have BA data, so will run using pBACIPS
KFM.macro.BA <- subset(KFM.macro, CA_MPA_Name_Short == "Scorpion SMR" | 
                        CA_MPA_Name_Short == "Harris Point SMR" | CA_MPA_Name_Short == "Gull Island SMR")
#These sites have CI data, so will do linear regression type analyses
KFM.macro.CI <- subset(KFM.macro, CA_MPA_Name_Short == "Anacapa Island SMR 2003" | CA_MPA_Name_Short == "South Point SMR" |
                       CA_MPA_Name_Short == "Santa Barbara Island SMR")

x <- unique(KFM.macro.BA$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)

for (j in x) {
  temp <- subset(KFM.macro.BA, KFM.macro.BA$CA_MPA_Name_Short == j, select = 
                   c(CA_MPA_Name_Short, mpa, reference, lnDiff, time.true, time.model, BA))
  temp = subset(temp, temp$BA == "Before")
  lmBefore <- lm(lnDiff ~ time.true, data = temp)
  pl <- summary(lmBefore)$coefficients[,4]
  LinearBefore <- rbind(LinearBefore, pl)
}

LinearBefore <- na.omit(LinearBefore) #Remove first line with no data
LinearBefore$site <- x #Join mpa info so we know which one we're looking at

#Harris Point is linear in before

KFM.macro.BA$lnmpa <- log(KFM.macro.BA$mpa)
KFM.macro.BA$lnreference <- log(KFM.macro.BA$reference)

x <- unique(KFM.macro.BA$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)
likelihood.poly <- data.frame(matrix(ncol = 5, nrow = 0)) #Initialize dataframe to put likelihood values from loop
pvals.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues
Norm.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues

for (j in x) {
  temp <- subset(KFM.macro.BA, KFM.macro.BA$CA_MPA_Name_Short == j, select = c(CA_MPA_Name_Short, lnmpa, 
                                                                       lnreference, lnDiff, time.true, time.model))
  ProgressiveChangeBACIPS(control=temp$lnreference,
                          impact=temp$lnmpa,
                          time.true=temp$time.true,
                          time.model=temp$time.model)
  ps <- summary(step.Model)[[1]][["Pr(>F)"]][1]
  pl <- summary(linear.Model)$coefficients[,4]
  psi <- summary(sigmoid.Model)$coefficients[,4]
  pa <- summary(asymptotic.Model)$coefficients[,4]
  pval <- c(ps, pl, psi, pa)
  likelihood.poly <- rbind(likelihood.poly, w)
  pvals.poly <- rbind(pvals.poly,pval)
  t <- shapiro.test(temp$lnDiff)
  Norm.poly <- rbind(Norm.poly, t$p.value)
}

#Rename column headers
colnames(likelihood.poly) <- c("Step","Linear","Asymptotic","Sigmoid")
colnames(pvals.poly) <- c("Step","Linear","Asymptotic M","Asymptotic B","Asymptotic L","Sigmoid M","Sigmoid B","Sigmoid L","Sigmoid K")
likelihood.poly <- cbind(likelihood.poly, x)


#KFM Macro Gull Island SMR######################
KFM.mfran <- filter(KFM.macro.BA, CA_MPA_Name_Short == 'Gull Island SMR') 
mod1 <-lm(data = KFM.mfran, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = KFM.mfran)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year

SumStats[224,] <- rbind('M. pyrifera', 'Gull Island SMR', mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Bio', 'Y', 'Y','pBACIPS', 'N')

##KFM Macro Scorpion SMR###############
KFM.macro.stipe <- filter(KFM.macro.BA, CA_MPA_Name_Short == 'Scorpion SMR') 

#Need these to fit sigmoid function above
time.model.of.impact=max(which(KFM.macro.stipe$time.model==0))
time.true = KFM.macro.stipe$time.true

sigmoid.Model<-mySIGfun(delta=KFM.macro.stipe$lnDiff,time.model=KFM.macro.stipe$time.model)
Scoef <- coef(sigmoid.Model) #These are the coefficients for the model
summary(sigmoid.Model)
res <- nlsResiduals(sigmoid.Model) # Extract residuals to test for normality and autocorrelation
test.nlsResiduals(res) # Test for autocorrelation

#Gives fitted value and standard error around estimates
interval <- data.frame(predFit(sigmoid.Model, newdata = KFM.macro.stipe, interval = "prediction"))
interval$se <- abs((interval$lwr - interval$fit)/1.96)
interval$time <- KFM.macro.stipe$time.model

KFM.macro.stipe.plot <- ggplot(interval, aes(x = time, y = fit)) + 
  geom_line(size=1) +
  geom_point(data = KFM.macro.stipe, aes(x = time.model, y = lnDiff), size = 1.5, color = "purple") +
  labs(y="Log(RR)", x="Time elaspsed (year)") +
  theme_classic() 
KFM.macro.stipe.plot 

##Calculate effect size
interval$stdev <- interval$se*sqrt(32)

pSD <- sqrt((((31)*(interval$stdev[1]^2))+
               ((31)*(interval$stdev[29]^2)))
            /(32))
pSE <- pSD*sqrt((1/(32)) + (1/(32)))
pCI <- pSE*1.96
mean <- interval$fit[29]-interval$fit[1] #We need to calculate the effect size at each year

SumStats[225,] <- rbind('M. pyrifera', 'Scorpion SMR', mean, pSE, pSD, pCI, 'Sigmoid', 'KFM', 'Bio', 'Y', 'Y','pBACIPS', 'N')


##KFM CI Macro#############################
x <- unique(KFM.macro.CI$CA_MPA_Name_Short)
l <- length(x)
for(i in 1:l) {
  idx <- which(KFM.macro.CI$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = KFM.macro.CI[idx,])
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(KFM.macro.CI[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('M. pyrifera', x[i], mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Bio', 'N', 'Y','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('M. pyrifera', x[i], mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Bio', 'N', 'N','CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('M. pyrifera', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'KFM', 'Bio', 'N', 'N','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('M. pyrifera', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'KFM', 'Bio', 'N', 'Y','CI', 'NA')
  }}


##Now on to Panulirus interruptus################################################################
KFM.lob <- subset(All.RR.sub.trans, source == "KFM" & y == "Panulirus interruptus")
KFM.lob$time.model <- NA #time.model sets all before to 0 and then all after sequential starting with 1
KFM.lob$time.true <- NA #time.true all time points are sequential
KFM.lob$CA_MPA_Name_Short <- as.character(KFM.lob$CA_MPA_Name_Short)
x <- unique(KFM.lob$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)
KFM.lob <- KFM.lob[order(as.numeric(KFM.lob$year),decreasing = FALSE), ]

#Create time.model and time.true columns
for(i in 1:l) {
  idx <- which(KFM.lob$CA_MPA_Name_Short == x[i])
  min.id <- min(KFM.lob$year[idx])
  max.id <- max(KFM.lob$year[idx])
  KFM.lob$time.model[idx] <-c(rep(0,length(which(KFM.lob$BA[idx]=="Before"))),seq(0:length(which(KFM.lob$BA[idx]=="After"))))
  KFM.lob$time.true[idx] <-seq(1:length(idx))
}

LinearBefore = data.frame(pInt = NA, pSlope = NA) #initialize dataframe to hold pvalues for linear regressions

#These sites have BA data, so will run using pBACIPS; Dropped Harris, no lobsters in before suggests. not good sampling to me
KFM.lob.BA <- subset(KFM.lob, CA_MPA_Name_Short == "Scorpion SMR" | CA_MPA_Name_Short == "Santa Barbara Island SMR" | 
                         CA_MPA_Name_Short == "Gull Island SMR")
#These sites have CI data, so will do linear regression type analyses
KFM.lob.CI <- subset(KFM.lob, CA_MPA_Name_Short == "Anacapa Island SMR 2003" | CA_MPA_Name_Short == "South Point SMR")

x <- unique(KFM.lob.BA$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)


for (j in x) {
  temp <- subset(KFM.lob.BA, KFM.lob.BA$CA_MPA_Name_Short == j, select = 
                   c(CA_MPA_Name_Short, mpa, reference, lnDiff, time.true, time.model, BA))
  temp = subset(temp, temp$BA == "Before")
  lmBefore <- lm(lnDiff ~ time.true, data = temp)
  pl <- summary(lmBefore)$coefficients[,4]
  LinearBefore <- rbind(LinearBefore, pl)
}

LinearBefore <- LinearBefore[-1,] #Remove first line with no data
LinearBefore$site <- x #Join mpa info so we know which one we're looking at

KFM.lob.BA$lnmpa <- log(KFM.lob.BA$mpa)
KFM.lob.BA$lnreference <- log(KFM.lob.BA$reference)

likelihood.poly <- data.frame(matrix(ncol = 5, nrow = 0)) #Initialize dataframe to put likelihood values from loop
pvals.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues
Norm.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues
KFM.lob.sub <- subset(KFM.lob.BA, CA_MPA_Name_Short != "Harris Point SMR" & CA_MPA_Name_Short != "Santa Barbara Island SMR") #Harris doesn't have any lobsters
x <- unique(KFM.lob.sub$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)

for (j in x) {
  temp <- subset(KFM.lob.sub, KFM.lob.sub$CA_MPA_Name_Short == j, select = c(CA_MPA_Name_Short, lnmpa, 
                                                                         lnreference, lnDiff, time.true, time.model))
  ProgressiveChangeBACIPS(control=temp$lnreference,
                          impact=temp$lnmpa,
                          time.true=temp$time.true,
                          time.model=temp$time.model)
  ps <- summary(step.Model)[[1]][["Pr(>F)"]][1]
  pl <- summary(linear.Model)$coefficients[,4]
  psi <- summary(sigmoid.Model)$coefficients[,4]
  pa <- summary(asymptotic.Model)$coefficients[,4]
  pval <- c(ps, pl, psi, pa)
  likelihood.poly <- rbind(likelihood.poly, w)
  pvals.poly <- rbind(pvals.poly,pval)
  t <- shapiro.test(temp$lnDiff)
  Norm.poly <- rbind(Norm.poly, t$p.value)
}

#Rename column headers
colnames(likelihood.poly) <- c("Step","Linear","Asymptotic","Sigmoid")
colnames(pvals.poly) <- c("Step","Linear","Asymptotic M","Asymptotic B","Asymptotic L","Sigmoid M","Sigmoid B","Sigmoid L","Sigmoid K")
likelihood.poly <- cbind(likelihood.poly, x)

#KFM Lobsters Gull Island SMR############################
KFM.macro.stipe <- filter(KFM.lob.sub, CA_MPA_Name_Short == 'Gull Island SMR') 


# Fit the sigmoid model
#Need these to fit sigmoid function above
time.model.of.impact=max(which(KFM.macro.stipe$time==0))
time.true = KFM.macro.stipe$time

sigmoid.Model<-mySIGfun(delta=KFM.macro.stipe$lnDiff,time.model=KFM.macro.stipe$time.model)
Scoef <- coef(sigmoid.Model) #These are the coefficients for the model
summary(sigmoid.Model)
res <- nlsResiduals(sigmoid.Model) # Extract residuals to test for normality and autocorrelation
test.nlsResiduals(res) # Test for autocorrelation

#Gives fitted value and standard error around estimates
interval <- data.frame(predFit(sigmoid.Model, newdata = KFM.macro.stipe, interval = "prediction"))
interval$se <- abs((interval$lwr - interval$fit)/1.96)
interval$time <- KFM.macro.stipe$time

KFM.macro.plot <- ggplot(interval, aes(x = time, y = fit)) + 
  geom_line(size=1) +
  geom_point(data = KFM.macro.stipe, aes(x = time, y = lnDiff), size = 1.5, color = "purple") +
  labs(y="Log(RR)", x="Time elaspsed (year)") +
  theme_classic() 
KFM.macro.plot 

##Calculate effect size
interval$stdev <- interval$se*sqrt(36)

pSD <- sqrt((((36+1)*(interval$stdev[1]^2))+
               ((36+1)*(interval$stdev[32]^2)))
            /(36+2))
pSE <- pSD*sqrt((1/(36+2)) + (1/(36+2)))
pCI <- pSE*1.96
mean <- interval$fit[32]-interval$fit[1] #We need to calculate the effect size at each year

SumStats[232,] <- rbind('P. interruptus', 'Gull Island SMR', mean, pSE, pSD, pCI, 'Sigmoid', 'KFM', 'Den', 'Y', 'Y','pBACIPS', 'N')



#KFM Lobsters Scorpion SMR###############################
KFM.lob.scorp <- filter(KFM.lob.sub, CA_MPA_Name_Short == 'Scorpion SMR') 

mod1 <-lm(data = KFM.lob.scorp, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = KFM.lob.scorp)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year

SumStats[233,] <- rbind('P. interruptus', 'Scorpion SMR', mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Den', 'Y', 'Y','pBACIPS', 'N')



########KFM Lobsters Santa Barbara Island SMR####################
#Now we'll fit and compete models for 
KFM.lobs <- subset(KFM.lob.BA, CA_MPA_Name_Short == "Santa Barbara Island SMR")

## Fit a linear model
linear.Model<-lm(data = KFM.lobs, lnDiff ~ time.model)
anova(linear.Model)
summary(linear.Model)

step.Model<-lm(data = KFM.lobs, lnDiff ~ BA)
anova(linear.Model)
summary(linear.Model)

time.model.of.impact=max(which(KFM.lobs$time.model==1))
time.true = KFM.lobs$time.true

#Fit the asymptotic model
#asymptotic.Model<-myASYfun(delta=KFM.lobs$diff,time.model=KFM.lobs$time.model)

# Fit the sigmoid model
#sigmoid.Model<-mySIGfun(delta=KFM.lobs$diff,time.model=KFM.lobs$time.model)

AIC.test<-AIC(linear.Model, step.Model) # could not fit Asymptotic or Sigmoid functions
AIC.test

##Calculate effect size for linear, since best fit
mod1 <-lm(data = KFM.lobs, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = KFM.lobs)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year
SumStats[234,] <- rbind('P. interruptus', 'Santa Barbara Island SMR', mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Den', 'Y', 'Y','pBACIPS', 'N')



##KFM CI Lobsters#################################################
KFM.lob.CI <- subset(KFM.lob.CI)

x <- unique(KFM.lob.CI$CA_MPA_Name_Short)
l <- length(x)
for(i in 1:l) {
  idx <- which(KFM.lob.CI$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = KFM.lob.CI[idx,])
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(KFM.lob.CI[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('P. interruptus', x[i], mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Den', 'N', 'Y','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('P. interruptus', x[i], mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Den', 'N', 'N','CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('P. interruptus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'KFM', 'Den', 'N', 'N','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('P. interruptus', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'KFM', 'Den', 'N', 'Y','CI', 'NA')
  }}


###Now Sheephead######################################################################################
Sheep.den <- subset(All.RR.sub.trans, y == "Semicossyphus pulcher" & resp == "Den" & source == "KFM")
Sheep.den$time.model <- NA #time.model sets all before to 0 and then all after sequential starting with 1
Sheep.den$time.true <- NA #time.true all time points are sequential
x <- unique(Sheep.den$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)
Sheep.den <- Sheep.den[order(as.numeric(Sheep.den$year),decreasing = FALSE), ]

#Create time.model and time.true columns
for(i in 1:l) {
  idx <- which(Sheep.den$CA_MPA_Name_Short == x[i])
  min.id <- min(Sheep.den$year[idx])
  max.id <- max(Sheep.den$year[idx])
  Sheep.den$time.model[idx] <-c(rep(0,length(which(Sheep.den$BA[idx]=="Before"))),seq(0:length(which(Sheep.den$BA[idx]=="After"))))
  Sheep.den$time.true[idx] <-seq(1:length(idx))
}

#These sites have BA data, so will run using pBACIPS
KFM.sheep.BA <- subset(Sheep.den, CA_MPA_Name_Short == "Scorpion SMR" | CA_MPA_Name_Short == "Santa Barbara Island SMR" | 
                         CA_MPA_Name_Short == "Harris Point SMR" | CA_MPA_Name_Short == "Gull Island SMR")
#These sites have CI data, so will do linear regression type analyses
KFM.sheep.CI <- subset(Sheep.den, CA_MPA_Name_Short == "Anacapa Island SMR 2003" | CA_MPA_Name_Short == "South Point SMR")


LinearBefore = data.frame(pInt = NA, pSlope = NA) #initialize dataframe to hold pvalues for linear regressions
x <- unique(KFM.sheep.BA$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)

for (j in x) {
  temp <- subset(KFM.sheep.BA, KFM.sheep.BA$CA_MPA_Name_Short == j, select = 
                   c(CA_MPA_Name_Short, mpa, reference, lnDiff, time.true, time.model, BA))
  temp = subset(temp, temp$BA == "Before")
  lmBefore <- lm(lnDiff ~ time.true, data = temp)
  pl <- summary(lmBefore)$coefficients[,4]
  LinearBefore <- rbind(LinearBefore, pl)
}

LinearBefore <- na.omit(LinearBefore) #Remove first line with no data
LinearBefore$site <- x #Join mpa info so we know which one we're looking at

#Need to assess models separately to figure out which models will and won't run

#KFM SPul Harris Point SMR#######################
KFM.lobs <- filter(KFM.sheep.BA, CA_MPA_Name_Short == 'Harris Point SMR') 

## Fit a step model
step.Model<- lm(data = KFM.lobs, lnDiff ~ BA)
anova(step.Model)
summary(step.Model)

## Fit a linear model
linear.Model<-lm(data = KFM.lobs, lnDiff ~ time.model)
anova(linear.Model)
summary(linear.Model)

AIC.test<-AIC(step.Model, linear.Model) # could not fit Asymptotic or Sigmoid functions
AIC.test

mod1 <-lm(data = KFM.lobs, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
qqPlot(mod1) #Plot residuals

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) ))
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year

SumStats[238,] <- rbind('S. pulcher', 'Harris Point SMR', mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Den', 'Y', 'Y','pBACIPS', 'N')


#KFM SPul Scorpion SMR############################
LTER.nap.SPUL.den <- filter(KFM.sheep.BA, CA_MPA_Name_Short == 'Scorpion SMR') 

## Fit a step model
step.Model<- lm(data = LTER.nap.SPUL.den, lnDiff ~ BA)
anova(step.Model)
summary(step.Model)

## Fit a linear model
linear.Model<-lm(data = LTER.nap.SPUL.den, lnDiff ~ time.model)
anova(linear.Model)
summary(linear.Model)

AIC.test<-AIC(step.Model, linear.Model) # could not fit Asymptotic or Sigmoid functions
AIC.test

mod1 <-lm(data = LTER.nap.SPUL.den, lnDiff ~ BA)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LTER.nap.SPUL.den)

##Calculate effect size
spurp.lter.nap <- tidy(emmeans(mod1, ~ BA))
spurp.lter.nap$stdev <- spurp.lter.nap$std.error*sqrt(spurp.lter.nap$df+2)
pSD <- sqrt((((spurp.lter.nap$df[1]+1)*(spurp.lter.nap$stdev[1]^2))+
               ((spurp.lter.nap$df[2]+1)*(spurp.lter.nap$stdev[2]^2)))
            /(spurp.lter.nap$df[1]+2))
pSE <- pSD*sqrt((1/(spurp.lter.nap$df[1]+2)) + (1/(spurp.lter.nap$df[2]+2)))
pCI <- pSE*1.96
mean <- spurp.lter.nap$estimate[1]-spurp.lter.nap$estimate[2] #We need to calculate the effect size at each year
SumStats[239,] <- rbind('S. pulcher', 'Scorpion SMR', mean, pSE, pSD, pCI, 'Step', 'KFM', 'Den', 'Y', 'Y','BACI', 'N')


##KFM SPul CI##############################
KFM.sheep.CI <- subset(KFM.sheep.CI)

x <- unique(KFM.sheep.CI$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)

for(i in 1:l) {
  idx <- which(KFM.sheep.CI$CA_MPA_Name_Short == x[i] )
  Lm.Ab <- lm(lnDiff ~ time, data = KFM.sheep.CI[idx,])
  before <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 0) )) 
  after <- tidy(emmeans(Lm.Ab, ~ time, at = list(time = 11) ))
  before$stdev <- before$std.error*sqrt(before$df+2)
  after$stdev <- after$std.error*sqrt(after$df+2)
  pSD <- sqrt((((before$df[1]+1)*(before$stdev[1]^2))+
                 ((after$df[1]+1)*(after$stdev[1]^2)))
              /(after$df[1]+2))
  pSE <- pSD*sqrt((1/(before$df[1]+2)) + (1/(after$df[1]+2)))
  pCI <- pSE*1.96
  mean <- after$estimate[1]-before$estimate[1] #We need to calculate the effect size at each year
  CP.mean <- summarySE(KFM.sheep.CI[idx,], measurevar = "lnDiff" )
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('S. pulcher', x[i], mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Den', 'N', 'Y','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('S. pulcher', x[i], mean, pSE, pSD, pCI, 'Linear', 'KFM', 'Den', 'N', 'N','CI', 'NA')
  }
  t <- nrow(SumStats)
  if (p <= 0.05) {
    SumStats[t+1,] <- rbind('S. pulcher', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'KFM', 'Den', 'N', 'N','CI', 'NA')
  } else {
    SumStats[t+1,] <- rbind('S. pulcher', x[i], CP.mean$lnDiff, CP.mean$se, CP.mean$sd, CP.mean$ci, 'Mean', 'KFM', 'Den', 'N', 'Y','CI', 'NA')
  }}
  

#######################################################################################################

#######################################################################################################
################LANDSAT analyses###################################################################################
fn <- read.csv("data/LANDSAT/MPA_Runs_new.csv") #Import TBell spreadsheet of annual maximum biomass

#Transforming to long format
bio.summarize.long <- gather(fn, Year, Biomass, -MPA,-status, -rep, -lat, -lon)
bio.summarize.long$year <-  rep(1984:2021, each = 60)
bio.summarize.long <- bio.summarize.long %>% arrange(MPA)

bio.ave.mpa <- bio.summarize.long %>%
  group_by(MPA, status, year) %>% 
  summarise_at(c("Biomass"), mean, na.rm = TRUE) %>%
  ungroup()

bio.sum.max.short <- bio.ave.mpa %>% 
  spread(status, Biomass)

bio.sum.max.short <- bio.sum.max.short[complete.cases(bio.sum.max.short), ]


#Back to long format
bio.sum.max.long <- gather(bio.sum.max.short, status, meanBio, -MPA,-year)
x <- unique(bio.sum.max.long$MPA) #Need to redo since sites have been removed
l = length(x)
bio.sum.max.long$Prop <- 0
bio.sum.max.long$PropCorr <- 0


#Turn maximum annual biomass into proportion of time series maximum biomass
for(i in 1:l) {
  idx <-which(x[i] == bio.sum.max.long$MPA)
  m <- max(bio.sum.max.long$meanBio[idx]) #Find the max of the mpa specific time series to use to create a proportion of the max
  bio.sum.max.long$Prop[idx] <- bio.sum.max.long$meanBio[idx]/m
  bio.sum.max.long$PropCorr[idx] <- bio.sum.max.long$Prop[idx] + 0.01 #Add the min value so that there are no zeroes in the dataset
}

bio.sum.max.long$taxon_name <- "Macrocystis pyrifera"
bio.sum.max.long.sub <- bio.sum.max.long[, colnames(bio.sum.max.long)[c(1,2,7,3,6)]]
bio.sum.max.LANDSAT.diff <- bio.sum.max.long.sub %>% 
  spread(status, PropCorr)
bio.sum.max.LANDSAT.diff$Diff <- bio.sum.max.LANDSAT.diff$mpa/bio.sum.max.LANDSAT.diff$reference #Calculate the response ratio
bio.sum.max.LANDSAT.diff$lnDiff <- log(bio.sum.max.LANDSAT.diff$Diff) #Log response ratio

bio.sum.max.LANDSAT.diff$BA <- NA
x <- unique(bio.sum.max.LANDSAT.diff$MPA) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)

for(i in 1:l) {
  idx <-which(x[i] == bio.sum.max.LANDSAT.diff$MPA)
  j <-which(x[i] == Site$CA_MPA_Name_Short)
  m <- length(idx)
  if (Site$MPA_Start[j] == 2012) {
    for (n in 1:m) {
      if (bio.sum.max.LANDSAT.diff$year[idx[n]] <=2012) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "Before"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2013) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2014) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2015) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2016) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2017) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2018) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2019) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2020) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2021) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2022) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2023) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2024) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      }}
    } else if (Site$MPA_Start[j] == 2003) {
    for (n in 1:m) {
      if (bio.sum.max.LANDSAT.diff$year[idx[n]] <=2003) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "Before"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2004) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2005) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2006) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2007) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2008) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2009) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2010) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2011) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2012) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2013) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2014) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2015) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2016) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2017) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2018) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2019) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2020) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2021) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2022) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2023) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      }}
    
    } else if (Site$MPA_Start[j] == 1978) {
    for (n in 1:m) {
      if (bio.sum.max.LANDSAT.diff$year[idx[n]] <= 1978) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "Before"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1982) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1983) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1984) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1985) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1987) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1988) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1989) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1990) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1991) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1992) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1993) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1994) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1995) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1996) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1997) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1998) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 1999) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2000) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2001) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2002) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2003) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2004) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2005) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2006) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2007) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2008) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2009) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2010) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2011) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2012) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2013) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2014) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2015) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2016) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2017) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2018) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2019) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2020) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2021) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2022) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      } else if (bio.sum.max.LANDSAT.diff$year[idx[n]] == 2023) {
        bio.sum.max.LANDSAT.diff$BA[idx[n]] = "After"
      }}}
  else {}
}


#Create binary column of whether mpa was protected or not
#Initialize Protection status
bio.sum.max.LANDSAT.diff$time.model <- NA #time.model sets all before to 0 and then all after sequential starting with 1
bio.sum.max.LANDSAT.diff$time.true <- NA #time.true all time points are sequential

x <- unique(bio.sum.max.LANDSAT.diff$MPA) #Need to redo since we're looking at total mpa by reference pairs
l = length(x)
bio.sum.max.LANDSAT.diff <- bio.sum.max.LANDSAT.diff[order(as.numeric(bio.sum.max.LANDSAT.diff$year),decreasing = FALSE), ]

#Create time.model and time.true columns
for(i in 1:l) {
  idx <- which(bio.sum.max.LANDSAT.diff$MPA == x[i])
  min.id <- min(bio.sum.max.LANDSAT.diff$year[idx])
  max.id <- max(bio.sum.max.LANDSAT.diff$year[idx])
  bio.sum.max.LANDSAT.diff$time.model[idx] <-c(rep(0,length(which(bio.sum.max.LANDSAT.diff$BA[idx]=="Before"))),seq(0:length(which(bio.sum.max.LANDSAT.diff$BA[idx]=="After"))))
  bio.sum.max.LANDSAT.diff$time.true[idx] <-seq(1:length(idx))
}

names(bio.sum.max.long)[1] <- "CA_MPA_Name_Short"
names(bio.sum.max.LANDSAT.diff)[1] <- "CA_MPA_Name_Short"

###############################################################################################################
########Run Individual Models##################################################################################
###############################################################################################################
LinearBefore = data.frame(pInt = NA, pSlope = NA) #initialize dataframe to hold pvalues for linear regressions
colnames(bio.sum.max.LANDSAT.diff)[1] <- "CA_MPA_Name_Short"

for (j in x) {
  temp <- subset(bio.sum.max.LANDSAT.diff, bio.sum.max.LANDSAT.diff$CA_MPA_Name_Short == j, select = 
                   c(CA_MPA_Name_Short, mpa, reference, lnDiff, time.true, time.model, BA))
  temp = subset(temp, temp$BA == "Before")
  lmBefore <- lm(lnDiff ~ time.true, data = temp)
  pl <- summary(lmBefore)$coefficients[,4]
  LinearBefore <- rbind(LinearBefore, pl)
}

LinearBefore <- na.omit(LinearBefore) #Remove first line with no data
LinearBefore$site <- x #Join mpa info so we know which one we're looking at



MPA_implement <- filter(Site, CA_MPA_Name_Short %in% x)
colnames(MPA_implement)[7] <- "CA_MPA_Nam"
x <- unique(bio.sum.max.LANDSAT.diff$MPA) #Need to redo since we're looking at total mpa by reference pairs
colnames(bio.sum.max.long)[1] <- "CA_MPA_Name_Short"


likelihood.poly <- data.frame(matrix(ncol = 5, nrow = 0)) #Initialize dataframe to put likelihood values from loop
pvals.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues
Norm.poly <- data.frame(matrix(ncol = 9, nrow = 0)) #Initialize dataframe for all pvalues

bio.sum.max.LANDSAT.diff$lnmpa <- log(bio.sum.max.LANDSAT.diff$mpa)
bio.sum.max.LANDSAT.diff$lnreference <- log(bio.sum.max.LANDSAT.diff$reference)
bio.sum.max.LANDSAT.diff.sub <- subset(bio.sum.max.LANDSAT.diff, CA_MPA_Name_Short != "Arrow Point to Lion Head Point SMCA" &
                                         CA_MPA_Name_Short != "Santa Barbara Island SMR" & CA_MPA_Name_Short != "Point Dume SMR" &
                                         CA_MPA_Name_Short != "Point Dume SMCA" & CA_MPA_Name_Short != "South La Jolla SMR" &
                                         CA_MPA_Name_Short != "South Point SMR" & CA_MPA_Name_Short != "Naples SMCA")
x <- unique(bio.sum.max.LANDSAT.diff.sub$CA_MPA_Name_Short) #Need to redo since we're looking at total mpa by reference pairs

for (j in x) {
  temp <- subset(bio.sum.max.LANDSAT.diff.sub, bio.sum.max.LANDSAT.diff.sub$CA_MPA_Name_Short == j, select = c(CA_MPA_Name_Short, lnmpa, 
                                                                         lnreference, lnDiff, time.true, time.model))
  ProgressiveChangeBACIPS(control=temp$lnreference,
                          impact=temp$lnmpa,
                          time.true=temp$time.true,
                          time.model=temp$time.model)
  ps <- summary(step.Model)[[1]][["Pr(>F)"]][1]
  pl <- summary(linear.Model)$coefficients[,4]
  psi <- summary(sigmoid.Model)$coefficients[,4]
  pa <- summary(asymptotic.Model)$coefficients[,4]
  pval <- c(ps, pl, psi, pa)
  likelihood.poly <- rbind(likelihood.poly, w)
  pvals.poly <- rbind(pvals.poly,pval)
  t <- shapiro.test(temp$lnDiff)
  Norm.poly <- rbind(Norm.poly, t$p.value)
}

#Rename column headers
colnames(likelihood.poly) <- c("Step","Linear","Asymptotic","Sigmoid")
colnames(pvals.poly) <- c("Step","Linear","Asymptotic M","Asymptotic B","Asymptotic L","Sigmoid M","Sigmoid B","Sigmoid L","Sigmoid K")
likelihood.poly <- cbind(likelihood.poly, x)

#Start with MPAs that didn't run through for loop########################################
#########################################################################################

#Arrow Point to Lion Head Point SMCA################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff, CA_MPA_Name_Short == 'Arrow Point to Lion Head Point SMCA') 

## Fit a step model
step.Model<- lm(data = LS.macro.SBI, lnDiff ~ BA)
anova(step.Model)
summary(step.Model)

## Fit a linear model
linear.Model<-lm(data = LS.macro.SBI, lnDiff ~ time.model)
anova(linear.Model)
summary(linear.Model)

AIC.test<-AIC(step.Model, linear.Model) # could not fit Asymptotic or Sigmoid functions
AIC.test

#step wins
mod1 <-lm(data = LS.macro.SBI, lnDiff ~ BA)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

##Calculate effect size
spurp.lter.nap <- tidy(emmeans(mod1, ~ BA))
spurp.lter.nap$stdev <- spurp.lter.nap$std.error*sqrt(spurp.lter.nap$df+2)
pSD <- sqrt((((spurp.lter.nap$df[1]+1)*(spurp.lter.nap$stdev[1]^2))+
               ((spurp.lter.nap$df[2]+1)*(spurp.lter.nap$stdev[2]^2)))
            /(spurp.lter.nap$df[1]+2))
pSE <- pSD*sqrt((1/(spurp.lter.nap$df[1]+2)) + (1/(spurp.lter.nap$df[2]+2)))
pCI <- pSE*1.96
mean <- spurp.lter.nap$estimate[1]-spurp.lter.nap$estimate[2] #We need to calculate the effect size at each year
w <- mean/pSD
SumStats[244,] <- rbind('M. pyrifera', 'Arrow Point to Lion Head Point SMCA', mean, pSE, pSD, pCI, 'Step', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')

#Point Dume SMR################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff, CA_MPA_Name_Short == 'Point Dume SMR') 

## Fit a step model
step.Model<- lm(data = LS.macro.SBI, lnDiff ~ BA)
anova(step.Model)
summary(step.Model)

## Fit a linear model
linear.Model<-lm(data = LS.macro.SBI, lnDiff ~ time.model)
anova(linear.Model)
summary(linear.Model)

AIC.test<-AIC(step.Model, linear.Model) # could not fit Asymptotic or Sigmoid functions
AIC.test

#step wins, barely, basically the same
mod1 <-lm(data = LS.macro.SBI, lnDiff ~ BA)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

##Calculate effect size
spurp.lter.nap <- tidy(emmeans(mod1, ~ BA))
spurp.lter.nap$stdev <- spurp.lter.nap$std.error*sqrt(spurp.lter.nap$df+2)
pSD <- sqrt((((spurp.lter.nap$df[1]+1)*(spurp.lter.nap$stdev[1]^2))+
               ((spurp.lter.nap$df[2]+1)*(spurp.lter.nap$stdev[2]^2)))
            /(spurp.lter.nap$df[1]+2))
pSE <- pSD*sqrt((1/(spurp.lter.nap$df[1]+2)) + (1/(spurp.lter.nap$df[2]+2)))
pCI <- pSE*1.96
mean <- spurp.lter.nap$estimate[1]-spurp.lter.nap$estimate[2] #We need to calculate the effect size at each year
w <- mean/pSD
SumStats[245,] <- rbind('M. pyrifera', 'Point Dume SMR', mean, pSE, pSD, pCI, 'Step', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')


#Point Dume SMCA################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff, CA_MPA_Name_Short == 'Point Dume SMCA') 

## Fit a step model
step.Model<- lm(data = LS.macro.SBI, lnDiff ~ BA)
anova(step.Model)
summary(step.Model)

## Fit a linear model
linear.Model<-lm(data = LS.macro.SBI, lnDiff ~ time.model)
anova(linear.Model)
summary(linear.Model)

AIC.test<-AIC(step.Model, linear.Model) # could not fit Asymptotic or Sigmoid functions
AIC.test

#step wins, barely, basically the same
mod1 <-lm(data = LS.macro.SBI, lnDiff ~ BA)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

##Calculate effect size
spurp.lter.nap <- tidy(emmeans(mod1, ~ BA))
spurp.lter.nap$stdev <- spurp.lter.nap$std.error*sqrt(spurp.lter.nap$df+2)
pSD <- sqrt((((spurp.lter.nap$df[1]+1)*(spurp.lter.nap$stdev[1]^2))+
               ((spurp.lter.nap$df[2]+1)*(spurp.lter.nap$stdev[2]^2)))
            /(spurp.lter.nap$df[1]+2))
pSE <- pSD*sqrt((1/(spurp.lter.nap$df[1]+2)) + (1/(spurp.lter.nap$df[2]+2)))
pCI <- pSE*1.96
mean <- spurp.lter.nap$estimate[1]-spurp.lter.nap$estimate[2] #We need to calculate the effect size at each year
w <- mean/pSD
SumStats[246,] <- rbind('M. pyrifera', 'Point Dume SMCA', mean, pSE, pSD, pCI, 'Step', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')

#South La Jolla SMR################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff, CA_MPA_Name_Short == 'South La Jolla SMR') 

## Fit a step model
step.Model<- lm(data = LS.macro.SBI, lnDiff ~ BA)
anova(step.Model)
summary(step.Model)

## Fit a linear model
linear.Model<-lm(data = LS.macro.SBI, lnDiff ~ time.model)
anova(linear.Model)
summary(linear.Model)

AIC.test<-AIC(step.Model, linear.Model) # could not fit Asymptotic or Sigmoid functions
AIC.test

#linear wins, barely, basically the same
mod1 <-lm(data = LS.macro.SBI, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year
w <- mean/pSD

SumStats[247,] <- rbind('M. pyrifera', 'South La Jolla SMR', mean, pSE, pSD, pCI, 'Linear', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')


#South Point SMR################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff, CA_MPA_Name_Short == 'South Point SMR') 

## Fit a step model
step.Model<- lm(data = LS.macro.SBI, lnDiff ~ BA)
anova(step.Model)
summary(step.Model)

## Fit a linear model
linear.Model<-lm(data = LS.macro.SBI, lnDiff ~ time.model)
anova(linear.Model)
summary(linear.Model)

AIC.test<-AIC(step.Model, linear.Model) # could not fit Asymptotic or Sigmoid functions
AIC.test

#step wins
mod1 <-lm(data = LS.macro.SBI, lnDiff ~ BA)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

##Calculate effect size
spurp.lter.nap <- tidy(emmeans(mod1, ~ BA))
spurp.lter.nap$stdev <- spurp.lter.nap$std.error*sqrt(spurp.lter.nap$df+2)
pSD <- sqrt((((spurp.lter.nap$df[1]+1)*(spurp.lter.nap$stdev[1]^2))+
               ((spurp.lter.nap$df[2]+1)*(spurp.lter.nap$stdev[2]^2)))
            /(spurp.lter.nap$df[1]+2))
pSE <- pSD*sqrt((1/(spurp.lter.nap$df[1]+2)) + (1/(spurp.lter.nap$df[2]+2)))
pCI <- pSE*1.96
mean <- spurp.lter.nap$estimate[1]-spurp.lter.nap$estimate[2] #We need to calculate the effect size at each year
SumStats[248,] <- rbind('M. pyrifera', 'South Point SMR', mean, pSE, pSD, pCI, 'Step', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')

#Abalone Cove SMCA######################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff.sub, CA_MPA_Name_Short == 'Abalone Cove SMCA') 

##Calculate effect size
mod1 <-lm(data = LS.macro.SBI, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year
w <- mean/pSD

SumStats[249,] <- rbind('M. pyrifera', 'Abalone Cove SMCA', mean, pSE, pSD, pCI, 'Linear', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')


#Point Vicente SMCA######################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff.sub, CA_MPA_Name_Short == 'Point Vicente SMCA') 

##Calculate effect size

mod1 <-lm(data = LS.macro.SBI, lnDiff ~ BA)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

spurp.lter.nap <- tidy(emmeans(mod1, ~ BA))
spurp.lter.nap$stdev <- spurp.lter.nap$std.error*sqrt(spurp.lter.nap$df+2)
pSD <- sqrt((((spurp.lter.nap$df[1]+1)*(spurp.lter.nap$stdev[1]^2))+
               ((spurp.lter.nap$df[2]+1)*(spurp.lter.nap$stdev[2]^2)))
            /(spurp.lter.nap$df[1]+2))
pSE <- pSD*sqrt((1/(spurp.lter.nap$df[1]+2)) + (1/(spurp.lter.nap$df[2]+2)))
pCI <- pSE*1.96
mean <- spurp.lter.nap$estimate[1]-spurp.lter.nap$estimate[2] #We need to calculate the effect size at each year
w <- mean/pSD
SumStats[nrow(SumStats)+1,] <- rbind('M. pyrifera', 'Point Vicente SMCA', mean, pSE, pSD, pCI, 'Step', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')

#Swamis SMCA#########################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff.sub, CA_MPA_Name_Short == 'Swamis SMCA') 

##Calculate effect size
mod1 <-lm(data = LS.macro.SBI, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year
w <- mean/pSD

SumStats[nrow(SumStats)+1,] <- rbind('M. pyrifera', 'Swamis SMCA', mean, pSE, pSD, pCI, 'Linear', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')

#Cabrillo SMR#########################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff.sub, CA_MPA_Name_Short == 'Cabrillo SMR') 

##Calculate effect size
mod1 <-lm(data = LS.macro.SBI, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year
w <- mean/pSD

SumStats[nrow(SumStats)+1,] <- rbind('M. pyrifera', 'Cabrillo SMR', mean, pSE, pSD, pCI, 'Linear', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')

#Campus Point SMCA######################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff.sub, CA_MPA_Name_Short == 'Campus Point SMCA') 

##Calculate effect size

mod1 <-lm(data = LS.macro.SBI, lnDiff ~ BA)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

spurp.lter.nap <- tidy(emmeans(mod1, ~ BA))
spurp.lter.nap$stdev <- spurp.lter.nap$std.error*sqrt(spurp.lter.nap$df+2)
pSD <- sqrt((((spurp.lter.nap$df[1]+1)*(spurp.lter.nap$stdev[1]^2))+
               ((spurp.lter.nap$df[2]+1)*(spurp.lter.nap$stdev[2]^2)))
            /(spurp.lter.nap$df[1]+2))
pSE <- pSD*sqrt((1/(spurp.lter.nap$df[1]+2)) + (1/(spurp.lter.nap$df[2]+2)))
pCI <- pSE*1.96
mean <- spurp.lter.nap$estimate[1]-spurp.lter.nap$estimate[2] #We need to calculate the effect size at each year
w <- mean/pSD
SumStats[nrow(SumStats)+1,] <- rbind('M. pyrifera', 'Campus Point SMCA', mean, pSE, pSD, pCI, 'Step', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')

#Carrington Pt SMR#########################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff.sub, CA_MPA_Name_Short == 'Carrington Pt SMR') 

##Calculate effect size
mod1 <-lm(data = LS.macro.SBI, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year
w <- mean/pSD

SumStats[nrow(SumStats)+1,] <- rbind('M. pyrifera', 'Carrington Pt SMR', mean, pSE, pSD, pCI, 'Linear', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')


#Gull Island SMR ########################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff.sub, CA_MPA_Name_Short == 'Gull Island SMR') 

mod1 <-lm(data = LS.macro.SBI, lnDiff ~ BA)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

##Calculate effect size
spurp.lter.nap <- tidy(emmeans(mod1, ~ BA))
spurp.lter.nap$stdev <- spurp.lter.nap$std.error*sqrt(spurp.lter.nap$df+2)
pSD <- sqrt((((spurp.lter.nap$df[1]+1)*(spurp.lter.nap$stdev[1]^2))+
               ((spurp.lter.nap$df[2]+1)*(spurp.lter.nap$stdev[2]^2)))
            /(spurp.lter.nap$df[1]+2))
pSE <- pSD*sqrt((1/(spurp.lter.nap$df[1]+2)) + (1/(spurp.lter.nap$df[2]+2)))
pCI <- pSE*1.96
mean <- spurp.lter.nap$estimate[1]-spurp.lter.nap$estimate[2] #We need to calculate the effect size at each year
w <- mean/pSD
SumStats[nrow(SumStats)+1,] <- rbind('M. pyrifera', 'Gull Island SMR', mean, pSE, pSD, pCI, 'Step', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')

#Harris Point SMR########################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff.sub, CA_MPA_Name_Short == 'Harris Point SMR') 

mod1 <-lm(data = LS.macro.SBI, lnDiff ~ BA)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

##Calculate effect size
spurp.lter.nap <- tidy(emmeans(mod1, ~ BA))
spurp.lter.nap$stdev <- spurp.lter.nap$std.error*sqrt(spurp.lter.nap$df+2)
pSD <- sqrt((((spurp.lter.nap$df[1]+1)*(spurp.lter.nap$stdev[1]^2))+
               ((spurp.lter.nap$df[2]+1)*(spurp.lter.nap$stdev[2]^2)))
            /(spurp.lter.nap$df[1]+2))
pSE <- pSD*sqrt((1/(spurp.lter.nap$df[1]+2)) + (1/(spurp.lter.nap$df[2]+2)))
pCI <- pSE*1.96
mean <- spurp.lter.nap$estimate[1]-spurp.lter.nap$estimate[2] #We need to calculate the effect size at each year
w <- mean/pSD
SumStats[nrow(SumStats)+1,] <- rbind('M. pyrifera', 'Harris Point SMR', mean, pSE, pSD, pCI, 'Step', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')

#Judith Rk SMR########################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff.sub, CA_MPA_Name_Short == 'Judith Rk SMR') 

mod1 <-lm(data = LS.macro.SBI, lnDiff ~ BA)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

##Calculate effect size
spurp.lter.nap <- tidy(emmeans(mod1, ~ BA))
spurp.lter.nap$stdev <- spurp.lter.nap$std.error*sqrt(spurp.lter.nap$df+2)
pSD <- sqrt((((spurp.lter.nap$df[1]+1)*(spurp.lter.nap$stdev[1]^2))+
               ((spurp.lter.nap$df[2]+1)*(spurp.lter.nap$stdev[2]^2)))
            /(spurp.lter.nap$df[1]+2))
pSE <- pSD*sqrt((1/(spurp.lter.nap$df[1]+2)) + (1/(spurp.lter.nap$df[2]+2)))
pCI <- pSE*1.96
mean <- spurp.lter.nap$estimate[1]-spurp.lter.nap$estimate[2] #We need to calculate the effect size at each year
w <- mean/pSD
SumStats[nrow(SumStats)+1,] <- rbind('M. pyrifera', 'Judith Rk SMR', mean, pSE, pSD, pCI, 'Step', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')

#Naples SMCA- Exclude since linear in before##########################

#Point Conception SMR#########################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff.sub, CA_MPA_Name_Short == 'Point Conception SMR') 

##Calculate effect size
mod1 <-lm(data = LS.macro.SBI, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year
w <- mean/pSD

SumStats[nrow(SumStats)+1,] <- rbind('M. pyrifera', 'Point Conception SMR', mean, pSE, pSD, pCI, 'Linear', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')

#San Miguel Island SC#########################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff.sub, CA_MPA_Name_Short == 'San Miguel Island SC') 

##Calculate effect size
mod1 <-lm(data = LS.macro.SBI, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year
w <- mean/pSD

SumStats[nrow(SumStats)+1,] <- rbind('M. pyrifera', 'San Miguel Island SC', mean, pSE, pSD, pCI, 'Linear', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')

#Skunk Pt SMR#########################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff.sub, CA_MPA_Name_Short == 'Skunk Pt SMR') 

#Need these to fit sigmoid function above
time.model.of.impact=max(which(LS.macro.SBI$time.model==0))
time.true = LS.macro.SBI$time.true

sigmoid.Model<-mySIGfun(delta=LS.macro.SBI$lnDiff,time.model=LS.macro.SBI$time.model)
Scoef <- coef(sigmoid.Model) #These are the coefficients for the model
summary(sigmoid.Model)
res <- nlsResiduals(sigmoid.Model) # Extract residuals to test for normality and autocorrelation
test.nlsResiduals(res) # Test for autocorrelation

#Gives fitted value and standard error around estimates
interval <- data.frame(predFit(sigmoid.Model, newdata = LS.macro.SBI, interval = "prediction"))
interval$se <- abs((interval$lwr - interval$fit)/1.96)
interval$time <- LS.macro.SBI$time.model

LS.macro.SBI.plot <- ggplot(interval, aes(x = time, y = fit)) + 
  geom_line(size=1) +
  geom_point(data = LS.macro.SBI, aes(x = time.model, y = lnDiff), size = 1.5, color = "purple") +
  labs(y="Log(RR)", x="Time elaspsed (year)") +
  theme_classic() 
LS.macro.SBI.plot 

##Calculate effect size
interval$stdev <- interval$se*sqrt(32)

pSD <- sqrt((((31)*(interval$stdev[1]^2))+
               ((31)*(interval$stdev[29]^2)))
            /(32))
pSE <- pSD*sqrt((1/(32)) + (1/(32)))
pCI <- pSE*1.96
mean <- interval$fit[29]-interval$fit[1] #We need to calculate the effect size at each year

SumStats[nrow(SumStats)+1,] <- rbind('M. pyrifera', 'Skunk Pt SMR', mean, pSE, pSD, pCI, 'Sigmoid', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')

#Farnsworth Onshore SMCA#########################################################################
LS.macro.SBI <- filter(bio.sum.max.LANDSAT.diff.sub, CA_MPA_Name_Short == 'Farnsworth Onshore SMCA') 

##Calculate effect size
mod1 <-lm(data = LS.macro.SBI, lnDiff ~ time.model)
mod1_an <- anova(mod1)
mod1_an
summary(mod1)
leveneTest(lnDiff ~ as.factor(BA), data = LS.macro.SBI)

##Calculate effect size
mpyr.lter.before <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 0) )) 
mpyr.lter.after <- tidy(emmeans(mod1, ~ time.model, at = list(time.model = 11) ))
mpyr.lter.before$stdev <- mpyr.lter.before$std.error*sqrt(mpyr.lter.before$df+2)
mpyr.lter.after$stdev <- mpyr.lter.after$std.error*sqrt(mpyr.lter.after$df+2)

pSD <- sqrt((((mpyr.lter.before$df[1]+1)*(mpyr.lter.before$stdev[1]^2))+
               ((mpyr.lter.after$df[1]+1)*(mpyr.lter.after$stdev[1]^2)))
            /(mpyr.lter.after$df[1]+2))
pSE <- pSD*sqrt((1/(mpyr.lter.before$df[1]+2)) + (1/(mpyr.lter.after$df[1]+2)))
pCI <- pSE*1.96
mean <- mpyr.lter.after$estimate[1]-mpyr.lter.before$estimate[1] #We need to calculate the effect size at each year
w <- mean/pSD

SumStats[nrow(SumStats)+1,] <- rbind('M. pyrifera', 'Farnsworth Onshore SMCA', mean, pSE, pSD, pCI, 'Linear', 'Landsat', 'Bio', 'Y', 'Y','pBACIPS', 'N')

#######################################################################################################

#######################################################################################################
#Updated Forest Plots####################################################################
pal <- wes_palette("Zissou1", 20, type = "continuous")

SumStats$Mean <- as.numeric(SumStats$Mean)
SumStats$SE <- as.numeric(SumStats$SE)
SumStats$CI <- as.numeric(SumStats$CI)
SumStats$SD <- as.numeric(SumStats$SD)

SumStats.sub <- SumStats[complete.cases(SumStats), ]


SumStats.sub$Taxa <- factor(SumStats.sub$Taxa, levels = c("S. purpuratus", "M. franciscanus", "M. pyrifera", "P. interruptus", "S. pulcher"), ordered = TRUE)
SumStats.sub$Resp <- factor(SumStats.sub$Resp, levels = c("Bio","Den"), ordered = TRUE)
SumStats.sub$Source <- factor(SumStats.sub$Source, levels = c("KFM","LTER","PISCO","Landsat"), ordered = TRUE)
SumStats.sub$Type <- factor(SumStats.sub$Type, levels = c("pBACIPS","BACI","CI"), ordered = TRUE)

SumStats.sub <- merge(SumStats.sub, Site, by.x = "MPA", by.y = "CA_MPA_Name_Short")

######Figure 3###########################################################################
SumStats.Final <- subset(SumStats.sub, Type.x == "pBACIPS" | Type.x == "CI" & LinearBefore == "NA" & Primary == "Y")

#Remove Painted Cave SMCA (lobster take allowed), San Miguel Island SC (weird location), Arrow Point (finfish fishing allowed),
#Judith Rk SMR overlaps San Miguel Island SC, Swamis potentially affected by landslides as per Dan Pondella (need to follow up)
SumStats.Final.edit <- subset(SumStats.Final, MPA != "Painted Cave SMCA" & MPA != "San Miguel Island SC"
                              & MPA != "Arrow Point to Lion Head Point SMCA" & MPA != "Judith Rk SMR"
                              & MPA != 'Point Conception SMR')
fp2 <- ggplot(data=SumStats.Final.edit, aes(x=MPA, y=Mean, ymin=Mean-CI, ymax=Mean+CI, shape = Source, color = Resp)) +
  geom_pointrange() + 
  geom_point(size = 4) +
  facet_wrap(~Taxa, ncol = 2, scales = "free") +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Marine Protected Area") + ylab("Effect Size (LnRR)") +
  theme_classic() +
  theme(strip.text = element_text(face = "italic")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16)) +
  theme(legend.position = c(0.9, 0.1)) +
  guides(shape = guide_legend(title = "Source"), color = guide_legend(title = "Response")) + 
  scale_color_manual(values = c(pal[7], pal[16]), labels = c("Density", "Biomass")) +
  scale_shape_manual(labels = c("NPS-KFM", "LTER", "PISCO", "LANDSAT"), values = c(15,16,17,18))
fp2


ggsave(plot = fp2, file = "plots/MPA_ES.png", 
       device = "png",  bg = "white",
       width = 28, height = 30, units = "cm", dpi = 300)

SumStats.Final$Resp <- factor(SumStats.Final$Resp, levels = c("Bio", "Den"), ordered = TRUE)
MeanES <- summarySE(data = SumStats.Final, measurevar = "Mean", groupvars = c("Taxa", "Resp"))
MeanES$Taxa <- factor(MeanES$Taxa, levels = c("S. purpuratus", "M. franciscanus", "M. pyrifera", "P. interruptus", "P. interruptus legal", "S. pulcher", "S. pulcher legal"), ordered = TRUE)
MeanES$Resp <- factor(MeanES$Resp, levels = c("Bio", "Den"), ordered = TRUE)
MeanES$cat <- paste(MeanES$Taxa, MeanES$Resp)

##############Figure 2######################################################################
KFM.purps <- subset(All.RR.sub.trans, CA_MPA_Name_Short == "Scorpion SMR" & source == "KFM" & y == "Strongylocentrotus purpuratus")
x <- unique(KFM.purps$CA_MPA_Name_Short) #Need to redo since sites have been removed
MPA_implement.sub <- filter(MPA_implement, Site %in% x)
kfm.ave.purps <- subset(kfm.ave.ave, y == "Strongylocentrotus purpuratus" & CA_MPA_Name_Short == "Scorpion SMR")

pal <- wes_palette("Zissou1", 20, type = "continuous")

kfm.purps.raw <- ggplot(data=kfm.ave.purps, aes(x=year, y=den, color = status, fill = status)) +
  geom_point(size = 1.5) + 
  geom_line(size = 1) +
  geom_vline(xintercept = 2003, linetype = "dashed") +
  ylab(bquote('Density (ind '~m^-2*')')) + 
  xlab("Years") + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  scale_color_manual(name = "Site", values = c(pal[10], pal[18])) +
  scale_fill_manual(name = "Site", values = c(pal[10], pal[18])) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position ="top", legend.key.width = unit(1,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank())
kfm.purps.raw

KFM.purps.long <- gather(KFM.purps, status, value, 4:5)
KFM.purps.long <- subset(KFM.purps.long, resp == "Den" )

KFM.purps.prop <- ggplot(data=KFM.purps.long, aes(x=year, y=value, color = status, fill = status)) +
  geom_point(size = 1.5) + 
  geom_line(size = 1) +
  geom_vline(xintercept = 2003, linetype = "dashed") +
  ylab("Proportion of maximum") + 
  xlab("Years") + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  scale_color_manual(name = "Site", values = c(pal[10], pal[18])) +
  scale_fill_manual(name = "Site", values = c(pal[10], pal[18])) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position ="top", legend.key.width = unit(1,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank())
KFM.purps.prop

subsetPurp <- subset(KFM.purps, BA == "After")
KFM.purps$BA <- factor(KFM.purps$BA, levels = c("Before", "After"), ordered = TRUE)

KFM.purps.diff <- ggplot(data=KFM.purps, aes(x=year, y=lnDiff, shape = BA)) +
  geom_point(size = 1.5, color = pal[7]) + 
  #geom_smooth(data = subsetPurp, aes(x = year, y = lnDiff), size = 1, method = "lm", fill = pal[7], color = pal[7]) +
  geom_vline(xintercept = 2003, linetype = "dashed") +
  ylab("Log(RR)") + 
  xlab("Years") + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position ="top", legend.key.width = unit(1,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank()) 
KFM.purps.diff

KFM.purps.diff2 <- ggplot(data=KFM.purps, aes(x=year, y=lnDiff, shape = BA)) +
  geom_point(size = 1.5, color = pal[7]) + 
  geom_vline(xintercept = 2003, linetype = "dashed") +
  geom_smooth(data = subsetPurp, aes(x = year, y = lnDiff), size = 1, method = "lm", color = pal[7], se = FALSE) +
  ylab("Log(RR)") + 
  xlab("Years") + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position ="top", legend.key.width = unit(1,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank()) 
KFM.purps.diff2

KFM.purp.fig <- ggarrange(kfm.purps.raw, KFM.purps.prop, KFM.purps.diff, KFM.purps.diff2, labels = c("(a)","(b)","(c)","(d)"),  
                          hjust = 0.05, ncol = 4, align = "hv")
KFM.purp.fig
ggsave(plot = KFM.purp.fig, file = "plots/Fig2_Purp_KFM_Example.png", 
       device = "png",  bg = "white",
       width = 40, height = 10, units = "cm", dpi = 300)

#######Fig.1########################################################################################

CP.purps <- subset(All.Resp.sub, CA_MPA_Name_Short == "Campus Point SMCA" & taxon_name == "Macrocystis pyrifera" & resp == "Bio")

Purps.CP <- ggplot(data=CP.purps, aes(x=year, y=value, color = status, fill = status, shape = source)) +
  geom_point(size = 2) + 
  geom_line(size = 1) +
  ylab(bquote('Biomass (g '~m^-2*')')) + 
  xlab("Years") + 
  geom_vline(xintercept = 2012, linetype = "dashed") +
  scale_color_manual(name = "Site", values = c(pal[10], pal[18])) +
  scale_fill_manual(name = "Site", values = c(pal[10], pal[18])) +
  scale_shape_manual(name = "Source", labels = c("LTER", "PISCO"), values = c(16,17)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.key.width = unit(0.5,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank(),
        axis.text=element_text(size=14), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
Purps.CP

PV.purps <- subset(All.Resp.sub, CA_MPA_Name_Short == "Point Vicente SMCA" & taxon_name == "Macrocystis pyrifera" & resp == "Bio")

Purps.PV <- ggplot(data=PV.purps, aes(x=year, y=value, color = status, fill = status, shape = source)) +
  geom_point(size = 2) + 
  geom_line(size = 1) +
  ylab(bquote('Biomass (g '~m^-2*')')) + 
  xlab("Years") + 
  geom_vline(xintercept = 2012, linetype = "dashed") +
  scale_color_manual(name = "Site", values = c(pal[10], pal[18])) +
  scale_fill_manual(name = "Site", values = c(pal[10], pal[18])) +
  scale_shape_manual(name = "Source", labels = c("PISCO"), values = c(17)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.key.width = unit(0.5,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank(),
        axis.text=element_text(size=14), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
Purps.PV

HI.purps <- subset(All.Resp.sub, CA_MPA_Name_Short == "Harris Point SMR" & taxon_name == "Macrocystis pyrifera" & resp == "Bio")

Purps.HI <- ggplot(data=HI.purps, aes(x=year, y=value, color = status, fill = status, shape = source)) +
  geom_point(size = 2) + 
  geom_line(size = 1) +
  ylab(bquote('Biomass (g '~m^-2*')')) + 
  xlab("Years") + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  scale_color_manual(name = "Site", values = c(pal[10], pal[18])) +
  scale_fill_manual(name = "Site", values = c(pal[10], pal[18])) +
  scale_shape_manual(name = "Source", labels = c("KFM", "PISCO"), values = c(15,17)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.key.width = unit(0.5,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank(),
        axis.text=element_text(size=14), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
Purps.HI

S.purps <- subset(All.Resp.sub, CA_MPA_Name_Short == "Santa Barbara Island SMR" & taxon_name == "Macrocystis pyrifera" & resp == "Bio")

Purps.S <- ggplot(data=S.purps, aes(x=year, y=value, color = status, fill = status, shape = source)) +
  geom_point(size = 2) + 
  geom_line(size = 1) +
  ylab(bquote('Biomass (g '~m^-2*')')) + 
  xlab("Years") + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  scale_color_manual(name = "Site", values = c(pal[10], pal[18])) +
  scale_fill_manual(name = "Site", values = c(pal[10], pal[18])) +
  scale_shape_manual(name = "Source", labels = c("KFM", "PISCO"), values = c(15,17)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.key.width = unit(0.5,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank(),
        axis.text=element_text(size=14), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
Purps.S

GI.purps <- subset(All.Resp.sub, CA_MPA_Name_Short == "Gull Island SMR" & taxon_name == "Macrocystis pyrifera" & resp == "Bio")

Purps.GI <- ggplot(data=GI.purps, aes(x=year, y=value, color = status, fill = status, shape = source)) +
  geom_point(size = 2) + 
  geom_line(size = 1) +
  ylab(bquote('Biomass (g '~m^-2*')')) + 
  xlab("Years") + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  scale_color_manual(name = "Site", values = c(pal[10], pal[18])) +
  scale_fill_manual(name = "Site", values = c(pal[10], pal[18])) +
  scale_shape_manual(name = "Source", labels = c("KFM", "PISCO"), values = c(15,17)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.key.width = unit(0.5,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank(),
        axis.text=element_text(size=14), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
Purps.GI

SP.purps <- subset(All.Resp.sub, CA_MPA_Name_Short == "South Point SMR" & taxon_name == "Macrocystis pyrifera" & resp == "Bio")

Purps.SP <- ggplot(data=SP.purps, aes(x=year, y=value, color = status, fill = status, shape = source)) +
  geom_point(size = 2) + 
  geom_line(size = 1) +
  ylab(bquote('Biomass (g '~m^-2*')')) + 
  xlab("Years") + 
  geom_vline(data = MPA_implement.sub, mapping = aes(xintercept = MPA_Start), linetype = "dashed") +
  scale_color_manual(name = "Site", values = c(pal[10], pal[18])) +
  scale_fill_manual(name = "Site", values = c(pal[10], pal[18])) +
  scale_shape_manual(name = "Source", labels = c("KFM", "PISCO"), values = c(15,17)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.key.width = unit(0.5,'cm'), legend.title.align=0.5,
        legend.box.background = element_rect(colour = "black"), legend.background = element_blank(),
        axis.text=element_text(size=14), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
Purps.SP


Purps.plots <- ggarrange(Purps.CP, Purps.PV, Purps.HI, Purps.S, Purps.GI, Purps.SP,
                          hjust = 0.05, ncol = 1, nrow = 6, align = "hv")
Purps.plots
ggsave(plot = Purps.plots, file = "plots/Fig1_plots.pdf", 
       device = "pdf",  bg = "white",
       width = 20, height = 60, units = "cm", dpi = 300)


############################################################################################
###########Whole Model w/data source nested in mpa as random effect#########################
############################################################################################
SumStats.Final.edit <- subset(SumStats.sub, Type.x == "pBACIPS" | Type.x == "CI" & LinearBefore == "NA" & Primary == "Y")

#Remove Painted Cave SMCA (lobster take allowed), San Miguel Island SC (weird location), Arrow Point (finfish fishing allowed),
#Judith Rk SMR overlaps San Miguel Island SC, Swamis potentially affected by landslides as per Dan Pondella (need to follow up)
SumStats.Final.edit <- subset(SumStats.Final.edit, MPA != "Painted Cave SMCA" & MPA != "San Miguel Island SC"
                              & MPA != "Arrow Point to Lion Head Point SMCA" & MPA != "Judith Rk SMR"  
                              & MPA != 'Point Conception SMR')

SumStats.Final.edit.den <- subset(SumStats.Final.edit, Resp == "Den")
SumStats.Final.edit.den$ES_ID <- seq.int(nrow(SumStats.Final.edit.den))
SumStats.Final.edit.den$vi <- SumStats.Final.edit.den$SD^2
SumStats.Final.edit.den$yi <- SumStats.Final.edit.den$Mean

SumStats.Final.edit.bio <- subset(SumStats.Final.edit, Resp == "Bio")
SumStats.Final.edit.bio$ES_ID <- seq.int(nrow(SumStats.Final.edit.bio))
SumStats.Final.edit.bio$vi <- SumStats.Final.edit.bio$SD^2
SumStats.Final.edit.bio$yi <- SumStats.Final.edit.bio$Mean

##
taxa.den.rma<-rma.mv(yi, vi, mods=~Taxa-1, random = ~1|MPA, test = "t", method="REML", data=SumStats.Final.edit.den) 
summary(taxa.den.rma)

# S3 method for rma.mv
cooks.taxa.den <- cooks.distance(taxa.den.rma)
SumStats.Final.den.cooks2 <- SumStats.Final.edit.den %>% cbind(cooks.taxa.den) %>% dplyr::rename(cooks.dist2 = cooks.taxa.den) %>%
  mutate(cooks.out2 = ifelse(cooks.dist2 > 4/62, "Y", "N"))

SumStats.Final.den.cooks <- subset(SumStats.Final.den.cooks2, cooks.out2 == "N")

taxa.den.rma<-rma.mv(yi, vi, mods=~Taxa -1, random = ~1|MPA, test = "t", method="REML", data=SumStats.Final.den.cooks) 
summary(taxa.den.rma)

# extract parameter estimates for table
taxa.mixed.est<-as.data.frame(cbind(coef(taxa.den.rma), taxa.den.rma$se, taxa.den.rma$zval, taxa.den.rma$pval, taxa.den.rma$ci.lb, taxa.den.rma$ci.ub)) %>% rownames_to_column()
colnames(taxa.mixed.est)<-c("mod.level", "estimate", "SE", "tvalue", "pvalue", "CI_low", "CI_up")
taxa.mixed.den.est <- taxa.mixed.est %>% mutate(mod.level = str_remove_all(mod.level, "Taxa"),
                                            Resp = "Den")#create column for x-axis spacing
taxa.mixed.den.est

taxa.bio.rma <- rma.mv(yi, vi, mods=~Taxa-1, random = ~1|MPA/Type.x, test = "t", method="REML", data=SumStats.Final.edit.bio) 
summary(taxa.bio.rma)

cooks.taxa.bio <- cooks.distance(taxa.bio.rma)
SumStats.Final.bio.cooks2 <- SumStats.Final.edit.bio %>% cbind(cooks.taxa.bio) %>% dplyr::rename(cooks.dist2 = cooks.taxa.bio) %>%
  mutate(cooks.out2 = ifelse(cooks.dist2 > 4/79, "Y", "N"))

SumStats.Final.bio.cooks <- subset(SumStats.Final.bio.cooks2, cooks.out2 == "N")

taxa.bio.rma<-rma.mv(yi, vi, mods=~Taxa-1, random = ~1|MPA/Type.x, test = "t", method="REML", data=SumStats.Final.bio.cooks) 
summary(taxa.bio.rma)

# extract parameter estimates for table
taxa.mixed.est<-as.data.frame(cbind(coef(taxa.bio.rma), taxa.bio.rma$se, taxa.bio.rma$zval, taxa.bio.rma$pval, taxa.bio.rma$ci.lb, taxa.bio.rma$ci.ub)) %>% rownames_to_column()
colnames(taxa.mixed.est)<-c("mod.level", "estimate", "SE", "tvalue", "pvalue", "CI_low", "CI_up")

taxa.mixed.bio.est <- taxa.mixed.est %>% mutate(mod.level = str_remove_all(mod.level, "Taxa"),
                                            Resp = "Bio")#create column for x-axis spacing
taxa.mixed.bio.est

taxa.mixed.est <- rbind(taxa.mixed.bio.est, taxa.mixed.den.est)
taxa.mixed.est$mod.level <- factor(taxa.mixed.est$mod.level, levels = c("S. purpuratus", "M. franciscanus", "M. pyrifera", "P. interruptus", "P. interruptus legal", "S. pulcher", "S. pulcher legal"), ordered = TRUE)

SumStats.cooks <- rbind(SumStats.Final.bio.cooks2,SumStats.Final.den.cooks2)
SumStats.cooks.outlier <- subset(SumStats.cooks, cooks.out2 == "Y")
SumStats.cooks.in <- subset(SumStats.cooks, cooks.out2 == "N")
SumStats.cooks.in$Resp <- as.character(SumStats.cooks.in$Resp)
SumStats.cooks.in$Resp[SumStats.cooks.in$Taxa=="M. pyrifera"] <- "Both"
SumStats.cooks.in$Resp <- factor(SumStats.cooks.in$Resp, levels = c("Den", "Bio"), ordered = TRUE)
SumStats.cooks.in.subLand <- subset(SumStats.cooks.in, Source == "Landsat")
SumStats.cooks.in.sub <- subset(SumStats.cooks.in, Source != "Landsat")
SumStats.cooks.in.sub$Resp <- factor(SumStats.cooks.in.sub$Resp, levels = c("Den", "Bio", "Both"), ordered = TRUE)
SumStats.cooks.in.subLand$Resp <- factor(SumStats.cooks.in.subLand$Resp, levels = c("Den", "Bio", "Both"), ordered = TRUE)
SumStats.cooks.in$Resp <- factor(SumStats.cooks.in$Resp, levels = c("Den", "Bio", "Both"), ordered = TRUE)

taxa.mixed.est$Resp[taxa.mixed.est$mod.level=="M. pyrifera"] <- "Bio"
SumStats.cooks.in$Resp[SumStats.cooks.in$Taxa=="M. pyrifera"] <- "Bio"

taxa.mixed.est$Resp <- factor(taxa.mixed.est$Resp, levels = c("Den", "Bio", "Both"), ordered = TRUE)
SumStats.cooks.in$Resp <- factor(SumStats.cooks.in$Resp, levels = c("Den", "Bio", "Both"), ordered = TRUE)
#####Figure 4#######################################################################################################
fp3 <- ggplot(data=taxa.mixed.est, aes(x=mod.level, y=estimate, color = Resp)) + 
  geom_point(position = position_dodge(0.5), size = 4,  shape = 15) + 
  geom_point(data = SumStats.cooks.in, aes(x=Taxa, y=Mean, color = Resp, alpha = 0.1), position = position_dodge(0.5), size = 1.5, shape = 15) +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_up), position = position_dodge(0.5), width = 0.1) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Taxa") + ylab("Effect Size (lnRR)") +
  theme_classic() +
  labs(color = "Response") +
  theme(legend.position = c(0.1, 0.85), ) +
  #theme(axis.text.x = element_text(face = "italic", angle = 45,  hjust = 1)) + 
  theme(axis.text.x = element_text(face = "italic")) + 
  scale_color_manual(values = c(pal[7], pal[16], pal[3]), labels = c("Density", "Biomass", "Combined")) +
  scale_alpha(guide = 'none') 
  #geom_text(data = label.df, label = "*", size = 10, color = "black")
fp3

ggsave(plot = fp3, file = "plots/Fig4_Mean_ESwide.png", 
       device = "png",  bg = "white",
       width = 15, height = 10, units = "cm", dpi = 300)



####Figure 5#########################################################################################################
All.RR.sub.trans$source <- factor(All.RR.sub.trans$source, levels = c("KFM","LTER","PISCO","Landsat"), ordered = TRUE)
All.RR.sub.trans$y <- factor(All.RR.sub.trans$y, levels = c("Strongylocentrotus purpuratus","Mesocentrotus franciscanus","Macrocystis pyrifera",
                                                            "Panulirus interruptus", "Semicossyphus pulcher"), ordered = TRUE)

All.RR.sub.trans.N <- subset(All.RR.sub.trans, CA_MPA_Name_Short == "Naples SMCA")
All.RR.sub.trans.N <- subset(All.RR.sub.trans.N, y != "Semicossyphus pulcher legal" & y != "Panulirus interruptus legal"
                              & y != "Panulirus interruptus sublegal" & y != "All urchins" & resp == "Den")

taxa <- c("Semicossyphus pulcher" = "#377eb8", "Macrocystis pyrifera" = "#4daf4a", "Panulirus interruptus" = "#ff7f00", "Strongylocentrotus purpuratus" = "#984ea3", "Mesocentrotus franciscanus" = "#e41a1c")
labels <- c("KFM" = 15, "LTER" = 16, "PISCO" = 17, "Landsat" = 18)

fp4 <- ggplot(data=All.RR.sub.trans.N, aes(x=year, y=lnDiff, color = y, shape = source)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) +
  geom_vline(xintercept=2012, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Year") + ylab("Effect Size (lnRR)") +
  theme_classic() +
  labs(color = "Taxa", shape = "Source") +
  theme(axis.text.x = element_text(face = "italic")) +
  scale_color_manual(values = taxa) +
  scale_shape_manual(values = labels) +
  theme(legend.position = c(0.9, 0.2), legend.text = element_text(size=6), legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'),     legend.spacing = unit(0.01, "cm")) 
fp4

All.RR.sub.trans.S <- subset(All.RR.sub.trans, CA_MPA_Name_Short == "Scorpion SMR")
All.RR.sub.trans.S <- subset(All.RR.sub.trans.S, y != "Semicossyphus pulcher legal" & y != "Panulirus interruptus legal"
                              & y != "Panulirus interruptus sublegal" & y != "All urchins" & resp == "Den")

fp5 <- ggplot(data=All.RR.sub.trans.S, aes(x=year, y=lnDiff, color = y, shape = source)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) +
  geom_vline(xintercept=2005, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Year") + ylab("Effect Size (lnRR)") +
  theme_classic() +
  labs(color = "Taxa", shape = "Source") +
  theme(axis.text.x = element_text(face = "italic")) +
  scale_color_manual(values = taxa) +
  scale_shape_manual(values = labels) +
  theme(legend.position = c(0.9, 0.2), legend.text = element_text(size=6), legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'),     legend.spacing = unit(0.01, "cm")) 
fp5

All.RR.sub.trans.AI <- subset(All.RR.sub.trans, CA_MPA_Name_Short == "Anacapa Island SMR 2003")
All.RR.sub.trans.AI <- subset(All.RR.sub.trans.AI, y != "Semicossyphus pulcher legal" & y != "Panulirus interruptus legal"
                             & y != "Panulirus interruptus sublegal" & y != "All urchins" & resp == "Den")

fp6 <- ggplot(data=All.RR.sub.trans.AI, aes(x=year, y=lnDiff, color = y, shape = source)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) +
  geom_vline(xintercept=2005, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Year") + ylab("Effect Size (lnRR)") +
  theme_classic() +
  labs(color = "Taxa", shape = "Source") +
  theme(axis.text.x = element_text(face = "italic")) +
  scale_color_manual(values = taxa) +
  scale_shape_manual(values = labels) +
  theme(legend.position = c(0.9, 0.2), legend.text = element_text(size=6), legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'),     legend.spacing = unit(0.01, "cm")) 
fp6

All.taxa <- ggarrange(fp4, fp5, fp6, labels = c("(a)","(b)","(c)"),  
                          hjust = 0.05, ncol = 1, align = "hv",
                          common.legend = TRUE, legend = "right")
All.taxa
ggsave(plot = All.taxa, file = "plots/Fig5_All_taxa.png", 
       device = "png",  bg = "white",
       width = 20, height = 30, units = "cm", dpi = 300)



























