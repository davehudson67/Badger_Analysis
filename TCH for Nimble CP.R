library(tidyverse)
library(data.table)
rm(list=ls())

CH<- read.table("all.diag.results.txt", header = TRUE, row.names = NULL)
CH$date <- as.character(CH$date)
CH$date<- as.Date(CH$date, format="%d/%m/%Y")

#Merge data frames
colnames(CH)[which(names(CH) == "tattoo")] <- "ID"
CH<- as.data.table(CH)
setindexv(CH, c('ID'))
CH<- CH[,c(1:6,13:19,31)]
summary(CH)
#CH<-CH[CH$captureyear!=1979]

#rm(gene.data.info)
#rm(y.gen.ind)

#Keep only badgers with known age/sex
CH<-CH[complete.cases(CH[,8])]
CH<-CH[complete.cases(CH[,10])]


#Create a trap month variable
CH$trap_month <- as.numeric(format(CH$date, "%m"))
CH$trap_month<- as.factor(CH$trap_month)
summary(CH$trap_month)

#Put in date order
CH<- CH[order(CH$date),]

#Change NAs
CH$statpak[is.na(CH$statpak)]<-3
CH$brock[is.na(CH$brock)]<-3
CH$IFNgamma[is.na(CH$IFNgamma)]<-3
CH$cult_SUM[is.na(CH$cult_SUM)]<--1
#Remove occassions when no tests
CH$no.tests<- if_else(CH$statpak==3 & CH$brock==3 & CH$cult_SUM==-1 & CH$IFNgamma==3 & CH$pm!="Yes",1,0 )
CH<- CH[CH$no.tests==0]
CH$no.tests<-NULL
summary(CH)

#Define what is infected
CH$infected<- as.factor(if_else((CH$statpak==1 | CH$brock==1 | CH$cult_SUM>0 | CH$IFNgamma==1),1,0))
CH$infected.as.cub<- as.factor(if_else((CH$statpak==1 | CH$brock==1 | CH$cult_SUM>0 | CH$IFNgamma==1) & CH$age=="CUB", 1, 0))

summary(CH$infected)
summary(CH$infected.as.cub)

#Adjust the variable so now its ever infected or never infected
CH<- arrange(CH,ID, desc(infected))
for (i in 1:nrow(CH)){
  ifelse(CH$ID[i]==(CH$ID[i-1]),CH$infected[i]<-(CH$infected[i-1]),CH$infected[i])
}
CH$infected<- as.factor(CH$infected)

#Copy infected as Cub across all capture history of each individual
CH<- arrange(CH, ID, desc(infected.as.cub))
for (i in 1:nrow(CH)){
  ifelse(CH$ID[i]==(CH$ID[i-1]),CH$infected.as.cub[i]<-(CH$infected.as.cub[i-1]),CH$infected.as.cub[i])
}
CH$infected.as.cub<- as.factor(CH$infected.as.cub)
summary(CH)

#Add if known death
CH<- arrange(CH, ID, desc(date))
CH$knowndeath <- ifelse(CH$pm=="Yes",1,0)
for (i in 1:nrow(CH)){
  ifelse(CH$ID[i]==(CH$ID[i-1]),CH$knowndeath[i]<-(CH$knowndeath[i-1]),CH$infected.as.cub[i])
}

#Create status variable (infected/uninfected throughout life)
CH$status<- ifelse(CH$sex=="MALE" & CH$infected==1,"Pm",0)
CH$status<- ifelse(CH$sex=="FEMALE" & CH$infected==1,"Pf",CH$status)
CH$status<- ifelse(CH$sex=="MALE" & CH$infected==0,"NPm",CH$status)
CH$status<- ifelse(CH$sex=="FEMALE" & CH$infected==0,"NPf",CH$status)
CH$status<- as.factor(CH$status)

#Create status variable (infected.as.cub/uninfected throughout life)
CH$status.cub<- ifelse(CH$sex=="MALE" & CH$infected.as.cub==1,"CPm",0)
CH$status.cub<- ifelse(CH$sex=="FEMALE" & CH$infected.as.cub==1,"CPf",CH$status.cub)
CH$status.cub<-as.factor(CH$status.cub)
summary(CH)

#See how many occasions each badger is caught
#counts<-CH %>% 
#  group_by(ID) %>%
#  summarise(no_rows = length(ID)) %>%
#  tbl_df %>% print(n=nrow(.))

#Merge count info with CH
#CH<- data.table(CH)
#counts<- data.table(counts)
#setindexv(CH, c('ID'))
#setindexv(counts, c('ID'))
#CH<- merge(CH, counts, by="ID")

#See if any counts are <=2
#CH$toofew<- as.factor(ifelse(CH$no_rows<=1 & CH$occasion==0,1,0))
#summary(CH)
#CH<-CH[CH$toofew==0]
#rm(counts)
#CH <- droplevels(CH)
#Keep only NP and CP badgers
#CH<-as.data.table(CH)
#CH<- CH[CH$status=="NPf"|CH$status=="NPm"|CH$status.cub=="CPf"|CH$status.cub=="CPm"]
summary(CH)

CH <- as.data.table(CH)
CH.CP<-as.data.frame(CH[CH$status.cub=="CPm" | CH$status.cub=="CPf"])
CH.NP<-as.data.frame(CH[CH$status=="NPm" | CH$status=="NPf"])
CH.CP <- droplevels(CH.CP)
summary(CH.CP)

#Split trap month into 4 different 'seasons', year split into 4 quarters
CH.CP$trap_season<- ifelse(CH.CP$trap_month==1| CH.CP$trap_month==2|CH.CP$trap_month==3, "1", 0)
CH.CP$trap_season<- ifelse(CH.CP$trap_month==6|CH.CP$trap_month==4|CH.CP$trap_month==5, "2", CH.CP$trap_season)
CH.CP$trap_season<- ifelse(CH.CP$trap_month==9|CH.CP$trap_month==7|CH.CP$trap_month==8, "3", CH.CP$trap_season)
CH.CP$trap_season<- ifelse(CH.CP$trap_month==12|CH.CP$trap_month==10|CH.CP$trap_month==11, "4", CH.CP$trap_season)

#Create occasion variable by combining capture year and trap season
CH.CP <- unite(CH.CP, "occasion", captureyear, trap_season, sep = ".", remove = FALSE)
CH.CP$occasion <- as.factor(CH.CP$occasion)
CH.CP$occasion <- as.numeric(CH.CP$occasion)

CH.CP$death <- ifelse(CH.CP$pm=="Yes", CH.CP$occasion, NA)
CH.CP<- as.data.table(CH.CP)
CH.CP<- arrange(CH.CP, ID, desc(date))
for (i in 1:nrow(CH.CP)){
  ifelse(CH.CP$ID[i]==(CH.CP$ID[i-1]),CH.CP$death[i]<-CH.CP$death[i-1], CH.CP$death[i])
}

summary(CH.CP)
#Get rid of redundant levels and variables
CH.CP <- droplevels(CH.CP)

#See how many occasions each badger is caught
#counts<-CH.CP %>% 
#  group_by(ID) %>%
#  summarise(no_rows = length(ID)) %>%
#  tbl_df %>% print(n=nrow(.))

#CH.CP<- data.table(CH.CP)
#counts<- data.table(counts)
#setindexv(CH.CP, c('ID'))
#setindexv(counts, c('ID'))
#CH.CP<- merge(CH.CP, counts, by="ID")

#See if any counts are <=2 and one of them is pm
#CH.CP$toofew<- as.factor(ifelse(CH.CP$no_rows.y<=1,1,0))
#summary(CH.CP)
#CH.CP<-CH.CP[CH.CP$toofew==0]
#rm(counts)
#CH.CP <- droplevels(CH.CP)

#Create age in quarter years variable
#Create first line for badger 1
CH.CP<- arrange(CH.CP, ID, occasion)
CH.CP$age_quart.yr<- 1
#Continue for loop for all other badgers
for (i in 1:nrow(CH.CP)){
  ifelse(CH.CP$ID[i]!=(CH.CP$ID[i-1]), CH.CP$age_quart.yr[i]<-1,
         CH.CP$age_quart.yr[i]<-((CH.CP$age_quart.yr[i-1])+(CH.CP$occasion[i]-CH.CP$occasion[i-1])))
}

#Remove duplicate badgers (Same Badger caught in the same season)
CH.CP<-as.data.table(CH.CP)
CH.CP$duplicate<- ifelse(duplicated(CH.CP[,c(1,7)]), CH.CP$duplicate<-1, CH.CP$duplicate<- 0)
CH.CP<-CH.CP[CH.CP$duplicate==0]
CH.CP$duplicate<- NULL
summary(CH.CP)
CH.CP <- droplevels(CH.CP)

#See how many occasions each badger is caught
#counts<-CH.CP %>% 
#  group_by(ID) %>%
#  summarise(no_rows = length(ID)) %>%
#  tbl_df %>% print(n=nrow(.))

#CH.CP<- data.table(CH.CP)
#counts<- data.table(counts)
#setindexv(CH.CP, c('ID'))
#setindexv(counts, c('ID'))
#CH.CP<- merge(CH.CP, counts, by="ID")

#See if any counts are <=2 and one of them is pm
#CH.CP$toofew<- as.factor(ifelse(CH.CP$no_rows.y<=2 & CH.CP$knowndeath==1,1,0))
#summary(CH.CP)
#CH.CP<-CH.CP[CH.CP$toofew==0]
#rm(counts)
CH.CP <- droplevels(CH.CP)

#CH.CP2<- length(CH.CP[,if(.N>2) .SD, by=ID])

#Create an age array of remaining badgers
CH.CP$occasion<-as.factor(CH.CP$occasion)
CP.age.array<- array(NA, dim=c(length(levels(CH.CP$ID)), length(levels(CH.CP$occasion))))
rownames(CP.age.array)<- levels(CH.CP$ID)
colnames(CP.age.array)<- levels(CH.CP$occasion)

#Total number of observations
n.obs<- length(levels(CH.CP$ID))*length(levels(CH.CP$occasion))
CH.CP$trap_season<-as.numeric(CH.CP$trap_season)

#for loop to create the age array
for(i in 1:n.obs){
  CP.age.array[CH.CP$ID[i], CH.CP$occasion[i]]<- CH.CP$age_quart.yr[i]+CH.CP$trap_season[i]-1
}

#Create a vector with occasion of first capture
CP.f <- apply(CP.age.array, 1, function(x){which(!is.na(x))[1]})

#Fill age array with ages of badgers when not observed
CP.age.array <- t(apply(CP.age.array, 1, function(x){
  n <- which(!is.na(x))[1]
  x[n:length(x)] <- x[n]:(x[n] + length(x) - n)
  x
}))

CH.CP$birth_yr<-as.numeric(CH.CP$occasion)-CH.CP$age_quart.yr

#Change NA in age array to 0 (Bugs doesnt like NAs)
#CP.age.array[is.na(CP.age.array)] <- 0

#State definition
#1 = test negative on all tests
#2 = test positive on one of any tests
#3 = positive culture from one site
#4 = positive culture from multiple sites
#5 = dead
CH.CP$disease_state <- 1
CH.CP$disease_state <- ifelse(CH.CP$statpak==1 | CH.CP$brock==1,2,CH.CP$disease_state)
CH.CP$disease_state <- ifelse(CH.CP$cult_SUM==1,3,CH.CP$disease_state)
CH.CP$disease_state <- ifelse(CH.CP$cult_SUM>1,4,CH.CP$disease_state)
CH.CP$disease_state <- ifelse(CH.CP$pm=="Yes",5, CH.CP$disease_state)
CH.CP$disease_state <- as.factor(CH.CP$disease_state)
#CH.CPeck disease state variable
summary(CH.CP$disease_state)

#Create a state array of remaining badgers
CP.state.array<- array(0, dim=c(length(levels(CH.CP$ID)), length(levels(CH.CP$occasion))))
rownames(CP.state.array)<- levels(CH.CP$ID)
colnames(CP.state.array)<- levels(CH.CP$occasion)

#for loop to create the state array
for(i in 1:n.obs){
  CP.state.array[CH.CP$ID[i], CH.CP$occasion[i]]<- CH.CP$disease_state[i]
}

#Recode 'not seen' - 0's to 5
#CP.rstatearray<- CP.state.array
#CP.rstatearray[CP.rstatearray==0] <- 5


#Setup CJS capture history, 1 for seen, 0 for unseen.
CP.CJS.array<- CP.state.array
CP.CJS.array[CP.CJS.array>=5] <- 0
CP.CJS.array[CP.CJS.array>=1] <- 1

#Covariate data
CP.CVdata<- data.table(CH.CP)
CP.CVdata<- arrange(CP.CVdata, ID, -as.numeric(occasion))
CP.CVdata<- distinct(CP.CVdata, CP.CVdata$ID, .keep_all = TRUE)

CP.dead <- CP.CVdata$death


#CP.CJS.array[1,]
