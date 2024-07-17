# Strain and species sharing across varying criteria for transmission

## load packages
library(readr)
library(readxl)
library(tidyverse)
library(reshape2)

## load data
fmt<-read_excel("FMT_strain_sharing.xlsx",sheet="FMT_strain_sharing")
donor_profiles<-read_excel("FMT_donor_profiles.xlsx",sheet="FMT_donor_profiles")

## 1. All strain sharing events
### calculate shared strains (proportion of species with 99.999% average nucleotide identity per dyad)
shared<-aggregate(list(fmt$popANI),by=list(fmt$subject1,fmt$subject2,fmt$type1,fmt$type2,fmt$same_triad),function(x){length(x[x>=0.99999])/length(x)})
colnames(shared)=c("subject1","subject2","type1","type2","same_triad","shared_strains")
### calculate number of shared species
shared$shared_species<-aggregate(list(fmt$popANI),by=list(fmt$subject1,fmt$subject2,fmt$type1,fmt$type2,fmt$same_triad),function(x){length(x)})[,6]

## 2. Strain sharing events that were absent from recipient prior to FMT
### identify instances of strain sharing between pre- and post-FMT samples from the same individual
recipients<-unique(c(unlist(fmt[fmt$type1=="after","subject1"],fmt[fmt$type2=="after","subject2"])))
for (rec in recipients){
    strains_in_before<-unlist(fmt[fmt$subject1==rec & fmt$subject2==rec & fmt$popANI>0.99999,"genome"])
    fmt[((fmt$subject1==rec & fmt$type1=="after") | (fmt$subject2==rec & fmt$type2=="after")) & fmt$genome%in%strains_in_before,"in_before"]<-"Y"
}
### exclude these strains from all comparisons involving this individual
fmt2<-fmt[fmt$in_before!="Y" | is.na(fmt$in_before),]
### calculate shared strains
shared2<-aggregate(list(fmt2$popANI),by=list(fmt2$subject1,fmt2$subject2,fmt2$type1,fmt2$type2,fmt2$same_triad),function(x){length(x[x>=0.99999])/length(x)})
colnames(shared2)=c("subject1","subject2","type1","type2","same_triad","shared_strains")
### calculate shared species
shared2$shared_species<-aggregate(list(fmt2$popANI),by=list(fmt2$subject1,fmt2$subject2,fmt2$type1,fmt2$type2,fmt2$same_triad),function(x){length(x)})[,6]

## 3. Species that were unique to a single donor prior to FMT
### identify donor-genome combinations that occur only once
donor_profiles$genome_donor<-paste(donor_profiles$genome,donor_profiles$donor)
donor_profiles<-donor_profiles[!duplicated(donor_profiles$genome_donor),]
freq_table<-sort(table(donor_profiles$genome))
private_sp<-names(freq_table[freq_table==1])
### include only these species
fmt3<-fmt2[fmt2$genome%in%private_sp,]
### calculate shared strains
shared3<-aggregate(list(fmt3$popANI),by=list(fmt3$subject1,fmt3$subject2,fmt3$type1,fmt3$type2,fmt3$same_triad),function(x){length(x[x>=0.99999])/length(x)})
colnames(shared3)=c("subject1","subject2","type1","type2","same_triad","shared_strains")
### calculate shared species
shared3$shared_species<-aggregate(list(fmt3$popANI),by=list(fmt3$subject1,fmt3$subject2,fmt3$type1,fmt3$type2,fmt3$same_triad),function(x){length(x)})[,6]

## 4. Species that were unique to a single recipient after FMT
### identify recipient-genome combinations that occur only once
recipient_profiles$genome_recipient<-paste(recipient_profiles$genome,recipient_profiles$name_R)
recipient_profiles<-recipient_profiles[!duplicated(recipient_profiles$genome_recipient),]
freq_table<-sort(table(recipient_profiles$genome))
private_sp4<-names(freq_table[freq_table==1])
### include only these species
fmt4<-fmt3[fmt3$genome%in%private_sp4,]
### calculate shraed strains
shared4<-aggregate(list(fmt4$popANI),by=list(fmt4$subject1,fmt4$subject2,fmt4$type1,fmt4$type2,fmt4$same_triad),function(x){length(x[x>=0.99999])/length(x)})
colnames(shared4)=c("subject1","subject2","type1","type2","same_triad","shared_strains")
### calculate shared species
shared4$shared_species<-aggregate(list(fmt4$popANI),by=list(fmt4$subject1,fmt4$subject2,fmt4$type1,fmt4$type2,fmt4$same_triad),function(x){length(x)})[,6]
