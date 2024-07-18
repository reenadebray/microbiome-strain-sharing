# Strain and species sharing across varying criteria for transmission

## load FMT data
library(readxl)
fmt<-read_excel("FMT_instrain.xlsx",sheet="FMT_strain_sharing")
donor_profiles<-read_excel("FMT_donor_profiles.xlsx",sheet="FMT_donor_profiles")
recipient_profiles<-read_excel("FMT_recipient_profiles.xlsx",sheet="FMT_recipient_profiles")

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

## Comparisons with 0 shared species drop out of these aggregate functions; identify them and add them back (with a value of zero)
donors<-unique(donor_profiles$donor)
recipients<-unique(recipient_profiles$recipient)
for (d in donors){
  for (r in recipients){
    if (paste(pair[1],pair[2])%in%c("D33 R9","R9 D3","DT001 R15","R15 DT001","DT001 R16","R16 DT001","DT002 R19","R19 DT002","DT002 R22","R22 DT002","DT019 R17","R17 DT019","DT019 R18","R18 DT019","DONOR3185 pR3_160304", "pR3_160304 DONOR3185")){same_triad="Y"}
    else{same_triad="N"}
    if (nrow(shared[(shared$subject1==d & shared$subject2==r & shared$type2=="after") | (shared$subject2==d & shared$subject1==r & shared$type1=="after"),])<1){shared[nrow(shared)+1,]<-c(d,r,"donor","after",same_triad,0,0)}
    if (nrow(shared2[(shared2$subject1==d & shared2$subject2==r & shared2$type2=="after") | (shared2$subject2==d & shared2$subject1==r & shared2$type1=="after"),])<1){shared2[nrow(shared2)+1,]<-c(d,r,"donor","after",same_triad,0,0)}
    if (nrow(shared3[(shared3$subject1==d & shared3$subject2==r & shared3$type2=="after") | (shared3$subject2==d & shared3$subject1==r & shared3$type1=="after"),])<1){shared3[nrow(shared3)+1,]<-c(d,r,"donor","after",same_triad,0,0)}
    if (nrow(shared4[(shared4$subject1==d & shared4$subject2==r & shared4$type2=="after") | (shared4$subject2==d & shared4$subject1==r & shared4$type1=="after"),])<1){shared4[nrow(shared4)+1,]<-c(d,r,"donor","after",same_triad,0,0)}
  }
}



## load Amboseli data
amboseli<-read_excel("Amboseli_instrain.xlsx",sheet="Amboseli_strain_sharing")
shared_amboseli<-aggregate(amboseli$popANI,by=list(amboseli$name1,amboseli$date1,amboseli$name2,amboseli$date2,amboseli$category),FUN=function(x){length(x[x>=0.99999])/length(x)})
colnames(shared_amboseli)=c("name1","date1","name2","date2","category","shared_strains")


