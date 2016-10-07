#Vellore environmental sample analysis;
setwd("~/stat/CMC/data/exposure data/")
re_dr<-read.csv("Drain_Sample_Database_092214_with_ec_conc.csv",header=TRUE)
hh_dat<-read.csv("Ind_data_102715.csv",header=TRUE)
pd_dat<-read.csv("phase2 public domain samples.csv",header=TRUE)

library(dplyr)
dat1<-re_dr %>% 
  select(SampleID,HouseID,ev_lat,ev_long,neighbor,ec_conc,ec_dil1,ec_dil2,ec_dil3,ec_dil4,
         ec_ecnt1,ec_ecnt2,ec_ecnt3,ec_ecnt4) %>% 
  mutate(study="re",SampleType="Drain",ec_dil1=1000/(10^ec_dil1),ec_dil2=1000/(10^ec_dil2),
         ec_dil3=1000/(10^ec_dil3),ec_dil4=1000/(10^ec_dil4))
dat2<-hh_dat %>% 
  select(SampleID,HouseID,ev_lat,ev_long,neighbor,SampleType,ec_conc,ec_dil1,ec_dil2,ec_dil3,
         ec_ecnt1,ec_ecnt2,ec_ecnt3) %>% 
  mutate(study="hh",ec_dil4=NA,ec_ecnt4=NA) %>%
  mutate(neighbor=ifelse(neighbor=="Old Town",1,ifelse(neighbor=="Cinna Allapuram",2,NA)))
dat3<-pd_dat %>% 
  select(sample_name,GPS_latitude,GPS_longitude,neighborhood,sample_type_name,ec_conc,eq_vol1,eq_vol2,
         count_num1,count_num2) %>% 
  mutate(study="pd",HouseID=NA,ec_dil3=NA,ec_ecnt3=NA,ec_dil4=NA,ec_ecnt4=NA) %>% 
  rename(SampleID=sample_name,ev_lat=GPS_latitude,ev_long=GPS_longitude,neighbor=neighborhood,
         SampleType=sample_type_name,ec_dil1=eq_vol1,ec_dil2=eq_vol2,ec_ecnt1=count_num1,ec_ecnt2=count_num2)

dat_all<-rbind(dat1,dat2,dat3)
dat_all<-cbind(1:nrow(dat_all),dat_all)
names(dat_all)[1]<-"ID"
dat_all$SampleTypeID<-0
dat_all$SampleTypeID[which(dat_all$SampleType=="Drain")]<-1
dat_all$SampleTypeID[which(dat_all$SampleType=="Bathing Water")]<-2
dat_all$SampleTypeID[which(dat_all$SampleType=="Child Hand Rinse")]<-3
dat_all$SampleTypeID[which(dat_all$SampleType=="Particulate")]<-4
dat_all$SampleTypeID[which(dat_all$SampleType=="Piped Water")]<-5
dat_all$SampleTypeID[which(dat_all$SampleType=="Produce")]<-6
dat_all$SampleTypeID[which(dat_all$SampleType=="Swabs")]<-7
dat_all$SampleTypeID[which(dat_all$SampleType=="Toy Feeding Spoon Rinse")]<-8
dat_all$ec_lnconc<-log10(dat_all$ec_conc)
dat_all$ec_lnconc[which(dat_all$ec_conc==0)]<--3
dat_all$ec_dil3[92]<-1

dat_all1<-dat_all

#calculate distance matrix
mat.dis<-matrix(NA,nrow=dim(dat_all1)[1],ncol=dim(dat_all1)[1])
for (i in 1:dim(dat_all1)[1]){
  for (j in 1:dim(dat_all1)[1]){
    mat.dis[i,j]<-sqrt((dat_all1$ev_long[i]-dat_all1$ev_long[j])^2+(dat_all1$ev_lat[i]-dat_all1$ev_lat[j])^2)
  }
}

gcd.slc <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
  return(d) # Distance in km
}

mat.dis.km<-matrix(NA,nrow=dim(dat_all1)[1],ncol=dim(dat_all1)[1])
for (i in 1:dim(dat_all1)[1]){
  for (j in 1:dim(dat_all1)[1]){
    mat.dis.km[i,j]<-gcd.slc(dat_all1$ev_long[i],dat_all1$ev_lat[i],dat_all1$ev_long[j],dat_all1$ev_lat[j])
    #change the NAs for correlation to 0;
    if (isTRUE(all.equal(sin(dat_all1$ev_lat[i])*sin(dat_all1$ev_lat[j]) + cos(dat_all1$ev_lat[i])*cos(dat_all1$ev_lat[j]) * cos(dat_all1$ev_long[j]-dat_all1$ev_long[i]),1))){mat.dis.km[i,j]<-0}
  }
}

trans_distance<-function(x){
  y<-ifelse(x>100,0,1/(x+1))
  return(y)
}

mat.corr<-trans_distance(mat.dis.km)

samtype<-dat_all1$SampleTypeID

mat.corr.samtype<-matrix(NA,nrow(dat_all1),nrow(dat_all1))

for (i in 1:nrow(dat_all1)){
  for (j in 1:nrow(dat_all1)){
    mat.corr.samtype[i,j]<-ifelse(samtype[i]!=samtype[j],0.4,ifelse(i==j,1,0.7)) #arbitary correlation number 0.4 between different sample types and 0.7 between same sample type but different samples..
  }
}

mat.correlation<-mat.corr*mat.corr.samtype
##############################################################################
n.smpl<-nrow(dat_all1)
num.comb<-length(unique(dat_all1$SampleType))

ind.env.smp <- array(NA,dim=c(nrow(dat_all1),3));
ind.env.smp[,1] <- as.numeric(dat_all1$ID)
ind.env.smp[,2] <- as.numeric(dat_all1$neighbor)
ind.env.smp[,3] <- as.numeric(dat_all1$SampleTypeID)

ec.eqv <- cbind(dat_all1$ec_dil1,dat_all1$ec_dil2,dat_all1$ec_dil3,dat_all1$ec_dil4);

ec.stdil <- array(NA,dim=c(nrow(ec.eqv),ncol(ec.eqv)));
for (i in 1:nrow(ec.eqv)){
    ec.stdil[i,] <- ifelse(dat_all1$SampleType[i]=="Particulate",2,
                            ifelse(dat_all1$SampleType[i]=="Swabs",14,
                                   ifelse(dat_all1$SampleType[i] %in% c("Child Hand Rinse","Produce","Tool Feeding Spoon Rinse"),500,100)
                            )
    )
}
ec.eqvol <- ec.eqv * ec.stdil;

ec.cnt <- cbind(dat_all1$ec_ecnt1,dat_all1$ec_ecnt2,dat_all1$ec_ecnt3,dat_all1$ec_ecnt4)

# count the number of replicate observations for each sample
mkrepl <- function(rawcnt,rawdil) {
  numsam <- length(rawcnt[,1]);
  replvec <- rep(NA,numsam);
  for(n in 1:numsam) {
    replvec[n] <- length(rawcnt[n,!is.na(rawcnt[n,])]);
  }
  return(round(replvec));
}

# generate dataset with three levels:
# n = sample number
# k = number aliquots tested (total of numaliq max)
# m = c(counted number, equivalent sample volume, indicator for censoring)
# if init =TRUE: same as previous function,
# but now with censored counts set at censorlimit+1
mkcnt.cens <- function(rawcnt,rawvol,numaliq,init=FALSE) {
  numsam <- length(rawcnt[,1]);
  replvec <- rep(NA,numaliq);
  cntdata <- array(NA,dim=c(numsam,numaliq,3));
  for(n in 1:numsam) { # loop through all samples
    for(k in 1:numaliq) { # initial fill of cntdata array with raw data
      cntdata[n,k,1] <- rawcnt[[n,k]];
      cntdata[n,k,2] <- rawvol[[n,k]];
      cntdata[n,k,3] <- NA; # ignore censoring now
    } # sort aliquot vectors: missing last
    cntdata[n,,] <- rbind(cntdata[n,!is.na(cntdata[n,,1]),],
                          cntdata[n,is.na(cntdata[n,,1]),]);
    for(k in 1:numaliq) {
      if(!is.na(cntdata[n,k,1])){ # count not missing?
        if(cntdata[n,k,1] > censorlimit){ # tntc?
          if(init){ # are we preparing init data?
            cntdata[n,k,1] <- censorlimit+1;
            cntdata[n,k,2] <- NA;
            cntdata[n,k,3] <- NA;
          }else{ # not preparing init data
            cntdata[n,k,1] <- NA;
            cntdata[n,k,3] <- TRUE;
          }
        }else{ # not tntc
          if(init){ # are we preparing init data?
            cntdata[n,k,1] <- NA;
            cntdata[n,k,2] <- NA;
            cntdata[n,k,3] <- NA;
          }else{ # not preparing init date
            cntdata[n,k,3] <- FALSE;
          }
        }
      }else{ # count missing!
        cntdata[n,k,2] <- NA; # then also no volume
      }
    }
  }
  return(cntdata);
}

censorlimit <- 200;
ec.repl <- mkrepl(ec.cnt,ec.eqvol);
ec.data <- mkcnt.cens(ec.cnt,ec.eqvol,4,FALSE);
ec.init <- mkcnt.cens(ec.cnt,ec.eqvol,4,TRUE);
#ec.init[notused,,] <- NA;
#I<-diag(n.smpl)

cntdata <- list("n.smpl"=n.smpl, "num.comb"=num.comb,
                "censorlimit"=censorlimit,
                "ind.env.smp"=ind.env.smp,
                "ec.repl"=ec.repl, "ec.data"=ec.data,
                "mat.corr"=mat.correlation);

cntinit <- list("ec.data"=ec.init);


