# WRITTEN BY: ELHUM A SHAMSHIRI (Elhum.Shamshiri@unige.ch)
# CODE SHOULD NOT BE DISTRIBUTED AND ANALYSIS SHOULD NOT BE CONDUCTED WITHOUT PRIOR CONSENT FROM THE AUTHOR
# If you would like to use this software for publication please contact

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

# Elhum A Shamshiri - please cite Shamshiri et al., 2017 (Hum. Brain Mapp.)

library(R.matlab)
library(signal)  
library(FIACH)
library(AnalyzeFMRI)

########################################################
## This first part makes the regressors that you need ##
########################################################
input<-"/Users/ellieshamshiri/Desktop/Regressing_IEDs/eegFmri/patients/Sub_01/session_1.cond.mat" # this is your conditions file (SPM compatible)
RT<-1/16 # this is your RT
time<-296 # this the number of timepoints you have/volumes
tr<-2.16 # this is your TR (in seconds)

spike.reg<-function(input=input_cond,RT,time,tr){
  spm_hrf<-function(RT,p=c(6, 16, 1, 1, 6, 0, 32)){
    fMRI_T<-16
    if(RT<1/16){
      FMRI_T<-RT
    }
    dt<- RT/fMRI_T
    u <-c(0:(p[7]/dt)) - p[6]/dt
    hrf <- dgamma(u,p[1]/p[3],dt/p[3]) - dgamma(u,p[2]/p[4],dt/p[4])/p[1]
    hrf <-hrf[(0:(p[7]/RT))*fMRI_T + 1];
    hrf <- hrf/sum(hrf)
    return(list(hrf=hrf,p=p))
  }
  basis.functions<-function(RT){
    can.hrf<- spm_hrf(RT)
    dp<-1
    p<-can.hrf$p
    p[6]<-p[6]+dp
    D<- (can.hrf$hrf - spm_hrf(RT,p)$hrf)/dp
    can.hrf$hrf<-cbind(can.hrf$hrf, D)
    p[6]<- p[6] - dp
    
    dp<-0.01;
    p[3]<-p[3] + dp
    D<- (can.hrf$hrf[,1] - spm_hrf(RT,p)$hrf)/dp
    can.hrf$hrf<-cbind(can.hrf$hrf, D)
    p[3]<- p[3] - dp
    return(can.hrf$hrf)
  }
  
  
  split.input<-unlist(strsplit(input,"/"))
  dir<-paste(unlist(split.input)[1:(length(split.input)-1)],collapse="/")
  file<-split.input[length(split.input)]
  file.base<-substr(file,1,nchar(file)-9)
  output.file<-paste(dir,"/spike_regs_",file.base,".txt",sep="")
  
  mat<-readMat(input)
  onsets<-mat$onsets
  durations<-mat$durations
  spike.names<-unlist(mat$names)
  
  spike.list<-list()
  for(i in 1:length(spike.names)){ #need this for different spike types
    spike.list[[i]]<-cbind(t(onsets[[i]]),durations[[i]])
  }
  
  basis<-basis.functions(RT)
  #f16<-list(lapply(spike.list[[1]],function(x){x*tr*1/RT}))
  f16<-list()
  for(i in 1:length(spike.names)){ 
  f16[[i]]<-lapply(spike.list[[i]],function(x){x*tr*1/RT}  )              #for now just for one spike type
  }
  
  stimulus<-rep(0,round(time*tr*1/RT))
  ends<-list()
  #ends<-list(list(f16[[1]]+t(f16[[2]])))
  for(i in 1:length(spike.names)){ 
  ends[[i]]<-t(f16[[i]][[1]])+f16[[i]][[2]]
  }
  
  stim.create<-function(stimulus,onsets,ends){
    for(i in 1:length(onsets)){
      stimulus[onsets[i]:ends[i]]<-1
    }
    return(stimulus)
  }
  
  stim.out<-list()
  for(i in 1:length(f16)){
    stim.out[[i]]<-stim.create(stimulus,f16[[i]][[1]],ends[[i]])
  }
  conv.len<-length(stimulus)+nrow(basis)-1
  
  micro.hrf<-list()
  for(i in 1:ncol(basis)){
    micro.hrf[[i]]<-matrix(nrow=conv.len,ncol=length(stim.out))
  }
  
  for(j in 1:ncol(basis)){
    for(i in 1:length(stim.out)){
      micro.hrf[[j]][,i]<-convolve(stim.out[[i]],rev(basis[,j]),type="o")      
    }
  }
  
  regs.micro<-do.call("cbind",micro.hrf)
  fmri.hrf<-apply(regs.micro,2,function(x){resample(x[1:(length(x)-nrow(basis)+1)],1,tr*1/RT)})
  plot.ts(fmri.hrf)
  if(nrow(fmri.hrf)!=time){stop("messed up the downsampling")}
  write.table(fmri.hrf,output.file,row.names=FALSE,col.names=FALSE)
}

# The output is named like so: "spike_regs_session_1.txt"

######################################################################
## This second part regressors out those regressors from the images ##
######################################################################
input_fmri <- list.files("/Users/ellieshamshiri/Desktop/Regressing_IEDs/eegFmri/patients/Sub_01/session_1/fmri/",pattern = "^filt",full.names = T)
files<-input_fmri # these are your preprocessed files
reg<- "/Users/ellieshamshiri/Desktop/Regressing_IEDs/eegFmri/patients/Sub_01/spike_regs_session_1.txt"
#basis <- #maybe you need to put this in?
freq<-128 # this is your high pass freq #800
tr<-2.16 # this is your tr

treat<-function(input,reg,freq,tr){
  #####################################
  ########## DATA READ ################
  #####################################
  data<-readNii(input, fourD = T)
  print("Data is Read")  
  #####################################
  ######### MATRIX CREATION  ##########
  #####################################                                                                               
  time<-dim(data)[4]                                                                            ## find how many time points
  mat<-fourDto2D(data,time)                                                                     ## convert to matrix
  med<-colMedian(mat)[[1]]
  print("Median image Created")
  #####################################
  ###### FILE MANAGEMENT ##############
  #####################################
  split.input<-strsplit(input,"/")                                                              ## split input into folders
  dir<-paste(paste(split.input[[1]][1:(length(split.input[[1]])-1)],collapse="/"),sep="")       ## recombine to form directory
  files<-lapply(split.input,function(x){x[length(x)]})                                          ## find how many folders in path
  file.base<-unlist(lapply(files,function(x){substr(x,1,nchar(x)-4)}))                          ## remove extension from file 
  output.file.name<-paste(dir,"/treated_",file.base,".nii",sep="")                              ## create output filename
  print("Output Files Created")                                                                 ## create folder
  #####################################
  ####### HEADER MANAGEMENT ###########
  #####################################
  header<-lapply(input,f.read.nifti.header)                                                     ## read in headers  
  for(i in 1:length(header)){                                                                   ## loop over all and ...  
    header[[i]]$file.name<-output.file.name[i]                                                  ## change filename
    #header[[i]]$magic<-"n+1"                                                                    ## cahnge magic(necessary if analyze input)
    #header[[i]]$vox.offset<-352                                                                 ## again necessary for analyze input
  }
  print("Headers Created")
  #####################################
  ####### CREATE DESIGN MATRIX ########
  #####################################
  print("Creating Design Matrix")
  x<-as.matrix(mat)
  N<-nrow(x)
  K<-round(2*(N*tr/(freq + 1)))
  n <- t(0:(N-1))
  C <- matrix(0,N,K)
  C[,1]=1/sqrt(N)
  
  for (i in 2:K){C[,i] = sqrt(2/N)*cos(pi*(2*n+1)*(i-1)/(2*N))}
  basis<-as.matrix(C[,2:ncol(C)])
  #basis <- read.table("/Volumes/Seagate Backup Plus Drive/EEG_fMRI/Patients/EmHi_20121109/MRI/session_1/rfmri/noise_basis6.txt")
  rgs<-read.table(reg)
  
  n.regs<-ncol(rgs)
  n.filter<-ncol(basis)
  design<-cbind(1,rgs,basis)
  print("Design Matrix Created")
  #####################################
  ####### TREATING TIME ###############
  #####################################
  print("Treatment Started")
  qr<-qr(design)
  qr.fit<-qr.fitted(qr,mat)
  resid<-mat-qr.fit
  treated.mat<-t(t(resid)+as.vector(med))
  treated.arr<-twoDto4D(treated.mat,dim=dim(data))
  print("Treatment Finished")
  #####################################
  ############ WRITING ################
  #####################################
  print("Writing")
  for(i in 1:length(input)){
    f.write.nifti(treated.arr[,,,i],                                                                   
                  output.file.name[i],                                        
                  size="float",                                                              
                  L=header[[i]],                                              
                  nii=TRUE)
  }
  ######################################
  ########## MEMORY CLEANUP ############
  ######################################
  print("Memory Cleanup") 
  rm(list=ls())
  garbage<-gc()
  print("End")
}

treat(files,reg,freq,tr)

# The output is appended from the original fMRI files and is called: "treated_fMRI.nii"
