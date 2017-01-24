files<- system2("find"," . -name *RData", T)
for(i in 1:length(files)){
  load(files[i])
}


allDist<- matrix(0,nrow=8,ncol=8)
M1Only <- ls()[grep("M1_M",ls())]
for(i in seq(1,length(M1Only), by=2)){
  print(paste(M1Only[i], M1Only[i+1]))
  if( get(M1Only[i])$mu_somatic == get(M1Only[i+1])$mu_somatic ){
    print(get(M1Only[i])$mu_somatic)
    index2<- as.numeric(sub("M1_M(.)_.*", "\\1", M1Only[i]))
    allDist[1,index2]<- get(M1Only[i])$mu_somatic
  }
}


MOthers <- ls()[grep("M[2-8]_M",ls())]
for(i in 1:length(MOthers)){
  print(get(MOthers[i])$mu_somatic)
  index<- strsplit(sub("M(.)_M(.)_.*", "\\1 \\2", MOthers[i]), " ")[[1]]
  index<- as.numeric(index)
  allDist[index[1], index[2]]<- get(MOthers[i])$mu_somatic
}

allDist<- allDist + (t(allDist))

as.dist(allDist)


#######################33


files<- system2("find"," . -name optimDngSomLib*RData", T)
for(i in 1:length(files)){
  load(files[i])
}

muSom<- matrix(0,nrow=8,ncol=8)
muLib<- matrix(0,nrow=8,ncol=8)
MAll <- ls()[grep("M[1-8]_M",ls())]
for(i in 1:length(MAll)){
  print(get(MAll[i])$mu_somatic)
  index<- strsplit(sub("M(.)_M(.)_.*", "\\1 \\2", MAll[i]), " ")[[1]]
  index<- as.numeric(index)
  muSom[index[1], index[2]]<- get(MAll[i])$mu_somatic
  muLib[index[1], index[2]]<- get(MAll[i])$mu_library
}

muSom<- muSom + (t(muSom))
muLib<- muLib + (t(muLib))

as.dist(allDist)



#########################################



files<- system2("find"," . -name optimDngLib*RData", T)
for(i in 1:length(files)){
  load(files[i])
}


muLibOnly<- vector(length=8)
MAll <- ls()[grep("M[1-8]_M",ls())]
for(i in 1:8){
  MEach <- ls()[grep( paste0("M",i,"_"),ls()) ]
  lib_each<- vector(length=length(MEach))
  for(j in 1:length(MEach)){
    lib_each[j]<- get(MEach[j])$mu_library
  }
  if( length(unique(lib_each)) == 1 ){
    muLibOnly[i] <- unique(lib_each)
  }
}


> muLibOnly
[1] 0.0016986623347644216 0.0019879176297104783 0.0019754392459755696
￼
[4] 0.0018543780177175853 0.0016242873572343715 0.0022419289930823826
[7] 0.0025842928070098792 0.0022734730783888516
> mean(muLibOnly)
[1] 0.0020300474329854423



#####################################


files<- system2("find"," . -name optimDngSomFixLib*RData", T)
for(i in 1:length(files)){
  load(files[i])
}

allSomFixLib<- matrix(0,nrow=8,ncol=8)
MAll <- ls()[grep("M[1-8]_M[1-8]",ls())]
for(i in 1:length(MAll)){
  
  index<- strsplit(sub("M(.)_M(.)_.*", "\\1 \\2", MAll[i]), " ")[[1]]
  index<- as.numeric(index)
  allSomFixLib[index[1], index[2]]<- get(MAll[i])$mu_somatic
  # muLib[index[1], index[2]]<- get(MAll[i])$mu_library
  
  
}


allSomFixLib<- allSomFixLib + (t(allSomFixLib))
