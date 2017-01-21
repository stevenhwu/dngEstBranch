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
