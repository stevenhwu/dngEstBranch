#args<- c("example/d1.tad", "1","2")
#args<- c("example/duo.tad", "1","2")

pid <- Sys.getpid()

args <- commandArgs(trailingOnly=TRUE)
print(args)
if(length(args)!=3){
  stop("3 parameters required: \"ad/tad, sample1, sample2\"")
}

dngLoglike <- "dng-loglike"
dngOtherArgs <- "-t 20"#paste0("--mu-somatic ", mu)
commandGetTotal <- " | tail -n 1 |  awk '{print $1}'"
commandGetHidden <- " | tail -n 1 |  awk '{print $2}'"
# log_likelihood	log_hidden	log_observed


args[2]<- match.arg(args[2], 1:8)
args[3]<- match.arg(args[3], 1:8)

tadFile <- args[1]
name <- c(paste0("M", args[2]), paste0("M", args[3]))
pedFile <- paste0("M", args[2], args[3], "_", pid, ".ped")

# name<- c("A","B")
# dngLoglike <- "./b_ninja/src/dng-loglike"
# pedFile <- "testDataSW/duo.ped"
# tadFile <- "testDataSW/duo.tad"


prefix<- paste0("\"1\t1\t0\t0\t1\t(", name[1], ":")

invLogit<- function(p){
  ep <- exp(p)
  return( ep/(1+ep))
}

callDngLikeli<- function(params){
  mu <- invLogit(params)

  # if(any(branch > 1e400 |  branch < 1e-30 | is.na(branch))){
  #   print(transformBranch)
  #   return(-Inf)
  # }
  # branch<- c(1,1)
  # pedigree<-paste0(prefix, branch[1], ",", name[2], ":", branch[2],"):0\"")
  # system2("echo", c(pedigree, " > ", pedFile))
  dngOtherArgs <- paste0("-t 20 --mu-somatic ", mu)
  args<- c("-p", pedFile, tadFile, dngOtherArgs, commandGetHidden)
  stdout<- system2(dngLoglike, args=args, stdout=T)
  ll <- as.numeric(stdout)
  return(ll)
}

callDngLikeliRelative <- function(params, baseLikeli){
  ll<-callDngLikeli(params)
  # baseLikeli<- params[3]
  print(ll)
  return(baseLikeli/ll)
}

# optim(mu, callDngLikeli, control=list(fnscale=-1), method="Brent", lower=1e-30, upper=1)
# mu<- runif(1)
# callDngLikeli(100L)
# callDngLikeli(branch)
options(digits=17)

# branch<- runif(2,0,1)
# #branch<- runif(2,-10,10)
# branch<- runif(2,1,10)
#
# result<- optim(branch, callDngLikeli, control=list(fnscale=-1
#    ,reltol=1e-12 , maxit=100)
# );
# result
# R1 - result$value

###############################

branch<- runif(2,0,2)
branch<- runif(2,-10,10)
# branch<- runif(2,exp(1),exp(2))
branch<- exp(runif(1))
#branch<- log(branch)


result<- optim(branch, callDngLikeliRelative,
  baseLikeli=baseLikeli,
  control=list(fnscale=-1
  # ,reltol=1e-12
  ,maxit=100
));
result
log(result$par)


R1 - log10(result$par)


R1 - result$value

# result$start<- branch
#
# result_id<- paste0(paste0(name,collapse="_"), "_", pid)
# outfile <- paste0("optimDngLoglikli_", result_id, ".RData" )
# print(result)
# assign(result_id, result)
# save(list=c("result_id", result_id), file=outfile)

pedigree<- vector(length=4)
baseLikeliAll<- vector(length=4)
pedigree[1]<-paste0(prefix, 1e-10, ",", name[2], ":", 1e-10,"):0\"")
pedigree[2]<-paste0(prefix, 10, ",", name[2], ":", 1e-10,"):0\"")
pedigree[3]<-paste0(prefix, 1e-10, ",", name[2], ":", 10,"):0\"")
pedigree[4]<-paste0(prefix, 10, ",", name[2], ":", 10,"):0\"")

for(i in 1:4){
  system2("echo", c(pedigree[i], " > ", pedFile))
  args<- c("-p", pedFile, tadFile, dngOtherArgs, commandGetHidden)
  stdout<- system2(dngLoglike, args=args, stdout=T)
  ll <- as.numeric(stdout)
  baseLikeliAll[i]<- ll
}

baseLikeliAll
baseLikeli<- min(baseLikeliAll)

########################################
pedigree[5]<-paste0(prefix, 1, ",", name[2], ":", 1,"):0\"")
system2("echo", c(pedigree[5], " > ", pedFile))
args<- c("-p", pedFile, tadFile, dngOtherArgs, commandGetHidden)
stdout<- system2(dngLoglike, args=args, stdout=T)


baseLikeli <- as.numeric(stdout)



mu_t<- rnorm(1)
result<- optim(mu_t, callDngLikeliRelative,
  baseLikeli=baseLikeli,
  control=list(fnscale=-1
  # ,reltol=1e-12
  ,maxit=100
));
result
invLogit(result$par)

r2<- optimize(callDngLikeliRelative, c(-100,100), maximum=T, baseLikeli=baseLikeli)
invLogit(r2$maximum)
