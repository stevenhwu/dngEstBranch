#args<- c("1","2")

pid <- Sys.getpid()
args <- commandArgs(trailingOnly=TRUE)


args[1]<- match.arg(args[1], 1:8)
args[2]<- match.arg(args[2], 1:8)

name<- c(paste0("M",args[1]), paste0("M",args[2]))
dngLoglike <- "dng-loglike"
pedFile <- paste0("M",args[1],args[2],".ped")
tadFile <- "bwa_2.ad"

otherArgs <- "-t 20"#paste0("--mu-somatic ", mu)
getFirst <- " | tail -n 1 |  awk '{print $1}'"


# name<- c("A","B")
# dngLoglike <- "./b_ninja/src/dng-loglike"
# pedFile <- "testDataSW/duo.ped"
# tadFile <- "testDataSW/duo.tad"


prefix<- paste0("\"1\t1\t0\t0\t1\t(", name[1], ":")

callDngLikeli<- function(branch){
  if(any(branch > 1e4 |  branch < 1e-30 )){
    return(-Inf)
  }
  pedigree<-paste0(prefix, branch[1], ",", name[2], ":", branch[2],"):0\"")
  system2("echo", c(pedigree, " > ", pedFile))

  args<- c("-p", pedFile, tadFile, otherArgs, getFirst)
  stdout<- system2(dngLoglike, args=args, stdout=T)
  ll <- as.numeric(stdout)
  return(ll)
}

# optim(mu, callDngLikeli, control=list(fnscale=-1), method="Brent", lower=1e-30, upper=1)
# mu<- runif(1)
# callDngLikeli(100L)
options(digits=16)
branch<- runif(2,0,1)
# callDngLikeli(branch)
result<- optim(branch, callDngLikeli, control=list(fnscale=-1
   ,reltol=1e-12 ,maxit=2)
);
result$start<- branch

result_id<- paste0(paste0(name,collapse="_"), "_", pid)
outfile <- paste0("optimDngLoglikli_", result_id, ".RData" )

assign(result_id, result)
save(list=c("result_id", result_id), file=outfile)
