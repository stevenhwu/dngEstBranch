options(digits=17)
#args<- c("example/d1.tad", 1, "argFile")
#args<- c("example/duo.tad", 1, "reedParams")

invLogit<- function(p){
  ep <- exp(p)
  return( ep/(1+ep))
}

callDngLikeli<- function(params, argFile){
  
  lib<- invLogit(params)

  dngOtherArgs <- paste("-t 10 --mu-library", lib,  " --arg-file ", argFile)
  # --mu-library 0.0020300474329854423")
  # --nuc-freqs=0.25,0.25,0.25,0.25")
  args<- c("-p", pedFile, tadFile, dngOtherArgs, commandGetHidden)
  stdout<- system2(dngLoglike, args=args, stdout=T)
  ll <- as.numeric(stdout)
  return(ll)
}

callDngLikeliRelative <- function(params, baseLikeli, argFile){
  ll<-callDngLikeli(params, argFile)
  return(baseLikeli/ll)
}

## End functions

pid <- Sys.getpid()

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3){
  stop("3 parameters required: \"ad/tad, sample1, argFile\"")
}
argFile <- args[3]

dngLoglike <- "dng-loglike"
dngOtherArgs <- paste("-t 10 --arg-file ", argFile)
commandGetTotal <- " | tail -n 1 |  awk '{print $1}'"
commandGetHidden <- " | tail -n 1 |  awk '{print $2}'"
# log_likelihood	log_hidden	log_observed

args[2]<- match.arg(args[2], 0:8)
# args[3]<- match.arg(args[3], 1:8)

tadFile <- args[1]
name <- paste0("M", args[2])
pedFile <- paste0(name, ".ped")

pedigree<- paste0("\"1\t1\t0\t0\t1\t", name, ":0\"")
# branch<- c(1,1)
# pedigree <- paste0(prefix, branch[1], ",", name[2], ":", branch[2],"):0\"")
system2("echo", c(pedigree, " > ", pedFile))
args<- c("-p", pedFile, tadFile, dngOtherArgs, commandGetHidden)
stdout<- system2(dngLoglike, args=args, stdout=T)
baseLikeli <- as.numeric(stdout)


params<- runif(1)
result<- optim(params, callDngLikeliRelative,
      baseLikeli=baseLikeli, argFile=argFile,
      method="Brent", lower=-50, upper=50,
      control=list(fnscale=-1, maxit=500),
);

mu <- invLogit(result$par)
names(mu)<- "library"
log_hidden <- baseLikeli/result$value

output<- list(mu_library=mu["library"], log_hidden=log_hidden, baseLikeli=baseLikeli, result=result)
print(output)

result_id<- paste0(paste0(name,collapse="_"), "_", pid)
outfile <- paste0("optimDngCustomLibOnly_", result_id, ".RData" )

assign(result_id, output)
save(list=c("result_id", result_id), file=outfile)
