options(digits=17)
#args<- c("example/d1.tad", "1","2")
#args<- c("example/duo.tad", "1","2")

invLogit<- function(p){
  ep <- exp(p)
  return( ep/(1+ep))
}

callDngLikeli<- function(params){
  som <- invLogit(params)

  dngOtherArgs <- paste("-t 20 --mu-somatic ", som, " --mu-library 0.0020300474329854423")
  # --nuc-freqs=0.25,0.25,0.25,0.25")
  args<- c("-p", pedFile, tadFile, dngOtherArgs, commandGetHidden)
  stdout<- system2(dngLoglike, args=args, stdout=T)
  ll <- as.numeric(stdout)
  return(ll)
}

callDngLikeliRelative <- function(params, baseLikeli){
  ll<-callDngLikeli(params)
  return(baseLikeli/ll)
}

## End functions

pid <- Sys.getpid()

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3){
  stop("3 parameters required: \"ad/tad, sample1, sample2\"")
}

dngLoglike <- "dng-loglike"
dngOtherArgs <- "-t 10 --mu-library 0.0020300474329854423"
commandGetTotal <- " | tail -n 1 |  awk '{print $1}'"
commandGetHidden <- " | tail -n 1 |  awk '{print $2}'"
# log_likelihood	log_hidden	log_observed

args[2]<- match.arg(args[2], 1:8)
args[3]<- match.arg(args[3], 1:8)

tadFile <- args[1]
name <- c(paste0("M", args[2]), paste0("M", args[3]))
pedFile <- paste0("M", args[2], args[3], ".ped")

prefix<- paste0("\"1\t1\t0\t0\t1\t(", name[1], ":")
branch<- c(1,1)
pedigree <- paste0(prefix, branch[1], ",", name[2], ":", branch[2],"):0\"")
system2("echo", c(pedigree, " > ", pedFile))
args<- c("-p", pedFile, tadFile, dngOtherArgs, commandGetHidden)
stdout<- system2(dngLoglike, args=args, stdout=T)
baseLikeli <- as.numeric(stdout)


params<- runif(1)
result<- optim(params, callDngLikeliRelative,  baseLikeli=baseLikeli,
              method="Brent", lower=-50, upper=50,
              control=list(fnscale=-1, maxit=500)
);

mu <- invLogit(result$par)
names(mu)<- c("somatic")
log_hidden <- baseLikeli/result$value

output<- list(mu_somatic=mu["somatic"], log_hidden=log_hidden,
 baseLikeli=baseLikeli, result=result)
print(output)

result_id<- paste0(paste0(name,collapse="_"), "_", pid)
outfile <- paste0("optimDngSomFixLib_", result_id, ".RData" )

assign(result_id, output)
save(list=c("result_id", result_id), file=outfile)
