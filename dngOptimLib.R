options(digits=17)
#args<- c("example/d1.tad", "1")


invLogit<- function(p){
  ep <- exp(p)
  return( ep/(1+ep))
}

callDngLikeli<- function(params){
  # som <- invLogit(params[1])
  lib<- invLogit(params[1])

  dngOtherArgs <- paste("-t 10 --mu-library", lib)
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
if(length(args)!=2){
  stop("2 parameters required: \"ad/tad, sample1\"")
}

dngLoglike <- "dng-loglike"
dngOtherArgs <- "-t 10"#paste0("--mu-somatic ", mu)
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


# result<- optimize(callDngLikeliRelative, c(-50,50), maximum=T,
#                 baseLikeli=baseLikeli)
#
# mu_somatic<- invLogit(result$maximum)
# log_hidden <- baseLikeli/result$objective

params<- runif(1)
result<- optim(params, callDngLikeliRelative,  baseLikeli=baseLikeli,
      method="Brent", lower=-50, upper=50,
      control=list(fnscale=-1, maxit=500),
);

mu <- invLogit(result$par)
names(mu)<- "library"
log_hidden <- baseLikeli/result$value

output<- list(mu_library=mu["library"], log_hidden=log_hidden, baseLikeli=baseLikeli, result=result)
print(output)

result_id<- paste0(paste0(name,collapse="_"), "_", pid)
outfile <- paste0("optimDngLib_", result_id, ".RData" )

assign(result_id, output)
save(list=c("result_id", result_id), file=outfile)
