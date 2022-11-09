# Script for running ABESS

# Configutation part #

# Input matrix folder after deduplication
input_dir = ""
# Output folder (make this directory)
output_dir = ""

######################


library(abess)

args=(commandArgs(TRUE))

if(length(args)!=2){
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

i <- args[1]
j <- args[2]


df <- read.table(paste0(input_dir, i, ".gt.domains.3.train.", j), sep=" ", header=F)
phenotypes <- read.table(paste0(input_dir, i, ".phen.domains.3.train.", j), header=F)
x <- as.matrix(df)
y <- phenotypes$V1


file_ex <- file(paste0(output_dir, i, ".domains.3.result.", j), open="w")

abess_fit <- abess(x, y, family = "binomial", normalize = 0, max.splicing.iter = 30, warm.start = FALSE, tune.type='gic', num.threads=24, seed=42)
writeLines(c('gic', toString(extract(abess_fit)$support.vars), toString(extract(abess_fit)$support.beta)), sep='\t', file_ex)
writeLines(c("\n"), sep="", file_ex)

close(file_ex)                
