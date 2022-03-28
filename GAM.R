
# parse options
library("optparse")
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input bias and rowsum", metavar="character"),
    make_option(c("-o", "--output"), type="character", default="residual.txt", 
              help="output residual file", metavar="character"),
    make_option(c("--log_folder"), type="character", default="log/", 
              help="log file, PDF file folder", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# check input and log folder
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("You need to specify a input file.", call.=FALSE)
}
dir.create(opt$log_folder, showWarnings = FALSE)
log_file = paste(opt$log_folder, '/GAM.log', sep = "")
pdf_file = paste(opt$log_folder, '/GAM.pdf', sep = "")

# record run info.
sink(log_file,append=T)
Sys.time()
paste("Input : ",opt$input, sep = "")
sink()

# GAM model
library(mgcv)
print("Loading input data ... ")
data <- read.table(opt$input,sep="	",header=T)
data <- data.frame(data)

print("Fitting GAM model ... ")
lm1 <- gam(rowsum ~ s(gcper) + s(mapscore) + s(cutSitesNumber),data=data, family=poisson(link=log))
sink(log_file, append=T)
summary(lm1)
sink()

print("Calculating residuals ...")
rs <- residuals.gam(lm1,type="working")
write.table(rs, file=opt$output, quote=F, row.names=F, col.names=F)
pdf(pdf_file)
par(mfrow=c(2,2))
plot(lm1)
dev.off()

print("Done!")
