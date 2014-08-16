args <- commandArgs(trailingOnly = TRUE)
print(args)

## set working direction
if(length(args) < 1){
	setwd('../work') 
}else{
	setwd(args[1])
}
## set feature file name
if (length(args) < 2){
	data_file = 'evaluate_raw.txt'
}else{
	data_file = args[2]
}
## set plot type
if (length(args) < 3){
	title = 'ROC Curves for Predicting Interface Residues'
#	title = 'Precision Curves for Predicting Interface Residues'
#	title = 'ROC Curves for Predicting Interacting Domains'
#	title = 'Precision Curves for Predicting Interacting Domains'
}else{
	title = args[3]
}

err_mode = 'vertical' ## such as none, vertical
err_type = 'none' ## such as none, stderror, stddev, boxplot

library('ROCR')

cat('Reading', data_file, '...\n')
## Ref: http://stackoverflow.com/questions/4106764/what-is-a-good-way-to-read-line-by-line-in-r
cin <- file(data_file, open='r')
name <- c()
data <- list()
n <- 0
while(length(readlines <- readLines(cin, n=2, warn=FALSE)) > 0) 
{
	lines <- strsplit(readlines, '\n')
	eles <- strsplit(lines[[1]], '\t')[[1]]
	cat(eles[1:2],'\t')
	the_type = eles[1]
	the_name = eles[2]
	if(length(grep(the_type, title))) 
	{
		real <- as.logical(eles[c(-1,-2)])
		cat(unique(real),'\t')
		eles <- strsplit(lines[[2]], '\t')[[1]]
		pred <- as.numeric(eles[c(-1,-2)])
		cat(length(pred), '\n')
		n <- n + 1
		name <- c(name, the_name)
		data[[n]] <- list(real=real, pred=pred)
	}else{
		cat('\n')
	}
}
cat('Finished :)\n')

cat('Start plotting ...\n')
log_file = paste(data_file, '.log', sep='')
write.table(args, file=log_file, sep='\t', append=FALSE, row.names=FALSE, col.names=FALSE)
pdf(paste(data_file, '.pdf', sep=''))

uname <- unique(name)
m <- length(uname)
colors <- rainbow(m)
linetype <- c(1:m)
aucs <- c()

for(i in 1:m)
{
	all_pred <- list()
	all_real <- list()
	for(j in 1:n)
	{
		if (uname[i] == name[j])
		{
			all_pred <- c(all_pred, list(data[[j]]$pred))
			all_real <- c(all_real, list(data[[j]]$real))
		}
	}
	pair <- prediction(all_pred, all_real)
	if(length(grep('Precision', title)))
		perf <- performance(pair, measure="prec", x.measure="rec")
	else
		perf <- performance(pair, measure="tpr", x.measure="fpr")
	auc <- as.numeric(performance(pair, "auc")@y.values)
	aucs <- c(aucs, mean(auc))
	write.table(c(title, uname[i], auc), file=log_file, sep='\t', append=TRUE, row.names=FALSE, col.names=FALSE)
	if(i == 1)
#		plot(perf, col=colors[i], lty=linetype[i], lwd=3, avg=err_mode, spread.estimate=err_type, main=title)
		plot(perf, col=colors[i], lty=linetype[i], lwd=3, avg=err_mode, spread.estimate=err_type, main=title, xlim=c(0,1), ylim=c(0,1))
	else
		plot(perf, col=colors[i], lty=linetype[i], lwd=3, avg=err_mode, spread.estimate=err_type, add=TRUE)
}

ord <- sort(aucs, decreasing=TRUE, index.return=TRUE)$ix

if(length(grep('ROC', title)))
	legend(0.7,0.6, uname[ord], col=colors[ord], lty=linetype[ord], lwd=3)
if(length(grep('Precision', title))) {
	if(length(grep('Residue', title)))
		legend(0.7,1.0, uname[ord], col=colors[ord], lty=linetype[ord], lwd=3)
	else
		legend(0.7,1.0, uname[ord], col=colors[ord], lty=linetype[ord], lwd=3)
}

dev.off()
cat('End successfully!\n')

