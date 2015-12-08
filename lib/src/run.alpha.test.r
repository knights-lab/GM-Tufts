# tests for differences in alpha diversity between groups at Timepoint 1 and 2 separately
run.alpha.test <- function(mapfile, alphafile, outputfile)
{
	map0 <- read.table(mapfile,sep='\t',head=T,row=1,comment='')
	alpha0 <- read.table(alphafile,sep='\t',head=T,row=1,comment='')

	# Paired t-tests: is diversity higher in one group after treatment, versus other group?
	# Stratify by Group, then do paired t.test on Pre_Post to see if significant difference
# 	map <- map0
# 	alpha <- alpha0
# 	samples_A <- rownames(map)[map$GROUP=="A"]
# 	samples_B <- rownames(map)[map$GROUP=="B"]
# 	
# 	paired.A <- t.test(alpha[samples_A, "PD_whole_tree"] ~ map[samples_A,"PRE_POST"], paired=T)
# 	paired.B <- t.test(alpha[samples_B, "PD_whole_tree"] ~ map[samples_B,"PRE_POST"], paired=T)
# 
# 	print(paired.A)
# 	print(paired.B)

	normality.pval <- shapiro.test(alpha0$PD_whole_tree)$p.value

	cat("Alpha Diversity between groups at Timepoint 1 and 2\n", file=outputfile)
	
	pdf(file="alpha-diversity-by-timepoint-beeswarm.pdf")
	for(timepoint in 1:2)
	{
		map <- map0[which(map0$PRE_POST==as.factor(timepoint)),]

		full <- merge(alpha0, map, by=0)	
		alpha <- alpha0[full$Row.names,]
		map <- map[full$Row.names,]

		if(normality.pval <= 0.05) {
			result <- wilcox.test(alpha$PD_whole_tree ~ map$GROUP)
		} else {
			result <- t.test(alpha$PD_whole_tree ~ map$GROUP)
		}
		cat("Timepoint ", timepoint, "\n", file=outputfile, append=T)
		cat("p-value", result$p.value, "\nGroup A and Group B Means", result$estimate, "\n", file=outputfile, append=T)

		# plot beeswarms of the alpha diversity to see if there are outliers
		beeswarm(alpha$PD_whole_tree ~ map$GROUP, main=paste("Alpha Diversity by Group, Timepoint", timepoint),
				ylab="PD Whole Tree")
		bxplot(alpha$PD_whole_tree ~ map$GROUP, add=T)

	}
	dev.off()
}