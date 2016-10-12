#results = data.frame of pvals, adj.pvals, and coeffs (return value of test.otu.features.r)
#x=otu_delta[samples_A,]
#y=map[samples_A,"clinical_delta"]
add.significant.plots <- function(taxa, clinical.var, pre.main, x, y, coeffs, pvals)
{
	short.taxa <- shorten.taxonomy(taxa)

	if(length(taxa) > 0)
	{
		for(i in 1:length(taxa))
		{
			coeff <- coeffs[i]
			pval <- pvals[i]
			plot(x=x[,taxa[i]], y=y, main=paste(pre.main, short.taxa[i]), xlab=short.taxa[i], ylab=clinical.var)
			abline(lm(y~x[,taxa[i]]), col="red")
			mtext(paste("Rho = ", coeff, ", adjusted pvalue = ", pval, sep=""), side=3)
		}
	}
}

# Tests for change in clinical_var ~ change in taxa * Group
test.delta.interactions <- function(otu_delta, map_delta, outputfile="interactions.results.txt", clinical.vars)
{
	cat("Delta Taxon ~ Delta Clinical * Group\n", file=outputfile)

	for(i in 1:length(clinical.vars))
	{	
		clinical.var <- clinical.vars[i]
		delta.clinical.var <- paste("delta", clinical.var,sep="_")

		# 1. change in clinical_var ~ change in taxa * Group
		ret1 <- test.otu.features(otu_delta, map_delta, response=delta.clinical.var, treatment="taxa", 
								controls=controls, treatment.interaction.var="Group")

		cat(clinical.var, file=outputfile, append=T)							
		cat("delta_clinical ~ delta_taxa * Group\n", file=outputfile, append=T)
		cat(names(ret1)[1], "\n\t", file=outputfile, append=T)
		write.table(ret1[[1]], file=outputfile, append=T, sep="\t", quote=F)
		cat("delta_clinical ~ delta_taxa * Group\n", names(ret1)[2], "\n\t",  file=outputfile, append=T)
		write.table(ret1[[2]], file=outputfile, append=T, sep="\t", quote=F)

	}
}

#change in clinical_var ~ change in taxa, stratified by Group
delta.taxon.v.delta.covariate.within.group <- function(otu_delta, map_delta, outputfile="delta-taxon-v-delta-covariate-within-group.txt", 
											clinical.vars, controls, plot.pvals=F, control.baseline=F, otu.baseline)
{
	corr.outputfile <- "delta-taxon-v-delta-covariate-within-group-correlation.txt"
	cat("Delta Taxon ~ Delta Clinical\n", file=outputfile)
	cat("Spearman's Correlation: Delta Taxon ~ Delta Clinical\n", file=corr.outputfile)
		
	samples_A <- rownames(map_delta)[map_delta$GROUP=="A"]
	samples_B <- rownames(map_delta)[map_delta$GROUP=="B"]
	participants_A <- as.factor(map_delta[samples_A,"PARTICIPANTID"])
	participants_B <- as.factor(map_delta[samples_B,"PARTICIPANTID"])
	all.adj.pvals <- NULL

	pdf(file="delta-taxon-v-delta-covariate-within-group.pdf")	
	no.sig.taxa <- TRUE
	n.clinical.vars <- length(clinical.vars)
	for(i in 1:n.clinical.vars)
	{	
		clinical.var <- clinical.vars[i]
		delta.clinical.var <- paste("delta", clinical.var,sep="_")	
	
		# if we want to control for the initial value of the clinical var and initial value of the taxa
		# NOTE: don't adjust for number of clinical variables
print(clinical.vars[i])
		if(control.baseline==T) {
print("--------------------1GROUP A COMPARISONS")
			ret2A <- test.otu.features2(otu=otu_delta[samples_A,], map=map_delta[samples_A,], response=delta.clinical.var, treatment="taxa", 
									controls=c(controls, clinical.var), otu.baseline=otu.baseline[samples_A,])	
print("--------------------2GROUP B COMPARISONS")
			ret2B <- test.otu.features2(otu=otu_delta[samples_B,], map=map_delta[samples_B,], response=delta.clinical.var, treatment="taxa", 
									controls=c(controls, clinical.var), otu.baseline=otu.baseline[samples_B,])	
		
		} else {
print("--------------------3GROUP A COMPARISONS")
			ret2A <- test.otu.features(otu=otu_delta[samples_A,], map=map_delta[samples_A,], response=delta.clinical.var, treatment="taxa", 
									controls=controls)	
print("--------------------4GROUP B COMPARISONS")
			ret2B <- test.otu.features(otu=otu_delta[samples_B,], map=map_delta[samples_B,], response=delta.clinical.var, treatment="taxa", 
									controls=controls)	
		}
	
		cat(clinical.var, "\n", file=outputfile, append=T)
		cat("delta_clinical ~ delta_taxa\nGroup A:\n\t", file=outputfile, append=T)
		write.table(ret2A[[1]], file=outputfile, append=T, sep="\t", quote=F)

		cat(clinical.var, "\n", file=outputfile, append=T)
		cat("delta_clinical ~ delta_taxa\nGroup B:\n\t",  file=outputfile, append=T)
		write.table(ret2B[[1]], file=outputfile, append=T, sep="\t", quote=F)

		ret2A.corr <- test.otu.features(otu=otu_delta[samples_A,], map=map_delta[samples_A,], response=delta.clinical.var, treatment="taxa", 
								controls=controls, use.corr=T)	
		ret2B.corr <- test.otu.features(otu_delta[samples_B,], map_delta[samples_B,], response=delta.clinical.var, treatment="taxa", 
								controls=controls, use.corr=T)	

		cat(clinical.var, "\n", file=corr.outputfile, append=T)
		cat("delta_clinical ~ delta_taxa\nGroup A:\n\t", file=corr.outputfile, append=T)
		write.table(ret2A.corr, file=corr.outputfile, append=T, sep="\t", quote=F)

		cat(clinical.var, "\n", file=corr.outputfile, append=T)
		cat("delta_clinical ~ delta_taxa\nGroup B:\n\t",  file=corr.outputfile, append=T)
		write.table(ret2B.corr, file=corr.outputfile, append=T, sep="\t", quote=F)

		# let's only print taxa that are significant in both LM and Spearman's
		sig.taxa.A <- intersect(rownames(ret2A.corr)[!is.na(ret2A.corr[,"adj.pvals",drop=F]) & ret2A.corr[,"adj.pvals",drop=F] < .25], rownames(ret2A[[1]][,"adj.pvals",drop=F])[!is.na(ret2A[[1]][,"adj.pvals"]) & ret2A[[1]][,"adj.pvals"] < .25])
		sig.taxa.B <- intersect(rownames(ret2B.corr)[!is.na(ret2B.corr[,"adj.pvals",drop=F]) & ret2B.corr[,"adj.pvals",drop=F] < .25], rownames(ret2B[[1]][,"adj.pvals",drop=F])[!is.na(ret2B[[1]][,"adj.pvals"]) & ret2B[[1]][,"adj.pvals"] < .25])

		# spit out significant taxa for GM
		if(length(sig.taxa.A)+length(sig.taxa.B) > 0)
		{
			if(no.sig.taxa==TRUE)
			{
				no.sig.taxa <- FALSE
				cat("Significant Delta Taxon ~ Delta Clinical within Groups\n", file=paste("significant-", outputfile,sep=""))
			}
			cat(clinical.var, "\n", file=paste("significant-", outputfile,sep=""), append=T)
		}
		if(length(sig.taxa.A) > 0)
		{
			cat("Group A\n", file=paste("significant-", outputfile,sep=""), append=T)
			write.table(ret2A[[1]][sig.taxa.A,], file=paste("significant-", outputfile,sep=""), sep="\t", quote=F, append=T)

			# For any taxa that are significant, let's plot them (simple scatter plot)
			add.significant.plots(taxa=sig.taxa.A, clinical.var=clinical.var, 
					x=otu_delta[samples_A,], y=map_delta[samples_A,delta.clinical.var],pre.main=paste("Group A, Delta",clinical.var,"~ Delta"),
					coeffs=ret2A.corr[sig.taxa.A,"coeffs"], pvals = ret2A.corr[sig.taxa.A,"adj.pvals"])
		}
		if(length(sig.taxa.B) > 0)
		{
			cat("Group B\n", file=paste("significant-", outputfile,sep=""), append=T)
			write.table(ret2B[[1]][sig.taxa.B,], file=paste("significant-", outputfile,sep=""), sep="\t", quote=F, append=T)
		
			add.significant.plots(taxa=sig.taxa.B, clinical.var=clinical.var, 
				x=otu_delta[samples_B,], y=map_delta[samples_B,delta.clinical.var],pre.main=paste("Group B, Delta",clinical.var,"~ Delta"),
				coeffs=ret2B.corr[sig.taxa.B,"coeffs"], pvals = ret2B.corr[sig.taxa.B,"adj.pvals"])
		}	
		
		
		# lets plot -log pvalues by group per clinical variable here
		temp <- cbind(ret2A[[1]]$adj.pvals, ret2B[[1]]$adj.pvals)
		colnames(temp) <- c(paste(clinical.var, "A",sep="-"), paste(clinical.var, "B",sep="-"))
		all.adj.pvals <- cbind(all.adj.pvals, temp)
		
	}
	dev.off()
		
	if(plot.pvals==T)
	{
		# make box and whiskers plot of all adjusted pvalues
		pdf("BOXPLOT-pvalues-v-covariate.pdf")
		par(mar=c(10.1,4.1,3.1,2.1))
		if(length(clinical.vars) > 30)
		{
			boxplot(-log10(all.adj.pvals),las=2,main="-Log(q-value) vs. delta clinical variables by group", cex.axis=.25)
		} else {
			boxplot(-log10(all.adj.pvals),las=2,main="-Log(q-value) vs. delta clinical variables by group")
		}
		abline(h=-log10(0.25),lty=2)
			dev.off()
	}
}

delta.taxon.v.delta.covariate.all.groups <- function(otu_delta, map_delta, outputfile="delta-taxon-v-delta-covariate-all-groups.txt", 
											clinical.vars, controls, plot.pvals=F, control.baseline=F, otu.baseline)
{
	corr.outputfile <- "delta-taxon-v-delta-covariate-all-groups-correlation.txt"
	cat("Delta Taxon ~ Delta Clinical\n", file=outputfile)
	cat("Spearman's Correlation: Delta Taxon ~ Delta Clinical\n", file=corr.outputfile)
		
	all.adj.pvals <- NULL

	pdf(file="delta-taxon-v-delta-covariate-all-groups.pdf")	

	no.sig.taxa <- TRUE
	n.clinical.vars <- length(clinical.vars)

	for(i in 1:length(clinical.vars))
	{	
		clinical.var <- clinical.vars[i]
		delta.clinical.var <- paste("delta", clinical.var,sep="_")

print(clinical.vars[i])
print("--------------------ALL GROUP COMPARISONS")
		if(control.baseline==T){
			ret <- test.otu.features2(otu=otu_delta, map=map_delta, response=delta.clinical.var, treatment="taxa", 
								controls=c(controls, clinical.var), otu.baseline=otu.baseline)	
		} else {
			ret <- test.otu.features(otu=otu_delta, map=map_delta, response=delta.clinical.var, treatment="taxa", 
								controls=controls)	
		}
		
		cat(clinical.var, "\n", file=outputfile, append=T)
		cat("delta_clinical ~ delta_taxa:\n\t", file=outputfile, append=T)
		write.table(ret[[1]], file=outputfile, append=T, sep="\t", quote=F)

		ret.corr <- test.otu.features(otu=otu_delta, map=map_delta, response=delta.clinical.var, treatment="taxa", 
								controls=controls, use.corr=T)	

		cat(clinical.var, "\n", file=corr.outputfile, append=T)
		cat("delta_clinical ~ delta_taxa:\n\t", file=corr.outputfile, append=T)
		write.table(ret.corr, file=corr.outputfile, append=T, sep="\t", quote=F)

		# let's only print taxa that are significant in both LM and Spearman's
		sig.taxa <- intersect(rownames(ret.corr[,"adj.pvals",drop=F])[!is.na(ret.corr[,"adj.pvals"]) & ret.corr[,"adj.pvals"] < .25], rownames(ret[[1]][,"adj.pvals",drop=F])[!is.na(ret[[1]][,"adj.pvals"]) & ret[[1]][,"adj.pvals"] < .25])

		# spit out significant taxa for GM
		if(length(sig.taxa) > 0)
		{
			if(no.sig.taxa==TRUE)
			{
				no.sig.taxa <- FALSE
				cat("Significant Delta Taxon ~ Delta Clinical for all Groups\n", file=paste("significant-", outputfile,sep=""))
			}
			
			cat(clinical.var, "\n", file=paste("significant-", outputfile,sep=""), append=T)

			write.table(ret[[1]][sig.taxa,], file=paste("significant-", outputfile,sep=""), sep="\t", quote=F, append=T)

			# For any taxa that are significant, let's plot them (simple scatter plot)
			add.significant.plots(taxa=sig.taxa, clinical.var=clinical.var, 
					x=otu_delta, y=map_delta[,delta.clinical.var],pre.main=paste("All groups, Delta",clinical.var,"~ Delta"),
					coeffs=ret.corr[sig.taxa,"coeffs"], pvals = ret.corr[sig.taxa,"adj.pvals"])
		}
		
		# lets plot -log pvalues by group per clinical variable here
		all.adj.pvals <- cbind(all.adj.pvals, ret[[1]]$adj.pvals)
		
	}
	dev.off()
		
	if(plot.pvals==T)
	{
		colnames(all.adj.pvals) <- clinical.vars
		pdf("BOXPLOT-pvalues-v-covariate-all-groups.pdf")
		# make box and whiskers plot of all adjusted pvalues
		par(mar=c(10.1,4.1,3.1,2.1))
		if(length(clinical.vars) > 30)
		{
			boxplot(-log10(all.adj.pvals),las=2,main="-Log(q-value) vs. delta clinical variables", cex.axis=.25)
		} else {
			boxplot(-log10(all.adj.pvals),las=2,main="-Log(q-value) vs. delta clinical variables") 
		}
		
		abline(h=-log10(0.25),lty=2)
		dev.off()
	}
}

# taxa at time 1 ~ clinical var at time 1 or time 2
taxon.v.covariate <- function(otu, map, clinical.vars, timepoint=1, controls)
{
	outputfile=paste("taxon-v-covariate-timepoint-", timepoint, ".txt", sep="")
	corr.outputfile=paste("taxon-v-covariate-timepoint-", timepoint, "-correlation.txt", sep="")
	
	pdf(file=paste("taxon-v-covariate-timepoint-", timepoint, ".pdf", sep=""))	

	no.sig.taxa <- TRUE
	n.clinical.vars <- length(clinical.vars)
	
	for(i in 1:length(clinical.vars))
	{	
		clinical.var <- clinical.vars[i]
		
		ret4 <- test.otu.features(otu, map, treatment=clinical.var, controls=controls)	

		cat(clinical.var, "\n", file=outputfile, append=T)
		cat("taxa ~ clinical_var\n\t", file=outputfile, append=T)
		write.table(ret4[[1]], file=outputfile, append=T, sep="\t", quote=F)

		ret4.corr <- test.otu.features(otu, map, response="taxa", treatment=clinical.var, controls=controls, use.corr=T)	

		cat(clinical.var, "\n", file=corr.outputfile, append=T)
		cat("taxa ~ clinical_var\n\t", file=corr.outputfile, append=T)
		write.table(ret4.corr, file=corr.outputfile, append=T, sep="\t", quote=F)

		# let's only print taxa that are significant in both LM and Spearman's
		# make sure that we ignore any that are NA 
		corr.names <- rownames(ret4.corr[,"adj.pvals",drop=F])[!is.na(ret4.corr[,"adj.pvals"]) & ret4.corr[,"adj.pvals"] < .25]
		lm.names <- rownames(ret4[[1]][,"adj.pvals",drop=F])[!is.na(ret4[[1]][,"adj.pvals"]) & ret4[[1]][,"adj.pvals"] < .25]
		sig.taxa <- intersect(corr.names,lm.names)

		# spit out significant taxa for GM
		if(length(sig.taxa) > 0)
		{
			if(no.sig.taxa==TRUE)
			{
				no.sig.taxa <- FALSE
				cat("Significant Taxon ~ Clinical at Timepoint\n", file=paste("significant-", outputfile,sep=""))
			} 
			
			cat(clinical.var, "\n", file=paste("significant-", outputfile,sep=""), append=T)
			write.table(ret4[[1]][sig.taxa,], file=paste("significant-", outputfile,sep=""), sep="\t", quote=F, append=T)
		
			add.significant.plots(taxa=sig.taxa, clinical.var=clinical.var, 
			x=otu, y=map[,clinical.var], pre.main=paste("Timepoint ", timepoint, ", ", clinical.var,"~"),
			coeffs=ret4.corr[sig.taxa,"coeffs"], pvals = ret4.corr[sig.taxa,"adj.pvals"])
		}
	}
	dev.off()
}

# Tests for alpha.div at time1, delta alpha div ~ delta clinical covariate AND again stratified by group
test.alpha.div <- function(samples1, samples2, map_delta, alphafile, outputfile="alpha.div.results.txt", clinical.vars, controls)
{
	alpha <- read.table(alphafile,sep='\t',head=T,row=1,comment='')
	
	# all groups combined
	pd_1 <- alpha[samples1,"PD_whole_tree"]
	pd_2 <- alpha[samples2,"PD_whole_tree"]

	alpha_div <- cbind(pd_1, delta_pd=(pd_2-pd_1)) # treat this like a 2 column OTU table
	

	cat("Alpha Diversity at Timepoint 1 and Delta Alpha Diversity associated with Delta Clinical Covariate\n", file=outputfile, append=F)

	n.clinical.vars <- length(clinical.vars)
	for(i in 1:length(clinical.vars))
	{	
		clinical.var <- clinical.vars[i]
		delta.clinical.var <- paste("delta", clinical.var,sep="_")

		ret3 <- test.otu.features(alpha_div, map_delta, treatment=delta.clinical.var, controls=controls)	
		cat(clinical.var, "\n", file=outputfile, append=T)
		rownames(ret3[[1]]) <- c("alpha_div1 ~ delta_clinical_var", "delta_alpha_div ~ delta_clinical_var")
		write.table(ret3[[1]], file=outputfile, append=T, sep="\t", quote=F)
	}	
	
}


run.clinical.tests <- function(mapfile, otufile, clinical.vars, alphafile, run.test=c(2,3,4), keep.pathways.file='data/pathways.to.keep.txt', controls=c("SEX","RACE","BMI"))
{
	library('beeswarm')

	# read in data files	
	minOTUInSamples = .001
	ret <- load.data(otufile=otufile, mapfile=mapfile, minOTUInSamples=minOTUInSamples, normalize=TRUE)
	otu <- ret$otu
	map <- ret$map
	kegg <- ret$kegg
	
	#remove Unassigned from entire OTU altogether
	otu <- otu[, colnames(otu) != "Unassigned;Other;Other;Other;Other;Other"]
	
	kegg_pathways <- NULL
	if(length(kegg)>0){
		# filter the kegg pathways and return the next level up for displaying
		kegg <- filter.pathways(kegg, keep.pathways.file)
		next.kegg <- get.next.kegg(kegg)
		names(next.kegg) <- names(kegg)
		kegg_pathways<-next.kegg

		# select only the kegg pathways of interest
		otu=otu[,names(kegg_pathways)]
	}
	
	# add AR to the mapping file (sum of ARA ARB ARC)
	map$AR <- rowSums(cbind(map$ARA,map$ARB,map$ARC), na.rm=T)
	
	participant_ids <- unique(map$PARTICIPANTID)
	ix_t1 <- which(map$PRE_POST=="1" & (map$PARTICIPANTID %in% participant_ids))
	ix_t2 <- which(map$PRE_POST=="2" & (map$PARTICIPANTID %in% participant_ids))
	otu_delta <- otu[ix_t2,] - otu[ix_t1,]
	rownames(otu_delta) <- rownames(map)[ix_t2] #participant_ids
	
#	write.table(otu_delta, file="otu_delta.txt", quote=F, sep="\t")
	
	clinical_deltas <- map[ix_t2, clinical.vars] - map[ix_t1, clinical.vars]
	colnames(clinical_deltas) <- paste("delta_",clinical.vars, sep="")
	map_delta <- cbind(map[ix_t2, ], clinical_deltas)

	if(1 %in% run.test){
		test.delta.interactions(otu_delta=otu_delta, map_delta=map_delta, clinical.vars=clinical.vars, controls=controls)
	} 
	if(2 %in% run.test) {
		delta.taxon.v.delta.covariate.within.group(otu_delta=otu_delta, map_delta=map_delta, clinical.vars=clinical.vars, plot.pvals=F, controls=controls)
	} 
	if (3 %in% run.test) {
		taxon.v.covariate(otu=otu[ix_t1,], map=map[ix_t1,], clinical.vars=clinical.vars, controls=controls) #at baseline
		taxon.v.covariate(otu=otu[ix_t2,], map=map[ix_t2,], clinical.vars=clinical.vars, timepoint=2, controls=controls)
	} 
	if (4 %in% run.test){	
		test.alpha.div(samples1=rownames(map)[ix_t1], samples2=rownames(map)[ix_t2], map_delta=map_delta, alphafile=alphafile,clinical.vars=clinical.vars, controls=controls, outputfile="alpha.div.results.all.groups.txt")
		
		participant_ids_A <- unique(map[map$GROUP=="A","PARTICIPANTID"]) # necessary for maintaining order between timepoints
		participant_ids_B <- unique(map[map$GROUP=="B","PARTICIPANTID"]) # necessary for maintaining order between timepoints
		samples_A1 <- rownames(map)[map$PRE_POST==1 & map$PARTICIPANTID %in% participant_ids_A]
		samples_A2 <- rownames(map)[map$PRE_POST==2 & map$PARTICIPANTID %in% participant_ids_A]
		samples_B1 <- rownames(map)[map$PRE_POST==1 & map$PARTICIPANTID %in% participant_ids_B]
		samples_B2 <- rownames(map)[map$PRE_POST==2 & map$PARTICIPANTID %in% participant_ids_B]

		test.alpha.div(samples1=samples_A1, samples2=samples_A2, map_delta=map_delta[samples_A2,], alphafile=alphafile,clinical.vars=clinical.vars, controls=controls, outputfile="alpha.div.results.GroupA.txt")

		test.alpha.div(samples1=samples_B1, samples2=samples_B2, map_delta=map_delta[samples_B2,], alphafile=alphafile,clinical.vars=clinical.vars, controls=controls, outputfile="alpha.div.results.GroupB.txt")

	}
	if(5 %in% run.test) {
		delta.taxon.v.delta.covariate.within.group(otu_delta=otu_delta, map_delta=map_delta, clinical.vars=clinical.vars, plot.pvals=T, controls=controls)
	} 
	if (6 %in% run.test) { # delta MB predicts delta clinical var
		predict.clinical.change(otu=otu_delta, map_delta=map_delta, clinical.vars=clinical.vars, n.iterations=10, otu.delta=TRUE)
	}
	if (7 %in% run.test) {
		delta.taxon.v.delta.covariate.all.groups(otu_delta=otu_delta, map_delta=map_delta, clinical.vars=clinical.vars, plot.pvals=T, controls=controls)
	}
	if (8 %in% run.test) { # initial MB predicts delta clinical var
		initial.otu <- otu[ix_t1,]
		rownames(initial.otu) <- rownames(otu[ix_t2,]) # just rename the rownames to be the 2nd timepoints because this is how the delta_map got renamed
		predict.clinical.change(otu=initial.otu, map_delta=map_delta, clinical.vars=clinical.vars, n.iterations=10, otu.delta=FALSE)
	}
	if(9 %in% run.test) {
		# this test runs the typical test, but instead controls for baseline clinical var and baseline taxa value
		rownames(otu_delta) <- rownames(map)[ix_t1] #participant_ids at timepoint 1 (not relevant since these are deltas but just important for merging to map purposes)
		map_delta <- cbind(map[ix_t1, ], clinical_deltas) # for use with controlling for baseline clinical vars 

print("--------------------0WITHIN GROUP COMPARISONS")
		delta.taxon.v.delta.covariate.within.group(otu_delta=otu_delta, map_delta=map_delta, clinical.vars=clinical.vars, plot.pvals=F, controls=controls, control.baseline=T, otu.baseline=otu[ix_t1,])
print("--------------------0ALL GROUP COMPARISONS")
		delta.taxon.v.delta.covariate.all.groups(otu_delta=otu_delta, map_delta=map_delta, clinical.vars=clinical.vars, plot.pvals=F, controls=controls, control.baseline=T, otu.baseline=otu[ix_t1,])
	} 
	if (10 %in% run.test) {
		taxon.v.covariate(otu=otu_delta, map=map[ix_t2,], clinical.vars=clinical.vars, timepoint=2, controls=controls)
	}
	if (11 %in% run.test) { # delta MB predicts measpredREE_2
		predict.clinical.change(otu=otu_delta, map_delta=map_delta, clinical.vars=clinical.vars, n.iterations=10, otu.delta=TRUE)
	}

}