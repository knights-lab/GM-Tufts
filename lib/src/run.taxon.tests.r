run.taxon.tests <- function(mapfile, otufile, outputfile, keep.pathways.file='data/pathways.to.keep.txt')
{
	minOTUInSamples = .001
	ret <- load.data(otufile=otufile, mapfile=mapfile, minOTUInSamples=minOTUInSamples, normalize=TRUE)
	otu <- ret$otu
	map <- ret$map
	kegg <- ret$kegg
	
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
	
	# Test 1 - Within subject comparisons (Timepoint 1 vs Timepoint 2, stratified by Group A and B)
		ix_A <- which(map$GROUP=="A")
		ix_B <- which(map$GROUP=="B")

		retA <- test.otu.features(otu[ix_A,], map[ix_A,], treatment="PRE_POST", controls=c("SEX","RACE","BMI"))
		cat("Timepoint 1 vs Timepoint 2, Group A\n", file=outputfile)
		write.table(retA[[1]], file=outputfile, append=T, sep="\t", quote=F)

		retB <- test.otu.features(otu[ix_B,], map[ix_B,], treatment="PRE_POST", controls=c("SEX","RACE","BMI"))
		cat("Timepoint 1 vs Timepoint 2, Group B\n", file=outputfile, append=T)
		write.table(retB[[1]], file=outputfile, append=T, sep="\t", quote=F)

	# Test 2a - Between subject comparison (Group A vs. Group B at TimePoint==1)
		ix_1 <- which(map$PRE_POST=="1")
		ret2a <- test.otu.features(otu[ix_1,], map[ix_1,], treatment="GROUP", controls=c("SEX","RACE","BMI"))
		cat("Group A vs. Group B at TimePoint 1\n", file=outputfile, append=T)
		write.table(ret2a[[1]], file=outputfile, append=T, sep="\t", quote=F)

	# Test 2 - Between subject comparison (Group A vs. Group B at TimePoint==2)
		ix_2 <- which(map$PRE_POST=="2")
		ret2 <- test.otu.features(otu[ix_2,], map[ix_2,], treatment="GROUP", controls=c("SEX","RACE","BMI"))
		cat("Group A vs. Group B at TimePoint 2\n", file=outputfile, append=T)
		write.table(ret2[[1]], file=outputfile, append=T, sep="\t", quote=F)

	# Test 3 - Between subject comparison (DELTA_Taxon, Group A vs. Group B, use TimePoint==2 for other vars)
		participant_ids <- unique(map$PARTICIPANTID)
		ix_t1 <- which(map$PRE_POST=="1" & (map$PARTICIPANTID %in% participant_ids))
		ix_t2 <- which(map$PRE_POST=="2" & (map$PARTICIPANTID %in% participant_ids))
		otu_delta <- otu[ix_t2,] - otu[ix_t1,]	
		ret3 <- test.otu.features(otu_delta, map[ix_t2,], treatment="GROUP", controls=c("SEX","RACE","BMI"))
		cat("Change in Taxon Abundances, Group A vs. Group B\n", file=outputfile, append=T)
		write.table(ret3[[1]], file=outputfile, append=T, sep="\t", quote=F)
		
		taxa <- rownames(ret3[[1]][,"adj.pvals",drop=F])[ret3[[1]][,"adj.pvals"] < .25]
		
		write.table(ret3[[1]][taxa,], file="significant-delta-taxon-v-group.txt", sep="\t", quote=F)
			
		participant_ids <- unique(map$PARTICIPANTID)
		ix_t1 <- which(map$PRE_POST=="1" & (map$PARTICIPANTID %in% participant_ids))
		ix_t2 <- which(map$PRE_POST=="2" & (map$PARTICIPANTID %in% participant_ids))
		otu_delta <- otu[ix_t2,] - otu[ix_t1,]
		rownames(otu_delta) <- rownames(map)[ix_t2] # rename to the sample IDs of timepoint 2

		samples_A <- rownames(map[ix_t2,])[map[ix_t2,]$GROUP=="A"]
		samples_B <- rownames(map[ix_t2,])[map[ix_t2,]$GROUP=="B"]

		short.taxa <- shorten.taxonomy(taxa)
		if(length(taxa) > 0)
		{	
			pdf(file="delta-taxon-v-group.pdf")
			for(i in 1:length(taxa))
			{
				beeswarm(x=list(GroupA=otu_delta[samples_A,taxa[i]], GroupB=otu_delta[samples_B,taxa[i]]), main=short.taxa[i],
						ylab="Delta Taxon Abundance")
				bxplot(x=list(GroupA=otu_delta[samples_A,taxa[i]], GroupB=otu_delta[samples_B,taxa[i]]),add=T)
			}
			dev.off()
		}
}
