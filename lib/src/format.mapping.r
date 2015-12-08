
# this is for preprocessing the General Mills Tufts mapping file only
format.mapping <- function(mapfile="data/GM_alldata_noIDs_formatted.txt", outputfile="data/GM_alldata_noIDs_fully_formatted.txt")
{
	map <- read.table(mapfile,sep='\t',head=T,row=1,comment='')
	colnames(map) <- toupper(colnames(map))

	ix_1 <- grep("_1$", colnames(map))
	ix_2 <- grep("_2$", colnames(map))
	names_1 <- sub("_1$", "", colnames(map)[ix_1])
	names_2 <- sub("_2$", "", colnames(map)[ix_2])
	
	matched_names <- names_2[which(names_2 %in% names_1)]
	leftover_names_1 <- names_1[which(!(names_1 %in% names_2))]

	map1 <- map[,ix_1]
	map2 <- map[,ix_2]
	
	colnames(map1) <- sub("_1$", "", colnames(map1))
	colnames(map2) <- sub("_2$", "", colnames(map2))
	
	new_map <- rbind(map1[,matched_names], map2[,matched_names])

	if(length(leftover_names_1) > 0)
	{
		map2.empty <- matrix(nrow=nrow(map2), ncol=length(leftover_names_1))
		colnames(map2.empty) <- leftover_names_1
	
		new_map <- cbind(new_map, rbind(map1[,leftover_names_1], map2.empty))
	}

	# now add a column for pre or post treatment "Pre_Post", make sure this is a factor
	new_map <- data.frame(new_map, Pre_Post=as.factor(c(rep(1, nrow(map1)), rep(2, nrow(map2)))))
	
	# now add back to the map file everything that didn't end in a _1 or _2 (repeat it twice)
	if(length(map[, -c(ix_1, ix_2)]) > 0)
	{
		new_map <- data.frame(new_map, rbind(map[, -c(ix_1, ix_2)], map[, -c(ix_1, ix_2)]))
	}
	
	participant.ids <- rownames(map)
	rnames <- c(paste("Sample", rownames(map), "1", sep="."), paste("Sample", rownames(map), "2", sep="."))
	rownames(new_map) <- rnames
	colnames(new_map) <- toupper(colnames(new_map))
	new_map <- data.frame(PARTICIPANTID=participant.ids, new_map)
	write.table(new_map, file=outputfile, quote=F, sep="\t")
}

extra.stuff<-function()
{
	# extraneous
	mapfile <- "data/GM_alldata_noIDs_fully_formatted.txt"
	map <- read.table(mapfile,sep='\t',head=T,row=1,comment='')
	# remove old values that exist, we'll replace them
	map <- map[, which(!(colnames(map) %in% c("ARA", "ARB", "ARC", "Timepoint")))]
	brandnew_map <- merge(map, new_map, by=0)
	write.table(brandnew_map, file="map_with_AR_SCFA.txt", quote=F, sep="\t")


	# add new columns before formatting
	map <- read.table("/Users/pvangay/Copy/UMN/KnightsLab/Gen Mills - Tufts/GM_alldata_noIDs_formatted.txt",sep='\t',head=T,row=1,comment='')
	map2 <- read.table("data/Updated AR and Stool Variables.txt",sep='\t',head=T,row=1,comment='')
	map3 <- read.table("data/Acet_to_Prop_Rat.txt",sep='\t',head=T,row=1,comment='')
	map4 <- read.table("data/BPSAVG.txt",sep='\t',head=T,row=1,comment='')

	# remove old values of ARs
	map <- map[, which(!(colnames(map) %in% c("ARA_1", "ARB_1", "ARC_1", "ARA_2", "ARB_2", "ARC_2")))]
	new.map <- merge(map, map2, by=0, all.x=T)
	rownames(new.map) <- new.map$Row.names
	# add Acet to Prop ratio
	new.map <- merge(new.map[,-1], map3[,c("Acettoprop_1","Acettoprop_2")], by=0, all.x=T)
	rownames(new.map) <- new.map$Row.names
	new.map <- merge(new.map[,-1], map4[,c("BPSAVG_1","BPSAVG_2","BPDAVG_1","BPDAVG_2")], by=0, all.x=T)
	rownames(new.map) <- new.map$Row.names
	write.table(new.map, file="temp map.txt", quote=F, sep="\t")
	# REMEMBER TO FIX THE ROWNAMES AFTER THIS!
	
}