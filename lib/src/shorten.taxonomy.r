"shorten.taxonomy" <- function(ids,delim=';'){
	ids <- gsub('[kpcofgs]__','',ids)
	newids <- ids
	ids <- strsplit(ids,delim)
	for(i in seq_along(ids)){
		n <- length(ids[[i]])
		j <- n
		while(ids[[i]][j] == 'Other' || ids[[i]][j] == '') j <- j - 1
		newids[i] <- ids[[i]][j]
	}
	return(newids)
}