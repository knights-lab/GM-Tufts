# Execute analysis for GM Tufts datasets

setwd("/Users/pvangay/Google Drive/UMN/KnightsLab/Gen Mills - Tufts/GitHub")

source("./lib/src/run.alpha.test.r")
source("./lib/src/run.taxon.tests.r")
source("./lib/src/run.clinical.tests.r")
source("./lib/src/predict.clinical.change.r")

source("./lib/src/collapse-features.r")
source('./lib/src/load.data.r')
source('./lib/src/test.otu.features.r')
source('./lib/src/get.next.kegg.r')	
source('./lib/src/filter.pathways.r')	
source('./lib/src/shorten.taxonomy.r')
source("./lib/src/rf.cross.validation.classification.r")


# this mapping file has been formatted to contain all new variables
mapfile <-  "./data/GM_alldata_noIDs_fully_formatted.txt"

otufile_L6 <- "./data/otu_table_mc2_w_tax_L6_s18.txt"
otufile_L2 <- "./data/otu_table_mc2_w_tax_L2_s18.txt"
alphafile <- "./data/alpha_subsampled_19490.txt" 
mg_file <- "./data/predicted.metagenomes.L3.txt"

# 1. Run all variations of testing for differential taxon
run.taxon.tests(mapfile, otufile_L2, outputfile="taxon_differences_L2.txt")
run.taxon.tests(mapfile, otufile_L6, outputfile="taxon_differences_L6.txt")

# 2. Look at differences in Alpha Diversity between groups at Time Point 1 and 2 (separately)
run.alpha.test(mapfile, alphafile, outputfile="alpha_diversity.txt")

# 3a. Use Initial Microbiome to Predict Change in Clinical Variables
clinical.vars <- c("FM", "CD4NAIVE", "STOOLWT", "DAILYSTOOLEN", "ARB", "TDIAM24HR", "CONA_5UG_PER_ML", "PLASMA_IL_6", "LPS_TNF_A", "ISO", "OGTT_AUCG", "ACETATE", "PROPIONATE", "BUTYRATE")	
run.clinical.tests(mapfile=mapfile, otufile=otufile_L6, clinical.vars=clinical.vars, alphafile=alphafile,run.test=c(6,8))
# 3b. Use Delta Microbiome to Predict Change in Clinical Variables
run.clinical.tests(mapfile=mapfile, otufile=otufile_L6, clinical.vars=clinical.vars, alphafile=alphafile,run.test=c(6))

# 4. Find Association and interactions between change in microbiome and change in clinical variables
clinical.vars2 <- c("FM", "CD4NAIVE", "STOOLWT", "DAILYSTOOLEN", "ARB", "CONA_5UG_PER_ML", "LPS_TNF_A", "ISO", "TOTALEM",
"CD4CD25", "TDIAM48HR", "AR", "ACETATE", "PROPIONATE", "BUTYRATE", "OGTT_AUCG", "ADJREE", "ACETTOPROP")
run.clinical.tests(mapfile=mapfile, otufile=otufile_L6, clinical.vars=clinical.vars2, alphafile=alphafile,run.test=c(3,4,5,7))

# 5. PICRUST: Run all variations of testing for differential taxon
run.taxon.tests(mapfile=mapfile, otufile=mg_file, outputfile="metagenome_differences.txt")

# 6. PICRUST: Find Association and interactions between change in microbiome and change in clinical variables
run.clinical.tests(mapfile=mapfile, otufile=mg_file, clinical.vars=clinical.vars2, alphafile=alphafile,run.test=c(3,4,5,7))


# clinical variables by lab
covariates.df <- read.table("./data/covariates_by_lab.txt",sep='\t',head=T,comment='')

	# NIL
	covariates <- toupper(covariates.df[,1])
	covariates <- covariates[covariates != ""]
	# METAGENOME
	run.clinical.tests(mapfile=mapfile, otufile=mg_file, clinical.vars=covariates, alphafile=alphafile,run.test=c(3,4,5))
	# TAXON
	run.clinical.tests(mapfile=mapfile, otufile=otufile_L6, clinical.vars=covariates, alphafile=alphafile,run.test=c(4))

	run.clinical.tests(mapfile=mapfile, otufile=otufile_L6, clinical.vars=covariates, alphafile=alphafile,run.test=c(3,4,5,7))
	run.clinical.tests(mapfile=mapfile, otufile=otufile_L6, clinical.vars=covariates, alphafile=alphafile,run.test=c(6,8))

	# EML
	covariates <- toupper(covariates.df[,2])
	covariates <- covariates[covariates != ""]
	# skip BMI for now - MG analysis
	run.clinical.tests(mapfile=mapfile, otufile=mg_file, clinical.vars=covariates[-2], alphafile=alphafile,run.test=c(3,4,5))
	run.clinical.tests(mapfile=mapfile, otufile=otufile_L6, clinical.vars=covariates[-2], alphafile=alphafile,run.test=c(3,4,5,7))
	run.clinical.tests(mapfile=mapfile, otufile=otufile_L6, clinical.vars=covariates, alphafile=alphafile,run.test=c(6,8))
	# BMI run separately
	run.clinical.tests(mapfile=mapfile, otufile=otufile_L6, clinical.vars=c("BMI","WT"), controls=c("SEX","RACE"), alphafile=alphafile,run.test=c(2,4,5,7))
	run.clinical.tests(mapfile=mapfile, otufile=mg_file, clinical.vars=c("BMI","WT"),controls=c("SEX","RACE"), alphafile=alphafile,run.test=c(2,4,5))

	# VBL - all samples 
	covariates <- toupper(covariates.df[,3])
	covariates <- covariates[covariates != ""]
	run.clinical.tests(mapfile=mapfile, otufile=mg_file, clinical.vars=covariates, alphafile=alphafile,run.test=c(3,4,5))
	run.clinical.tests(mapfile=mapfile, otufile=otufile_L6, clinical.vars=covariates, alphafile=alphafile,run.test=c(6,8))
	run.clinical.tests(mapfile=mapfile, otufile=otufile_L6, clinical.vars=covariates, alphafile=alphafile,run.test=c(3,4,5,7))
	# VBL - partial samples 
	run.clinical.tests(mapfile="./data/GM_alldata_noIDs_fully_formatted_VBL_PartialSamples.txt", otufile=mg_file, clinical.vars=covariates, alphafile=alphafile,run.test=c(3,4,5))
	run.clinical.tests(mapfile="./data/GM_alldata_noIDs_fully_formatted_VBL_PartialSamples.txt", otufile=otufile_L6, clinical.vars=covariates, alphafile=alphafile,run.test=c(4))
	run.clinical.tests(mapfile="./data/GM_alldata_noIDs_fully_formatted_VBL_PartialSamples.txt", otufile=otufile_L6, clinical.vars=covariates, alphafile=alphafile,run.test=c(3,4,5,7))
	run.clinical.tests(mapfile="./data/GM_alldata_noIDs_fully_formatted_VBL_PartialSamples.txt", otufile=otufile_L6, clinical.vars=covariates, alphafile=alphafile,run.test=c(6,8))

	# TMC 
	covariates <- toupper(covariates.df[,4])
	covariates <- covariates[covariates != ""]
	run.clinical.tests(mapfile=mapfile, otufile=mg_file, clinical.vars=covariates, alphafile=alphafile,run.test=c(3,4,5))
	run.clinical.tests(mapfile=mapfile, otufile=otufile_L6, clinical.vars=covariates, alphafile=alphafile,run.test=c(3,4,5,7))
	run.clinical.tests(mapfile=mapfile, otufile=otufile_L6, clinical.vars=covariates, alphafile=alphafile,run.test=c(6,8))
	
	#TKL 
	covariates <- toupper(covariates.df[,5])
	covariates <- covariates[covariates != ""]
	run.clinical.tests(mapfile=mapfile, otufile=mg_file, clinical.vars=covariates, alphafile=alphafile,run.test=c(3,4,5))
	run.clinical.tests(mapfile=mapfile, otufile=otufile_L6, clinical.vars=covariates, alphafile=alphafile,run.test=c(3,4,5,7))
	run.clinical.tests(mapfile=mapfile, otufile=otufile_L6, clinical.vars=covariates, alphafile=alphafile,run.test=c(6,8))










