README for /metagenome output/ folder

This folder contains the tests run for analysis using genus level taxon tables.

## Compare diversity and clinical variables

./alpha.div.results.txt
pvalues, adjusted pvalues, and coefficients of associations between (a) alpha diversity and change in clinical variables, and (b) change in alpha diversity and change in clinical variables

## Compare change in taxa between groups

./delta-taxon-v-group.pdf
beeswarm plots of any significant associations between change in taxa between groups
./delta-taxon-v-group.txt
pvalues, adjusted pvalues, and coefficients of associations between change in taxa between groups
./significant-delta-taxon-v-group.txt
pvalues, adjusted pvalues, and coefficients of significant associations between change in taxa between groups

## Compare taxon and clinical covariate at baseline

./taxon-v-covariate-baseline.pdf
scatter plots of significant associations between taxon and clinical covariate at baseline as identified by LM and spearman's correlation
./significant-taxon-v-covariate-baseline.txt
significant pvalues, adjusted pvalues, and coefficients when testing for associations between taxon and clinical covariates at baseline 
./taxon-v-covariate-baseline.txt
all pvalues, adjusted pvalues, and coefficients when testing for associations between taxon and clinical covariates at baseline 

## Compare change in taxon against change in covariate stratified by group

./delta-taxon-v-delta-covariate-all-groups.pdf
scatter plots of significant associations between change in relative abundance of taxa and change in clinical covariates within each treatment group as identified by LM and spearman's correlation

./significant-delta-taxon-v-delta-covariate-all-groups.txt
pvalues, adjusted pvalues, and coefficients of significant associations between change in relative abundance of taxa and change in clinical covariates within each treatment group as identified by LM and spearman's correlation

./delta-taxon-v-delta-covariate-all-groups.txt
pvalues, adjusted pvalues, and coefficients of all associations between change in relative abundance of taxa and change in clinical covariates within each treatment group as identified by LM and spearman's correlation

./delta-taxon-v-delta-covariate-all-groups-correlation.txt
spearman's correlation adjusted pvalues for change in relative abundance of taxa and change in clinical covariates within each treatment group

./BOXPLOT-pvalues-v-covariate-all-groups.pdf
Boxplot of all pvalues of whether specific changes in taxa are significantly associated with changes in clinical covariates, across both treatment groups combined. Points above the dashed line denote bugs that were significantly associated using a linear model but may not appear in significant-delta-taxon-v-delta-covariate-all-groups.txt because they weren't also flagged as significantly correlated. Therefore some of these bugs may be outliers.