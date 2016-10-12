README for /metagenome output/ folder

This folder contains the tests run for analysis using PICRUST metagenome predictions (as opposed to taxa relative abundances).

## Compare change in taxon between groups

./significant-delta-mg-v-group.txt
pvalues, adjusted pvalues, and coefficients of significant associations between change in relative abundance of pathways in Group A vs. Group B using LM and spearman's correlation

./delta-mg-v-group.txt
all pvalues, adjusted pvalues, and coefficients for testing for linear association between change in relative abundance of pathways between Group A and Group B

## Compare taxon and covariates at baseline

./significant-mg-v-covariate-baseline.txt
pvalues, adjusted pvalues, and coefficients of significant associations between pathways and clinical covariates at baseline as identified by LM and spearman's correlation

./mg-v-covariate-baseline.txt
all LM pvalues, adjusted pvalues, and coefficients when testing for associations between pathways and clinical covariates at baseline 

## Compare taxon between timepoints stratified by group

./mg-v-group.txt
pvalues, adjusted pvalues, and coefficients of LM associations between relative abundance of pathways at timepoint 1 and timepoint 2, stratified by treatment group. 

## Compare change in taxon against change in covariate stratified by group

./delta-mg-v-delta-covariate-within-group.pdf
scatter plots of significant associations between change in relative abundance of pathways and change in clinical covariates within each treatment group as identified by LM and spearman's correlation

./significant-delta-mg-v-delta-covariate-within-group.txt
pvalues, adjusted pvalues, and coefficients of significant associations between change in relative abundance of pathways and change in clinical covariates within each treatment group as identified by LM and spearman's correlation

./delta-mg-v-delta-covariate-within-group-correlation.txt
spearman's correlation adjusted pvalues for change in relative abundance of pathways and change in clinical covariates within each treatment group

./delta-mg-v-delta-covariate-within-group.txt
pvalues, adjusted pvalues, and coefficients of all associations between change in relative abundance of pathways and cahnge in clinical covariates within each treatment group as identified by LM and spearman's correlation

./BOXPLOT-pvalues-v-covariate.pdf
Boxplot of all pvalues of whether specific changes in pathways are significantly associated with changes in clinical covariates, across both treatment groups combined. Points above the dashed line denote pathways that were significantly associated using a linear model but may not appear in significant-delta-mg-v-delta-covariate-all-groups.txt because they weren't also flagged as significantly correlated. Therefore some of these pathways may be outliers.
