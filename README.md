# tufts-fiber
Analysis of tufts whole-grain fiber study 201

# NOTE
Additional readmes with detailed descriptions of analyses and results can be found in corresponding /results folders

# QIIME processing steps

## preprocessing
cd raw
## join paired ends
multiple_join_paired_ends.py -i run1 -o run1-join --read1_indicator ".R1." --read2_indicator ".R2."

cd run1-join
## rename joined fastq files with their sample IDs
for folder in *; do mv $folder/fastqjoin.join.fastq $folder.fastq; done
## convert fastq to fasta
for file in *.fastq; do echo $file; convert_fastaqual_fastq.py -f $file -c fastq_to_fastaqual; done
## combine all fasta's and add sample IDs to fasta sequence headers
time add_qiime_labels.py -i . -m ../../sample-map-run1-final.txt -c FileName -o .
cd ..

## same preprocessing for run2
cd run2-join
for folder in *; do mv $folder/fastqjoin.join.fastq $folder.fastq; done
for file in *.fastq; do echo $file; convert_fastaqual_fastq.py -f $file -c fastq_to_fastaqual; done
time add_qiime_labels.py -i . -m ../../sample-map-run2-final.txt -c FileName -o .
cd ..

cd ..
## merge run1 and run2
cat raw/run1-join/combined_seqs.fna raw/run2-join/combined_seqs.fna > merged.fna


## pick open-ref otus
## with params file contents:
## assign_taxonomy:assignment_method	sortmerna
time pick_open_reference_otus.py -i merged.fna -r ~/public/ref/gg_13_8_otus/rep_set/97_otus.fasta -o otus-smr -p params.txt -m sortmerna_sumaclust -v 


## post-processing
cd otus-smr

## drop OTUs in < 10% of samples
filter_otus_from_otu_table.py -i otu_table_mc2_w_tax.biom -s 18 -o otu_table_mc2_w_tax_s18.biom

## get OTU table stats
biom summarize-table -i otu_table_mc2_w_tax.biom -o stats.txt
biom summarize-table -i otu_table_mc2_w_tax_s18.biom -o stats_s18.txt

## get taxon summaries
summarize_taxa.py -i otu_table_mc2_w_tax.biom -L 2,3,4,5,6,7 -o taxa

## drop taxa in < 10% of samples
for file in taxa/*.biom; do echo $file; filter_otus_from_otu_table.py -i $file -s 18 -o taxa/`basename $file .biom`_s18.biom; done
for file in taxa/*_s18.biom; do echo $file; biom convert -i $file -o taxa/`basename $file .biom`.txt --to-tsv --table-type "Taxon table"; done

## calculate alpha and beta diversity metrics
beta_diversity.py -i otu_table_mc2_w_tax.biom -t rep_set.tre -o beta
alpha_diversity.py -i otu_table_mc2_w_tax.biom -t rep_set.tre -o alpha.txt
cd ..

tar cvzf otus-smr.tgz otus-smr/alpha.txt otus-smr/beta/ otus-smr/log_20150512213758.txt otus-smr/otu_table_mc2_w_tax.biom otus-smr/otu_table_mc2_w_tax_s18.biom otus-smr/rep_set.fna otus-smr/rep_set.tre otus-smr/stats.txt otus-smr/taxa/*tax_L*s18*

