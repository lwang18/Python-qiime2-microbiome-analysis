# Python-qiime2-microbiome-analysis
Bacterial 16S &amp; fungal ITS analysis pipelines with Qiime2

1. Installation
```
*Install Miniconda*
conda update conda
conda install wget

*Install QIIME 2 within a conda environment*
wget https://data.qiime2.org/distro/core/qiime2-2019.10-py36-osx-conda.yml
conda env create -n qiime2-2019.10 --file qiime2-2019.10-py36-osx-conda.yml

*OPTIONAL CLEANUP*
rm qiime2-2019.10-py36-osx-conda.yml

*Activate the conda environment*
conda activate qiime2-2019.10
conda deactivate
```

Convert table.qza into csv file
```
qiime taxa collapse \
  --i-table table-with-phyla-no-mitochondria-no-chloroplast.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 7 \
  --output-dir level7/
 
qiime tools export \
--input-path level7/collapsed_table.qza \
--output-path level7/

biom convert -i level7/feature-table.biom -o level7/feature-table.tsv  --to-tsv --header-key taxonomy

$ <feature-table_level5tax.csv awk -F"\t" 'BEGIN{print "OTUs,Domain,Phylum,Class,Order,Family,Genus"}{gsub(" ","_",$0);gsub("\"","",$0);print $1","$3","$6","$9","$12","$15","$18}' > feature-table_level5_SplitTax.csv

```

Calypso needs the following files:
- metadata
- BIOM: feature-table.biom (OTU table)
- taxonomy.tsv

```
*1. Convert table.qza into BIOM file *
qiime tools export \
--input-path "0. Qiime2 files"/table-with-phyla-no-mitochondria-no-chloroplast.qza \
--output-path "0. Qiime2 files"/
**Output: feature-table.biom **

*2. Convert taxonomy.qza into taxonomy.tsv *
qiime tools export \
--input-path "0. Qiime2 files"/taxonomy.qza \
--output-path "0. Qiime2 files"/
**Output: taxonomy.tsv **
```
  
Qiime ITS

```
#### Step 1
Mkdir emp-paired-end-sequences

#Gzip both these files
gzip Undetermined_S0_L001_I1_001.fastq
gzip Undetermined_S0_L001_I1_001.fastq

#Mapping file ’sample-metadata.tsv’
qiime tools import \
  --type EMPSingleEndSequences \
  --input-path emp-single-end-sequences \
  --output-path emp-single-end-sequences.qza
```

```
#### Step 2
#rev comp barcodes are needed

qiime demux emp-single \
  --i-seqs emp-single-end-sequences.qza \
  --m-barcodes-file 3rd_its_sample-metadata.tsv \
  --m-barcodes-column BarcodeSequence \
  --p-rev-comp-barcodes \
  --o-per-sample-sequences demux.qza \
  --o-error-correction-details demux-details.qza

qiime demux summarize \
  --i-data demux2.qza \
  --o-visualization demux2.qzv
qiime tools view demux2.qzv


#Interactive tool
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 251 \
  --o-representative-sequences rep-seqs2.qza \
  --o-table table2.qza \
  --o-denoising-stats stats2.qza
```

```
qiime metadata tabulate \
  --m-input-file stats2.qza \
  --o-visualization stats2.qzv
qiime tools view stats2.qzv
qiime feature-table summarize \
  --i-table table2.qza \
  --o-visualization table2.qzv \
  --m-sample-metadata-file sample-metadata.tsv
qiime tools view table2.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs2.qza \
  --o-visualization rep-seqs.qzv
qiime tools view rep-seqs.qzv

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs2.qza \
  --o-alignment aligned-rep-seqs2.qza \
  --o-masked-alignment masked-aligned-rep-seqs2.qza \
  --o-tree unrooted-tree2.qza \
  --o-rooted-tree rooted-tree2.qza
```

Ref: https://github.com/gregcaporaso/2017.06.23-q2-fungal-tutorial
```
wget https://unite.ut.ee/sh_files/sh_qiime_release_20.11.2016.zip
unzip sh_qiime_release_20.11.2016.zip
gunzip sh_qiime_release_s_04.02.2020.gz

qiime tools import \
 --type FeatureData[Sequence] \
 --input-path sh_qiime_release_s_04.02.2020\
 --output-path unite97-sh_qiime_release_s_04.02.2020.qza
qiime tools import \
 --type FeatureData[Sequence] \
 --input-path sh_qiime_release_20.11.2016/sh_refs_qiime_ver7_99_20.11.2016.fasta \
 --output-path unite-ver7-99-seqs-20.11.2016.qza

qiime tools import \
 --type FeatureData[Sequence] \
 --input-path sh_qiime_release_20.11.2016/sh_refs_qiime_ver7_99_20.11.2016.fasta \
 --output-path unite-ver7-99-seqs-20.11.2016.qza

qiime tools import \
 --type FeatureData[Taxonomy] \
 --input-path sh_qiime_release_20.11.2016/sh_taxonomy_qiime_ver7_99_20.11.2016.txt \
 --output-path unite-ver7-99-tax-20.11.2016.qza \
 --input-format HeaderlessTSVTaxonomyFormat

qiime feature-classifier fit-classifier-naive-bayes \
 --i-reference-reads unite-ver7-99-seqs-20.11.2016.qza \
 --i-reference-taxonomy unite-ver7-99-tax-20.11.2016.qza \
 --o-classifier unite-ver7-99-classifier-20.11.2016.qza

qiime feature-classifier classify-sklearn \
  --i-classifier unite-ver7-99-classifier-20.11.2016.qza \
  --i-reads rep-seqs2.qza \
  --o-classification taxonomy.qza
#shows confidence in taxons found in all data
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy2.qzv
qiime tools view taxonomy2.qzv

qiime taxa barplot \
  --i-table table2.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata_HETEDIT.tsv \
  --o-visualization taxa-bar-plots.qzv
qiime tools view taxa-bar-plots.qzv

```

Phylogeny for ITS 
```
qiime diversity core-metrics \
  --i-table table2.qza \
  --p-sampling-depth 100 \
  --m-metadata-file sample-metadata_HETEDIT.tsv \
  --output-dir core-metrics-results_depth2
```

Alpha diversity
```
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results_depth2/shannon_vector.qza \
  --m-metadata-file sample-metadata_HETEDIT.tsv\
  --o-visualization core-metrics-results_depth2/shannon-significance.qzv
qiime tools view core-metrics-results_depth2/shannon-significance.qzv
```

Beta diversity
```
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_depth2/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample-metadata_HETEDIT.tsv\
  --m-metadata-column genotype_location \
  --o-visualization core-metrics-results_depth2/bray_curtis_geno_location.qzv \
  --p-pairwise
qiime tools view core-metrics-results_depth2/bray_curtis_geno_location.qzv 

qiime emperor plot \
  --i-pcoa core-metrics-resultsnophylo500/bray_curtis_pcoa_results.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-resultsnophylo500/bray_curtis_healthlocation.qzv

qiime tools view core-metrics-resultsnophylo500/bray_curtis_healthlocation.qzv
```
