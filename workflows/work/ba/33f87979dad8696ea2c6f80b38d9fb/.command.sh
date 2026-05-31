#!/bin/bash -ue
biom convert \
       -i asv_table.tsv \
       -o table.biom \
       --table-type="OTU table" \
       --to-hdf5

   qiime tools import \
       --type 'FeatureTable[Frequency]' \
       --input-path table.biom \
       --input-format BIOMV210Format \
       --output-path table.qza

   qiime tools import \
       --type 'FeatureData[Sequence]' \
       --input-path rep-seqs.fasta \
       --output-path rep-seqs.qza

   qiime phylogeny align-to-tree-mafft-fasttree \
       --i-sequences rep-seqs.qza \
       --o-alignment aligned-rep-seqs.qza \
       --o-masked-alignment masked-aligned-rep-seqs.qza \
       --o-tree unrooted-tree.qza \
       --o-rooted-tree rooted-tree.qza

   qiime feature-table filter-samples \
       --i-table table.qza \
--m-metadata-file sample-metadata.tsv \
--p-where "[healthy] IN ('yes', 'no')" \
--o-filtered-table table_filtered.qza
