## PICRUSt 2 analysis
```
#conert ASV table to biom then txt
qiime tools export --input-path SERDP_NOVI_rarefied_table_4600.qza --output-path asv_table
biom convert -i asv_table/feature-table.biom -o asv_table.txt --to-tsv

#run picrust2 pipeline
picrust2_pipeline.py -s LDM_RepSeqs.fasta -i asv_table.txt -o SERDP_picrust_out

#This is the set of poorly aligned input sequences to be excluded: 33e87c97cad7c2d4a64c7072fd1bbc53, ae8df71631d245aac23521343cde4ec2
```
