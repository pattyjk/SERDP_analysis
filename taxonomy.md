## Process taxonomy

Summaries in QIIME2
```
qiime taxa collapse --i-taxonomy Taxonomy.qza --o-collapsed-table tax_table.txt --i-table SERDP_NOVI_rarefied_table_4600.qza --p-level 7
qiime tools export --input-path tax_table.txt.qza --output-path tax_table
```

R analysis
```
#packaging
library(ggplot2)
library(plyr)
library(tidyr)
library(reshape2)
library(stringi)

#read ASV table
asv.tbl<-read.delim("/Users/patty/Desktop/SERDP/SERDP/asv_table.txt", row.names=1, header=T)

#read metadata for sequenced samples only
serdp.meta<-read.delim("/Users/patty/Desktop/serdp_meta_only_seq.txt", header=T, na.strings ="N/A")

#read taxonomy
taxo<-read.delim("/Users/patty/Desktop/tax_table.txt", header=T, row.names=1)

#convert from abundance to relative abundance
taxo<-sweep(taxo, 2, colSums(taxo), '/')

#get row sums (e.g. of taxonomy) & order by abundance
taxo$sum<-rowSums(taxo)
taxo<-taxo[order(taxo$sum, decreasing=T) , ]

#get top twent taxa
taxo<-taxo[1:20,]


#remove sum column
taxo<-as.data.frame(taxo[,-grep('sum', names(taxo))])

#get 'others' category (things not in top20)
others<-1-colSums(taxo)
taxo<-rbind(taxo, others)
rownames(taxo)[21]<-"Others;Others;Others;Others;Others;Others;Others"

#add taxonomy back
taxo$taxonomy<-row.names(taxo)

#melt data
taxo_m<-melt(taxo)
names(taxo_m)<-c("Taxonomy", "SampleID", "Rel_abun")

#change to percent (out of 100)
taxo_m$Rel_abun<-taxo_m$Rel_abun*100

#split taxonomy
taxo_m_split<-separate(taxo_m, Taxonomy, sep=";", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#remove added 'X' from sample ID
taxo_m_split$SampleID<-gsub("X","",as.character(taxo_m_split$SampleID))

#bind metadata to taxonomy
taxo_m_split<-merge(taxo_m_split, serdp.meta, by="SampleID")

#read in the best pallette
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")

#remove other sample type ('NA')

#plot data as stacked barplot
ggplot(taxo_m_split, aes(SampleID, Rel_abun, fill=Genus.x))+
  geom_bar(stat='identity')+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=pal)+
  guides(fill=guide_legend(ncol=1))+
  xlab("")+
  ylab("Relative Abundance")+
  theme_bw()+
  facet_wrap(~Location, scales = 'free')+
  theme(text = element_text(size=14), axis.text.x = element_blank())

```
