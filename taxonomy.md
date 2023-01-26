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

#get top twenty taxa
taxo2<-taxo[1:20,]


#remove sum column
taxo<-as.data.frame(taxo[,-grep('sum', names(taxo))])
taxo2<-as.data.frame(taxo2[,-grep('sum', names(taxo2))])

#get 'others' category (things not in top20)
others<-1-colSums(taxo2)
taxo2<-rbind(taxo2, others)
rownames(taxo2)[21]<-"Others;Others;Others;Others;Others;Others;Others"

#add taxonomy back
taxo$taxonomy<-row.names(taxo)
taxo2$taxonomy<-row.names(taxo2)

#melt data
taxo_m<-melt(taxo)
names(taxo_m)<-c("Taxonomy", "SampleID", "Rel_abun")
taxo_m2<-melt(taxo2)
names(taxo_m2)<-c("Taxonomy", "SampleID", "Rel_abun")

#change to percent (out of 100)
taxo_m$Rel_abun<-taxo_m$Rel_abun*100
taxo_m2$Rel_abun<-taxo_m2$Rel_abun*100

#split taxonomy
taxo_m_split<-separate(taxo_m, Taxonomy, sep=";", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
taxo_m_split2<-separate(taxo_m2, Taxonomy, sep=";", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))


#remove added 'X' from sample ID
taxo_m_split$SampleID<-gsub("X","",as.character(taxo_m_split$SampleID))
taxo_m_split2$SampleID<-gsub("X","",as.character(taxo_m_split2$SampleID))

#bind metadata to taxonomy
taxo_m_split<-merge(taxo_m_split, serdp.meta, by="SampleID")
taxo_m_split2<-merge(taxo_m_split2, serdp.meta, by="SampleID")

#make a combined taxonomy column
taxo_m_split$comb_tax <- paste(taxo_m_split$Kingdom,"_",taxo_m_split$Phylum,"_",taxo_m_split$Class, taxo_m_split$Order.x, "_", taxo_m_split$Family, "_", taxo_m_split$Genus.x)
taxo_m_split2$comb_tax <- paste(taxo_m_split2$Kingdom,"_",taxo_m_split2$Phylum,"_",taxo_m_split2$Class, taxo_m_split2$Order.x, "_", taxo_m_split2$Family, "_", taxo_m_split2$Genus.x)


#get means/sd/se of taxonomic groups by site
taxo_sum<-ddply(taxo_m_split, c("Location", "comb_tax"), summarize, mean=mean(Rel_abun), sd=sd(Rel_abun), n=length(Rel_abun), se=sd/n)

#identify things with >1% relative abundance
one_abun<-which(taxo_sum$mean>1)

#make table only things >1% relative abundance
taxo_sum<-taxo_sum[one_abun,]

#plot data
ggplot(taxo_sum, aes(comb_tax, mean))+
  geom_point(stat="identity", width = 1, color = "black")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)+
  facet_wrap(~Location)+
  xlab("")+
  title("Taxa with >1% relatie abundance")+
  ylab("Mean Relative abundance w/SE")+
  theme_bw()+
  coord_flip()+
  theme(
    plot.title = element_text(size=10), text = element_text(size=10), 
    axis.title.x = element_text(size=10),
    axis.title.y = element_text( size=10)
     
  )

#read in the best pallette
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")


#plot  top twenty taxa as stacked barplot
ggplot(taxo_m_split2, aes(SampleID, Rel_abun, fill=comb_tax))+
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
