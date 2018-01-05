gft<-read.table("/Users/Hintermann/Desktop/4C_wp_tb_brain/4C_tbw_plots/FilteredTranscriptsOfMus_musculus.GRCm38.91_ExonsOnly_UCSC.gtf", sep="\t", fill = TRUE,skip=0,header=F, quote="\"", stringsAsFactors=F)

gft2<-subset(gft,grepl('chr[12]', gft$V1))
gft3<-subset(gft2,grepl('gene_name [^G]', gft2$V9))

write.table(gft3,file="/Users/Hintermann/Desktop/4C_wp_tb_brain/4C_tbw_plots/FilteredTranscriptsOfMus_musculus.GRCm38.91_ExonsOnly_UCSC_T5.gtf", sep="\t",quote =F,row.names = F,col.names = F)
