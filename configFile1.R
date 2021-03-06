### Other script needed ###
pathForFunctionsBrowser<-"/home/ldelisle/Dropbox/scripts/toShare/BrowserFromAnoukAdaptedByLD/plot.functions_v04.R"


### Plotting area ###
chrToPlot<-"chr2" #Put the UCSC name.
wtStart<-74000000
wtEnd<-75000000  
genome<-"Inv_10lox11d13_lacZ_v2" #Name of the genome on which you want to plot (This need to be the genome in which the analysis have been done). Put wt if it is the wild-type genome else it needs to be a mutant in the brFile below.

### Global plotting parameters ###
leftMargin<-0.3 #In inches If you change it, the size of the pdf will be adjusted so the annotations will be the same size (default=0.3).
rightMargin<-0.02 #In inches If you change it, the size of the pdf will be adjusted so the annotations will be the same size (default=0.02).
removeExtremeValuesOfXAxis<-F #If you want not to plot the first and the last values of the x axis to be sure to fit in the margins.


###To plot on a mutant genome###
pathForFunctionsShift<-"/home/ldelisle/Dropbox/scripts/toShare/BrowserFromAnoukAdaptedByLD/20171204_shiftAnnotFunctions_compatibleInv.R"

### Mutant information ###
brFile<-"/home/ldelisle/Documents/sequencing/invertedGenomes/br_inv_final.txt" #used only if you want to plot on a mutant genome
chrMutant<-"chr2" #Which chromosome is concerned by the br file. Put the UCSC name.

### Annotations ###
##Genes
#If you already created the mutant Gtf file put it here
mutantGtfFile<-"/home/ldelisle/Documents/sequencing/testsBrowser/Inv_11-12_d13Z_vSplit_shifted_fromFilteredTranscriptsOfMus_musculus.GRCm38.90.gz_ExonsOnly_UCSC_pluslacZ_181217_corrected.gtf.gz" 
#Else put the wt gtf file and the annotations in the plotting area will be shifted.
#wtGtfFile<- 
biotypesToKeep<-c("protein_coding") #If you want all biotypes put NA, else put the list of biotypes you want for example, c("lincRNA","miRNA","protein_coding")
#All possible biotypes are c( "TEC","snRNA","protein_coding","processed_pseudogene","antisense","sense_intronic","lincRNA","processed_transcript","miRNA","snoRNA","misc_RNA","transcribed_unprocessed_pseudogene","unprocessed_pseudogene","unitary_pseudogene","sense_overlapping","rRNA","transcribed_processed_pseudogene","ribozyme","scaRNA","pseudogene","polymorphic_pseudogene","macro_lncRNA","bidirectional_promoter_lncRNA","3prime_overlapping_ncRNA","transcribed_unitary_pseudogene","TR_V_gene","TR_V_pseudogene","TR_D_gene","TR_J_gene","TR_C_gene","TR_J_pseudogene","IG_LV_gene","IG_V_gene","IG_V_pseudogene","IG_J_gene","IG_C_gene","sRNA","scRNA","IG_C_pseudogene","IG_D_gene","IG_D_pseudogene","IG_pseudogene","Mt_tRNA","Mt_rRNA")
genesToExclude<-NA #If you do not have genes to exclude from the plot, put NA, else put the list, for example: c("ENSMUSG00000099521","ENSMUSG00000086077") to remove Gm28309 and Gm14396...
col.genes.fwd<-"blue" #Put the color of genes fwd
col.genes.rev<-"red" #Put the color of genes rev
plot.axis.coordinates.Genes<-T #Do you want to plot the x axis below the genes:(T for yes and F for no)
coeffToLabGenes<-1 #If you want a bigger writting for genes names increase this number (or decrease for smaller)
col.lacZ<-"seagreen" #If you want a special color for the exon named "lacZ", else put NA
simplifyHoxGenesNames<-T #If you want to have d1 instead of Hoxd1:(T for yes and F for no)

##Annotations
#mutantAnnotationBedFile<-"" #If you already created the mutant annotation file put it here.
wtAnnotationBedFile<-"/home/ldelisle/Documents/sequencing/testsBrowser/annotations_mm10_simplified.txt" #Else, put the wt one and the annotations in the plotting area will be shifted.
col.annotations<-"black" #Put the color of the annotations.
coeffToLabAnnot<-1 #If you want a bigger writting for annotation names increase this number (or decrease for smaller)
annotToExclude<-NA #If you do not have annotation to exclude from the plot, put NA, else put the list, for example: c("GT1","GT2")

##If trying to plot on mutant with wt annotations
splitIfOverlap<-T #Specify if you want to split the annotation which totally cover a transgene or a deletion or an invertion (for inversion it would impact only for genes).

### Tracks to plot ###
trackRNAseqMetanephroFwd<-list(file="/home/ldelisle/Documents/sequencing/testsBrowser/accepted_hits_unique_coverage_positive_strand.bedgraph.gz",#Put the absolute path to the bedgraph (it can be gzip)
                     nb.reads=54727620, #Put the number of reads (see details in the presentation). Only used if normalize=T.
                     color="blue", #Put the color you want to use for this track.
                     plot.axis.coordinates=T, #Do you want to see the x axis with the coordinates (T for yes and F for no)
                     normalize=T, #Do you want to normalize to the number of million mapped reads (T for yes and F for no)
                     coeffToLabAxes=1.3, #If you want a bigger writting for axis increase this number (or decrease for smaller)
                     recompute.ylim=T, #This means that you want to have an adaptative limit for the y axis (T for yes and F for no).
                     ylim=NA, #If you want a fixed limit for y axis put the extremities here, separated by a comma for example c(0,10). Else put NA.
                     axis.interval=NA) #If you want to control the tick-marks, put here the values for example c(1,5,10) or if it is always the same interval you can write seq(0,10,3) which means from 0 to 10 every 3. In this case it will be 0,3,6 and 9.
#Check that each track has a closing parenthesis after axis.interval.
                     
trackRNAseqMetanephroRev<-list(file="/home/ldelisle/Documents/sequencing/testsBrowser/accepted_hits_unique_coverage_negative_strand.bedgraph.gz",#Put the absolute path to the bedgraph (it can be gzip)
                     nb.reads=54727620, #Put the number of reads (see details in the presentation). Only used if normalize=T.
                     color="red", #Put the color you want to use for this track.
                     plot.axis.coordinates=T, #Do you want to see the x axis with the coordinates (T for yes and F for no)
                     normalize=T, #Do you want to normalize to the number of million mapped reads (T for yes and F for no)
                     coeffToLabAxes=1, #If you want a bigger writting for axis increase this number (or decrease for smaller)
                     recompute.ylim=T, #This means that you want to have an adaptative limit for the y axis (T for yes and F for no).
                     ylim=NA, #If you want a fixed limit for y axis put the extremities here, separated by a comma for example c(0,10). Else put NA.
                     axis.interval=NA) #If you want to control the tick-marks, put here the values for example c(1,5,10) or if it is always the same interval you can write seq(0,10,3) which means from 0 to 10 every 3. In this case it will be 0,3,6 and 9.
#Check that each track has a closing parenthesis after axis.interval.

### Output ###
#Output location
outputPath="/home/ldelisle/Documents/sequencing/testsBrowser/figures/example" #The output file (without extension)
addExtension=T #Do you want to add to the output name the genome name and the WT coordinates.
#Output format
usePng<-F #If you want to use png replace F by T
pngRes<-96 #This is the resolution of the png file.

