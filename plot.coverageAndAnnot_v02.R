options(scipen=999)
options(stringsAsFactors=F)
rm(list=ls())
if(length(commandArgs(TRUE))>0){
  f<-commandArgs(TRUE)[1]
} else{
  #Ask for the config file
  #cat("Choose the configFile.\n")
  f <-"/Users/Hintermann/Desktop/4C_wp_tb_brain/4C_tbw_plots/configFile1_4Ctbw.R" #file.choose()
}
### Check the necessary options ###
#The config file
if(!file.exists(f)){
  stop("This file does not exist.")
}
source(f)

#The functions to plot
if(!exists("pathForFunctionsBrowser")){
  stop("No path for functions.")
} else if(!file.exists(pathForFunctionsBrowser)){
  stop("This file does not exist: ",pathForFunctionsBrowser)
}
source(pathForFunctionsBrowser)

#The region to plot
if(!exists("chrToPlot")){
  stop("The config file do not have chrToPlot definition.")
}
if(!exists("wtStart")){
  stop("The config file do not have wtStart definition.")
}
if(!is.numeric(wtStart)){
  stop("The wtStart is not numeric.")
}
if(!exists("wtEnd")){
  stop("The config file do not have wtEnd definition.")
}
if(!is.numeric(wtEnd)){
  stop("The wtEnd is not numeric.")
}
regionToPlot<-data.frame(chr=chrToPlot,start=wtStart,end=wtEnd)

#The genome
if(!exists("genome")){
  cat("The genome is not defined. Assuming it is the wild-type genome.\n")
  genome<-"wt"
}
#I need to define a fake function shiftDFFromBR
shiftDFFromBR<-function(annotationData,genome,br,colChr,colStart,colEnd,colStrand=0,verbose=T,chromoWithTg="chr2",splitIfOverlap=F){
  if(genome=="mm10"|tolower(genome)=="wt"){
    return(annotationData)
  } else {
    stop("The definition of shiftDFFromBR is not adapted.")    
  }
}

if(tolower(genome) != "wt"){
  #This is a mutant genome.
  #Check if it is the same chr which is affected by mutation and to plot
  if(!exists("chrMutant")){
    stop("The genome to plot on is not wt and the chrMutant is not defined in the configFile.")
  }
  if(chrToPlot==chrMutant){
    genomeToPlot<-genome
    #Needs functions
    if(!exists("pathForFunctionsShift")){
      stop("The genome to plot on is not wt and the pathForFunctionsShift is not defined in the configFile.")
    }
    if(!file.exists(pathForFunctionsShift)){
      stop("This file does not exist: ",pathForFunctionsShift)
    }
    source(pathForFunctionsShift)
    if(!exists("brFile")){
      stop("The genome to plot on is not wt and the brFile is not defined in the configFile.")
    }
    if(!file.exists(brFile)){
      stop("This file does not exist: ",brFile)
    }
    brDF<-read.delim(brFile)
    if(!genome%in%brDF[,1]){
      stop("The genome to plot on is not wt and the genome ",genome," is not part of the genome names (first column) of the br file.")
    }
    regionToPlot<-shiftDFFromBR(regionToPlot,genomeToPlot,brDF,1,2,3,colStrand=0,verbose=F,chromoWithTg=chrMutant,splitIfOverlap=F)
    if(nrow(regionToPlot)!=1){
      cat("The region to plot is incompatible with the mutant genome :\n")
      regionToPlot<-shiftDFFromBR(regionToPlot,genomeToPlot,brDF,1,2,3,colStrand=0,verbose=T,chromoWithTg=chrMutant,splitIfOverlap=F)
      stop("No more regions to plot.")
    }
  } else {
    genomeToPlot<-"wt"
    brDF<-NULL
    splitIfOverlap<-F
  }
} else {
  genomeToPlot<-"wt"
  brDF<-NULL
  chrMutant<-NULL
  splitIfOverlap<-F
}

### Genes ###
if(exists("mutantGtfFile")){
  if(!file.exists(mutantGtfFile)){
    cat("This file for mutantGtf does not exist: ",mutantGtfFile,".\n")
  } else {
    genesToPlot<-read.gtf.from.file(mutantGtfFile, chrToPlot)
  }
}
if(!exists("genesToPlot")){
  #It failed
  if(exists("wtGtfFile")){
    if(!file.exists(wtGtfFile)){
      cat("This file for wtGtf does not exist: ",wtGtfFile,".\n")
    } else {
      wtGenesToPlot<-read.gtf.from.file(wtGtfFile, chrToPlot)
      if(!exists("splitIfOverlap")){
        cat("splitIfOverlap has not been defined whereas needed. It is set to T.")
        splitIfOverlap<-T
      }
      if(!is.logical(splitIfOverlap)){
        cat("splitIfOverlap has been defined as",splitIfOverlap," whereas T or F expected. It is set to T.")
        splitIfOverlap<-T
      }
      genesToPlot<-shiftDFFromBR(wtGenesToPlot,genomeToPlot,brDF,which(colnames(wtGenesToPlot)=="seqid"),which(colnames(wtGenesToPlot)=="start"),which(colnames(wtGenesToPlot)=="end"),colStrand=which(colnames(wtGenesToPlot)=="strand"),verbose=F,chromoWithTg=chrMutant,splitIfOverlap=splitIfOverlap)
    }
  }
}

### Annotations ###
if(exists("mutantAnnotationBedFile")){
 if(!file.exists(mutantAnnotationBedFile)){
    cat("This bed file for mutantAnnotations does not exist: ",mutantAnnotationBedFile,".\n")
  } else {
    annotToPlot<-read.bed.from.file(mutantAnnotationBedFile)
  }
}
if(!exists("annotToPlot")){
  #It failed
  if(exists("wtAnnotationBedFile")){
    if(!file.exists(wtAnnotationBedFile)){
      cat("This file for wtAnnotations does not exist: ",wtAnnotationBedFile,".\n")
    } else {
      wtAnnotToPlot<-read.bed.from.file(wtAnnotationBedFile)
      if(!exists("splitIfOverlap")){
        cat("splitIfOverlap has not been defined whereas needed. It is set to T.")
        splitIfOverlap<-T
      }
      if(!is.logical(splitIfOverlap)){
        cat("splitIfOverlap has been defined as",splitIfOverlap," whereas T or F expected. It is set to T.")
        splitIfOverlap<-T
      }
      annotToPlot<-shiftDFFromBR(wtAnnotToPlot,genomeToPlot,brDF,1,2,3,colStrand=0,verbose=F,chromoWithTg=chrMutant,splitIfOverlap=splitIfOverlap)
    }
  }
}

### Tracks ###
tracks<-setdiff(ls()[grepl("track",tolower(ls()))],c("tracks","tracksToPlot","check.trackOptions"))
cat(paste("There are",length(tracks),"tracks identified.\n"))
tracksToPlot<-lapply(tracks,function(n){eval(parse(text=n))})
nameTracks<- sapply(tracks,function(f){strsplit(f,"track")[[1]][2]})

### Output ###
if(exists("usePng")){
  if(!is.logical(usePng)){
    cat("usePng is not boolean and is not taken into account.\n The output will be pdf.\n")
    usePng<-F
  }
} else{
  usePng<-F
}

if(usePng){
  if(! exists("pngRes")){
    cat("You need to specify a resolution for png files.(pngRes<-96 or whatever).\n By default it will be 96")
    pngRes<-96
  } else{
    if(! is.numeric(pngRes)){
      cat("The resolution for png files need to be a number, for example pngRes<-96 or whatever.\n By default it will be 96")
      pngRes<-96
    }
  }
}


#plotting parameters
if(!exists("leftMargin")){
  cat("The left Margin has not been defined. 0.3 will be used.\n")
  leftMargin<-0.3
} else {
  if(!is.numeric(leftMargin)){
    cat("leftMargin is not a number. 0.3 will be used.\n")
    leftMargin<-0.3
  }
}

if(!exists("rightMargin")){
  cat("The right Margin has not been defined. 0.02 will be used.\n")
  rightMargin<-0.02
} else {
  if(!is.numeric(leftMargin)){
    cat("rightMargin is not a number. 0.02 will be used.\n")
    rightMargin<-0.02
  }
}

if(!exists("removeExtremeValuesOfXAxis")){
  removeExtremeValuesOfXAxis<-F
} else {
  if(!is.logical(removeExtremeValuesOfXAxis)){
    cat("removeExtremeValuesOfXAxis is not logical. It will be set to F.")
    removeExtremeValuesOfXAxis<-F
  }
}

if(!exists("outputPath")){
  outputPath<-paste0(dirname(f),"/figures/",gsub(" ","_",Sys.time()))
  cat("The plots will be in :")
  cat(outputPath)
  cat("\n")
  dir.create(paste0(dirname(f),"/figures"),showWarnings=F)
} else{
  if(!dir.exists(dirname(outputPath))){
    dir.create(dirname(outputPath))
  }
}

if(!exists("addExtension")){
  addExtension<-F
} else {
  if(!is.logical(addExtension)){
    cat("The addExtension specified is not logical. It will be set to F.")
    addExtension<-F
  } else {
    if(addExtension){
      outputPath<-paste0(outputPath,"_",genome,"_",chrToPlot,":",wtStart,":",wtEnd)
    }
  }
}

#Evaluate the height needed:
heightOfPlot<-0
if(exists("genesToPlot")){
  heightOfPlot<-heightOfPlot+2
}
if(exists("annotToPlot")){
  heightOfPlot<-heightOfPlot+2
}
heightOfPlot<-heightOfPlot+2*length(tracksToPlot)

if(heightOfPlot==0){
  stop("There is nothing to plot.")
}

widthOfFile<-8+leftMargin-0.3+rightMargin-0.02

if (!usePng){
  pdf(paste0(outputPath,".pdf"),width=widthOfFile, height=heightOfPlot,title=basename(f))
} else {
  png(paste0(outputPath,".png"),width=widthOfFile, height=heightOfPlot,units="in",res=pngRes)
}
par(mfrow=c(heightOfPlot/2,1))
if(exists("genesToPlot")){
  ##I need to find another way to check the color
  genesPara<-check.genesOptions()
  plot.genes(genesToPlot, regionToPlot[1,1], regionToPlot[1,2], regionToPlot[1,3], biotypes=genesPara$biotypesToKeep, plot.axis=genesPara$plot.axis.coordinates.Genes,
             excluded=genesPara$genesToExclude, col.fwd=genesPara$col.genes.fwd, col.rev=genesPara$col.genes.rev, cex.names=genesPara$coeffToLabGenes,
             col.lacZ=genesPara$col.lacZ, simplifyHoxGenesNames=genesPara$simplifyHoxGenesNames,
             removeExtremeValuesOfXAxis=removeExtremeValuesOfXAxis,leftMargin=leftMargin,rightMargin=rightMargin)
}
if(exists("annotToPlot")){
  annotPara<-check.annotOptions()
  plot.annotations(annotToPlot, regionToPlot[1,1], regionToPlot[1,2], regionToPlot[1,3], col=annotPara$col.annotations, cex.label=annotPara$coeffToLabAnnot, annotToExclude=annotPara$annotToExclude,
  leftMargin=leftMargin,rightMargin=rightMargin)
}
if(length(tracksToPlot)>0){
  for (i in 1:length(tracksToPlot)){
    track<-tracksToPlot[[i]]
    cat("Checking",nameTracks[i],"\n")
    checkedTrack<-check.trackOptions(track)
    if(!all(is.na(checkedTrack))){
      coverage<-read.coverage.from.file(checkedTrack$file, regionToPlot[1,1], nbReads=checkedTrack$nb.reads, normalize=checkedTrack$normalize)
      plot.coverage(coverage, regionToPlot[1,1], regionToPlot[1,2], regionToPlot[1,3], plot.axis=checkedTrack$plot.axis.coordinates, 
                    recompute.ylim=checkedTrack$recompute.ylim, ylim=checkedTrack$ylim, col=checkedTrack$color, cex.axislab=checkedTrack$coeffToLabAxes,
                    axis.interval=checkedTrack$axis.interval,title=nameTracks[i],
                    removeExtremeValuesOfXAxis=removeExtremeValuesOfXAxis,leftMargin=leftMargin,rightMargin=rightMargin)
    }
  } 
}
dev.off()
