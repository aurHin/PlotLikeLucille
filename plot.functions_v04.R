if(!"rtracklayer"%in%installed.packages()){
  source("https://bioconductor.org/biocLite.R")
  biocLite("rtracklayer")
}
library(rtracklayer)
library(tools)
options(stringsAsFactors=F)

isWindows<-grepl("windows",tolower(sessionInfo()$running))
if(isWindows){
  cat("You are running on a Windows machine.\nThe access to the coverage and gtf will be done genome wide, it will be a bit longer.\nIf it fails because of lack of memory, please subset it in galaxy (for example only the chromosome you want).\n")
}
#####################################################################################
#####################################################################################

isValidColor <- function(colorname){
  if(is.numeric(colorname)){
    return(TRUE)
  }
  if(colorname %in% colors()){
    return(TRUE)
  }
  if(is.character(colorname)){
    if(nchar(colorname)==7 || nchar(colorname)==9){
      if(substr(colorname,1,1)=="#"){
        #I should do other checks
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

#####################################################################################

simplifyExonsKeepingInfo <- function(gtfLike){
  grgtfLike<-makeGRangesFromDataFrame(gtfLike,keep.extra.columns = T)
  grsimplified<-reduce(grgtfLike)
  origin<-as.data.frame(findOverlaps(grgtfLike,grsimplified))
  mcols(grsimplified)<-aggregate(mcols(grgtfLike),by=list(origin$subjectHits),function(x){paste(unique(x),collapse=",")})
  simplified<-as.data.frame(grsimplified)
  colnames(simplified)[colnames(simplified)=="seqnames"]<-"seqid"
  return(simplified[,colnames(gtfLike)])
}

#####################################################################################

read.coverage.from.file <- function(filePath, chrToPlot, nbReads=NA, normalize=FALSE){
  cat("Reading coverage for",filePath,chrToPlot,".\n")
  if(!isWindows){
    shortFile<-paste0(filePath,chrToPlot,".gz")
    if (file.exists(shortFile)) {
      shortExists<-T
    } else {
      shortExists<-F
      cat("A file with only",chrToPlot,"coverage will be created at",shortFile,".\n")
    }
    if (file_ext(filePath)=="gz"){
      commandToOpen<-"gunzip -c "
    } else {
      commandToOpen<-"cat "
    }
    if(!shortExists){
      system(paste0(commandToOpen,filePath," | grep ^",chrToPlot,"$'\t' | gzip > ",shortFile))
      coverage<-tryCatch(read.delim(gzfile(shortFile),h=F),error = function(e){NA})
      if(all(is.na(coverage))){
        cat("There is no coverage for",chrToPlot,".\nWill try with",gsub("chr","",chrToPlot),".\n")
        system(paste0(commandToOpen,filePath," | grep ^",gsub("chr","",chrToPlot),"$'\t' | awk -v OFS='\t' '{print \"chr\"$0}' | gzip > ",shortFile))
        coverage<-tryCatch(read.delim(gzfile(shortFile),h=F),error = function(e){NA})
      }
    } else {
      coverage<-tryCatch(read.delim(gzfile(shortFile),h=F),error = function(e){NA})
    }  
  } else {
    if (file_ext(filePath)=="gz"){
      firstChars<-substr(readLines(gzfile(filePath),n=1),1,5)
      if(firstChars=="track"){
        coverage<-read.table(gzfile(filePath),sep="\t",stringsAsFactors=F,comment.char="#",skip=1)
      } else {
        coverage<-read.table(gzfile(filePath),sep="\t",stringsAsFactors=F,comment.char="#")
      }
    } else {
      firstChars<-substr(readLines(filePath,n=1),1,5)
      if(firstChars=="track"){
        coverage<-read.table(filePath,sep="\t",stringsAsFactors=F,comment.char="#",skip=1)
      } else {
        coverage<-read.table(filePath,sep="\t",stringsAsFactors=F,comment.char="#")
      }
    }
    tempCov<-subset(coverage,coverage[,1]==chrToPlot)
    if(nrow(tempCov)==0){
      cat("There is no coverage for",chrToPlot,".\nWill try with",gsub("chr","",chrToPlot),".\n")
      tempCov<-subset(coverage,coverage[,1]==gsub("chr","",chrToPlot))
      if(nrow(tempCov)>0){
        tempCov[,1]<-paste0("chr",tempCov[,1])
      } else {
        tempCov<-NA #To be compatible with the non windows version
      }
    }
    coverage<-tempCov
  }  
  if(normalize==TRUE && ! all(is.na(coverage)) && !is.na(nbReads)){
    cat("Normalizing read coverage with respect to the number of mapped reads.\n")
    M=nbReads/1e6
    cat("Dividing by M =",M, "million mapped reads.\n")
    coverage[,4]=coverage[,4]/M
  }
  return(coverage)
}

#####################################################################################

read.gtf.from.file <- function(filePath, chrToPlot){
  cat("Reading gtf for",filePath,chrToPlot,".\n")
  if(!isWindows){
    shortFile<-paste0(filePath,chrToPlot,".gz")
    if (file.exists(shortFile)) {
      shortExists<-T
    } else {
      shortExists<-F
      cat("A file with only",chrToPlot,"genes and only exons will be created at",shortFile,".\n")
    }
    if (file_ext(filePath)=="gz"){
      commandToOpen<-"gunzip -c "
    } else {
      commandToOpen<-"cat "
    }
    if(!shortExists){
      system(paste0(commandToOpen,filePath," | awk -v OFS='\t' '$1==\"",chrToPlot,"\" && $3==\"exon\"{print}' | gzip > ",shortFile))
      gtf<-readGFF(shortFile)
      if(nrow(gtf)==0){
        cat("There is no genes for",chrToPlot,".\nWill try with",gsub("chr","",chrToPlot),".\n")
        system(paste0(commandToOpen,filePath," | awk -v OFS='\t' '$1==\"",gsub("chr","",chrToPlot),"\" && $3==\"exon\"{print \"chr\"$0}' | gzip > ",shortFile))
        gtf<-readGFF(shortFile)
      }
    } else {
      gtf<-readGFF(shortFile)
    }
  } else {
    gtf<-readGFF(filePath)
    tempGtf<-subset(gtf,gtf$seqid==chrToPlot)
    if(nrow(tempGtf)==0){
      cat("There is no genes for",chrToPlot,".\nWill try with",gsub("chr","",chrToPlot),".\n")
      tempGtf<-subset(gtf,seqid==gsub("chr","",chrToPlot))
      if(nrow(tempGtf)>0){
        tempGtf$seqid<-paste0("chr",tempGtf$seqid)
      } else {
        tempGtf<-NA #To be compatible with the non windows version
      }
    }
    gtf<-tempGtf  
  }
  return(gtf)
}

#####################################################################################

read.bed.from.file <- function(filePath){
  if (file_ext(filePath)=="gz"){
    firstChars<-substr(readLines(gzfile(filePath),n=1),1,5)
    if(firstChars=="track"){
      annotationData<-read.table(gzfile(filePath),sep="\t",stringsAsFactors=F,comment.char="#",skip=1)
    } else {
      annotationData<-read.table(gzfile(filePath),sep="\t",stringsAsFactors=F,comment.char="#")
    }
  } else {
    firstChars<-substr(readLines(filePath,n=1),1,5)
    if(firstChars=="track"){
      annotationData<-read.table(filePath,sep="\t",stringsAsFactors=F,comment.char="#",skip=1)
    } else {
      annotationData<-read.table(filePath,sep="\t",stringsAsFactors=F,comment.char="#")
    }
  }
  return(annotationData)
}

#####################################################################################

plot.annotations<-function(annotations, chrToPlot, realStart, realEnd, col="black", cex.label=1, annotToExclude=NA,leftMargin=0.3,rightMargin=0.02){
  #Plot annotations which are already in the good genome
  #Normally this is already done.
  transformed.annot=annotations[which(annotations[,1]%in%c(chrToPlot)),]
  
  if(!all(is.na(annotToExclude))){
    transformed.annot=transformed.annot[which(!transformed.annot[,4]%in%annotToExclude),]
  }
  
  par(mai=c(0.05,leftMargin,0.05,rightMargin))

  xlim=c(realStart,realEnd)
  plot(1,type="n",xlab="",ylab="",axes=F,xlim=xlim,ylim=c(0.2,0.8),xaxs="i")
  
  height=0.06
  small=0.07
  
  ypos=0.4
  
  abline(h=ypos+height/2,lty=1, lwd=0.5, col="gray40")
  
  for(i in 1:nrow(transformed.annot)){
    thisstart=transformed.annot[i,2]
    thisend=transformed.annot[i,3]
    name=transformed.annot[i,4]
    #Different from Anouk (I put border=NA)
    rect(thisstart,ypos,thisend,ypos+height,col=col,border=NA)
    text(x=(thisstart+thisend)/2, y=ypos+2*height, labels=name, adj=c(0.5,1), cex=cex.label)
  }
}

#####################################################################################

plot.genes <- function(gtf, chrToPlot, realStart, realEnd, biotypes=NA, plot.axis=FALSE,  excluded=NA, col.fwd="gray30", col.rev="gray50", cex.names=1, col.lacZ=NA, simplifyHoxGenesNames=FALSE,removeExtremeValuesOfXAxis=F,leftMargin=0.3,rightMargin=0.02){

  if(plot.axis){
    #Increase the bottom margin
    par(mai=c(0.5,leftMargin,0.05,rightMargin))
  } else{
    par(mai=c(0.05,leftMargin,0.05,rightMargin))
  }
  
  ylim=c(0.2,1.1)
  
  #Restrict the gtf to the start and end
  exonsToPlot<-subset(gtf,seqid==chrToPlot & start<=realEnd & end>=realStart)
  
  if(!all(is.na(biotypes))){
    exonsToPlot<-subset(exonsToPlot,gene_biotype%in%biotypes)
  }    
  
  if(!all(is.na(excluded))){
    exonsToPlot<-exonsToPlot[!exonsToPlot$gene_id%in%excluded,]
  }  
  
  exonsToPlot<-exonsToPlot[order(exonsToPlot$start),]
  

  xlim=c(realStart, realEnd)
  plot(1,type="n",xlab="",ylab="",axes=F,xlim=xlim,ylim=ylim,xaxs="i")

  genes=unique(exonsToPlot$gene_id)
 
  height=0.06
  #Not like Anouk
  #small=0.07
  lwdForLine=0.7
  
  #tinyx=(realEnd-realStart)/1000
  
  ypos=c(0.25,0.75)
  names(ypos)=c("-","+")

  col=c(col.fwd, col.rev)
  names(col)=c("+", "-")

  abline(h=ypos+height/2,lty=1, lwd=0.5, col="gray40")

  if(plot.axis){
    xaxis=pretty(c(realStart,realEnd))
    if(removeExtremeValuesOfXAxis){
      xaxis<-xaxis[xaxis<=realEnd & xaxis>=realStart]
      xaxis<-xaxis[2:(length(xaxis)-1)]
    }
    xlabels=paste(round(xaxis/1000, digits=0),"kb")
    axis(side=1, at=xaxis, labels=xlabels, cex.axis=0.6)
  }
  
  arrow.width=diff(xlim)/50
  
  for(g in genes){
    this.coords=subset(exonsToPlot,gene_id==g)
    name=this.coords$gene_name[1]
    #I simplify the exons:
    #this.coords<-simplifyExonsKeepingInfo(this.coords)
    #Finally I do not keep the other info:
    this.coords<-as.data.frame(reduce(makeGRangesFromDataFrame(this.coords)))
    
    nbex=nrow(this.coords)
    
    if(simplifyHoxGenesNames){
      name<-gsub("^Hox","",name)
    }

    minx=min(this.coords$start)
    maxx=max(this.coords$end)
    
    for(i in 1:nbex){
      #Draw a rectangle for each exon
      thisstart=this.coords$start[i]
      thisend=this.coords$end[i]
      strand=this.coords$strand[i]
      this.ypos=ypos[as.character(strand)]
      this.col=col[as.character(strand)]
      #Different from Anouk border=NA
      #if(this.coords[i,"exon_id"]=="lacZ" & !is.na(col.lacZ)){
      #  rect(thisstart,this.ypos,thisend,this.ypos+height,col=col.lacZ,border=NA)
      #} else {
        rect(thisstart,this.ypos,thisend,this.ypos+height,col=this.col,border=NA)
      #}
      thisy=this.ypos+height
      ymid=thisy+height/2
      
      if(i>1){
        #Link with the previous exon
        #If they were on the same strand:
        if(this.coords$strand[i]==this.coords$strand[i-1]){
          prevstart=this.coords$start[i-1]
          prevend=this.coords$end[i-1]
          
          midpoint=mean(c(prevend,thisstart))
          segments(prevend,thisy,midpoint,ymid,col=this.col,lwd=lwdForLine)
          segments(midpoint,ymid,thisstart,thisy,col=this.col,lwd=lwdForLine)
        }
      }
    }
    
    strand=this.coords$strand[1]
    this.ypos=ypos[as.character(strand)]
    this.col=col[as.character(strand)]
    if(strand=="+"){
      # #Put the name at the left of the gene
      # text(x=minx-tinyx,y=this.ypos+4*small,labels=name,cex=0.6*cex.names, adj=0, font=3) ## names of the genes aligned left
      # #Put the arrow
      # segments(minx,this.ypos+small, minx, this.ypos+height+small*1.3, col=this.col)
      # arrows(minx, this.ypos+height+small*1.3, minx+arrow.width, this.ypos+height+small*1.3, length=0.05, col=this.col)
      #New version
      #Put the name at the left of the gene
      text(x=minx,y=this.ypos+3*height,labels=name,cex=0.6*cex.names, adj=c(0,0), font=3) ## names of the genes aligned left and bottom
      #Put the arrow
      segments(minx,this.ypos+height, minx, this.ypos+2*height, lwd=lwdForLine, col=this.col)
      arrows(minx, this.ypos+2*height, minx+arrow.width, this.ypos+2*height, length=0.05, col=this.col,lwd=lwdForLine)
    }

    if(strand=="-"){
      # #Put the name at the right
      # text(x=maxx-tinyx,y=this.ypos+4*small,labels=name,cex=0.6*cex.names, adj=1, font=3) ## names of the genes aligned right
      # #Put the arrow
      # segments(maxx,this.ypos+small, maxx, this.ypos+height+small*1.3, col=this.col)
      # arrows(maxx, this.ypos+height+small*1.3, maxx-arrow.width, this.ypos+height+small*1.3, length=0.05, col=this.col)
      #New version
      #Put the name at the right
      text(x=maxx,y=this.ypos+3*height,labels=name,cex=0.6*cex.names, adj=c(1,0), font=3) ## names of the genes aligned right and bottom
      #Put the arrow
      segments(maxx,this.ypos+height, maxx, this.ypos+2*height, lwd=lwdForLine, col=this.col)
      arrows(maxx, this.ypos+2*height, maxx-arrow.width, this.ypos+2*height, length=0.05, col=this.col,lwd=lwdForLine)    }
  }
  #At the end plot the lacZ in a different color
  if(!is.na(col.lacZ)){
    lacZExon<-subset(gtf,exon_id=="lacZ")
    strand=lacZExon$strand[1]
    this.ypos=ypos[strand]
    thisstart=lacZExon$start[1]
    thisend=lacZExon$end[1]
    rect(thisstart,this.ypos,thisend,this.ypos+height,col=col.lacZ,border=NA)
  }
}
  
#####################################################################################
  
plot.coverage <- function(cov, chrToPlot, realStart, realEnd, plot.axis=TRUE, recompute.ylim=TRUE, ylim=NA, col="black", pretty.axis=TRUE, cex.axislab=1, axis.interval=NA, title="",removeExtremeValuesOfXAxis=F,leftMargin=0.3,rightMargin=0.02){

  if(all(is.na(cov))){
    nbcov<-0
  } else {
    #Restricted to the plot area
    cov=cov[which(cov[,1]==chrToPlot & cov[,2]<=realEnd & cov[,3]>=realStart),]
    #Remove the NA values
    cov<-cov[!is.na(cov[,4]),]
    #Remove the 0 values
    cov<-cov[cov[,4]!=0,]
    nbcov=nrow(cov)
  }
  
  par(mai=c(0.4,leftMargin,0.3,rightMargin))

  if(recompute.ylim) {
    if(nbcov>0){
      #There is something to plot
      ylim=c(0,max(cov[,4]))
    } else {
      ylim=c(0,1)
    }
  }
  
  plot(1,type="n",xlab="",ylab="",axes=F,xlim=c(realStart,realEnd), ylim=ylim, xaxs="i", main=title)
  ## manage Y axis

  if(recompute.ylim){
    axis(side=2, cex.axis=0.6*cex.axislab, mgp=c(3,0.5,0)) ## default axis
  } else{
    if(all(is.na(axis.interval))){
      axis(side=2, cex.axis=0.6*cex.axislab, mgp=c(3,0.5,0)) 
    } else{
      axis(side=2, cex.axis=0.6*cex.axislab, mgp=c(3,0.5,0), at=axis.interval)
    }
  }
  
  ## for the x axis 
  if(plot.axis){
    xaxis=pretty(c(realStart,realEnd))    
    if(removeExtremeValuesOfXAxis){
      xaxis<-xaxis[xaxis<=realEnd & xaxis>=realStart]
      xaxis<-xaxis[2:(length(xaxis)-1)]
    }
    xlabels=paste(round(xaxis/1000, digits=0),"kb")
    axis(side=1, at=xaxis, labels=xlabels, cex.axis=0.6*cex.axislab, mgp=c(3,0.35,0))
  }
  
  if(nbcov>0){
  ## plot the coverage
  #rect(cov$V2+1,0,cov$V3,cov$V4,col=col,border=NA or border=col,lwd=1 or 0.5)
  #rect(cov[,2]+1,0,cov[,3],cov[,4],col=col,border=NA,lwd=0.5)
    #I checked with UCSC and they are plotting from column2 to column3.
    rect(cov[,2],0,cov[,3],cov[,4],col=col,border=NA,lwd=0.5)
  }

}

#####################################################################################

check.trackOptions <- function(track){
  isPerfect<-T
  if(class(track)!="list"){
    cat("The track is not a list, check the configFile example.\n")
    return(NA)
  }
  if(!"file"%in%names(track)){
    cat("The track definition do not contain a file path. Check the configFile example.\n")
    return(NA)
  }
  if(!file.exists(track$file)){
    cat("The file provided in the track definition do not exists.\n")
    return(NA)
  }
  #There is a file
  if(!"normalize"%in%names(track)){
    cat("The track will not be normalized.\n")
    track$normalize<-F
    track$nb.reads<-NA
    isPerfect<-F
  } else if(!is.logical(track$normalize)) {
    cat("The normalize value should be T or F, here got ",track$normalize," the track will not be normalized.\n")
    track$normalize<-F
    track$nb.reads<-NA
    isPerfect<-F    
  } else if(track$normalize) {
    if(!is.numeric(track$nb.reads)) {
      cat("You asked for a normalization but you did not provided a correct nb of reads. Please provide a number in nb.reads\nThe track will not be normalized.\n")
      track$normalize<-F
      track$nb.reads<-NA
      isPerfect<-F
    }
  }
  #Now normalize and nb.reads are compatible and can be used in plot.coverage
  if(!"color"%in%names(track)){
    cat("No color have been specified. It will be black.\n")
    track$color<-"black"
    isPerfect<-F
  } else {
    if(!isValidColor(track$color)){
      cat("The color specified cannot be interpreted by R. Go to http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf or use rgb to create a compatible color.\nThe track will be black.\n")
      track$color<-"black"
      isPerfect<-F
    }
  }
  #The color is correct
  if(!"plot.axis.coordinates"%in%names(track)){
    cat("The x axis will not be plotted.\n")
    track$plot.axis.coordinates<-F
    isPerfect<-F
  } else if(!is.logical(track$plot.axis.coordinates)) {
    cat("The plot.axis.coordinates value should be T or F, here got ",track$plot.axis.coordinates," the x axis will not be plotted.\n")
    track$plot.axis.coordinates<-F
    isPerfect<-F
  }
  #The plot.axis.coordinates is correct
  if(!"coeffToLabAxes"%in%names(track)){
    track$coeffToLabAxes<-1
  } else {
    if(!is.numeric(track$coeffToLabAxes)){
      cat("The coefficient specified for the label of the axes is not a number. 1 will be used.\n")
      track$coeffToLabAxes<-1
      isPerfect<-F
    }
  }
  #The coeffToLabAxes is correct
  if(!"recompute.ylim"%in%names(track)){
    track$recompute.ylim<-F
    warningMessage<-"recompute.ylim was not specified"
  } else {
    if(!is.logical(track$recompute.ylim)){
      warningMessage<-paste("recompute.ylim was",track$recompute.ylim,"whereas T or F was expected")
      track$recompute.ylim<-F
    } else {
      warningMessage<-paste("recompute.ylim was set to F")
    }
  }
  #recompute.ylim is defined.
  if(!"ylim"%in%names(track)){
    if(track$recompute.ylim){
      track$ylim<-NA
    } else {
      #ylim is not defined and recompute.ylim is FALSE
      cat(warningMessage,"but ylim is not defined. The y limits will be computed automatically (recompute.ylim=T).\n")
      track$ylim<-NA
      track$recompute.ylim<-T
    }
  } else {
    if(!track$recompute.ylim){
      if(all(is.na(track$ylim))){
        #ylim is NA and recompute.ylim is FALSE
        cat(warningMessage,"but ylim is set as NA. The y limits will be computed automatically (recompute.ylim=T).\n")
        track$ylim<-NA
        track$recompute.ylim<-T 
        isPerfect<-F
      } else if(length(track$ylim)!=2) {
        cat(warningMessage,"but ylim not constituted of 2 values. It should have been defined like c(0,10). The y limits will be computed automatically (recompute.ylim=T).\n")
        track$ylim<-NA
        track$recompute.ylim<-T
        isPerfect<-F
      } else if(!all(is.numeric(track$ylim))){
        cat(warningMessage,"but ylim values are not numeric. It should have been defined like c(0,10). The y limits will be computed automatically (recompute.ylim=T).\n")
        track$ylim<-NA
        track$recompute.ylim<-T
        isPerfect<-F
      }
      #ylim is correctly defined.
    } else {
      #ylim is defined and recompute.ylim is set to T
      if(!all(is.na(track$ylim))){
        cat("ylim is defined but recompute.ylim is set to T. ylim is ignored.\n")
        track$ylim<-NA
        isPerfect<-F
      }
    }
  }
  #ylim and recompute.ylim are correctly defined.
  if(!"axis.interval"%in%names(track)){
    track$axis.interval<-NA
  } else {
    if(!all(is.na(track$axis.interval))){
      if(track$recompute.ylim){
        cat("axis.interval is defined but recompute.ylim is set to T. axis.interval values will not be used.\n")
        track$axis.interval<-NA
        isPerfect<-F
      } else {
        #There are ylim defined.
        if(!all(is.numeric(track$axis.interval))){
          cat("axis.interval are not numbers. They will be ignored.\n")
          track$axis.interval<-NA
          isPerfect<-F          
        }
      }
    }
  }
  #axis.interval is correctly defined.
  if(isPerfect){
    cat("The track was perfectly defined.\n")
  }
  return(track)
}

#####################################################################################

check.genesOptions <- function(){
  genesOptions<-list()
  if(!exists("biotypesToKeep")){
    genesOptions$biotypesToKeep<-NA
  } else {
    if(all(is.na(biotypesToKeep))){
      genesOptions$biotypesToKeep<-NA
    } else if(length(unique(genesToPlot$gene_id[genesToPlot$gene_biotype%in%biotypesToKeep]))==0){
      cat("There is no gene in the biotype specified. Do not use this filter.\n")
      genesOptions$biotypesToKeep<-NA
    } else {
      genesOptions$biotypesToKeep<-biotypesToKeep
    }
  }
  if(!exists("genesToExclude")){
    genesOptions$genesToExclude<-NA
  } else {
    genesOptions$genesToExclude<-genesToExclude
  }
  if(!exists("col.genes.fwd")){
    cat("The color for fwd genes is not specified. They will be blue.\n")
    genesOptions$col.genes.fwd<-"blue"
  } else {
    if(!isValidColor(col.genes.fwd)){
      cat("The color specified for fwd genes (",col.genes.fwd,") cannot be interpreted by R. Go to http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf or use rgb to create a compatible color.\nThey will be blue.\n")
      genesOptions$col.genes.fwd<-"blue"
    } else {
      genesOptions$col.genes.fwd<-col.genes.fwd
    }
  }
  if(!exists("col.genes.rev")){
    cat("The color for rev genes is not specified. They will be red.\n")
    genesOptions$col.genes.rev<-"red"
  } else {
    if(!isValidColor(col.genes.rev)){
      cat("The color specified for rev genes cannot be interpreted by R. Go to http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf or use rgb to create a compatible color.\nThey will be red.\n")
      genesOptions$col.genes.rev<-"red"
    } else {
      genesOptions$col.genes.rev<-col.genes.rev
    }
  }
  if(!exists("plot.axis.coordinates.Genes")){
    cat("The x axis will not be plotted.\n")
    genesOptions$plot.axis.coordinates.Genes<-F
  } else if(!is.logical(plot.axis.coordinates.Genes)) {
    cat("The plot.axis.coordinates.Genes value should be T or F, here got ",plot.axis.coordinates.Genes," the x axis will not be plotted.\n")
    genesOptions$plot.axis.coordinates.Genes<-F
  } else {
    genesOptions$plot.axis.coordinates.Genes<-plot.axis.coordinates.Genes
  }
  if(!exists("coeffToLabGenes")){
    genesOptions$coeffToLabGenes<-1
  } else {
    if(!is.numeric(coeffToLabGenes)){
      cat("The coefficient specified for the label of the genes is not a number. 1 will be used.\n")
      genesOptions$coeffToLabGenes<-1
    } else {
      genesOptions$coeffToLabGenes<-coeffToLabGenes
    }
  }
  if(!exists("col.lacZ")){
    genesOptions$col.lacZ<-NA
  } else if (is.na(col.lacZ)){
    genesOptions$col.lacZ<-NA
  } else if (!isValidColor(col.lacZ)){
    cat("The color specified for lacZ cannot be interpreted by R. Go to http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf or use rgb to create a compatible color.\nIt will not be different from other exons.\n")
    genesOptions$col.lacZ<-NA
  } else {
    genesOptions$col.lacZ<-col.lacZ
  }
  if(!exists("simplifyHoxGenesNames")){
    genesOptions$simplifyHoxGenesNames<-F
  } else if(!is.logical(simplifyHoxGenesNames)) {
    cat("The simplifyHoxGenesNames value should be T or F, here got ",simplifyHoxGenesNames," the Hox genes will not be simplified.\n")
    genesOptions$simplifyHoxGenesNames<-F
  } else {
    genesOptions$simplifyHoxGenesNames<-simplifyHoxGenesNames
  }
  return(genesOptions)
}

#####################################################################################

check.annotOptions <- function(){
  annotOptions<-list()
  if(!exists("col.annotations")){
    cat("The color for annotations is not specified. They will be black.\n")
    annotOptions$col.annotations<-"black"
  } else {
    if(!isValidColor(col.annotations)){
      cat("The color specified for annotations cannot be interpreted by R. Go to http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf or use rgb to create a compatible color.\nThey will be blue.\n")
      annotOptions$col.annotations<-"black"
    } else {
      annotOptions$col.annotations<-col.annotations
    }
  }
  if(!exists("coeffToLabAnnot")){
    annotOptions$coeffToLabAnnot<-1
  } else {
    if(!is.numeric(coeffToLabAnnot)){
      cat("The coefficient specified for the label of the annotations is not a number. 1 will be used.\n")
      annotOptions$coeffToLabAnnot<-1
    } else {
      annotOptions$coeffToLabAnnot<-coeffToLabAnnot
    }
  }
  if(!exists("annotToExclude")){
    annotOptions$annotToExclude<-NA
  } else {
    annotOptions$annotToExclude<-annotToExclude
  }
  return(annotOptions)
}

#####################################################################################
#####################################################################################
