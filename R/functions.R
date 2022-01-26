#' Calculate SlideWindow mean 
#'
#' This function takes computes the mean over a sequence 
#'
#' @param X vector
#' @param N length
#' @return A matrix of the infile
#' @export
#' 
SlideWindowMean <- function(X, N)
{
  sapply(1:length(X), function(Z)
  {
    if((Z+N) <= length(X))
    {
      mean(X[Z:(Z+N)]) 
    }
    else
    {
      mean(X[Z:length(X)])
    }
  })
}



#' Read sample
#'
#' This function Read a intensity file
#'
#' @param Rawfile path to file
#' @param skip number of lines to skip
#' @param chr if limit to any specific chromosome 
#' @return A matrix of the infile
#' @export
#' 
ReadSample <- function(RawFile="Test.txt", skip=0, LCR=NULL, PFB=NULL, chr=NA, SNPList=NULL, start=NULL, stop=NULL)
{
  suppressPackageStartupMessages(library(data.table))
  Sample <- fread(RawFile, head=T, sep="\t", skip=skip, verbose=FALSE)
  Sample <- as.data.frame(Sample)
  colnames(Sample) <- gsub(" ", ".", colnames(Sample))
  colnames(Sample)[colnames(Sample) %in% "Name"] <- "SNP.Name"
  colnames(Sample)[colnames(Sample) %in% "Chromosome"] <- "Chr"
  colnames(Sample)[colnames(Sample) %in% "Allele1.-.Top"] <- "Allele1"
  colnames(Sample)[colnames(Sample) %in% "Allele2.-.Top"] <- "Allele2"
  colnames(Sample)[colnames(Sample) %in% "BAF"] <- "B.Allele.Freq"
  colnames(Sample)[colnames(Sample) %in% "LRR"] <- "Log.R.Ratio"
  colnames(Sample)[grep("Log.R.Ratio", colnames(Sample))] <- "Log.R.Ratio" # Remove text from PennCNV format.
  colnames(Sample)[grep("B.Allele.Freq", colnames(Sample))] <- "B.Allele.Freq"
  colnames(Sample)[grep("Chr", colnames(Sample))] <- "Chr"
  colnames(Sample)[grep("Position", colnames(Sample))] <- "Position"
  
  if(sum(colnames(Sample) %in% c("Log.R.Ratio", "B.Allele.Freq", "Chr")) == 0)
  {
    stop("\n ERROR: Could not find header.\nPlease check if you need to skip some lines using skip variable: default = 10.\nExpected header:SNP.Name\tChr\tPosition\tB.Allele.Freq\tLog.R.Ratio.\n")
  }
  
  # Windows problem
  Sample$B.Allele.Freq <- as.numeric(gsub("\\r(?!\\n)","", Sample$B.Allele.Freq, perl=T))
  Sample$Log.R.Ratio <- as.numeric(gsub("\\r(?!\\n)","", Sample$Log.R.Ratio, perl=T))
  #CNV <- CNV[,c("SNP.Name","Chr", "Position", "Log.R.Ratio", "B.Allele.Freq", "Allele1", "Allele2")] # SNP.Name
  
  
  # Genotype together (deCODE)
  #if(!is.null(Sample$Genotype))
  #{
  #	Allele <- sapply(Sample$Genotype, function(X){ data.frame(Allele1= unlist(strsplit(X, ""))[1], Allele2=unlist(strsplit(X, ""))[2], stringsAsFactors=F) })
  #	Allele <- MatrixOrList2df(Allele)
  #	Sample$Allele1 <- Allele$Allele1
  #	Sample$Allele2 <- Allele$Allele2
  #}
  
  # SNP-positionFile
  if(!is.null(SNPList)){
    suppressPackageStartupMessages(library(data.table))
    SNPlist <- fread(SNPList, head=T, skip=skip, verbose=FALSE, stringsAsFactors=FALSE, colClasses =c("Chr"="factor"))
    SNPlist <- as.data.frame(SNPlist)
    colnames(SNPlist)[colnames(SNPlist) %in% "Name"] <- "SNP.Name"
    colnames(SNPlist)[colnames(SNPlist) %in% "Chromosome"] <- "Chr" 
    colnames(Sample)[colnames(Sample) %in% "BAF"] <- "B.Allele.Freq"
    colnames(Sample)[colnames(Sample) %in% "LRR"] <- "Log.R.Ratio"
    Sample$Chr <- NULL # remove previous values of Chr if present
    Sample$Position <- NULL # remove previous values of Position if present
    total <- merge(Sample, SNPlist,by="SNP.Name")
    Sample <- total    
  }
  
  if(sum(colnames(Sample) %in% c("Log.R.Ratio", "B.Allele.Freq", "Chr")) == 0)
  {
    stop("\n ERROR: Could not find header.\nPlease check if you need to skip some lines using skip variable: default = 10.\nExpected header:SNP.Name\tChr\tPosition\tB.Allele.Freq\tLog.R.Ratio.\n")
  }
  
  # removing chr from chromosome name (deCODE)
  Sample$Chr <- gsub("chr", "", Sample$Chr)
  
  
  # PFB
  if(is.null(PFB)){ Sample$PFB <- rep(0.5, nrow(Sample)) }else{ Sample$PFB <- PFB }
  
  # Subsetting 
  Sample <- subset(Sample, !Chr %in% c("XY", "0")) # "MT", "X", "Y",
  Sample <- subset(Sample, !is.na(Sample$B.Allele.Freq)) # Removes SNPs without BAF-value 
  Sample <- subset(Sample, !is.na(Sample$Log.R.Ratio)) #  Removes SNPs without LRR-value 
  
  # chr specific. Example chr="22"
  if(!is.na(chr)){ Sample <- subset(Sample, Chr %in% chr) }
  
  if(!is.null(start)){ Sample <- subset(Sample, Position >= start & Position <= stop) }
  
  if(!is.null(LCR))
  {
    Sample <- subset(Sample, !SNP.Name %in% LCR) 
  }
  
  Sample$LRR <- Sample$Log.R.Ratio # CNV$LRR is the
  
  Sample <- Sample[order(Sample$Chr, Sample$Position),]
  
  return(Sample)
}



#'Compute LRRs and BAF values
#'
#'  Given the path of a sample  computes  mean and sd for LRRs and  heterozigosity for BAF values per chromosome 
#'  
#' @import data.table
#' @return data frame with all computed values
#' @export 
#' 
computeValues<-function(path){
  print(path)
  
  start <- Sys.time()
  Sample<-ReadSample(RawFile = path,skip=10)
  Sample.ID<-Sample$Sample.ID[1]
  meanS<-setDT(Sample)[ , .(mean_LRR = mean(LRR)), by = Chr]
  sdS<-setDT(Sample)[ , .(mean_LRR = sd(LRR)), by = Chr]
  BAF<-setDT(Sample)[,  .(100*(sum(B.Allele.Freq > 0.4 & B.Allele.Freq < 0.6)/length(B.Allele.Freq))), by = Chr]
  
  meanS2<-transpose(meanS)
  meanS3<-meanS2[2,]
  colnames(meanS3)<-paste0('m',as.character(meanS2[1,]))
  
  
  sdS2<-transpose(sdS)
  sdS3<-sdS2[2,]
  colnames(sdS3)<-paste0('s',as.character(sdS2[1,]))
  head(sdS3)
  
  BAF2<-transpose(BAF)
  BAF3<-BAF2[2,]
  colnames(BAF3)<-paste0('baf',as.character(BAF2[1,]))
  head(BAF3)
  r0<-cbind(data.frame(Path=path,Sample.ID=Sample.ID),meanS3,sdS3,BAF3)
  end <- Sys.time()
  t5<-as.numeric (end - start, units = "mins") # or secs ..
  print(paste0('Computing time ',t5, 'mins'))
  return(r0)
}



#' Run in Parallele  computeValues
#' 
#' Applies the computeValues function to run in multiple cores
#' 
#' @import foreach
#' @import data.table
#' @import parallel
#' @import doParallel
#' @import R.utils
#' @return data frame with all computed values
#' @export 
#' 
runPar<-function(paths,cores=10){
  
  
  
  start <- Sys.time()
  cl <- parallel::makeCluster(cores[1]-1) #not to overload your computer
  doParallel::registerDoParallel(cl)
  computeValues<-function(path){
    print(path)
    
    Sample<-ReadSample(RawFile = path,skip=10)
    Sample.ID<-Sample$Sample.ID[1]
    meanS<-setDT(Sample)[ , .(mean_LRR = mean(LRR)), by = Chr]
    sdS<-setDT(Sample)[ , .(mean_LRR = sd(LRR)), by = Chr]
    BAF<-setDT(Sample)[,  .(100*(sum(B.Allele.Freq > 0.4 & B.Allele.Freq < 0.6)/length(B.Allele.Freq))), by = Chr]
    
    meanS2<-transpose(meanS)
    meanS3<-meanS2[2,]
    colnames(meanS3)<-paste0('m',as.character(meanS2[1,]))
    
    
    sdS2<-transpose(sdS)
    sdS3<-sdS2[2,]
    colnames(sdS3)<-paste0('s',as.character(sdS2[1,]))
    head(sdS3)
    
    BAF2<-transpose(BAF)
    BAF3<-BAF2[2,]
    colnames(BAF3)<-paste0('baf',as.character(BAF2[1,]))
    head(BAF3)
    r0<-cbind(data.frame(Path=path,Sample.ID=Sample.ID),meanS3,sdS3,BAF3)
    return(r0)
  }
  
  ReadSample <- function(RawFile="Test.txt", skip=0, LCR=NULL, PFB=NULL, chr=NA, SNPList=NULL, start=NULL, stop=NULL)
  {
    suppressPackageStartupMessages(library(data.table))
    Sample <- fread(RawFile, head=T, sep="\t", skip=skip, verbose=FALSE)
    Sample <- as.data.frame(Sample)
    colnames(Sample) <- gsub(" ", ".", colnames(Sample))
    colnames(Sample)[colnames(Sample) %in% "Name"] <- "SNP.Name"
    colnames(Sample)[colnames(Sample) %in% "Chromosome"] <- "Chr"
    colnames(Sample)[colnames(Sample) %in% "Allele1.-.Top"] <- "Allele1"
    colnames(Sample)[colnames(Sample) %in% "Allele2.-.Top"] <- "Allele2"
    colnames(Sample)[colnames(Sample) %in% "BAF"] <- "B.Allele.Freq"
    colnames(Sample)[colnames(Sample) %in% "LRR"] <- "Log.R.Ratio"
    colnames(Sample)[grep("Log.R.Ratio", colnames(Sample))] <- "Log.R.Ratio" # Remove text from PennCNV format.
    colnames(Sample)[grep("B.Allele.Freq", colnames(Sample))] <- "B.Allele.Freq"
    colnames(Sample)[grep("Chr", colnames(Sample))] <- "Chr"
    colnames(Sample)[grep("Position", colnames(Sample))] <- "Position"
    
    if(sum(colnames(Sample) %in% c("Log.R.Ratio", "B.Allele.Freq", "Chr")) == 0)
    {
      stop("\n ERROR: Could not find header.\nPlease check if you need to skip some lines using skip variable: default = 10.\nExpected header:SNP.Name\tChr\tPosition\tB.Allele.Freq\tLog.R.Ratio.\n")
    }
    
    # Windows problem
    Sample$B.Allele.Freq <- as.numeric(gsub("\\r(?!\\n)","", Sample$B.Allele.Freq, perl=T))
    Sample$Log.R.Ratio <- as.numeric(gsub("\\r(?!\\n)","", Sample$Log.R.Ratio, perl=T))
    #CNV <- CNV[,c("SNP.Name","Chr", "Position", "Log.R.Ratio", "B.Allele.Freq", "Allele1", "Allele2")] # SNP.Name
    
    
    # Genotype together (deCODE)
    #if(!is.null(Sample$Genotype))
    #{
    #	Allele <- sapply(Sample$Genotype, function(X){ data.frame(Allele1= unlist(strsplit(X, ""))[1], Allele2=unlist(strsplit(X, ""))[2], stringsAsFactors=F) })
    #	Allele <- MatrixOrList2df(Allele)
    #	Sample$Allele1 <- Allele$Allele1
    #	Sample$Allele2 <- Allele$Allele2
    #}
    
    # SNP-positionFile
    if(!is.null(SNPList)){
      suppressPackageStartupMessages(library(data.table))
      SNPlist <- fread(SNPList, head=T, skip=skip, verbose=FALSE, stringsAsFactors=FALSE, colClasses =c("Chr"="factor"))
      SNPlist <- as.data.frame(SNPlist)
      colnames(SNPlist)[colnames(SNPlist) %in% "Name"] <- "SNP.Name"
      colnames(SNPlist)[colnames(SNPlist) %in% "Chromosome"] <- "Chr" 
      colnames(Sample)[colnames(Sample) %in% "BAF"] <- "B.Allele.Freq"
      colnames(Sample)[colnames(Sample) %in% "LRR"] <- "Log.R.Ratio"
      Sample$Chr <- NULL # remove previous values of Chr if present
      Sample$Position <- NULL # remove previous values of Position if present
      total <- merge(Sample, SNPlist,by="SNP.Name")
      Sample <- total    
    }
    
    if(sum(colnames(Sample) %in% c("Log.R.Ratio", "B.Allele.Freq", "Chr")) == 0)
    {
      stop("\n ERROR: Could not find header.\nPlease check if you need to skip some lines using skip variable: default = 10.\nExpected header:SNP.Name\tChr\tPosition\tB.Allele.Freq\tLog.R.Ratio.\n")
    }
    
    # removing chr from chromosome name (deCODE)
    Sample$Chr <- gsub("chr", "", Sample$Chr)
    
    
    # PFB
    if(is.null(PFB)){ Sample$PFB <- rep(0.5, nrow(Sample)) }else{ Sample$PFB <- PFB }
    
    # Subsetting 
    Sample <- subset(Sample, !Chr %in% c("XY", "0")) # "MT", "X", "Y",
    Sample <- subset(Sample, !is.na(Sample$B.Allele.Freq)) # Removes SNPs without BAF-value 
    Sample <- subset(Sample, !is.na(Sample$Log.R.Ratio)) #  Removes SNPs without LRR-value 
    
    # chr specific. Example chr="22"
    if(!is.na(chr)){ Sample <- subset(Sample, Chr %in% chr) }
    
    if(!is.null(start)){ Sample <- subset(Sample, Position >= start & Position <= stop) }
    
    if(!is.null(LCR))
    {
      Sample <- subset(Sample, !SNP.Name %in% LCR) 
    }
    
    Sample$LRR <- Sample$Log.R.Ratio # CNV$LRR is the
    
    Sample <- Sample[order(Sample$Chr, Sample$Position),]
    
    return(Sample)
  }
  
  
  finalMatrix <- foreach::foreach(i=paths, .combine=rbind) %dopar% {
    tempMatrix = computeValues(i) #calling a function
    #do other things if you want
    tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
  }
  
  end <- Sys.time()
  t5<-as.numeric (end - start, units = "mins") # or secs ..
  print(paste0('Computing time ',t5, 'mins'))
  return(finalMatrix)
}



#' Create and optimize Clusters for Males
#' 
#' Given a data frame with mean LRR and BAF values computes SCA clusters for Males
#' 
#' @param dataFrame data frame with computed LRR & BAF values
#' @param minPtsRange vector with range of values to optimize the minPts parameter 
#' @import dbscan
#' @return data frame with all computed values
#' @export 
#' 
optimizeClusterMales<- function(dataFrame,minPtsRange=c(5,8,10,12,15,20,25,30,35,40,50,100),interClusProb=0.5){
  
  cls<-c()
  for (x in minPtsRange){
    print(paste0(as.character(length(minPtsRange)-which(minPtsRange==x)),' Optimizations left...'))
    cl1 <- hdbscan(dataFrame, minPts = x)
    cls<-c(cls,list(cl1))
  }
  
  
  clsQC1<-c()
  clsQC2<-data.frame()
  
  for (x in 1:length(minPtsRange)){
    if (sum(c('1','2','3') %in% names(table(cls[[x]]$cluster)))==3 & length(table(cls[[x]]$cluster))<5){ 
      clsQC1<-c(clsQC1,list(cls[[x]]))
      clsQC2<-rbind(clsQC2,data.frame(minPts=minPtsRange[x],prob=mean(cls[[x]]$membership_prob)))
    }
    best <- clsQC2[which.max(clsQC2$prob),]
    bestPos<-which(minPtsRange==best$minPts)}
  
  f1<-cls[[bestPos]]
  dataFrame$cluster<-f1$cluster
  dataFrame$prob<-f1$membership_prob
  sample2 <- subset(dataFrame,dataFrame$cluster!=0)
  
  chrX <- data.frame(aggregate(sample2[,'mX'], list(sample2$cluster), mean))
  names(chrX)<-c('group','mX')
  XXY<-as.integer(chrX$group[which.max(chrX$mX)])
  
  chrY <- data.frame(aggregate(sample2[,'mY'], list(sample2$cluster), mean))
  names(chrY)<-c('group','mY')
  XYY<-as.integer(chrY$group[which.max(chrY$mY)])
  
  clusters<-c(1,2,3)
  clusters
  clusters <- clusters[-which(clusters ==XXY)]
  clusters <- clusters[-which(clusters ==XYY)]
  XY<-clusters
  dataFrame$clusKario<-0
  dataFrame$clusKario[dataFrame$cluster==XY]<-1
  dataFrame$clusKario[dataFrame$cluster==XXY]<-2
  dataFrame$clusKario[dataFrame$cluster==XYY]<-3
  
  
  
  ## Compute samples between clusters
  dataFrame$clusKario2<-dataFrame$clusKario
  dataFrame$id<-1:dim(dataFrame)[1]
  sample2<-subset(dataFrame,dataFrame$clusKario==1) # subset XY Cluster
  
  lm3<-(lm(mY ~ mX, data = sample2)) # Compute regression line 
  b<-as.numeric(lm3$coefficients[1])
  mx<-as.numeric(lm3$coefficients[2])
  
  # Calculate samples over & under the line
  predY<-(sample2$mX*mx)+(b)
  dev<-sample2$mY-predY
  sample2$xyyDev<-dev
  
  #XYY
  mydataInterXYY<-subset(sample2,sample2$prob<interClusProb & sample2$xyyDev>0)
  #mydataInterXYY_indx<-which(sample2$prob<interClusProb & sample2$xyyDev>0)
  #dim(mydataInterXYY)
  #XXY
  mydataInterXXY<-subset(sample2,sample2$prob<interClusProb & sample2$xyyDev<0)
  #mydataInterXXY_indx<-which(sample2$prob<interClusProb & sample2$xyyDev<0)
  
  #dim(mydataInterXXY)
  
  # Plot results
  #plot(sample1$mX,sample1$mY,pch=20,col=optimizeClusterMales$clusKario+1)
  #abline(c(b,mx))
  #points(mydataInterXYY$mX,mydataInterXYY$mY,col='yellow',pch=20)
  #points(mydataInterXXY$mX,mydataInterXXY$mY,col='yellow',pch=20)
  dataFrame$clusKario3<-dataFrame$clusKario2
  
  dataFrame$clusKario3[dataFrame$id %in% mydataInterXYY$id]<-(-1)
  dataFrame$clusKario3[dataFrame$id %in% mydataInterXXY$id]<-(-1)
  
  
  
  obj1<-list(f1,dataFrame,mydataInterXYY,mydataInterXXY)
  return(obj1)
  
}



#' Plot Male Clusters
#' 
#' Given the optimizeClusterMales output is plots the results
#' 
#' @param dataFrame data frame with computed LRR & BAF values
#' @param optimizeClusterMales best cluster. Output  from optimizeClusterMales method
#' @return Plot with male Clusters
#' @export 
#' 
plotMaleClusters<- function(dataFrame,optimizeClusterMales){
  
  #optimizeClusterMales$clusKarioCol[optimizeClusterMales$clusKario==0]<-'#adadad'
  #optimizeClusterMales$clusKarioCol[optimizeClusterMales$clusKario==1]<-'#4F76C3'
  #optimizeClusterMales$clusKarioCol[optimizeClusterMales$clusKario==2]<-'#D80913'
  #optimizeClusterMales$clusKarioCol[optimizeClusterMales$clusKario==3]<-'#37C27B'
  colors<-c('#E9D400','#C31200','#1C0087' ,'#009421','grey')

  optimizeClusterMales$clusKarioCol[optimizeClusterMales$clusKario3==0]<-colors[1]
  optimizeClusterMales$clusKarioCol[optimizeClusterMales$clusKario3==1]<-colors[3]  
  optimizeClusterMales$clusKarioCol[optimizeClusterMales$clusKario3==2]<- colors[2]
  optimizeClusterMales$clusKarioCol[optimizeClusterMales$clusKario3==3]<-colors[4] # XYY
  optimizeClusterMales$clusKarioCol[optimizeClusterMales$clusKario3==-1]<-colors[5] # XYY
  
  
  plot(dataFrame$mX,dataFrame$mY,pch=20,col=optimizeClusterMales$clusKarioCol, xlab="Chromosome X", ylab='Chromosome Y')
  rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4], col = "#FFFFFA") # Color  FFFDE2
  points(dataFrame$mX,dataFrame$mY,pch=20,col=optimizeClusterMales$clusKarioCol, xlab="Chromosome X", ylab='Chromosome Y')
  #points(f1[[3]]$mX,f1[[3]]$mY,col=colors[5],pch=20)
  #points(f1[[4]]$mX,f1[[4]]$mY,col=colors[5],pch=20)
  
}


#' Create and optimize Clusters for Females
#' 
#' Given a data frame with mean LRR and BAF values computes SCA clusters for Females
#' 
#' @param dataFrame data frame with computed LRR & BAF values
#' @param minPtsRange vector with range of values to optimize the minPts parameter 
#' @import dbscan
#' @return data frame with all computed values
#' @export 
#' 
optimizeClusterFemales<- function(dataFrame,minPtsRange=c(5,8,10,12,15,20,25,30,35,40,50,100),interClusProb=0.5){
  
  cls<-c()
  for (x in minPtsRange){
    print(paste0(as.character(length(minPtsRange)-which(minPtsRange==x)),' Optimizations left...'))
    cl1 <- hdbscan(dataFrame, minPts = x)
    cls<-c(cls,list(cl1))
  }
  
  clsQC1<-c()
  clsQC2<-data.frame()
  
  for (x in 1:length(minPtsRange)){
    if (sum(c('1','2','3') %in% names(table(cls[[x]]$cluster)))==3 & length(table(cls[[x]]$cluster))<5){ 
      clsQC1<-c(clsQC1,list(cls[[x]]))
      clsQC2<-rbind(clsQC2,data.frame(minPts=minPtsRange[x],prob=mean(cls[[x]]$membership_prob)))
    }
    best <- clsQC2[which.max(clsQC2$prob),]
    bestPos<-which(minPtsRange==best$minPts)}
  
  f1<-cls[[bestPos]]
  dataFrame$cluster<-f1$cluster
  dataFrame$prob<-f1$membership_prob
  sample2 <- subset(dataFrame,dataFrame$cluster!=0)
  
  chrX <- data.frame(aggregate(sample2[,'mX'], list(sample2$cluster), mean))
  names(chrX)<-c('group','mX')
  XXX<-as.integer(chrX$group[which.max(chrX$mX)])
  
  X<-as.integer(chrX$group[which.min(chrX$mX)])
  
  clusters<-c(1,2,3)
  clusters
  clusters <- clusters[-which(clusters ==XXX)]
  clusters <- clusters[-which(clusters ==X)]
  XX<-clusters
  dataFrame$clusKario<-0
  dataFrame$clusKario[dataFrame$cluster==XX]<-1
  dataFrame$clusKario[dataFrame$cluster==XXX]<-2
  dataFrame$clusKario[dataFrame$cluster==X]<-3
  
  
  
  ## Compute samples between clusters
  dataFrame$clusKario2<-dataFrame$clusKario
  dataFrame$id<-1:dim(dataFrame)[1]
  sample2<-subset(dataFrame,dataFrame$clusKario==1) # subset XX Cluster
  
  lm3<-(lm(mA ~ mX, data = sample2)) # Compute regression line 
  b<-as.numeric(lm3$coefficients[1])
  mx<-as.numeric(lm3$coefficients[2])
  
  # Calculate samples over & under the line
  predA<-(sample2$mX*mx)+(b)
  dev<-sample2$mA-predA
  sample2$xxDev<-dev
  
  #XXX
  mydataInterXXX<-subset(sample2,sample2$prob<interClusProb & sample2$xxDev>0)
  dim(mydataInterXXX)
  
  #X
  mydataInterX<-subset(sample2,sample2$prob<interClusProb & sample2$xxDev<0)
  
  
  dataFrame$clusKario3<-dataFrame$clusKario2
  
  
  dataFrame$clusKario3[dataFrame$id %in% mydataInterXXX$id]<-(-1)
  dataFrame$clusKario3[dataFrame$id %in%   mydataInterX$id]<-(-1)
  
  

  
  obj1<-list(f1,dataFrame,mydataInterXXX,mydataInterX)
  return(obj1)
  
}



#' Plot Female Clusters
#' 
#' Given the optimizeClusterFemales output is plots the results
#' 
#' @param dataFrame data frame with computed LRR & BAF values
#' @param optimizeClusterFemales best cluster. Output  from optimizeClusterFemales method
#' @return Plot with male Clusters
#' @export 
#' 
plotFemaleClusters<- function(dataFrame,optimizeClusterFemales){
  
  colors<-c('#E9D400','#C31200','#1C0087' ,'#009421','grey')
  optimizeClusterFemales$clusKarioCol[optimizeClusterFemales$clusKario3==0]<-colors[1]
  optimizeClusterFemales$clusKarioCol[optimizeClusterFemales$clusKario3==1]<-colors[3]
  optimizeClusterFemales$clusKarioCol[optimizeClusterFemales$clusKario3==2]<-colors[2]
  optimizeClusterFemales$clusKarioCol[optimizeClusterFemales$clusKario3==3]<-colors[4]
  optimizeClusterFemales$clusKarioCol[optimizeClusterFemales$clusKario3==-1]<-colors[5]
  
  
  plot(dataFrame$mX,dataFrame$mA,pch=20,col=optimizeClusterFemales$clusKarioCol, xlab="Chromosome X", ylab='Autosomal Chromosomes')
  rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4], col = "#FFFFFA") # Color  FFFDE2
  points(dataFrame$mX,dataFrame$mA,pch=20,col=optimizeClusterFemales$clusKarioCol)
  #points(f1[[3]]$mX,f1[[3]]$mA,col='#FFBC56',pch=20)
  #points(f1[[4]]$mX,f1[[4]]$mA,col='#FFBC56',pch=20)
  
  
}



#' Plot individual samples X and Y chromosomes
#' 
#' Given the optimizeClusterFemales output is plots the results
#' 
#' @param Sample data frame with computed LRR & BAF values
#' @param Name path/name where to save plot
#' @param Type Output file format one from: 'png' (default), 'pdf', 'svg'
#' @return Plot with X and Y chromosomes
#' @export 
#' 

PlotXandYChr <- function(Sample, Name="Test.png", Type='png',Save=FALSE)
{
  
  library(gridExtra)
  library(ggplot2)
  #library(ggbio)
  library(RColorBrewer)
  
  tmp <- subset(Sample, Chr %in% "X")
  tmp2 <- subset(Sample, Chr %in% "Y")
  
  tmp <- tmp[order(tmp$Position),]
  tmp2 <- tmp2[order(tmp2$Position),]
  
  # Fixing LRR 
  tmp$Log.R.Ratio[abs(tmp$Log.R.Ratio) > 10] <- 0
  tmp2$Log.R.Ratio[abs(tmp2$Log.R.Ratio) > 10] <- 0
  
  Colors <- c('#f93800','#283350')
  Mean <- SlideWindowMean(tmp$Log.R.Ratio, 50)
  tmp$Mean <- Mean
  
  Mean <- SlideWindowMean(tmp2$Log.R.Ratio, 25)
  tmp2$Mean <- Mean
  
  ## Chromosome X ##
  
  # Limits
  Min <- median(tmp$Log.R.Ratio) - (sd(tmp$Log.R.Ratio)*4)
  Max <- median(tmp$Log.R.Ratio) + (sd(tmp$Log.R.Ratio)*4)
  
  # LRR
  p1 <- ggplot(tmp, aes(Position, Log.R.Ratio)) + labs(y = "LRR", x=NULL) 
  p1 <- p1 + geom_point(alpha=0.05, size=1, colour=Colors[1]) + ylim(-1, 1)	
  p1 <- p1 + geom_line(data=tmp, aes(x=Position, y = Mean), size = 0.5, alpha=0.5, colour="black") 
  # BAF
  p2 <- ggplot(tmp, aes(Position, y = B.Allele.Freq)) + labs(y = "BAF") + labs(y = "BAF", x= NULL) + ggtitle('Chromosome X')
  p2 <- p2 + geom_point(aes(col="B.Allele.Freq"), alpha=0.3, size=1, colour=Colors[2])  
  
  
  ## Chromosome Y ##
  chr <- unique(tmp2$Chr)
  
  # Limits
  Min <- median(tmp2$Log.R.Ratio) - (sd(tmp2$Log.R.Ratio)*4)
  Max <- median(tmp2$Log.R.Ratio) + (sd(tmp2$Log.R.Ratio)*4)
  
  # LRR
  p5 <- ggplot(tmp2, aes(Position, Log.R.Ratio)) + labs(y = "LRR", x= NULL)
  p5 <- p5 + geom_point(alpha=0.05, size=1, colour=Colors[1]) + ylim(Min, Max)		
  p5 <- p5 + geom_line(data=tmp2, aes(x=Position, y = Mean), size = 0.5, alpha=0.5, colour="black") 
  
  # BAF
  p6 <- ggplot(tmp2, aes(Position, y = B.Allele.Freq))   + labs(y = "BAF", x= NULL)+ ggtitle('Chromosome Y') 
  p6 <- p6 + geom_point(aes(col="B.Allele.Freq"), alpha=0.5, size=1, colour=Colors[2]) 
  
  
  # Sample Name
  ID <- unique(Sample$Sample.ID)[1]
  Title <- paste("Sample: ", ID, sep="", collapse="")
  
  if (Save==TRUE & Type=='pdf'){ 
    print('pdf')
    pdf(paste0(Name,'.pdf'))
    grid.arrange(p2,p1,p6,p5,ncol=1,newpage = F)  
    dev.off() # Close the file
  }
  
  if (Save==TRUE & Type=='png'){ 
    print('png')
    png(paste0(Name,'.png'))
    grid.arrange(p2,p1,p6,p5,ncol=1,newpage = F)  
    dev.off() # Close the file
  }
  
  if (Save==TRUE & Type=='svg'){ 
    print('svg')
    svg(paste0(Name,'.svg'))
    grid.arrange(p2,p1,p6,p5,ncol=1,newpage = F)  
    dev.off() # Close the file
  }
  
  
  if (Save==FALSE){return(grid.arrange(p2,p1,p6,p5,ncol=1,newpage = F))}
}



############## END ################



#' Run in Parallele 
#' @import dbscan
#' @return data frame with all computed values
#' @export 
#' 
ClusterSamplesOld <- function(File){
  df1<-File # All samples Log R Ratio values
  
  # Calculate mean LRR across chromosomes
  df1$meanALLautosChr<-(df1$m1+df1$m2+df1$m3++df1$m4+df1$m5+df1$m6+df1$m7+df1$m8+df1$m9+df1$m10
                        +df1$m11+df1$m12+df1$m13+df1$m14+df1$m15+df1$m16+df1$m17+df1$m18
                        +df1$m19+df1$m20+df1$m21+df1$m22)/22
  
  df1$sdALLautosChr<-(df1$s1+df1$s2+df1$s3+df1$s4+df1$s5+df1$s6+df1$s7+df1$s8+df1$s9+df1$s10
                      +df1$s11+df1$s12+df1$s13+df1$s14+df1$s15+df1$s16+df1$s17+df1$s18
                      +df1$s19+df1$s20+df1$s21+df1$s22)/22
  
  mydata<-data.frame(xp=df1$mX,yp=df1$mY)
  cl1 <- hdbscan(mydata, minPts = 8)
  mydata$cluster<-cl1$cluster
  mydata$prob<-cl1$membership_prob
  return(mydata)
  
  
}

