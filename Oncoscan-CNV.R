setwd("~/workspace/Malawi/results/tsv/CopyNumber")
#source("https://bioconductor.org/biocLite.R")

#biocLite("rCGH")
#biocLite("DNAcopy")

library(rCGH)
library(DNAcopy)

#############################################################
#### reading files from Affymetrix Power Tool using rCGH 
#############################################################

path = "."
file.names <- dir(path, pattern =".tsv")
for(i in 1:length(file.names)){
  sampleName <- gsub("\\..*","",file.names[i])
  print(paste0("* Running analysis for the :" , sampleName))
  #outfile <- read.csv(file.names[i], sep = "\t", header = TRUE, skip = 14)
  cgh <- readAffyOncoScan(file.names[i], sampleName = sampleName ,labName = "Vikki")
  cgh <- adjustSignal(cgh, nCores=4)
  cgh <- segmentCGH(cgh, UndoSD = 1, nCores=4)
  cgh <- EMnormalize(cgh, verbose = TRUE)
  segTable <- getSegTable(cgh)
  #bygene <- byGeneTable(segTable)
  #out.file <- rbind(out.file, file)
  # plotDensity(cgh)
  # multiplot(cgh, symbol = c("dlg2", "gpc5"))
  # recenter(cgh) <- 2
  # pdf("test.file")
  # plotProfile(cgh, symbol = c("dlg2", "gpc5"))
  # dev.off()
  write.csv(segTable, paste0(sampleName, ".csv", sep = "")) 
}

#############################################################
#### merging cvs files
#############################################################
# pasting file name in the first column 
# path = "."
# file.names <- dir(path, pattern =".csv")
# for(i in 1:length(file.names)){
#   sampleName <- gsub("\\..*","",file.names[i])
#   print(paste0("* Running analysis for the :" , sampleName))
#   file <- read.csv(file.names[i], header =TRUE)
#   file <- file[-1,]
#   samplesN = rep(sampleName,length(rownames(file)))
#   file <- cbind(samplesN,file)
#   write.csv(file, paste0(sampleName, ".csv", sep = "")) 
# }
setwd("~/workspace/Malawi/results/tsv/CopyNumber/SegFile/")
filenames <- list.files(path = ".", pattern = ".csv")
merge <- do.call("rbind", lapply(filenames, read.csv, header = TRUE))
names(merge)
segfile <- (merge)[,2:7]
write.table(segfile,"132OncoscanSample_UndoSD=1_PlusEMnormalised.seg", sep = "\t", row.names = FALSE)

###################################################################
#### reading files from Affymetrix Power Tool using DNACopy Package 
###################################################################
setwd("~/workspace/Malawi/results/tsv/CopyNumber")
path = "."
file.names <- dir(path, pattern =".tsv")
for(i in 1:length(file.names)){
  data <- read.csv(file.names[i], header = TRUE, sep = "\t")
  sampleName <- gsub("\\..*","",file.names[i])
  print(paste0("* Running analysis for the :" , sampleName))
  CNA.object <- CNA(cbind(data$Log2Ratio), data$Chromosome,data$Position, data.type="logratio",sampleid=sampleName)
  smoothed.CNA.object <- smooth.CNA(CNA.object)
  segment.smoothed.CNA.object <- segment(smoothed.CNA.object, min.width = 50, undo.splits = "sdundo", verbose=2)
  segTable <- segment.smoothed.CNA.object$output
  write.csv(segTable, paste0(sampleName, ".csv", sep = ""), row.names = FALSE, quote = FALSE)
}

#########

setwd("~/workspace/Malawi/results/tsv/CopyNumber/SegFile/")
filenames <- list.files(path = ".", pattern = ".csv")
merge <- do.call("rbind", lapply(filenames, read.csv, header = TRUE))
segTable <- merge
colnames(segTable)[2:4] <- c ("chromosome","start", "end")
segfile <- smoothSeg(segTable, min.diff=0.25)
write.table(segfile,"MalawiSDNULL.seg", sep = "\t", row.names = FALSE, quote = FALSE)


###################################################################
#### modified package of DNA COPY number functions
###################################################################

changepoints <- function(genomdat, data.type="logratio", alpha=0.01, weights=
                           NULL, sbdry, sbn, nperm=10000, p.method="hybrid", 
                         min.width=2, kmax=25, nmin=200, trimmed.SD=NULL, 
                         undo.splits="none", undo.prune=0.05, undo.SD=3,
                         verbose=1, ngrid=100, tol=1e-6)
{
  n <- length(genomdat)
  if (missing(trimmed.SD)) trimmed.SD <- mad(diff(genomdat))/sqrt(2)
  #   start with the whole 
  seg.end <- c(0,n)
  k <- length(seg.end)
  change.loc <- NULL
  weighted <- ifelse(is.null(weights), FALSE, TRUE) 
  while (k > 1)
  {
    current.n <- seg.end[k]-seg.end[k-1]
    if (verbose>=3) cat(".... current segment:",seg.end[k-1]+1,"-",seg.end[k],"\n")
    if(current.n >= 2*min.width) {
      current.genomdat <- genomdat[(seg.end[k-1]+1):seg.end[k]]
      #   check whether hybrid method needs to be used
      hybrid <- FALSE
      delta <- 0
      if ((p.method=="hybrid") & (nmin < current.n)) {
        hybrid <- TRUE
        delta <- (kmax+1)/current.n
      }
      #   call the changepoint routine
      if (weighted) {
        #   get the weights for the current set of probes
        current.wts <- weights[(seg.end[k-1]+1):seg.end[k]]
        current.rwts <- sqrt(current.wts)
        current.cwts <- cumsum(current.wts)/sqrt(sum(current.wts))
        #   if all values of current.genomdat are the same don't segment
        if (isTRUE(all.equal(diff(range(current.genomdat)), 0))) {
          zzz <- list()
          zzz$ncpt <- 0
        } else {
          #   centering the current data will save a lot of computations later
          current.avg <- sum(current.genomdat*current.wts)/sum(current.wts)
          current.genomdat <- current.genomdat - current.avg
          #   need total sum of squares too
          current.tss <- sum(current.wts*(current.genomdat^2))
          zzz <- .Fortran("wfindcpt",
                          n=as.integer(current.n),
                          x=as.double(current.genomdat),
                          tss=as.double(current.tss),
                          wts=as.double(current.wts),
                          rwts=as.double(current.rwts),
                          cwts=as.double(current.cwts),
                          px=double(current.n),
                          sx=double(current.n),
                          nperm=as.integer(nperm),
                          cpval=as.double(alpha),
                          ncpt=integer(1),
                          icpt=integer(2),
                          hybrid=as.logical(hybrid),
                          al0=as.integer(min.width),
                          hk=as.integer(kmax),
                          mncwt=double(kmax),
                          delta=as.double(delta),
                          ngrid=as.integer(ngrid),
                          sbn=as.integer(sbn),
                          sbdry=as.integer(sbdry),
                          tol= as.double(tol),
                          PACKAGE="DNAcopy")
        }
      } else { 
        #   if all values of current.genomdat are the same don't segment
        if (isTRUE(all.equal(diff(range(current.genomdat)), 0))) {
          zzz <- list()
          zzz$ncpt <- 0
        } else {
          #   centering the current data will save a lot of computations later
          current.avg <- mean(current.genomdat)
          current.genomdat <- current.genomdat - current.avg
          #   need total sum of squares too
          current.tss <- sum(current.genomdat^2)
          zzz <- .Fortran("fndcpt",
                          n=as.integer(current.n),
                          x=as.double(current.genomdat),
                          tss=as.double(current.tss),
                          px=double(current.n),
                          sx=double(current.n),
                          nperm=as.integer(nperm),
                          cpval=as.double(alpha),
                          ncpt=integer(1),
                          icpt=integer(2),
                          ibin=as.logical(data.type=="binary"),
                          hybrid=as.logical(hybrid),
                          al0=as.integer(min.width),
                          hk=as.integer(kmax),
                          delta=as.double(delta),
                          ngrid=as.integer(ngrid),
                          sbn=as.integer(sbn),
                          sbdry=as.integer(sbdry),
                          tol= as.double(tol),
                          PACKAGE="DNAcopy")
        }
      }
    } else {
      zzz <- list()
      zzz$ncpt <- 0
    }
    if(zzz$ncpt==0) change.loc <- c(change.loc,seg.end[k])
    seg.end <- switch(1+zzz$ncpt,seg.end[-k],
                      c(seg.end[1:(k-1)],seg.end[k-1]+zzz$icpt[1],seg.end[k]),
                      c(seg.end[1:(k-1)],seg.end[k-1]+zzz$icpt,seg.end[k]))
    k <- length(seg.end)
    if(verbose>=3) cat(".... segments to go:",seg.end,"\n")
  }
  seg.ends <- rev(change.loc)
  nseg <- length(seg.ends)
  lseg <- diff(c(0,seg.ends))
  if (nseg > 1) {
    if (undo.splits == "prune") {
      lseg <- changepoints.prune(genomdat, lseg, undo.prune)
    }
    if (undo.splits == "sdundo") {
      lseg <- changepoints.sdundo(genomdat, lseg, trimmed.SD, undo.SD)
    }
  }
  segmeans <- 0*lseg
  ll <- uu <- 0
  for(i in 1:length(lseg)) {
    uu <- uu + lseg[i]
    if (weighted) {
      segmeans[i] <- sum(genomdat[(ll+1):uu]*weights[(ll+1):uu])/sum(weights[(ll+1):uu])
    } else {
      segmeans[i] <- mean(genomdat[(ll+1):uu])
    }
    ll <- uu
  }
  list("lseg" = lseg, "segmeans" = segmeans)
}

changepoints.prune <- function(genomdat, lseg, change.cutoff=0.05) {
  n <- length(genomdat)
  nseg <- length(lseg)
  ncpt <- nseg-1
  zzz <- .Fortran("prune",
                  as.integer(n),
                  as.double(genomdat),
                  as.integer(nseg),
                  as.integer(lseg),
                  as.double(change.cutoff),
                  double(nseg),
                  as.integer(ncpt),
                  loc=integer(ncpt),
                  integer(2*ncpt),
                  pncpt=integer(1), PACKAGE="DNAcopy")
  pruned.ncpt <- zzz$pncpt
  pruned.cpts <- cumsum(lseg)[zzz$loc[1:pruned.ncpt]]
  pruned.lseg <- diff(c(0,pruned.cpts,n))
  pruned.lseg
}

changepoints.sdundo <- function(genomdat, lseg, trimmed.SD, change.SD=3) {
  change.SD <- trimmed.SD*change.SD
  cpt.loc <- cumsum(lseg)
  sdundo <- TRUE
  while(sdundo) {
    k <- length(cpt.loc)
    if (k>1) {
      segments0 <- cbind(c(1,1+cpt.loc[-k]),cpt.loc)
      segmed <- apply(segments0, 1, function(i,x) {median(x[i[1]:i[2]])}, genomdat)
      adsegmed <- abs(diff(segmed))
      if (min(adsegmed) < change.SD) {
        i <- which(adsegmed == min(adsegmed))
        cpt.loc <- cpt.loc[-i]
      } else {
        sdundo <- FALSE
      }
    } else {
      sdundo <- FALSE
    }
  }
  lseg.sdundo <- diff(c(0,cpt.loc))
  lseg.sdundo
}

trimmed.variance <- function(genomdat, trim=0.025)
{
  n <- length(genomdat)
  n.keep <- round((1-2*trim)*(n-1))
  inflfact(trim)*sum((sort(abs(diff(genomdat)))[1:n.keep])^2 / (2*n.keep))
}

inflfact <- function(trim)
{
  a <- qnorm(1-trim)
  x <- seq(-a,a,length.out=10001)
  x1 <- (x[-10001] + x[-1])/2
  1/(sum(x1^2*dnorm(x1)/(1-2*trim))*(2*a/10000))
}


########################################################
##### changing the mim.width upto 50
########################################################
segment <- function(x, weights=NULL, alpha=0.01, nperm=10000, p.method= 
                      c("hybrid","perm"), min.width=2, kmax=25, nmin=200, 
                    eta=0.05, sbdry=NULL, trim = 0.025, undo.splits=
                      c("none","prune", "sdundo"), undo.prune=0.05, undo.SD=3,
                    verbose=1)
{
  if (!inherits(x, 'CNA')) stop("First arg must be a copy number array object")
  call <- match.call()
  if (min.width < 2 | min.width > 100) stop("minimum segment width should be between 2 and 100") ##### changed ---> ("minimum segment width should be between 2 and 100")
  if (nmin < 4*kmax) stop("nmin should be >= 4*kmax")
  if (missing(sbdry)) {
    if (nperm==10000 & alpha==0.01 & eta==0.05) {
      if (!exists("default.DNAcopy.bdry")) data(default.DNAcopy.bdry, package="DNAcopy",envir=environment())
      sbdry <- get("default.DNAcopy.bdry", envir=environment())
    } else {
      max.ones <- floor(nperm*alpha) + 1
      sbdry <- getbdry(eta, nperm, max.ones)
    }
  }
  weighted <- ifelse(missing(weights), FALSE, TRUE)
  #   rudimentary error checking for weights
  if (weighted) {
    if (length(weights) != nrow(x)) stop("length of weights should be the same as the number of probes")
    if (min(weights) <= 0) stop("all weights should be positive")
  }
  sbn <- length(sbdry)
  nsample <- ncol(x)-2
  sampleid <- colnames(x)[-(1:2)]
  uchrom <- unique(x$chrom)
  data.type <- attr(x, "data.type")
  p.method <- match.arg(p.method)
  undo.splits <- match.arg(undo.splits)
  segres <- list()
  segres$data <- x
  allsegs <- list()
  allsegs$ID <- NULL
  allsegs$chrom <- NULL
  allsegs$loc.start <- NULL
  allsegs$loc.end <- NULL
  allsegs$num.mark <- NULL
  allsegs$seg.mean <- NULL
  segRows <- list()
  segRows$startRow <- NULL
  segRows$endRow <- NULL
  for (isamp in 1:nsample) {
    if (verbose>=1) cat(paste("Analyzing:", sampleid[isamp],"\n"))
    genomdati <- x[,isamp+2]
    ina <- which(is.finite(genomdati))
    genomdati <- genomdati[ina]
    trimmed.SD <- sqrt(trimmed.variance(genomdati, trim))
    chromi <- x$chrom[ina]
    #      maploci <- x$maploc[ina]
    if (weighted) {
      wghts <- weights[ina]
    } else {
      wghts <- NULL
    }
    sample.lsegs <- NULL
    sample.segmeans <- NULL
    for (ic in uchrom) {
      if (verbose>=2) cat(paste("  current chromosome:", ic, "\n"))
      segci <- changepoints(genomdati[chromi==ic], data.type, alpha, wghts,
                            sbdry, sbn, nperm, p.method, min.width, kmax,
                            nmin, trimmed.SD, undo.splits, undo.prune,
                            undo.SD, verbose)
      sample.lsegs <- c(sample.lsegs, segci$lseg)
      sample.segmeans <- c(sample.segmeans, segci$segmeans)
    }
    sample.nseg <- length(sample.lsegs)
    sample.segs.start <- ina[cumsum(c(1,sample.lsegs[-sample.nseg]))]
    sample.segs.end <- ina[cumsum(sample.lsegs)]
    allsegs$ID <- c(allsegs$ID, rep(isamp,sample.nseg))
    allsegs$chrom <- c(allsegs$chrom, x$chrom[sample.segs.end])
    allsegs$loc.start <- c(allsegs$loc.start, x$maploc[sample.segs.start])
    allsegs$loc.end <- c(allsegs$loc.end, x$maploc[sample.segs.end])
    allsegs$num.mark <- c(allsegs$num.mark, sample.lsegs)
    allsegs$seg.mean <- c(allsegs$seg.mean, sample.segmeans)
    segRows$startRow <- c(segRows$startRow, sample.segs.start)
    segRows$endRow <- c(segRows$endRow, sample.segs.end)
  }
  allsegs$ID <- sampleid[allsegs$ID]
  allsegs$seg.mean <- round(allsegs$seg.mean, 4)
  allsegs <- as.data.frame(allsegs)
  allsegs$ID <- as.character(allsegs$ID)
  segres$output <- allsegs
  segres$segRows <- as.data.frame(segRows)
  segres$call <- call
  if (weighted) segres$weights <- weights
  class(segres) <- "DNAcopy"
  segres
}

########################################################################
smoothSeg <- function(seg, ...){ # require(multicore)
  #get vector of sample names 
  allSamples <- as.character(unique(seg$ID))
  #extract each sample in turn and smooth the segments for each sample within the main seg file 
  seg <- lapply(allSamples, function(sample, segs=seg)
  { sampleSeg <- segs[segs$ID==sample,] 
    sampleSeg <- mergeSample(sampleSeg, ...) }) #lapply returns a list, rbind the results back together as a data frame 
  seg <- do.call(rbind, seg) }

