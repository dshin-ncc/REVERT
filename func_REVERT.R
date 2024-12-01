

#############################################################################################
## Function for identifying switching genes. 
## input: cds data for a selected subtrajectory
## output: switching genes 
#############################################################################################
library(GeneSwitches)
library(SingleCellExperiment)

func_switchgenes = function(cds_subset){

  ordered.meta <- dplyr::arrange(pData(cds_subset), Pseudotime)
  monocle.exp.data <- Biobase::exprs(cds_subset)[,rownames(ordered.meta)] # genes x cells

  sce.gs <- SingleCellExperiment(assays = list(counts = monocle.exp.data))
  sce.gs <- scuttle::logNormCounts(sce.gs)
  assay(sce.gs, "expdata") <- logcounts(sce.gs)
  colData(sce.gs)$Pseudotime <- cds_subset$Pseudotime

  pca <- prcomp(t(assays(sce.gs)$expdata), scale.=FALSE)
  rd_PCA <- pca$x[,1:2]
  reducedDims(sce.gs) <- SimpleList(PCA = rd_PCA)
  h <- hist(as.matrix(assays(sce.gs)$logcounts), breaks = 200, plot = FALSE)
  {plot(h, freq = FALSE, xlim = c(0,2), ylim = c(0,1), main = "Histogram of gene expression",
        xlab = "Gene expression", col = "darkgoldenrod2", border = "grey")
    abline(v=0.2, col="blue")}

  bn_cutoff <- 0.2
  sce.gs <- binarize_exp(sce.gs, fix_cutoff = TRUE, binarize_cutoff = bn_cutoff)

  sce.gs <- find_switch_logistic_fastglm(sce.gs, downsample = FALSE, show_warning = FALSE)
  sg_allgenes <- filter_switchgenes(sce.gs, allgenes = TRUE, topnum = 15)
  sg_total <- filter_switchgenes(sce.gs, allgenes = TRUE)

  return(list(sg_top15 = sg_allgenes, sg_total = sg_total))

}







#############################################################################################
## Function for basic visualization, PCA, UMAP, and temporal gene expression, of pseudo-time
## meta should have pseudo-time info. meta$Pseudotime
## output includes correlation value of each gene with pseudo-time, highly correlated gene list.
#############################################################################################
func_pseudotime <- function(exp, meta, normalized = FALSE, ifpseudoga = FALSE){
  library(SingleCellExperiment)
  library(pseudoga)
  library(scran)
  library(scater)
  library(parallel)
  
  # randomize the order of cells
  # idx.rand <- sample(dim(exp)[2])
  # exp <- exp[,idx.rand]
  # meta <- meta[idx.rand,]
  
  if (normalized == FALSE){
    #building sce
    sce = SingleCellExperiment(list(counts = exp), colData = meta)
    clust <- quickCluster(sce) 
    sce <- computeSumFactors(sce, cluster=clust, min.mean=0.1)
    sce <- logNormCounts(sce)
    
    #making a fresh sce
    sce.fresh = SingleCellExperiment(assay = list(normalized = as.matrix(assay(sce, 2))), colData = colData(sce)) #only taking the second assay of the sce, which is normalized
    
  }else{
    #making a fresh sce
    sce = SingleCellExperiment(list(counts = exp), colData = meta)
    clust <- quickCluster(sce) 
    sce <- computeSumFactors(sce, cluster=clust, min.mean=0.1)
    sce.fresh = SingleCellExperiment(assay = list(normalized = as.matrix(assay(sce, 1))), colData = colData(sce)) #
  }

  if (ifpseudoga == TRUE){
        library(tictoc)
  tic()
  # sce.fresh = pseudoga_parallel(sce.fresh, type = "normalized") #default setting.
  sce.fresh = pseudoga(sce.fresh, type = "normalized") #default setting.
  toc() #fast

  }

  
  #quick visualization confirm
  library(PCAtools)
  library(BiocSingular)
  #pca.ck = pca(mat = ck.mat, rank = 5)
  pca.reg = pca(mat = assay(sce.fresh), rank = 2, BSPARAM = FastAutoParam()) #really only need 2
  pca.coords = pca.reg$rotated[,1:2]
  
  reducedDim(sce.fresh, "pca_lognorm") = pca.coords
  
  #plotting
  library(gridExtra)
  png(filename = paste("pseudotime_PCA_plot_",Sys.Date(),".png",sep=""), width = 2048, height = 1024, units = "px", res = 300)
  g1<-print(plotReducedDim(sce.fresh, dimred = "pca_lognorm", colour_by = "Pseudotime"))
  g2<-print(plotReducedDim(sce.fresh, dimred = "pca_lognorm", colour_by = paste(names(meta)[2])))
  grid.arrange(g1,g2,nrow=1)
  dev.off()
  
  png(filename = paste("pseudotime_UMAP_plot_",Sys.Date(),".png",sep=""), width = 2048, height = 1024, units = "px", res = 300)
  sce.fresh <- runUMAP(sce.fresh,dimred="pca_lognorm")
  p1=print(plotReducedDim(sce.fresh, dimred = "UMAP", colour_by = "Pseudotime"))
  p2=print(plotReducedDim(sce.fresh, dimred = "UMAP", colour_by = paste(names(meta)[2])))
  grid.arrange(p1,p2,nrow=1)
  dev.off()
  
  ############ Genes with highest correlation with the pseudotime
  
  cors1 <- NULL
  data <- assay(sce.fresh,1)
  cors1 <- cor(t(data), colData(sce.fresh)$Pseudotime, method="spearman")
  # cors1 <- cor(t(data), colData(sce.fresh)$Pseudotime, method="pearson")
  abscors1 <- abs(cors1)
  geneord1 <- order(abscors1, decreasing = TRUE)
  pseudogene1<-rownames(data)[geneord1]
  corr.v <- cors1[geneord1]
  # head(pseudogene1)


  ##################### Temporal changes of 25 genes at top or near elbow points 
  library(smerc)
  abscorr <- abs(corr.v)
  abscorr <- abscorr[!is.na(abscorr)] # remove NA for the genes that have zero expression over all cells
  opt.pt <- elbow_point(1:length(abscorr), abscorr) # opt.pt$x, opt.pt$y

  png(filename = paste("Corr_genes_pseudotime.png",sep=""),
      width = 1024, height = 1024, units = "px", res = 300)
  plot(abscorr, main = "Correlation according to ranking", xlab = "ranking", ylab = "correaltion")
  points(opt.pt$x, opt.pt$y, pch=16, col="red", cex=1.3)
  text(opt.pt$x, opt.pt$y, labels=paste("Elbow pt.(",opt.pt$x,", ",round(opt.pt$y, digits=2),")"), pos = 4, cex=1 )
  dev.off()
  
  rnkidx <- rank(colData(sce.fresh)$Pseudotime,ties.method="random")
  library(gridExtra)
  png(filename = paste("Highly_ranked_genes_25.png",sep=""), 
      width = 4096, height = 4096, units = "px", res = 300)
  par(mfrow=c(5,5))
  for (i in 1:25) {
    plot(rnkidx,data[pseudogene1[i],],col="red",main=pseudogene1[i],
         xlab="Sequence",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)
    }
  dev.off()
  
  png(filename = paste("Genes_near_elbow_25.png",sep=""), 
      width = 4096, height = 4096, units = "px", res = 300)
  par(mfrow=c(5,5))
  for (i in (opt.pt$x - 12):(opt.pt$x + 12)) {
    plot(rnkidx,data[pseudogene1[i],],col="red",main=pseudogene1[i],
         xlab="Sequence",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)
    }
  dev.off()

  return(list(sce.fresh = sce.fresh, exp.data = data, pseudogenes = pseudogene1, corr.v = corr.v ))
}







#############################################################################################
## Function for smoothing and binarizing the given expression data
## Input expression data should be ouput of 'func_pseudotime', pseudotime.res.
## moving.width = width of moving window for smoothing
## moving.full.length = TRUE : the length of smoothing data is the same as that of original
#############################################################################################
func_smoothing_binarization <- function(pseudotime.res, moving.width, moving.full.length = TRUE){
  
    
  pga.sce <- pseudotime.res$sce.fresh
  pga.data <- pseudotime.res$exp.data
  pga.genes <- pseudotime.res$pseudogenes
  pga.corr <- pseudotime.res$corr.v
  
  library(smerc)
  abscorr <- abs(pga.corr)
  abscorr <- abscorr[!is.na(abscorr)] # remove NA for the genes that have zero expression over all cells
  abscorr2 <- abscorr[1:floor(length(abscorr)/2)] # consider 1/2 part from the top rank. otherwise, elbow_point may give the second elbow point.
  opt.pt <- elbow_point(1:length(abscorr2), abscorr2) # opt.pt$x, opt.pt$y

    
  ##### smoothing #####
  ptime.seq <- pga.sce$Pseudotime
  names(ptime.seq) <- pga.sce$sample_ID
  
  ord.ptime <- order(ptime.seq, decreasing = FALSE)
  ord.cell.name <- names(ptime.seq)[ord.ptime]
  
  all.order.exp <- pga.data[pga.genes, ord.cell.name]
  window.width <- moving.width
  movavg.order.exp<-t(edgeR::movingAverageByCol(t(all.order.exp), width = window.width, 
                                         full.length = moving.full.length))
  colnames(movavg.order.exp) <- ord.cell.name # 'movingAverageByCol' has a bug where a part of name are lost (half length of moving window?)
  
  t.time <- ptime.seq[ord.ptime[window.width:length(ord.ptime)]]
  png(filename = paste("smoothing_25highly_ranked_genes.png",sep=""), 
      width = 4096, height = 4096, units = "px", res = 300)
  par(mfrow=c(5,5))
  for (i in 1:25) {
    plot(movavg.order.exp[i,],col="red",main=pga.genes[i],
         xlab="Sequence",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)
    # plot(colData(pga.sce)$Pseudotime,pga.data[pga.genes[i],],col="red",main=paste("Expression of ",pga.genes[i],sep=""),xlab="Pseudotime",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)
  }
  dev.off()
  
  png(filename = paste("smoothing_25_genes_near_elbow.png",sep=""), 
      width = 4096, height = 4096, units = "px", res = 300)
  par(mfrow=c(5,5))
  for (i in (opt.pt$x - 12):(opt.pt$x + 12)) {
    plot(movavg.order.exp[i,],col="red",main=pga.genes[i],
         xlab="Sequence",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)
    # plot(colData(pga.sce)$Pseudotime,pga.data[pga.genes[i],],col="red",main=paste("Expression of ",pga.genes[i],sep=""),xlab="Pseudotime",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)
  }
  dev.off()
  
   ############################# Binarization #################################
  library(BoolNet)
  
  bin.order.exp <- binarizeTimeSeries(all.order.exp)
  bin.order.exp1 <- data.frame(t(bin.order.exp$binarizedMeasurements))
  bin.order.exp2 <- cbind(rownames(bin.order.exp1), bin.order.exp1)
  colnames(bin.order.exp2)[1] <- "sample_ID"

  binplot0 <- as.matrix(t(bin.order.exp1))

  png(filename = paste("binarized_count_exp_25highly_ranked_genes.png",sep=""), 
      width = 4096, height = 4096, units = "px", res = 300)
  par(mfrow=c(5,5))
  for (i in 1:25) {
    plot(binplot0[i,],col="red",main=pga.genes[i],
         xlab="Sequence",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)
    # plot(colData(pga.sce)$Pseudotime,pga.data[pga.genes[i],],col="red",main=paste("Expression of ",pga.genes[i],sep=""),xlab="Pseudotime",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)
  }
  dev.off()
  
  
  bin.movavg <- binarizeTimeSeries(movavg.order.exp)
  bin.movavg.exp <- data.frame(t(bin.movavg$binarizedMeasurements))
  bin.movavg.exp2 <- cbind(rownames(bin.movavg.exp), bin.movavg.exp)
  colnames(bin.movavg.exp2)[1] <- "sample_ID"
  
  binplot <- as.matrix(t(bin.movavg.exp))

  png(filename = paste("binarized_25highly_ranked_genes.png",sep=""), 
      width = 4096, height = 4096, units = "px", res = 300)
  par(mfrow=c(5,5))
  for (i in 1:25) {
    plot(binplot[i,],col="red",main=pga.genes[i],
         xlab="Sequence",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)
    # plot(colData(pga.sce)$Pseudotime,pga.data[pga.genes[i],],col="red",main=paste("Expression of ",pga.genes[i],sep=""),xlab="Pseudotime",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)
  }
  dev.off()
  
  png(filename = paste("binarized_25_genes_near_elbow.png",sep=""), 
      width = 4096, height = 4096, units = "px", res = 300)
  par(mfrow=c(5,5))
  for (i in (opt.pt$x - 12):(opt.pt$x + 12)) {
    plot(binplot[i,],col="red",main=pga.genes[i],
         xlab="Sequence",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)
    # plot(colData(pga.sce)$Pseudotime,pga.data[pga.genes[i],],col="red",main=paste("Expression of ",pga.genes[i],sep=""),xlab="Pseudotime",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)
  }
  dev.off()
  
  save(pseudotime.res, movavg.order.exp, bin.movavg.exp2, file = "pseudotime_smooth_binary_result.RData")
  return(list(pseudotime.res = pseudotime.res, smoothing.data = movavg.order.exp, binary.cnt.res = bin.order.exp2, binary.res = bin.movavg.exp2, opt.pt=opt.pt ))
}







#############################################################################################
## Function for converting scns output rules to logic forms
## 1. excute 'rule2logic' in the folder where scns results are saved
## 2. replace several characters with different forms needed for Boolnet
## 3. the final logic is the 'Union' of all possible logics
## 4. move the final results to the input directory(directory.name) with input names(in.suffix)
## 
#############################################################################################
func_rule2logic <- function(logic.res.folder, directory.name, in.suffix){

  
  ### converting logic inference results to simple form of logics
  ################################################################
  setwd(logic.res.folder)
  filePath<-paste0(logic.res.folder,'/')
  files <- list.files(path=filePath, pattern='*.txt')
  files

  
  output_file_whole <- c()
  agreement_level_whole <- c()
  for(file in files){
    # read input file.
    raw_file <- utils::read.table(paste(filePath, file, sep=""), stringsAsFactors = F, sep='\t', quote="\"'")
    
    # distinguish the target name from the input file name like 'ADH1B_boolean_rules_5.txt'.
    target <- strsplit(file, '_')[[1]][1]
    
    # extract agreement level from the first line of input file.
    agreement_level <- as.numeric(strsplit(raw_file[1,], ' = ')[[1]][2])
    
    # find where the best rules are located.
    best_rule_idx <- which(raw_file == 'The best rules were:')
    other_rule_idx <- which(raw_file == 'Other rules found were:')
    
    best_rule <- raw_file[(best_rule_idx+1):(other_rule_idx-1),]
    
    # process each rules into logics using 'rule2logic' function.
    bool_logic_file <- c()
    for(rule in best_rule){
      bool_logic_file <- rbind(bool_logic_file, paste(target, '=' , rule2logic(rule), sep=' '))
    }
    
    # stack each agreement levels and output logics.
    agreement_level_whole <- rbind(agreement_level_whole, c(target, agreement_level))
    output_file_whole <- rbind(output_file_whole, bool_logic_file)
  }
  # save the output file.
  logic.folder <- paste(filePath,"rule2logic",sep="")
  if(!dir.exists(logic.folder)) dir.create(logic.folder)
  setwd(logic.folder)

  filenames <- paste(logic.folder,"/whole_logic.txt",sep="")
  filenamesagree <- paste(logic.folder,"/agreement_level_whole.txt",sep="")

  write.table(agreement_level_whole, file = filenamesagree, sep="\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(output_file_whole, file = filenames, quote = FALSE, col.names = FALSE, row.names = FALSE)


  # utils::write.table(agreement_level_whole, 'agreement_level_whole.txt', col.names = F, row.names = F, quote = F, sep = "\t")
  # utils::write.table(output_file_whole, 'whole_logic.txt', col.names = F, row.names = F, quote = F)

  ### Saving multiple logics to union of them for each gene
  ################################################################

  # filenames <- paste(logic.folder,"/whole_logic.txt",sep="")
  # filenamesagree <- paste(logic.folder,"/agreement_level_whole.txt",sep="")
  filenamesunion <- paste(logic.folder,"/whole_logic_union.txt",sep="")
  
  # Replace "=" and "~" with "," and "!", respectively
  x <- readLines(filenames)
  y <- gsub( " =", ",", x )
  y <- gsub( "~", "!", y )
  cat(y, file=filenames, sep="\n")
  
  ## Review output
  cat(readLines(filenames), sep="\n")
  
  ## connect multiple logics for a gene with union 
  logics <- readLines(filenames)
  agr <- read.table(filenamesagree, header = FALSE)
  genes.logic <- agr[,1]
  
  yy <- array()
  for( ii in 1:length(genes.logic) ){ 
    temp.gene <- paste(genes.logic[ii],", ",sep="")
    log.idx <- grep(temp.gene, logics)
    
    if(length(log.idx) > 1){
      temp.split <- strsplit(logics[log.idx], split = temp.gene)
      union.logic <- paste("(",paste(unlist(lapply(1:length(log.idx), function(x) temp.split[[x]][2])), collapse = ")|("), ")", sep = "")
      yy[ii] <- paste(temp.gene, union.logic, sep = "")
    } else {
      yy[ii] <- grep(temp.gene, logics, value =T)
    }
  }
  cat(yy, file = filenamesunion, sep = "\n")
  
  system(paste("(echo targets,factors && cat ",filenames,") > ",logic.folder,"/logic.txt",sep=""))
  system(paste("(echo targets,factors && cat ",filenamesunion,") > ",logic.folder,"/logic_union.txt",sep=""))
  
  system(paste("rm ",filenames,sep=""))
  system(paste("rm ",filenamesunion,sep=""))
  
  ## move the output files and rename 
  topath <- directory.name
  renamefile <- paste("logic_",in.suffix,".txt",sep = "")
  renamefile1 <- paste("logic_union_",in.suffix,".txt",sep = "")
  renamefile2 <- paste("agreement_level_",in.suffix,".txt",sep = "")
  
  file.rename( from = file.path(logic.folder, "logic.txt") ,
               to = file.path(topath, renamefile) )
  file.rename( from = file.path(logic.folder, "logic_union.txt") ,
               to = file.path(topath, renamefile1) )
  file.rename( from = file.path(logic.folder, "agreement_level_whole.txt") ,
               to = file.path(topath, renamefile2) )
  
  system(paste("rm ",filePath,"*.txt",sep=""))
  
  # barplot agreement level
  agr <- read.table(paste(topath,"/",renamefile2,sep=""), header = FALSE)
  agr2 <- agr[,2]
  names(agr2)<-agr[,1]
  
  library(scater)
  png(filename = paste(topath,"/agreement_",in.suffix,".png",sep=""), width = 1024, height = 800,  res = 100)
  barplot(agr2)
  dev.off()
}






rule2logic <- function(rule_line){
  # extract only logics out of other properties. logics should look like "[(logic)]".
  logic <- regmatches(rule_line, gregexpr("(?<=\\[).*?(?=\\])", rule_line, perl=T))[[1]]
  
  # separate 'logic' into elements.
  logic_parse <- regmatches(logic, gregexpr("(?<=\\().*?(?=\\))", logic, perl=T))[[1]]
  
  # convert 'logic_parse' into matrix
  logic_mat <- t(sapply(logic_parse, function(x){strsplit(x, ', ')[[1]]}))
  rownames(logic_mat) <- NULL;
  
  # find activator indices(has 'a') and repressor indices(has 'r') out of 'logic_parse'
  a_idx <- which(stri_count_fixed(logic_mat[,1], 'a') != 0)
  r_idx <- which(stri_count_fixed(logic_mat[,1], 'r') != 0)
  
  # conversion matrix between 'or', 'and' and '|', '&', respectively.
  or_and <- matrix(c("or","and",' | ',' & '), nrow=2)
  
  # make logic for activators #
  # if there is two activators, 'a0' is one of 'and' or 'or'.
  # if there is one or less activators, 'a0' is a node name.
  if(length(a_idx) %in% c(0,1)){
    a_line <- logic_mat[a_idx,2]
  }else{
    oper_ind<- which(logic_mat[,2] %in% or_and[,1])
    oper_ind <- oper_ind[oper_ind %in% a_idx]
    a_or_and <- or_and[or_and[,1] %in% logic_mat[oper_ind,2],2]
    # a_line <- paste('(', paste(logic_mat[a_idx[-1],2], collapse = a_or_and), ')', sep='')
    # a_line <- paste('(', paste(logic_mat[a_idx[-oper_ind],2], collapse = a_or_and), ')', sep='')
    a_line <- paste('(', paste(logic_mat[a_idx[a_idx!=oper_ind],2], collapse = a_or_and), ')', sep='')
  }
  
  # make logic for repressors #
  # if there is two repressors, 'r0' is one of 'and' or 'or'.
  # if there is one or less repressors, 'r0' is a node name.
  if(length(r_idx) %in% c(0,1)){
    r_line <- logic_mat[r_idx,2]
  }else{
    oper_ind<- which(logic_mat[,2] %in% or_and[,1])
    oper_ind <- oper_ind[oper_ind %in% r_idx]
    r_or_and <- or_and[or_and[,1] %in% logic_mat[oper_ind,2],2]
    r_line <- paste('(', paste(logic_mat[r_idx[-which(r_idx == oper_ind)],2], collapse=r_or_and), ')', sep='')
  }
  
  # combine activator line and repressor line #
  if(length(r_line)==0){
    ar_line <- a_line
  }else if(length(a_line)==0){
    ar_line <- paste('~', r_line, sep='')
  }else{
    ar_line <- paste(a_line, r_line, sep=' & ~')
  }
  
  return(ar_line)
}

















#############################################################################################
## Function for logic inference
## The structure of the different network models is determined by the width of the moving window for smoothing and the DEGs considered.
## To determine the network structure, we extract the network filtered by DEGs from the complex GRN obtained from pyscenic.
## Then filter again to links with high correlation between two genes (ex,abs value >0.8).
## Finally, we get the strongly connected component to extract the core structure.
##
## Input: Result directory, maximum number of cores, moving width list, expression data, network edges, marker.genes, switching.genes 
## Output: the below list will be saved in the directory for each smoothing width.
## avg.agr: agreement level that represents how simulated data fit well with real data
## top.att: top 5 attractors   
## top.basin: the corresponding basin size of the top attractors
## basinsizelist: the list of basin size for all the attractors
## logics inferred, figures of pruned networks, dynamic patterns of genes, and so on will be saved.
##
#############################################################################################

func_logic_infer <- function(res.dir, maxCore, moving.width.list, corr.cutoff, sig_gene_names, pseudotime.res, edges_all, marker.genes = NULL, switching.genes = NULL, elbowTF = NULL){

  Nloop <- length(sig_gene_names)*length(moving.width.list)
  pb = txtProgressBar(min = 0, max = Nloop, style = 3) 
  stepi = 0

  for (imov in 1:length(moving.width.list)){
    moving.width <- moving.width.list[imov]
    
    sub.folder <- paste(res.dir,"/smoothing_",moving.width,sep="")
    dir.create(sub.folder)
    setwd(sub.folder)
    
    ##### Smoothing expression data with moving window and binarizing it
    #####################################################################
    pseudotime.bin.res <- func_smoothing_binarization(pseudotime.res, moving.width, moving.full.length = TRUE)
    save(pseudotime.bin.res, file = "pseudotime_smoothing_bin_results.RData")
    # load("pseudotime_smoothing_bin_results.RData")
    elbow.genes <- pseudotime.res$pseudogenes[1:pseudotime.bin.res$opt.pt$x]

      
    glist<-rep(0, times = length(sig_gene_names)) # number of genes in each DEG case
    
    ##### Determine NETWORK STRUCTURES 
    ##### GRN from Pyscenic --> filtering with DEGs --> filtering with g-g correlation --> get SCC
    ##############################################################################################
    for (ideg in 1:length(sig_gene_names)) {
      
      # genes.oi<-unique(c(switching.genes$geneID, elbow.genes, sig_gene_names[[ideg]])) # genes of interest: switching genes + elbow genes + DEGs
      # genes.oi<-sig_gene_names[[ideg]] # DEGs

      if(elbowTF == TRUE){
        genes.oi<-unique(c(sig_gene_names[[ideg]], switching.genes, elbow.genes))
        }else{
          genes.oi<-unique(c(sig_gene_names[[ideg]], switching.genes))
        }
      
      ## Filtering with DEGs and higly correlated genes 
      overlap.genes <- intersect(rownames(pseudotime.bin.res$smoothing.data), genes.oi)
      ttcorr <- cor(t(pseudotime.bin.res$smoothing.data[overlap.genes,]), method = "spearman") # correlation matrix

      corr.sign <- (ttcorr > corr.cutoff) + (-1*(ttcorr < -corr.cutoff)) # matrix of correlation sign

      ##### Pruning the GRN through correlation matrix
      ################################################
      message("pruning the GRN through correaltion matrix")
      edges.grn.hrgs <- data.frame(source = edges_all$source, target = edges_all$target,
                                   direction = sapply(1:dim(edges_all)[1], function(x)
                                     ifelse(is.element(edges_all[x,1], rownames(corr.sign)) & is.element(edges_all[x,2], rownames(corr.sign)),
                                            corr.sign[edges_all[x,1], edges_all[x,2]],0)))
      
      save(edges.grn.hrgs, sig_gene_names, switching.genes,elbow.genes, file=paste("EDGES_GRN_", ideg,"DEG.RData"))

      
      ##### Extract SCC from the pruned GRN
      #####################################
      xxk<-dplyr::filter(edges.grn.hrgs, direction !=0) # 
      gxxk <- graph_from_data_frame(xxk)
      gxxkc <- components(gxxk, mode = c("strong"))
      nodess.total = which(gxxkc$membership == which.max(gxxkc$csize)) # maximal scc
      # nodess.total = which(gxxkc$membership %in% which(gxxkc$csize != 1)) # all sccs
      gxxkscc <- induced_subgraph(gxxk, nodess.total)
      
      marker.nodes <- intersect(marker.genes, names(gxxkc$membership)[nodess.total]) # to show marker genes with different color
      V(gxxkscc)$color = "orange"
      V(gxxkscc)$frame.color <- "white"
      V(gxxkscc)[marker.nodes]$color = "#076fa2"
      
      png(filename = paste("SCC_net_",ideg,"DEG.png",sep=""), width = 1920, height = 1920, units = "px", res = 300)
      plot(gxxkscc, edge.arrow.size = 0.5, edge.curved = .1,  vertex.size = igraph::degree(gxxkscc), 
                   edge.color = ifelse(E(gxxkscc)$direction == -1, "red", "skyblue"),
                   vertex.label.color = "black", vertex.label.dist=1, vertex.label.family = "Helvetica",
                   vertex.label.cex = 0.8)
      dev.off()
      
      ##### plot temporal changes of SCC network nodes from smoothing data
      ####################################################################
      smooth.data <- pseudotime.bin.res$smoothing.data
      scc.nodes <- names(V(gxxkscc))
      library(gridExtra)
      png(filename = paste("pseudotime_scc_genes_", ideg,"DEG.png",sep=""), 
          width = 4096, height = 4096, units = "px", res = 300)
      par(mfrow=c(5,5))
      for (i in 1:min(length(scc.nodes), 25)) {
        plot(smooth.data[scc.nodes[i],],col="red",main=scc.nodes[i],
             xlab="Sequence",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)
      }  
      dev.off()
      
      ##### plot temporal changes of SCC network nodes from binarized data
      ####################################################################
      bin.cnt.data <- pseudotime.bin.res$binary.cnt.res
      scc.nodes <- names(V(gxxkscc))
      library(gridExtra)
      png(filename = paste("pseudotime_bin_cnt_scc_genes_", ideg,"DEG.png",sep=""), 
          width = 4096, height = 4096, units = "px", res = 300)
      par(mfrow=c(5,5))
      for (i in 1:min(length(scc.nodes), 25)) {
        plot(bin.cnt.data[,scc.nodes[i]],col="red",main=scc.nodes[i],
             xlab="Sequence",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)
        
      }
      dev.off()
      
      
      ##### save nodes and edges data of the scc network
      ##################################################
      g.scc <- gxxkscc
      subnet.name <- names(V(g.scc))
      subnet.edges <- data.frame(source = as_edgelist(g.scc)[,1], target = as_edgelist(g.scc)[,2], direction = E(g.scc)$direction)
      
      ## Every genes should have at least one activation input, otherwise, there will be an error in SCNS
      ## Add self-activation to the genes that have no activation input.
      for (igene in 1:length(subnet.name)){
        tempgene <- subnet.name[igene]
        if(sum(dplyr::filter(subnet.edges, target == tempgene)$direction == "1") == 0) {
          subnet.edges <- rbind(subnet.edges, c(tempgene, tempgene, "1"))
        }
      }
      
      save.subnet.bin.cnt.mat <- pseudotime.bin.res$binary.cnt.res[,c("sample_ID", subnet.name)]
      write.table( save.subnet.bin.cnt.mat, 
                   file = paste(sub.folder,"/bin_cnt_exp_",ideg,"_DEG_g",length(subnet.name),".txt",sep=""),
                   sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
      
      save.subnet.bin.mat <- pseudotime.bin.res$binary.res[,c("sample_ID", subnet.name)]
      write.table( save.subnet.bin.mat, 
                   file = paste(sub.folder,"/bin_exp_",ideg,"_DEG_g",length(subnet.name),".txt",sep=""),
                   sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
      
      edges_subnet <- data.frame(from.gene = subnet.edges$source, relation = subnet.edges$direction, to.gene = subnet.edges$target)
      write.table( edges_subnet, 
                   file = paste(sub.folder,"/net_edges_",ideg,"_DEG_g",length(subnet.name),".txt",sep=""),
                   sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
      
      glist[ideg] <- length(subnet.name)
      
    }
    
    
    ##### Determine BOOLEAN LOGICs of the network (Pseudotime-network-inference)
    ####################################################

    for (ideg in 1:length(glist)) {
      
      ##### input data for sncs: expression data and network structure
      #################################################################
      smooth.exp.data.path <- paste(sub.folder,"/bin_exp_",ideg,"_DEG_g",glist[ideg],".txt",sep="")
      edge.data.path <- paste(sub.folder,"/net_edges_",ideg,"_DEG_g",glist[ideg],".txt",sep="")
      
      
      alltstep.list <- moving.width.list
      ## tstep (input-output time step) larger than smoothing width is reasonable 
      tstep.list <- alltstep.list[alltstep.list >= moving.width] 
      modewidth <- 1 # another hyper-parameter to be selected in SCNS. just set to 1
      
          
      library(BoolNet)
      library(igraph)
      avg.agr <- rep(0, times = length(tstep.list))
      top.att <- list()
      top.basin <- list()
      
      for(i.step in 1:length(tstep.list)) {

        tstep <- tstep.list[i.step]

        ##### RUN Logic inference function based on Pseudotime-network-inference
        #####       (https://github.com/fionahamey/Pseudotime-network-inference)
        ##### The outputs of Pseudotime-network-inference will be saved to ROOT.FOLDER/scns_results/.
        ##### put rule2logic.R file into the folder 'scns_results'.
        #####################################################################
              
        logic.res.folder <- paste(res.dir, "/logic_inf_results",sep="")  
        if(!dir.exists(logic.res.folder)) dir.create(logic.res.folder) # create a temporal folder for logic inference results ------
        setwd(res.dir)
        
        pythonexe <- "/home/dkshin/anaconda3/envs/pyscenic_env/bin/python3" # the path of python3 installed
        logicpy <- "/mnt/gluster_server/dkshin/logic_inf/net_logic_inference.py"      
        ncores <- ifelse(glist[ideg]<maxCore, glist[ideg], maxCore)
        system(paste(pythonexe,logicpy, tstep, modewidth, ncores, smooth.exp.data.path, edge.data.path, sep=" "))
        
        
        ##### run rule2logic: convert Pseudotime-network-inference output rules to logic forms 
        ##### modify the folder paths in rule2logic as your environments before running.
        ################################################################################
        # setwd(logic.res.folder)      
        
        in.suffix <- paste(ideg,"DEG_step",tstep,"_mode",modewidth,sep="")
        func_rule2logic(logic.res.folder,sub.folder, in.suffix)
        
        ## Plot network model from the inferred logics
        setwd(sub.folder)
        
        orig.net <- read.table(paste(sub.folder,"/net_edges_",ideg,"_DEG_g",glist[ideg],".txt",sep=""), header = TRUE )
        testnet<- loadNetwork(paste("logic_union_",in.suffix,".txt", sep=""))
        logic.net <- plotNetworkWiring(testnet)

        logic.net.edge <- get.edgelist(logic.net)
        logic.net.sign <- data.frame(source = logic.net.edge[,1], target = logic.net.edge[,2], 
                                  direction = sapply(1:dim(logic.net.edge)[1], function(x) 
                                    dplyr::filter(orig.net, from.gene == logic.net.edge[x,1] & to.gene == logic.net.edge[x,2])$relation))

        g.logic <- graph_from_data_frame(logic.net.sign)
        V(g.logic)$color = "orange"
        V(g.logic)$frame.color = "white"
        png(filename = paste("logic_net_",in.suffix,".png",sep=""), width = 1920, height = 1920, units = "px", res = 300)
        plot(g.logic, edge.arrow.size = 0.5, edge.curved = .1, vertex.size = igraph::degree(g.logic),
             edge.color = ifelse(E(g.logic)$direction == -1, "red", "skyblue"),
             vertex.label.color = "black", vertex.label.dist=1, vertex.label.family = "Helvetica",
             vertex.label.cex = 0.8)
        dev.off()

        
        ##### Get attractors and their basins
        ######################################
        
        ## For a large network, the number of random initial states can be set to 100000 
        temp.gene.n <- length(testnet$genes)
        if(temp.gene.n >= 18) {
          testatt = getAttractors(testnet, method = "random", startStates = 100000)
        }else{
          testatt = getAttractors(testnet)
        }
        
        ## agreement
        agreem <- read.table(paste("agreement_level_",in.suffix,".txt", sep=""), header = FALSE)
        avg.agr[i.step] <- mean(agreem[,2])
        
        ## basin size distibution
        basinsizelist <- sapply(1:length(testatt$attractors), function(x) testatt$attractors[[x]]$basinSize)
        # hist(basinsizelist)
        png(filename = paste("basin_histogram_",in.suffix,".png",sep=""), width = 1920, height = 1920, units = "px", res = 300)
        hist(basinsizelist/sum(basinsizelist)*100, 10) # %
        dev.off()
        
        ## sort attractors according to the basin size
        temp.max <- ifelse(length(basinsizelist) < 5, length(basinsizelist), 5)
        maxbasin <- order(basinsizelist, decreasing = TRUE)[1:temp.max]
        maxatt <- sapply(1:length(maxbasin), function(x) getAttractorSequence(testatt, maxbasin[x]))
        top.att[[i.step]] <- maxatt
        top.basin[[i.step]] <- basinsizelist[maxbasin]
        
        ## plot attractors 
        png(filename = paste("top_attractors_",in.suffix,".png",sep=""), width = 1920, height = 1920, units = "px", res = 300)
        par(mar = c(5, 7, 5, 5))
        plotAttractors(testatt, maxbasin[1:temp.max],title = "Top 5 attractors", mode = "table", 
               onColor = "#4daf4a",offColor = "#d73027", borderColor="white",eps=0.1, allInOnePlot = FALSE, drawLegend = FALSE)
        dev.off()
        
        
        
      } # tstep loop

      setwd(sub.folder)
      save(avg.agr, top.att, top.basin, basinsizelist, file = paste("save_attractor_analysis_smoothing_",moving.width,"_DEG_",ideg,".RData",sep=""))
      
      stepi = stepi+1
      setTxtProgressBar(pb, stepi)
      
    } # deg loop
    
    
    
  } # smoothing loop
  close(pb)



}



















library(ggplot2)
library(tidyr)
library(tidyverse)
library(gridExtra)
library(BoolNet)
library(igraph)
library(monocle)
library(ggExtra)
library(foreach)
library(doParallel)


#############################################################################################
## Function for summarizing the results of logic inference
## Input: Result directory, moving width list, Nnc 
## Output: the below list will be saved.
## agree.mat: import the agreement level, output of scns
## net.size: network size of the inferred network. 
## att.mat: average gene activity of the top attractor, or ratio of "1" of the top attractor
## basin.mat: basin size of the top attractor
## distance.n.mat: distance to the normal attractor from the closest attractor among top 5 big attractors.
## distance.c.mat: distance to the cancer attractor from the closest attractor among top 5 big attractors.
## normal/cancer attracter defined by the binarization of the average activities of 
## the first/last 1000 cells in the inferred pseudo-time axis. 
## 
#############################################################################################
func_attractor_analysis <- function(res.dir, moving.width.list, ratio.att, Ndeg.set){
  setwd(res.dir)
  ### euclidean distance of two vectors
  euclidean <- function(a, b) sqrt(sum((a - b)^2))/sqrt(length(a))

  alltstep.list <- moving.width.list

  ## initialize
  agree.mat <- list()
  net.size <- list()
  att.mat <- list()
  basin.mat <- list()
  n.c.att <- list()
  distance.n.mat <- list()
  distance.c.mat <- list()

  Nloop <- Ndeg.set*length(moving.width.list)
  pb = txtProgressBar(min = 0, max = Nloop, style = 3) 
  stepi = 0

  for (ideg in 1:Ndeg.set){
    
    agree.mat[[ideg]] <- matrix(NA, nrow = length(alltstep.list), ncol = length(moving.width.list))
    net.size[[ideg]] <- matrix(NA, nrow = length(alltstep.list), ncol = length(moving.width.list))
    att.mat[[ideg]] <- matrix(NA, nrow = length(alltstep.list), ncol = length(moving.width.list))
    basin.mat[[ideg]] <- matrix(NA, nrow = length(alltstep.list), ncol = length(moving.width.list))
    distance.n.mat[[ideg]] <- matrix(NA, nrow = length(alltstep.list), ncol = length(moving.width.list))
    distance.c.mat[[ideg]] <- matrix(NA, nrow = length(alltstep.list), ncol = length(moving.width.list))
    
    smooth.nc.att <- list()
    
    for (ismooth in 1:length(moving.width.list)) {
      moving.width <- moving.width.list[ismooth]
      sub.folder <- paste(res.dir,"/smoothing_",moving.width,sep="")
      
      load(paste(sub.folder, "/pseudotime_smoothing_bin_results.RData",sep=""))  
      load(paste(sub.folder,"/save_attractor_analysis_smoothing_",moving.width,"_DEG_",ideg,".RData",sep=""))
      # import avg.agr, top.att, top.basin, and basinsizelist
      
      agree.mat[[ideg]][(length(alltstep.list)-length(avg.agr)+1):length(alltstep.list),ismooth] <- avg.agr
      
      init.tstep <- which(alltstep.list == moving.width)
      step.nc.att <- list()
      for (istep in init.tstep:length(alltstep.list)){
        
        if (is.null(dim(top.att[[istep-init.tstep+1]]))){ ## if network size == 1, top.att[[x]] is a list, not a matrix.
          tt <- top.att[[istep-init.tstep+1]][1]
        }else{
          tt <- top.att[[istep-init.tstep+1]][,1]
        }
        
        n.node <- length(tt) # network size
        att.first <- unlist(lapply(1:n.node, function(x) mean(tt[x][[1]]) )) # for cyclic attractors, use mean activities.
        att.first.mean <- mean(att.first) # the top attractor
        
        bb <- top.basin[[istep-init.tstep+1]]
        basin.first <- bb[1]/sum(basinsizelist) # top basin
        
        net.size[[ideg]][istep, ismooth] <- n.node
        att.mat[[ideg]][istep, ismooth] <- att.first.mean
        basin.mat[[ideg]][istep, ismooth] <- basin.first
        
        if (length(tt) < 2) { ## we don't care about the case of (network size < 2)
          distance.n.mat[[ideg]][istep,ismooth] <- NA
          distance.c.mat[[ideg]][istep,ismooth] <- NA
        }else{
          attlist <- top.att[[istep-init.tstep+1]]
          smooth.data <-pseudotime.bin.res$binary.res 
          Nnc <- ceiling(dim(smooth.data)[1]*ratio.att)
          normal.att <- colMeans(smooth.data[1:Nnc,rownames(attlist)])
          cancer.att <- colMeans(smooth.data[(dim(smooth.data)[1]-Nnc+1):dim(smooth.data)[1],rownames(attlist)])
          step.nc.att[[istep]] <- data.frame(normal.att = normal.att, cancer.att = cancer.att)
          distance.n.mat[[ideg]][istep, ismooth] <- min(sapply(1:dim(attlist)[2], function(x) euclidean(unlist(attlist[,x]), (normal.att>=0.5) )))
          distance.c.mat[[ideg]][istep, ismooth] <- min(sapply(1:dim(attlist)[2], function(x) euclidean(unlist(attlist[,x]), (cancer.att>=0.5) )))
          
        } # if
        
      } # istep loop
      smooth.nc.att[[ismooth]] <- step.nc.att
      
      stepi = stepi+1
      setTxtProgressBar(pb, stepi)
      

    } # ismooth loop
    n.c.att[[ideg]] <- smooth.nc.att
  } # ideg loop

  close(pb)

  save(agree.mat, net.size, att.mat, basin.mat, n.c.att, distance.n.mat, distance.c.mat, file="Summarized_results_of_attractors.Rdata")

  ############################################################################################
  ## Plot the results of Logic inference 
  ## 앞에서 정리한 특징들을 plot하는 단계. 
  ## 6개의 결과 그림이 한 파일에 저장됨
  ## (agree.mat, net.size, att.mat, basin.mat, distance.n.mat, distance.c.mat) 
  ## DEG 별로 위의 그림이 다른 파일에 저장됨.
  ############################################################################################


  plot.obj <- list(agree.mat, net.size, att.mat, basin.mat, distance.n.mat, distance.c.mat)
  plot.obj.string <- c("Agreement level", "Network size", "Mean activity of top attractor", "Basin size of top attractor", 
                       "Distance to normal attractor", "Distance to cancer attractor")
  smooth.char <- as.character(moving.width.list)
  step.char <- smooth.char


  for (ideg in 1:Ndeg.set) {
    DEGnum <- ideg
    
    png(filename = paste("All_PLOTS_",ideg,"DEG.png",sep=""),
        width = 4608, height = 3072, units = "px", res = 300)
    
    g.plot <- list()
    ccut <- list()
    for (iplot in 1:length(plot.obj)) {

      plot.mat <- plot.obj[[iplot]]
      temp.mat <- data.frame(plot.mat[[DEGnum]])
      colnames(temp.mat) <- smooth.char
      temp.mat <- cbind(step.char, temp.mat)
      colnames(temp.mat)[1] <- "Tstep"
      
      ## matrix to pair (Tstep, smoothwindow, level)
      temp.mat2 <- temp.mat %>%
        pivot_longer(
          cols = -Tstep,
          names_to = "smoothwindow",
          values_to = "level"
        )
      
      plot.mat.factor <- data.frame(temp.mat2)
      
      re.step <- plot.mat.factor$Tstep
      re.smooth <- plot.mat.factor$smoothwindow
      
      plot.mat.factor$Tstep <- factor(plot.mat.factor$Tstep, levels = unique(re.step))
      plot.mat.factor$smoothwindow <- factor(plot.mat.factor$smoothwindow, levels = unique(re.smooth))
      
      ## cut off for text color
      d.level <- max(plot.mat.factor$level,na.rm = TRUE)-min(plot.mat.factor$level,na.rm = TRUE)
      color.cut <- d.level*0.8+min(plot.mat.factor$level,na.rm = TRUE)
     
      ccut[[iplot]]<-color.cut
      gg <- ggplot(plot.mat.factor, aes(Tstep, fct_rev(smoothwindow),
                                  fill = level)) +
        geom_tile(color = "white",lwd = 1.5,linetype = 1)+
        geom_text(aes(label = format(round(level, 2), nsmall = 2), color = abs(level) < color.cut)) + # agree.mat:0.69, net.size:15, att.mat:0.75
        # geom_text(aes(label = format(round(level, 2), nsmall = 2))) + # agree.mat:0.69, net.size:15, att.mat:0.75
        coord_fixed(expand = FALSE) +
        scale_color_manual(values = c("white", "black"),
                           guide = "none", na.value = "white")+
        scale_fill_distiller(
          palette = "PuOr", na.value = "white",
          direction = 1 #, limits = c(0,1)
        ) +
        labs(title =plot.obj.string[iplot], x = "Step", y = "Smoothing window width")

      g.plot[[iplot]] <- gg
      
    }
    grid.arrange(g.plot[[1]],g.plot[[3]],g.plot[[5]],g.plot[[2]],g.plot[[4]],g.plot[[6]], nrow=2, ncol=3)
    dev.off()
    
  }


}












#################################### EOF ###############################################
