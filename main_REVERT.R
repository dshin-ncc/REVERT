'''
Pipeline 'REVERT' (REVERse Transition) 
created by Dongkwan Shin
on 2022.09.18 

Main code of REVERT
original version for D.Shin
'''

rm(list = ls())

# Set the directory where the `main_REVERT.R` file is located
root.folder <- "/YourFolder"
setwd(root.folder)

library(Seurat)
library(pheatmap)
library(dendextend)
library(monocle)
library(smerc)
source(paste(root.folder,"/func_REVERT.R", sep=""))

########### Create a folder to save all the results
###################################################
foldername <- "CRC_organoid_results"
res.dir <- paste(root.folder, "/Result/", foldername, sep="")
dir.create(res.dir)

# Detect available CPU cores
maxCore <- 40 # Maximum number of CPU cores for logical inference step

corr.cutoff <- 0.8 # Correlation cutoff for gene-gene correlation in GRN



'''#########################################################################################################
A. Data Preparation 
This step involves loading expression data and filtering genes based on HVGs and relevant marker genes.
Users can preprocess the data using their own methods. Approximately 3000–7000 HVGs are considered.
Important CRC-related super enhancer genes (e.g., `H3K27ac_res_annotation.rda`) and marker genes from MSigDB are used.
Metadata files include labels for normal and tumor samples in the `nmtu` column.
'''

##### Load expression data 
###########################
load("Organoid_seurat_TSdata.rda")
TS.seurat <- TSdata
# Source: normal, cancer

##### (i) Highly variable genes 
################################
n.hvg <- 7000 # Number of highly variable genes
TS.seurat <- FindVariableFeatures(TS.seurat, selection.method = "vst", nfeatures = n.hvg)
top.hvg <- head(VariableFeatures(TS.seurat), n.hvg)

##### (ii) Super enhancer genes
################################
pval <- 0.05 # P-value for selecting super enhancer genes
load("H3K27ac_res_annotation.rda")
super.e <- H3K27ac_res_annotation

se.sig <- dplyr::filter(super.e, pvalue < pval)
se.set <- data.frame(symbol = se.sig$SYMBOL, log2fc = se.sig$log2FoldChange)
uniq.se <- names(table(se.set$symbol))

# Filter unique or duplicated genes with consistent log2foldchange
temp.set <- sapply(uniq.se, function(x) 
  if((dim(dplyr::filter(se.set, symbol == x))[1] == 1) | prod(sign(dplyr::filter(se.set, symbol == x)$log2fc))){
    c(x, mean(dplyr::filter(se.set, symbol == x)$log2fc))
  })

uniq.se.set <- data.frame(symbol = temp.set[1,], log2fc = temp.set[2,]) 
se.genes <- uniq.se.set$symbol

##### (iii) Colon marker genes from MSigDB
###########################################
marker.up <- unique(read.table("colon_cancer_up.txt", header = FALSE))[,1]
marker.dn <- unique(read.table("colon_cancer_dn.txt", header = FALSE))[,1]
marker.cna <- unique(read.table("colon_cancer_cna.txt", header = FALSE))[,1]

marker.genes <- union(marker.up, marker.dn)

##### Final genes of interest
###############################
crc.genes <- union(se.genes, marker.genes)
goi <- union(top.hvg, crc.genes)

sub.genes <- intersect(rownames(TS.seurat), goi) 
sub.genes <- sub.genes[-grep("^MT-", sub.genes)] # Remove mitochondrial genes

exp.mat <- data.frame(GetAssayData(TS.seurat, slot="counts")[sub.genes,]) # Genes x cells 
meta <- data.frame(sample_ID = colnames(exp.mat), nmtu = TS.seurat$source)
rownames(meta) <- colnames(exp.mat)








'''#########################################################################################################
B. Monocle - Extract a trajectory from normal to cancer 
This section infers pseudo-time trajectories using Monocle2 and extracts the most interesting sub-trajectories.
Results may vary depending on hyperparameters such as `number_cells_expressed`, `qval`, and `max_components`.
Three sets of DEGs (differentially expressed genes) are selected based on different p-value thresholds.
These DEG sets are used in the next step (SCENIC) to extract the core network structure from the constructed GRN.
'''

##### Monocle2
################################
exp <- as.matrix(exp.mat[rowSums(exp.mat) != 0, ]) # Gene x cell matrix

# Calculate the number of cells expressing each gene
number_cells_expressed <- rowSums(ifelse(exp == 0, 0, 1))

# Create sample sheet
sample_sheet <- data.frame(row.names = meta$sample_ID, nmtu = meta$nmtu)
pd <- new("AnnotatedDataFrame", data = sample_sheet)

# Create gene annotation
gene_annotation <- data.frame(
  row.names = rownames(exp),
  gene_short_name = rownames(exp),
  number_cells_expressed = number_cells_expressed
)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

# Initialize a CellDataSet
cds <- newCellDataSet(
  as(exp, "sparseMatrix"), 
  phenoData = pd, 
  featureData = fd, 
  expressionFamily = negbinomial.size()
)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# Define hyperparameters
Ncellexp <- 60  # Minimum number of cells expressing each gene
Vqval <- 0.01  # Q-value threshold for DEG selection
Nmaxcomp <- 2  # Maximum number of components for dimensionality reduction

# Filter genes based on expression
expressed_genes <- row.names(subset(fData(cds), number_cells_expressed >= Ncellexp))

# Perform differential expression analysis
diff_test_res <- differentialGeneTest(cds[expressed_genes, ], fullModelFormulaStr = "~nmtu")
ordering_genes <- row.names(subset(diff_test_res, qval < Vqval))

# Set ordering genes for trajectory analysis
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)

# Perform dimensionality reduction and cell ordering
cds <- reduceDimension(cds, max_components = Nmaxcomp, method = 'DDRTree', auto_param_selection = TRUE)
cds <- orderCells(cds)

# Plot cell trajectory
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "nmtu")
plot_cell_trajectory(cds, color_by = "State") + facet_wrap(~State, nrow = 1)

##### Choose a sub-trajectory
################################
# Select specific states (e.g., State 1&3)
# cds_s <- cds[, c(which(cds$State == 1))]
cds_s <- cds[,c(which(cds$State == 1), which(cds$State == 3))]
                   
# Define genes to test
to_be_tested <- row.names(fData(cds_s))
cds_subset <- cds_s[to_be_tested, ]

##### Select DEGs by q-value
######################################
# qval.list <- c(0.00001, 0.001, 0.05)
qval.list <- c(0.001, 0.05)
diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_diff <- diff_test_res[, c("gene_short_name", "pval", "qval")]

# Generate lists of significant DEGs for each q-value threshold
sig_gene_names <- sapply(qval.list, function(x) row.names(subset(sig_diff, qval < x)))

####### Expression and metadata for a sub-trajectory
setwd(res.dir)
save(cds_subset, file = "save_cds_subset.RData")

sub.exp <- Biobase::exprs(cds_subset)
sub.meta <- meta[colnames(sub.exp), ]


'''#######################################################################
B-1. GeneSwitches - Check transition genes along the pseudo-time trajectory  
This step identifies switching genes along the pseudo-time trajectory using the GeneSwitches algorithm.  
Switching genes represent genes that change significantly as cells transition along pseudo-time.  
The DEGs identified in step B and the switching genes are combined to create a final set of genes of interest.  
Switching genes often overlap with DEGs identified by Monocle.  
Using only DEGs without switching genes may also yield similar results.
'''

###### GeneSwitches 
##################################
# Identify switching genes along the pseudo-time trajectory
switch.genes <- func_switchgenes(cds_subset) # Outputs include sg_top15 and sg_total










'''#########################################################################################################
C. pySCENIC - Construct a GRN through the imputation of expression data  
This step applies pySCENIC to filtered expression data to construct a large-scale GRN (Gene Regulatory Network).  
Imputed data generated using MAGIC is used for smoothing expression data in the next steps.  
The expression data from the sub-trajectory is provided as input to pySCENIC, which is executed in Python.  
Results are saved in the `/pyscenic` folder within the input data directory.  
The output includes regulons and their target genes.  
The GRN structure is represented as `reg.per.gene`, a binary matrix where rows are target genes and columns are regulons.  
Additionally, `reg.per.cell` represents regulon activity per cell. Some column names may require corrections to match formatting.  
'''

### Imputation through MAGIC
library(Rmagic)

# Prepare data for MAGIC imputation
magic.exp <- t(sub.exp)
keep_cols <- colSums(magic.exp > 0) > 10 # Retain genes expressed in at least 10 cells
magic.exp <- magic.exp[, keep_cols] # Cell x genes

# Normalize and transform data
magic.exp <- library.size.normalize(magic.exp)
magic.exp <- sqrt(magic.exp)

# Perform MAGIC imputation
magic.exp <- magic(magic.exp)
save(magic.exp, file = "magic_exp_data.RData")

##### Read pySCENIC results 
###########################
library(stringr) # For string manipulation

# Read regulon activity per cell
reg.per.cell <- read.table(
  paste(root.folder, "/Result/", foldername, "/pyscenic/regulon_per_cell.csv", sep = ""), 
  header = TRUE
)
rownames(reg.per.cell) <- rownames(sub.meta)

# Fix column names by removing trailing characters (e.g., "...")
temp.list <- colnames(reg.per.cell)
colnames(reg.per.cell) <- str_sub(temp.list, 1, nchar(temp.list) - 3)

# Transpose to get regulon x cell matrix
regulon.exp <- as.matrix(t(reg.per.cell))

# Read regulon-target gene matrix
reg.per.gene <- read.table(
  paste(root.folder, "/Result/", foldername, "/pyscenic/regulon_per_gene.csv", sep = ""), 
  header = TRUE
)
temp.list <- colnames(reg.per.gene)
colnames(reg.per.gene) <- str_sub(temp.list, 1, nchar(temp.list) - 3)

##### Make GRN from pySCENIC results: reg.per.gene --> edge list 
#################################################################
library(igraph)

# Transpose regulon-target matrix and extract edges
trpg <- t(reg.per.gene)
xy.ind <- which(trpg == 1, arr.ind = TRUE)

# Create edge list
edges_all <- data.frame(
  source = rownames(trpg)[xy.ind[, 1]], 
  target = colnames(trpg)[xy.ind[, 2]], 
  direction = rep(1, dim(xy.ind)[1])
)
# gall <- graph_from_data_frame(edges_all)

'''#######################################################################
C-1. Visualize pseudo-time and find genes highly correlated with pseudo-time  
Visualize pseudo-time in PCA/UMAP space and calculate correlations between genes and pseudo-time.
'''

### Data preparation for pseudo-time analysis
sub.meta.ptime <- sub.meta
sub.meta.ptime$Pseudotime <- cds_subset$Pseudotime

# Analyze pseudo-time correlations
pseudotime.res <- func_pseudotime(
  t(magic.exp$result), 
  sub.meta.ptime, 
  normalized = FALSE, 
  ifpseudoga = FALSE
)
save(pseudotime.res, file = "pseudotime_result.RData")










'''#########################################################################################################
D. Logic inference - Construct a dynamical GRN  
This section infers Boolean logical rules for dynamic GRN modeling.  
Different network structures are inferred based on varying hyperparameters such as smoothing window width,  
the set of DEGs, and high-correlation links (e.g., |correlation| > 0.8).  
Strongly connected components are extracted to identify core structures of the network.  
The binarized expression data (after smoothing) is required for logic inference.  
Different logic rules may be inferred depending on the time steps, typically set to be larger than the smoothing width.  
Three DEG sets (from step B) based on p-value thresholds are considered, and the results are analyzed  
for different time steps (smoothing window width).  
'''

# Define the list of smoothing window widths
moving.width.list <- c(1, 50, 100, 300, 500) # Smoothing window widths for expression data

# Define switching genes obtained from the GeneSwitches algorithm
switching.genes <- switch.genes$sg_total$geneID

# Set a flag for elbow point calculation
elbowtruefalse <- TRUE

# Perform logic inference
func_logic_infer(
  res.dir,                # Directory for saving results
  maxCore,                # Maximum number of CPU cores to use
  moving.width.list,      # List of smoothing window widths
  corr.cutoff,            # Correlation cutoff for high-correlation links
  sig_gene_names,         # Significant genes based on q-value thresholds
  pseudotime.res,         # Pseudo-time analysis results
  edges_all,              # Edge list from pySCENIC GRN
  crc.genes,              # CRC-related genes
  switching.genes,        # Switching genes from GeneSwitches
  elbowTF = elbowtruefalse # Elbow calculation flag
)










'''#########################################################################################################
E. Summarize the results of attractor simulation 
This section summarizes the characteristics of logical models inferred using various hyperparameters.  
The following metrics are calculated for each model:
- **agree.mat**: Agreement level matrix, indicating the level of agreement between inferred models.
- **net.size**: The size of the inferred network (number of nodes or edges).
- **att.mat**: Average activity of genes in the top attractor or the ratio of "1" states in the top attractor.
- **basin.mat**: Basin size of the top attractor.
- **distance.n.mat**: Distance from the closest attractor to the normal attractor among the top 5 largest attractors.
- **distance.c.mat**: Distance from the closest attractor to the cancer attractor among the top 5 largest attractors.

Normal and cancer attractors are defined by binarizing the average activities of the first and last 1000 cells  
in the inferred pseudo-time axis, respectively.
'''

# Define hyperparameters
moving.width.list <- c(1, 50, 100, 300, 500) # Smoothing window widths for attractor simulation
ratio.att <- 0.2 # Fraction of cells at the endpoints used to define normal and cancer attractors
Ndeg.set <- 3    # Number of DEG sets to analyze

# Perform attractor analysis
func_attractor_analysis(
  res.dir,            # Directory for saving results
  moving.width.list,  # Smoothing window widths
  ratio.att,          # Ratio for defining normal and cancer attractors
  Ndeg.set            # Number of DEG sets
)






#############################################################################################
#############################################################################################

