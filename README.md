# REVERT
<img width="1062" alt="image" src="https://github.com/user-attachments/assets/ed62196f-8520-431d-b9e1-ac980bc6696e">

REVERT aims to reconstruct the dynamic network model of GRN that can represent cellular dynamics in the tumor transition state and to identify potential intervention targets for cancer reversion through in silico perturbation simulations with attractor landscape analysis. Dynamic network modeling is crucial for understanding and predicting the behavior of a system in response to perturbations under untested circumstances, enabling the identification of an optimal reversion switch that can drive the cellular system toward a desired normal cell state upon the attractor landscape. 

This document explains each section in main_REVERT.R.


## Section A: Data Preparation

Overview

This section prepares the input data for downstream analysis. It involves:
	1.	Loading expression data.
	2.	Filtering genes based on criteria like Highly Variable Genes (HVGs) and CRC-specific marker genes.
	3.	Generating a matrix of selected genes and corresponding metadata.

Key Steps

	•	Load Expression Data: Load the raw expression data into an R object.
	•	Filter HVGs: Identify and select genes with high variability across samples.
	•	Incorporate CRC-Specific Genes:
	•	Use super enhancer genes and CRC-related markers from MSigDB.
	•	Combine these with HVGs to generate a comprehensive gene set of interest.
	•	Generate Expression Matrix: Create a filtered expression matrix and corresponding metadata for downstream analysis.



## Section B: Monocle - Extract a Trajectory from Normal to Cancer

Overview

This section uses Monocle2 to infer pseudo-time trajectories, representing transitions from normal to cancer states. It identifies differentially expressed genes (DEGs) along the trajectory.

Key Steps

	•	Pseudo-Time Inference:
	•	Use Monocle2 to compute the pseudo-time trajectory of cells.
	•	Dimensionality reduction is performed using DDRTree.
	•	Cells are ordered based on pseudo-time.
	•	Select DEGs:
	•	Identify genes that are significantly differentially expressed along the pseudo-time trajectory.
	•	Use q-value thresholds to filter DEGs for further analysis.
	•	Sub-Trajectory Extraction:
	•	Extract specific sub-trajectories (e.g., specific cell states) for more focused analyses.

Section B-1: GeneSwitches - Identify Switching Genes

Overview

This supplementary analysis identifies switching genes along the pseudo-time trajectory. Switching genes are critical for understanding the dynamic molecular changes during transitions.

Key Steps

	•	Use the GeneSwitches algorithm to identify genes that switch between different expression states along pseudo-time.
	•	Combine switching genes with DEGs identified in Section B to create a comprehensive list of genes of interest.

Biological Importance

Switching genes often overlap with DEGs but add granularity to understanding the trajectory by pinpointing specific transitions.



## Section C: pySCENIC - Construct a GRN

Overview

This section constructs a large-scale Gene Regulatory Network (GRN) using pySCENIC. It also smooths expression data using MAGIC for downstream analyses.

Key Steps

	1.	Imputation with MAGIC:
	•	Smooth expression data to reduce noise and prepare it for GRN construction.
	2.	pySCENIC Integration:
	•	Use pySCENIC (executed in Python) to compute:
	•	Regulon activity per cell (reg.per.cell).
	•	Regulon-target relationships (reg.per.gene).
	•	Clean and format results for R-based analyses.
	3.	GRN Construction:
	•	Convert pySCENIC results into an edge list for GRN visualization and modeling.
	4.	Pseudo-Time Analysis:
	•	Visualize pseudo-time in PCA/UMAP space and calculate correlations between genes and pseudo-time.



## Section D: Logic Inference - Construct a Dynamical GRN

Overview

This section performs Boolean logic inference to model dynamic GRNs. It identifies core regulatory pathways by analyzing interactions and extracting key structures.

Key Steps

	1.	Define Smoothing Parameters:
	•	Specify window widths for smoothing expression data.
	2.	Logical Inference:
	•	Use Boolean logic to infer regulatory rules for the GRN.
	•	Focus on highly correlated links and DEGs identified in earlier steps.
	3.	Core Network Extraction:
	•	Identify strongly connected components in the network to extract key regulatory pathways.

Output

Dynamic GRNs that model state transitions in the biological system.



## Section E: Attractor Analysis - Summarize Simulation Results

Overview

This section evaluates the stability and dynamics of the GRNs using attractor simulations. Attractors represent stable states of the network, such as normal and cancer states.

Key Steps

	1.	Calculate Attractor Metrics:
	•	Agree.mat: Agreement between different network models.
	•	Net.size: Network size (nodes or edges).
	•	Att.mat: Average gene activity in the top attractor.
	•	Basin.mat: Basin size for the top attractor.
	•	Distance Metrics: Distances from attractors to normal and cancer states.
	2.	Define Normal and Cancer Attractors:
	•	Use the first and last cells in the pseudo-time trajectory to define normal and cancer states.

Biological Importance

Attractors provide insights into the stability of normal and cancer states, as well as the transitions between them.

