# iscience 2023
## Identity and Nature of Neural Stem Cells in the Adult Human Subventricular Zone
### Abstract
The existence of neural stem cells (NSCs) in adult human brain neurogenic regions
remains unresolved. To address this, we created a cell atlas of the adult human
subventricular zone (SVZ) derived from fresh neurosurgical samples using single-cell
transcriptomics. We discovered 2 adult radial glia (RG)-like populations, aRG1 and
aRG2. aRG1 shared features with fetal early RG (eRG) and aRG2 were
transcriptomically similar to fetal outer RG (oRG). We also captured early neuronal and
oligodendrocytic NSC states. We found that the biological programs driven by their
transcriptomes support their roles as early-lineage NSCs. Finally, we show that these
NSCs have the potential to transition between states and along lineage trajectories.
These data reveal that multipotent NSCs reside in the adult human SVZ.
## Main contents

### Seuerat_clustering
  1. Running seurat clustering on SVZ samples and selecting the NSC-cluster.
  2. Selecting progenitor cells from 2 fetal human brain datasets.
  3. Merging the selected cells from 3 datasets.
  4. Some visualizations
#### data
   1- results of FindAllMarkers and some input files used in the main script.

### Velocity
   1. Basic preprocessing: filtering genes and normalizing data
   2. Compute the first- and second-order moments (needed for velocity estimation)
   3. Estimation of velocities
   4. Calculating probabilities of one cell transitioning into another cell, using cosine correlation
   5. PAGA graph abstraction (Wolf et al 2019) provides a graph-like map of the data topology with weighted edges corresponding to the connectivity between two clusters. In scvelo, PAGA is extended by velocity-inferred directionality.
   6. visualization


 
