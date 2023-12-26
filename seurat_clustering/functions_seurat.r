RunSeurat <- function(counts,project="project",meta.data,nPCA,resolution,nFeature_RNA_min=200,nFeature_RNA_max=3500,
	percent.mt_max=8,nfeatures_var=3000,vars.to.regress=NULL,PLOT=FALSE){

obj <- CreateSeuratObject(counts = counts, project = project , min.cells = 3,
min.features = 200, meta.data = meta.data)

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")

# Filter cells
obj <- subset(obj, subset = nFeature_RNA > nFeature_RNA_min & nFeature_RNA < nFeature_RNA_max & percent.mt < percent.mt_max)

# Normalization and dimentsional reduction
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures_var)
all.genes <- rownames(obj)
if(is.null(vars.to.regress)){
obj <- ScaleData(obj, features = all.genes)}else{
  obj <- ScaleData(obj, features = all.genes,vars.to.regress=vars.to.regress)  
}
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
ElbowPlot(obj)
obj=JackStraw(obj);obj=ScoreJackStraw(obj,dims = 1:nPCA)
JackStrawPlot(obj,reduction = "pca",dims = 1:nPCA)

# clustering and visualization
obj <- FindNeighbors(obj, dims = 1:nPCA,k.param = 10,reduction = "pca")
obj <- FindClusters(obj, resolution = resolution)
obj <- RunTSNE(obj, dims = 1:nPCA)
obj <- RunUMAP(obj, dims = 1:nPCA)
DimPlot(obj, reduction = "umap",pt.size=pt.size,label.size = label.size)

# Cell Cycle Scores
cc.genes <- readLines(con = "data/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43];g2m.genes <- cc.genes[44:97]
obj <- CellCycleScoring(obj,s.features = s.genes,g2m.features = g2m.genes)

}



pheatmap_w <- function(x,...){
  pheatmap(x,display_numbers = TRUE,number_format = "%0.f",cluster_cols = F,cluster_rows = F,color = c("white"),legend = FALSE,...)
}

subset_SVZ <- function(obj){
  ss <- subset(obj,subset = annotation %in% c("aRG1","aRG2","opNSC","npNSC"))
}

cols.fp=c("royalblue","skyblue1","orange","red")
####






