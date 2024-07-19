#读取seurat对象
seurat_obj <- readRDS('~/Seurat/seurat_obj.rds')

#创建WGCNA对象
seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "variable",
    wgcna_name = "WGCNA" 
)

#创建metacells
seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("seurat_clusters", "orig.ident"),
    k = 25, 
    reduction = "harmony", 
    slot='counts',
    max_shared = 10,
    ident.group = 'seurat_clusters'
)

seurat_obj <- NormalizeMetacells(seurat_obj)

metacell_obj <- GetMetacellObject(seurat_obj)

seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = c("0", "1","2","3","4","5","6", "7", "8", "9" ,"10"), 
    group.by='seurat_clusters',
    assay = 'RNA', 
    use_metacells = TRUE, 
    slot = 'data'
)

#选取软阈值
seurat_obj <- TestSoftPowers(
  seurat_obj,
  powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2)))

plot_list <- PlotSoftPowers(seurat_obj,
                            point_size = 5,
                            text_size = 3)
                            # selected_power = NULL

wrap_plots(plot_list, ncol=2)    

#计算拓扑重叠矩阵
seurat_obj <- ConstructNetwork(
    seurat_obj, 
    soft_power=6,
    setDatExpr=FALSE,
    corType = "pearson",
    networkType = "signed",
    TOMType = "signed",
    detectCutHeight = 0.995,
    minModuleSize = 50,
    mergeCutHeight = 0.2,
    numericLabels = FALSE,
    tom_outdir = "SG-TOM",
    tom_name = 'SG' # name of the topoligical overlap matrix written to disk
)
PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')

seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# 计算模块特征值
seurat_obj <- ModuleEigengenes(
    seurat_obj,
    scale.model.use = "linear", #  choices are "linear", "poisson", or "negbinom"
    assay = NULL, # 默认:DefaultAssay(seurat_obj)
    pc_dim = 1,
    group.by.vars="orig.ident" # 根据样本去批次化 harmonize
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

cur_traits <- c('nCount_RNA', 'nFeature_RNA', 'percent.mito')

str(seurat_obj@meta.data[,cur_traits])

#比较模块与性状间的相关性
seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = 'seurat_clusters', 
    corFnc = "bicor", # to obtain Pearson correlation
    corOptions = "use='p'", # to obtain Pearson correlation
    harmonized = TRUE,
    assay = NULL,
    slot = "data", # default to normalized 'data' slot
    group_name = c("0", "1","2","3","4","5","6", "7", "8", "9" ,"10") 
)


hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']