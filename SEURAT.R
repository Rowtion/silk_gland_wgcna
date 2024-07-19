#合并样本
dir='~/samples_dir/'
samples=list.files( dir  )
sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro) 
  sce=CreateSeuratObject(counts =  Read10X(file.path(dir,pro) ) ,
                         project =  pro )
  
  return(sce)
})

library(stringr)

names(sceList) =  samples

sce.all <- merge(sceList[[1]], y= sceList[ -1 ] ,
                 add.cell.ids =  samples) 

#SEURAT处理
seurat_ob = sce.all

#筛选线粒体基因
mito.genes <- c("COX3", "ND3", "ND5", "ND4", "ND4L", "ND6", "CYTB", "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6")
raw.counts = GetAssayData(seurat_ob, slot = "counts")
percent.mito <- Matrix::colSums(raw.counts[mito.genes, ])/Matrix::colSums(raw.counts)
seurat_ob <- AddMetaData(object = seurat_ob, metadata = percent.mito, col.name = "percent.mito")

#筛选细胞
UMIs_per_cell = Matrix::colSums(raw.counts)
genes_per_cell = Matrix::colSums(raw.counts>0) 
df = data.frame(UMIs_per_cell=UMIs_per_cell, genes_per_cell=genes_per_cell)
df = df[order(df$UMIs_per_cell),] 
m <- rlm(genes_per_cell~UMIs_per_cell,data=df) 
pb <- data.frame(predict(m, interval='prediction',level = 1-1e-3, type="response"))

outliers <- rownames(df)[df$genes_per_cell > pb$upr | df$genes_per_cell < pb$lwr]
seurat_ob = subset(seurat_ob, cells = colnames(raw.counts)[!colnames(raw.counts) %in% outliers])

seurat_by_sample =SplitObject(seurat_ob, split.by = "orig.ident" )

for( idx in 1:length(seurat_by_sample) ){
    object = seurat_by_sample[[idx]]
    upper_bound <- 10^(mean(log10(object@meta.data[,"nFeature_RNA"][object@meta.data[,"nFeature_RNA"]>0])) + 2*sd(log10(object@meta.data[,"nFeature_RNA"][object@meta.data[,"nFeature_RNA"]>0])))
    lower_bound <- 10^(mean(log10(object@meta.data[,"nFeature_RNA"][object@meta.data[,"nFeature_RNA"]>0])) - 2*sd(log10(object@meta.data[,"nFeature_RNA"][object@meta.data[,"nFeature_RNA"]>0])))
object <- subset(object, subset = nFeature_RNA >= lower_bound & nFeature_RNA <= upper_bound)
    seurat_by_sample[[idx]] = object
}
for( idx in 1:length(seurat_by_sample) ){
    object = seurat_by_sample[[idx]]
    upper_bound <- 10^(mean(log10(object@meta.data[,"nCount_RNA"][object@meta.data[,"nCount_RNA"]>0])) + 2*sd(log10(object@meta.data[,"nCount_RNA"][object@meta.data[,"nCount_RNA"]>0])))
    lower_bound <- 10^(mean(log10(object@meta.data[,"nCount_RNA"][object@meta.data[,"nCount_RNA"]>0])) - 2*sd(log10(object@meta.data[,"nCount_RNA"][object@meta.data[,"nCount_RNA"]>0])))
object <- subset(object, subset = nCount_RNA >= lower_bound & nCount_RNA <= upper_bound)
    seurat_by_sample[[idx]] = object
}
for( idx in 1:length(seurat_by_sample) ){
    object = seurat_by_sample[[idx]]
object <- subset(object, subset = percent.mito >= -Inf & percent.mito <= 0.3)
    seurat_by_sample[[idx]] = object
}
merged_seurat = seurat_by_sample[[1]]   
for( idx in 2:length(seurat_by_sample) ){
    merged_seurat = merge(x = merged_seurat, y = seurat_by_sample[[idx]]
)
}
seurat_ob = merged_seurat

#标准化
seurat_ob <- NormalizeData(object = seurat_ob,normalization.method = "LogNormalize",scale.factor = 10000)

#识别高变基因
seurat_ob = FindVariableFeatures(object= seurat_ob, loess.span = 0.3,
                    clip.max = "auto", mean.function = "FastExpMean",
                    dispersion.function = "FastLogVMR", num.bin = 20,
                    nfeature = 4000, binning.method = "equal_width" )

top10 <- head(VariableFeatures(seurat_ob), 10)

#Harmony降维去批次效应
seurat_ob <- ScaleData(object = seurat_ob, features = rownames(seurat_ob),
                    vars.to.regress = c("nCount_RNA","percent.mito"), verbose = T )

seurat_ob <- RunPCA(seurat_ob, features = VariableFeatures(object = seurat_ob))

seurat_ob <- RunHarmony(seurat_ob, "orig.ident")

#聚类
seurat_ob = FindNeighbors( seurat_ob, reduction = "harmony", dims = 1:10,
                                features = VariableFeatures(seurat_ob),
                               nn.eps = 0, force.recalc = T, verbose = F)
seurat_ob <- FindClusters(object = seurat_ob, resolution = 0.4, algorithm = 1, verbose = F)
seurat_ob = RunUMAP(seurat_ob,dims = 1:10,verbose = F,seed.use=1,reduction = "harmony", n.components = 2)

#识别markergene
global_DEGs = FindAllMarkers(object = seurat_ob,only.pos = T,test.use = "bimod", logfc.threshold = 0, min.pct = 0.25)
global_DEGs = global_DEGs %>% mutate( gene_diff = round(global_DEGs$pct.1 / global_DEGs$pct.2, 3)) %>% select( gene, everything())

#提取top10 markergene
topn_markers  = global_DEGs %>% group_by(cluster) %>% 
            arrange(p_val,desc(gene_diff)) %>%
            top_n(10,gene_diff)

saveRDS(seurat_ob, file = "~/Seurat/seurat_obj.rds")

