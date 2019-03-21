suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("--loom", nargs='+',
                    help = "The path to loom file [*.loom] or the path of samples which shoud contain <SAMPLE>/velocyto/*.loom")
parser$add_argument("--seurat",
                    help = "RDS file from Seurat.")
parser$add_argument("--subset",
                    help = "Choose sub idents of seurat.")
parser$add_argument("--cell_total_thresh", default=1e3,
                    help = "Filter out cell total expression less than threshold.")
parser$add_argument("--pagoda2_opts",
                    default="modelType='plain',trim=10,log.scale=T,do.par=T,gam.k=10,nPcs=100,n.odgenes=3e3,maxit=300,k=16,center=T,distance='cosine',method=multilevel.community,name='multilevel',embeddingType='tSNE',perplexity=50,verbose=F,pType='PCA'",
                    help = "Options for pagoda2")
parser$add_argument("--estimate_opts",
                    default="deltaT=1,kCells=16,fit.quantile=0.2",
                    help = "Options for gene.relative.velocity.estimates")
parser$add_argument("--velocity_opts",
                    default="n=1000,scale='linear',cex=0.8,arrow.scale=3,show.grid.flow=T,min.grid.cell.mass=0.5,grid.n=50,arrow.lwd=1,do.par=F,cell.border.alpha = 0,main='Cell Velocity'",
                    help = "Options for show.velocity.on.embedding.cor")
parser$add_argument("--height",
                    default=10,
                    help = "Height of plot (inch)")
parser$add_argument("--width",
                    default=10,
                    help = "Width of plot (inch)")
parser$add_argument("--n.cores",
                    default=parallel::detectCores(),
                    help = "How many cores to use.")
parser$add_argument("--dump_RDS",
                    default=FALSE, action='store_true',
                    help = "Dump RDS file.")
parser$add_argument("--output",
                    help = "output prefix")

args <- parser$parse_args()

load_looms <- function(loom_files){
  looms <- lapply(loom_files, function(loom_file){
    read.loom.matrices(loom_file)
  })
  
  names(looms) <- tools::file_path_sans_ext(basename(loom_files))
  looms
}

merge_sample_looms <- function(looms){
  types <- c('spliced', 'unspliced', 'ambiguous')
  merged_looms <- lapply(types, function(type){
    m <- do.call(cbind, lapply(looms, function(x){x[[type]]}))
    colnames(m) <- gsub('x$', '', gsub(':', '_', colnames(m)))
    m
  })
  names(merged_looms) <- types
  merged_looms
}


suppressPackageStartupMessages(library(velocyto.R))

if(length(args$loom) == 0){
  stop("loom file path not specified.")
}

message("Load loom files")
looms <- load_looms(args$loom)
merged_looms<-merge_sample_looms(looms)

if(length(args$seurat) == 0){
  # pagoda2
  suppressPackageStartupMessages(library(pagoda2))
  emat <- merged_looms$spliced
  emat <- emat[,colSums(emat) >= args$cell_total_thresh]
  idx<-!duplicated(row.names(emat))
  emat <- emat[idx,]
  
  p2_opts <- eval(parse(text=sprintf("list(%s)", args$pagoda2_opts)))
  p2_opts$x <- emat
  p2_opts$n.cores <- args$n.cores
  p2_opts$plot <- F
  p2_opts$type <- p2_opts$pType
  
  r <- do.call(Pagoda2$new, p2_opts[c('x', 'modelType', 'trim', 'log.scale', 'n.cores')])
  do.call(r$adjustVariance, p2_opts[c('plot', 'do.par', 'gam.k')])
  do.call(r$calculatePcaReduction, p2_opts[c('nPcs', 'n.odgenes', 'maxit')])
  do.call(r$makeKnnGraph, p2_opts[c('k', 'type', 'center', 'distance')])
  do.call(r$getKnnClusters, p2_opts[c('method', 'type', 'name')])
  do.call(r$getEmbedding, p2_opts[c('type', 'embeddingType', 'perplexity', 'verbose')])
  rm(p2_opts)
  
  cell_ids <- rownames(r$counts)
  cluster.label <- r$clusters$PCA$multilevel # take the cluster factor that was calculated by p2
  cell.colors <- pagoda2:::fac2col(cluster.label)
  # take embedding form p2
  emb <- r$embeddings$PCA$tSNE
  cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))
}else{
  # seruat
  message("Load data from Seurat RDS file")
  suppressPackageStartupMessages(library(Seurat))
  seurate_data <- readRDS(file=args$seurat)
  if(length(args$subset) != 0){
    sub_ident <- eval(parse(text=sprintf('c(%s)', args$subset)))
    seurate_data <- SubsetData(seurate_data, ident.use = sub_ident, subset.raw = T)
  }
  cluster.label <- seurate_data@ident
  emb <- seurate_data@dr$tsne@cell.embeddings
  cell_ids <- rownames(seurate_data@meta.data)
  cell.dist <- as.dist(1-armaCor(t(emb)))
  
  g <- TSNEPlot(object=seurate_data, do.label = TRUE, pt.size = 1.5)
  cell.colors <- as.list(ggplot_build(g)$data[[1]]$colour)
  names(cell.colors) <- rownames(emb)
  rm(seurate_data)
}

emat <- merged_looms$spliced[, cell_ids]
nmat <- merged_looms$unspliced[, cell_ids]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.2)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)

if(!all(labels(cell.dist) == colnames(emat))){
  stop("cell id not identical.")
}

# Main velocity estimation
message("Estimate velocity")
estimate_opts <- eval(parse(text=sprintf("list(%s)", args$estimate_opts)))
velocity_opts <- eval(parse(text=sprintf("list(%s)", args$velocity_opts)))

estimate_opts$emat <- emat
estimate_opts$nmat <- nmat
estimate_opts$cell.dist <- cell.dist
estimate_opts$n.cores <- args$n.cores
velocity_opts$vel <- do.call(gene.relative.velocity.estimates, args = estimate_opts)

message("Show velocity")
velocity_opts$emb <- emb
velocity_opts$cell.colors <- ac(cell.colors, alpha=0.5)
velocity_opts$n.cores <- args$n.cores
if(args$dump_RDS) saveRDS(velocity_opts, paste0(args$output, '.velocity.RDS'))

pdf(paste0(args$output, '.velocity.pdf'), height=args$height, width=args$width)
p1 <- do.call(show.velocity.on.embedding.cor, args=velocity_opts)
dev.off()
message("Success.")