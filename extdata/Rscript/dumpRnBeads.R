suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("RnBeads"))
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("--RData", required=T,
                    help = "The path to RnBeads RData file")
parser$add_argument("--regions", nargs='*',
                    help = "Choose from [sites, tiling, genes, promoters, cpgislands]")
parser$add_argument("--outdir", required=T,
                    help = "The path to save RnBeads meth matrix and annotation")
args <- parser$parse_args()

message('Loading RData')
load(args$RData)

if(is.null(args$regions)) args$regions <- c('sites', names(object@regions))

for(region in args$regions){
  message("Proccessing ", region)
  m <- as.data.frame(meth(object, type=region, row.names=TRUE))
  id <- row.names(m)
  m <- cbind(id, id, m)
  names(m)[1:2] <- ''

  write.table(
    m,
    sep='\t', row.names = F,
    file=file.path(args$outdir, paste0(region, '.meth.tsv'))
  )
  
  annotation(object, type=region) %>%
    write.table(file=file.path(args$outdir, paste0(region, '.meth.anno.tsv')), sep='\t')
}