library(magrittr)
library(GenomicRanges)
library(EnsDb.Hsapiens.v75)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)

load_tfbsConsSites <- function(tfbsConsSitesFile, tfbsConsFactorsFile){
  tfbs_anno <- gzfile(tfbsConsFactorsFile) %>%
    read.delim(col.names = c('feature', 'ac', 'species', 'group', 'id'))

  tfbs <- gzfile(tfbsConsSitesFile) %>%
    read.delim(
      col.names=c('bin', 'chromosome', 'start', 'end', 'feature', 'score', 'strand', 'zScore') ) %>%
    dplyr::left_join(tfbs_anno, by='feature') %>%
    dplyr::filter(species == 'human') %>%
    GRanges()
}

load_RegulatoryFeatureActivity <- function(regulatory_activity_gff){
  lung_regulatory <- rtracklayer::import.gff3(regulatory_activity_gff)
  idx <- mcols(lung_regulatory)$activity != 'NA'
  lung_regulatory <- lung_regulatory[idx,]
  lung_regulatory <- keepSeqlevels(lung_regulatory, seqlevels(lung_regulatory)[1:24])
  seqlevels(lung_regulatory) <- paste0('chr', seqlevels(lung_regulatory))
  mcols(lung_regulatory)$id <- mcols(lung_regulatory)$type
  mcols(lung_regulatory)$group <- mcols(lung_regulatory)$activity
  lung_regulatory
}

load_humantfs_PWMs <- function(PWM_dir){
  pwm_files <- dir(PWM_dir)
  PWMs <- lapply(pwm_files, function(x){
    read.delim(file.path(PWM_dir, x)) %>%
      dplyr::select(-Pos) %>%
      as.matrix() %>%
      t()
  })
  names(PWMs) <- gsub('.txt', '', pwm_files)
  PWMs
}


load_humantfs_TF_Motif <- function(Human_TF_MotifList_txt_file){
  pick_best <-function(x, identity, is_best){
    if(length(x) == 1){
      x
    }else if(all(is.na(is_best))){
      x[head(which.max(identity), n=1)]
    }else{
      x[head(which(is_best == TRUE), n=1)]
    }
  }

  read.delim(Human_TF_MotifList_txt_file) %>%
    dplyr::rename(group=HGNC.symbol, Ensembl=Ensembl.ID, id=Motif.ID) %>%
    dplyr::mutate(
      identity = as.numeric(gsub('.*\\(|%.*', '', Motif.evidence)),
      identity = ifelse(Motif.evidence == 'Direct', 100, identity)
    ) %>%
    dplyr::group_by(CIS.BP.ID) %>%
    dplyr::summarise(
      id=pick_best(id, identity, Best.Motif.s....Figure.2A.),
      group=pick_best(group, identity, Best.Motif.s....Figure.2A.),
      Ensembl=pick_best(Ensembl, identity, Best.Motif.s....Figure.2A.) )
}


load_rnb <- function(rnb_RData){
  load(rnb_RData)
  CpG <- annotation(object, type='sites', add.names=TRUE) %>%
    dplyr::select(Chromosome, Start, End, Strand) %>%
    GRanges()
  mcols(CpG) <- meth(object, row.names=F)
  CpG
}


load_sample_groups <- function(sample_csv){
  samples_info <- read.csv(sample_csv) %>%
    dplyr::select(TCGA_ID, Type)
  sample_groups <- samples_info$Type
  names(sample_groups) <- samples_info$TCGA_ID
  sample_groups
}


load_site_diff_meth <- function(diffMethTable_site_csv_gz){
  gzfile(diffMethTable_site_csv_gz) %>%
    read.csv() %>%
    dplyr::select(cgid, diffmeth.p.adj.fdr, combinedRank) %>%
    dplyr::mutate(
      log10FDR = -log10(diffmeth.p.adj.fdr),
      log10Rank = -log10(combinedRank / max(combinedRank, na.rm=T))
      ) %>%
    tibble::column_to_rownames('cgid')
}


load_genomic_range <- function(x){
  read.csv(x, row.names = 1) %>%
    GRanges()
}

load_ranges <- function(meth_dir, pattern='.csv'){
  gr <- meth_dir %>%
    dir(pattern = pattern, full.names = T) %>%
    lapply(load_genomic_range)
  names(gr) <- tools::file_path_sans_ext(dir(meth_dir, pattern='.csv'))
  gr
}


extend <- function(x, width){
  start(x) <- start(x) - width
  end(x) <- end(x) + width
  x
}

prepare_gene_model <- function(gr, edb, byID=F, self=F){
  if(self) return(gr)

  if(byID){
    gene_model <- transcripts(
      edb, columns=c('symbol'),
      filter = GeneIdFilter(as.character(gr$id)) )
    return(gene_model)
  }

  chr <- as.character(seqnames(gr))
  repeat{
    gene_model <- try(getGeneRegionTrackForGviz(
      edb, chromosome = chr,
      start = start(gr), end = end(gr)) )
    if(class(gene_model) == 'GRanges'){
      break
    }else{
      gr <- extend(gr, 1000)
      message('Extend 1k on both side. ', width(gr))
    }
  }
  gene_model
}

prepare_query_region <- function(region_list){
  region <- Reduce(c, region_list)
  GRanges(seqnames = seqnames(region),
          IRanges(start = min(start(region)),
                  end = max(end(region))))
}

get_CpG_tracks <- function(CpG_region, track_list=c('beta', 'log10FDR', 'log10Rank'), groups = NULL){
  tracks <- list()
  strand(CpG_region) <- '*'

  if('beta' %in% track_list){
    tracks[['beta']] <- DataTrack(
      CpG_region,
      name = 'beta',
      type = c('confint', 'a', 'p'),
      groups = groups)
  }

  if('log10FDR' %in% track_list){
    mcols(CpG_region) <- diff_meth[names(CpG_region), 'log10FDR']
    tracks[['log10FDR']] <- DataTrack(CpG_region, name = '-log10(P)', type='p')
  }

  if('log10Rank' %in% track_list){
    mcols(CpG_region) <- diff_meth[names(CpG_region), 'log10Rank']
    tracks[['log10Rank']] <- DataTrack(CpG_region, name = '-log10(Rank)', type='p')
  }

  tracks[['CpG']] <- AnnotationTrack(
    CpG_region,
    showFeatureId = F,
    name='CpG',
    fill = 'black',
    stacking = "dense",
    shape = 'box')

  tracks
}

get_regulatory_track <- function(regulatory_region, name){
  colormap <- c(ACTIVE = 'darkgreen', INACTIVE = 'darkgrey',
                POISED = 'darkblue', REPRESSED = 'darkred')
  AnnotationTrack(
    regulatory_region, name=name,
    featureAnnotation = "id",
    groupAnnotation = 'group',
    stacking = 'squish',
    mergeGroups = F,
    showOverplotting = F,
    showFeatureId = TRUE,
    cex.feature = 0.5,
    fill = colormap[regulatory_region$group],
    shape = "box")
}

prepare_humantf <- function(motif_region, min.score="90%", cores=parallel::detectCores()){
  query_region_seq <- getSeq(genome, motif_region)[[1]]
  cl <- makeCluster(cores)
  clusterEvalQ(cl, {
    library(Biostrings)
    library(BSgenome.Hsapiens.UCSC.hg19)
  })
  clusterExport(cl, c('PWMs', 'query_region_seq', 'motif_region', 'TF_Motif'),
                envir = environment())
  tfbs_predict <- parallel::parLapply(cl, names(PWMs), function(x){
    hits <- matchPWM(PWMs[[x]], query_region_seq, min.score=min.score)
    if(length(hits) > 0){
      hits <- GRanges(seqnames(motif_region), ranges(hits))
      end(hits) <- end(hits) + start(motif_region)
      start(hits) <- start(hits) + start(motif_region)
      idx <- TF_Motif$CIS.BP.ID == x
      mcols(hits) <- TF_Motif[idx,]
      hits
    }else{
      NULL
    }
  })
  stopCluster(cl)
  tfbs_predict <- tfbs_predict[!sapply(tfbs_predict, is.null)]
  Reduce(c, tfbs_predict)
}

plotRegion <- function(
  gr,
  track_list=c('gene', 'tfbs', 'humantf', 'activity', 'beta', 'log10FDR', 'log10Rank'),
  byID = F, self = F, activity.organ = c('lung', 'brain'),
  humantf.CpG_only = T, humantf.min.score = "90%", humantf.cores = parallel::detectCores(), ... ){

  tracks <- list()
  regions <- list()

  chr <- as.character(seqnames(gr))
  regions[['hightlight']] <- highlight <- gr
  regions[['query']] <- gr

  if(width(gr) <= 500){
    regions[['query']] <- gr <- extend(gr, 1000)
  }

  itrack <- IdeogramTrack(genome = 'hg19', chromosome = chr)
  tracks[['GenomeAxis']] <- GenomeAxisTrack()

  if('gene' %in% track_list){
    regions[['gene']] <- gene_model <- prepare_gene_model(gr, edb, byID, self)
    tracks[['gene']] <- GeneRegionTrack(gene_model, name='gene model')
  }

  query_region <- prepare_query_region(regions)

  if(any(c('beta', 'log10FDR', 'log10Rank') %in% track_list)){
    regions[['CpG']] <- CpG_region <- subsetByOverlaps(CpG, query_region)
    groups <- sample_groups[colnames(mcols(CpG_region))]
    tracks <- c(tracks, get_CpG_tracks(regions[['CpG']], track_list, groups))
  }

  if('tfbs'%in% track_list){
    regions[['tfbs']] <- TFBS_region <- subsetByOverlaps(tfbsConsSites, CpG_region)
    if(length(TFBS_region) > 0){
      tracks[['tfbs']] <- AnnotationTrack(
        TFBS_region, name='TFBS',
        groupAnnotation = 'group',
        mergeGroups = F,
        showOverplotting = TRUE,
        showFeatureId = FALSE,
        shape = "box")
    }
  }

  if('activity' %in% track_list){
    regulatory_region <- lapply(
      regulatory_activity[activity.organ], function(x){
        subsetByOverlaps(x, CpG_region)})
    idx <- sapply(regulatory_region, length) > 0
    regions[['regulatory_activity']] <- regulatory_region <- regulatory_region[idx]

    regulatory_tracks <- lapply(
      names(regulatory_region), function(x){
        get_regulatory_track(regulatory_region[[x]], x)
        })
    tracks <- c(tracks, regulatory_tracks)
  }

  if('humantf' %in% track_list){
    regions[['humantf']] <- predict_humantf_region <- prepare_humantf(
      motif_region = extend(highlight, 20),
      min.score = humantf.min.score,
      cores = humantf.cores)
    if(humantf.CpG_only){
      predict_humantf_region <- subsetByOverlaps(predict_humantf_region, CpG_region)
    }
    if(length(predict_humantf_region) > 0){
      tracks[['humantf']] <- predict_humantf_track <- AnnotationTrack(
        predict_humantf_region, name='humantf(predict)',
        featureAnnotation = "id",
        groupAnnotation = 'group',
        stacking = 'squish',
        mergeGroups = F,
        showOverplotting = T,
        showFeatureId = TRUE,
        cex.feature = 0.5,
        shape = "box")
    }
  }

  extend_size <- as.integer(0.01 * max(width(highlight)))
  tracks <- tracks[c("GenomeAxis", track_list, 'CpG')]
  idx <- !sapply(tracks, is.null)
  hltrack <- HighlightTrack(
    trackList = tracks[idx],
    start = min(start(highlight)) - extend_size,
    end = max(end(highlight)) + extend_size,
    chromosome = chr
    )

  main <- ifelse(is.null(gr$symbol) || is.na(gr$symbol), as.character(gr$id), as.character(gr$symbol))
  plotTracks(
    list(itrack, hltrack),
    transcriptAnnotation = "transcripts",
    #featureAnnotation = "id",
    fontcolor.feature = "lightgrey",
    extend.left = 0.05, extend.right = 0.05,
    col = NULL, main = main, cex.main=1.2,
    ...)

  mcols(regions[['CpG']]) <- diff_meth[names(CpG_region), ]
  idx <- sapply(regions, length) > 0
  invisible(regions[idx])
}
