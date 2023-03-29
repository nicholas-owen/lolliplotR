
#' Create lolliplot gene plot
#'
#' @param gene_name `character` vector of the HGNC symbol gene of interest
#' @param variant_filename provided `csv` formatted variant file or `dataframe` of variants
#' @param variant_input_type `file` or `dataframe` depending on input provided
#' @param varPalette
#' @param color_by_data
#' @param min_intron_size
#' @param exon_height `numeric` representing the draw height of `exon` features
#'
#' @return output to gene symbol prefixed `SVG` formatted plots, and printed to local graphics device
#' @export
#'
#' @examples
create_gene_lolliplot<-function(gene_name=GOI, variant_filename=varFile, variant_input_type=varInputType, varPalette=varPalette, color_by_data="published", exon_height=exonHeight, min_intron_size=minIntronSize){
  stopifnot("`gene_name` must be a character." = is.character(gene_name))
  #stopifnot("`minIntronSize` must be a numeric" = is.numeric(minIntronSize))

  if (variant_input_type=="file"){
    #read variant file
    dfVariants<-read_variant_data(varFile = varFile, do_count=FALSE)
  } else {
    message("Taking variant input from dataframe..")
    dfVariants<-variant_filename
  }
  goiInfo<-get_EnsemblIDs(gene_name)
  goiInfo$strand <- ifelse(goiInfo$strand=="1", "+", "-")
  goiInfo$height<-exonHeight
  goiInfo$names <- paste0("exon ",goiInfo$rank)
  colPal<-RColorBrewer::brewer.pal(dim(goiInfo)[1], exonPalette)
  goiInfo$fill<-colorRampPalette(colPal)(dim(goiInfo)[1])
  goiGR<-GenomicRanges::makeGRangesFromDataFrame(goiInfo,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field=c("seqnames", "seqname",
                                                   "chromosome", "chrom",
                                                   "chr", "chromosome_name",
                                                   "seqid"),
                                  start.field="exon_chrom_start",
                                  end.field="exon_chrom_end",
                                  strand.field="strand",
                                  starts.in.df.are.0based=FALSE)
  names(goiGR)<-goiInfo$names
  # rescale introns
  goiGR<-rescale_introns(goiGR, min_intron_size)
  variantGR<-rescale_variants(dfVariants, goiGR, min_intron_size)
  variantGR$Type<-as.factor(variantGR$Type)
  # legends <- list(labels=c("published", "novel"),
  #                fill=c("black", "white"),
  #                color=c("gray80", "gray80"))

  #recolor variants by published/freq/vartype
  # varColPal<-RColorBrewer::brewer.pal(dim(goiInfo)[1], varPalette)
  # goiInfo$fill<-colorRampPalette(varColPal)(dim(goiInfo)[1])


  if (color_by_data=="published"){
    varColPal<-RColorBrewer::brewer.pal(2, varPalette)
    variantGR$color <- ifelse(dfVariants$published=="1", "white", varColPal[3])
    #message("Colors",varColPal[3] )
    legends <- list(labels=c("published", "novel"),
                    fill= c("white", varColPal[3]), # goiInfo$fill,
                    color=c("gray80", "gray80"))

  }
  if (color_by_data=="type"){
    varColPal<-RColorBrewer::brewer.pal(nlevels(variantGR$Type), varPalette)
    variantGR$color <- varColPal[variantGR$Type]
    legends <- list(labels=levels(variantGR$Type),
                    fill=varColPal,
                    color=rep("gray80", nlevels(variantGR$Type)))
  }

  pnas.height.inches<- 225/25.4
  pnas.width.single.column.inches<-87/25.4
  pnas.width.double.column.inches<-178/25.4
  pnas.width.onehalf.column.inches<-114/25.4
  output_filename<-paste0(goiInfo$external_gene_name[1], "_",
                          goiInfo$ensembl_transcript_id[1], "_transcript_varmap.svg")
  message("Saving to: ", output_filename)
  svg(output_filename, width = pnas.width.double.column.inches*3, height = 2*pnas.height.inches/1.1)
  trackViewer::lolliplot(SNP.gr = variantGR, features = goiGR, ylab = "Variant Frequency", legend = legends)
  dev.off()
  trackViewer::lolliplot(SNP.gr = variantGR, features = goiGR, ylab = "Variant Frequency", legend = legends)
}


