
#' Create lolliplot peptide plot
#'
#' @param gene_name `character` vector of the HGNC symbol gene of interest
#' @param variant_filename provided `csv` formatted variant file or `dataframe` of variants
#' @param variant_input_type `file` or `dataframe` depending on input provided
#' @param varPalette
#' @param color_by_data
#' @param exon_height `numeric` representing the draw height of `exon` features
#'
#' @return output to gene symbol prefixed `SVG` formatted plots, and printed to local graphics device
#' @export
#'
#' @examples
#'
create_peptide_lolliplot<-function(gene_name=GOI, variant_filename=varFile, variant_input_type=varInputType, varPalette=varPalette, color_by_data="published", exon_height=exonHeight){
  stopifnot("`gene_name` must be a character." = is.character(gene_name))
  if (variant_input_type=="file"){
    #read variant file
    dfVariants<-read_variant_data(varFile = varFile, do_count=FALSE)
  } else {
    message("Taking variant input from dataframe..")
    dfVariants<-variant_filename
  }
  goiInfo<-get_Ensembl_peptide(gene_name)
  goiInfo$height<-exonHeight
  goiInfo$names <- as.factor(goiInfo$interpro_short_description)
  colPal<-RColorBrewer::brewer.pal(nlevels(goiInfo$names), exonPalette)
  # scales::show_col(colPal)
  goiInfo$fill<-colorRampPalette(colPal)(nlevels(goiInfo$names))[goiInfo$names]
  goiGR<-GenomicRanges::makeGRangesFromDataFrame(goiInfo,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=TRUE,
                                  seqinfo=NULL,
                                  seqnames.field=c("seqnames", "seqname",
                                                   "chromosome", "chrom",
                                                   "chr", "chromosome_name",
                                                   "seqid"),
                                  start.field="interpro_start",
                                  end.field="interpro_end",
                                  #strand.field="strand",
                                  starts.in.df.are.0based=FALSE)
  names(goiGR)<-goiInfo$names
  goiGR$alpha<-0.5
  variantGR <- GenomicRanges::GRanges(unique(seqnames(goiGR)), IRanges(dfVariants$aa_pos, width=1, score=dfVariants$freq,SNPsideID=ifelse(dfVariants$phenotype=="1", "top", "bottom"),names=dfVariants$aa_change, Type=dfVariants$Type))
  variantGR$Type<-as.factor(variantGR$Type)

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

  xaxis <- c(1, goiInfo$interpro_end[1])

  pnas.height.inches<- 225/25.4
  pnas.width.single.column.inches<-87/25.4
  pnas.width.double.column.inches<-178/25.4
  pnas.width.onehalf.column.inches<-114/25.4
  output_filename<-paste0(goiInfo$external_gene_name[1], "_",
                          goiInfo$ensembl_peptide_id[1], "_peptide_varmap.svg")
  message("Saving to: ", output_filename)
  svg(output_filename, width = pnas.width.double.column.inches*3, height = 2*pnas.height.inches/1.1)
  trackViewer::lolliplot(SNP.gr = variantGR, features = goiGR, ylab = "Variant Frequency", legend = legends, xaxis = xaxis)
  dev.off()
  trackViewer::lolliplot(SNP.gr = variantGR, features = goiGR, ylab = "Variant Frequency", legend = legends, xaxis = xaxis)
}

