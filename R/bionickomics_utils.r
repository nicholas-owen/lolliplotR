
#' Title
#'
#' @param gene_of_interest
#'
#' @return
#' @export
#'
#' @examples
get_ClinVar_variants<-function(gene_of_interest){
  clinvarfile<-download_ClinVar()
  dfClinVar<-import_ClinVar(clinvarfile, GeneOfInterest = gene_of_interest)
  dfVars<-tidy_ClinVar(dfClinVar)
  return(dfVars)
}









#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
tidy_ClinVar<-function(df){
  #check format is correct

  #remove any with - positions
  df<-df[df$PositionVCF>0,]

  #split NAME for cDNA and AA changes
  df$variant<-gsub("^.*?:","",df$Name)
  df$variant<-gsub(" .*$", "", df$variant)
  df$cDNA_pos<-gsub("^.*?:","",df$Name)
  df$cDNA_pos<-gsub("c.", replacement = "", x=df$variant)
  df$cDNA_pos<-str_first_number(df$cDNA_pos)
  df$aa_pos<-gsub("^.*?:","",df$Name)
  df$aa_pos<-gsub("^.*\\(", "", df$aa_pos)
  df$aa_pos<-gsub(")", "", df$aa_pos)
  df$aa_change<-df$aa_pos
  df$aa_pos<-gsub("p.","", df$aa_pos)
  df$aa_pos<-str_first_number(df$aa_pos)
  df$freq<-df$NumberSubmitters
  df$published<-1
  df$phenotype<-1
  df$phenotypelist<-df$PhenotypeList
  # NOTE have GRCh37 and GRCh38 builds in there so likely double needed amounts!
  # for now keep only GRCh38!
  df<-df[df$Assembly=="GRCh38",]
  df$Type<-gsub("single nucleotide variant", "SNV", df$Type)
  df$Type<-gsub("Deletion", "Del", df$Type)
  df$Type<-gsub("Duplication", "Dup", df$Type)
  df$Type<-gsub("Insertion", "Ins", df$Type)
  df$Type<-gsub("Microsatellite", "Ms", df$Type)
  df<- df %>%
    dplyr::select(variant, cDNA_pos, aa_change, aa_pos,Type, freq, published, phenotype , phenotypelist)
  return(df)
}







#' Title
#'
#' @return
#' @export
#'
#' @examples
download_ClinVar<-function(){

  dir<-paste0(getwd(), "/")
  file <- sprintf("%sclinvar_%s.txt.gz", dir, Sys.Date())
  if (file.exists(file)){
    message("Latest ClinVar database found..")
  } else {
    message("Downloading latest ClinVar database: ", file)
    download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
                  destfile = file, method = "curl")
  }
  message("ClinVar file:", file)
  return(file)
}

#' Title
#'
#' @param ClinVarFile
#' @param GeneOfInterest
#' @param pathogenic
#'
#' @return
#' @export
#'
#' @examples
import_ClinVar<-function(ClinVarFile, GeneOfInterest, pathogenic=TRUE){
  # import file
  message("Importing ClinVar data..")
  dfClinVar<-read_tsv(ClinVarFile, show_col_types = FALSE)
  # subset on GOI
  dfClinVar<-dfClinVar[dfClinVar$GeneSymbol==GeneOfInterest,]
  if (pathogenic==TRUE){
    dfClinVar<-dfClinVar[dfClinVar$ClinicalSignificance=="Pathogenic",]
  }

  return(dfClinVar)
}



#' Title
#'
#' @param gene_of_interest
#' @param pathogenic
#'
#' @return
#' @export
#'
#' @examples
get_variants_ClinVar<-function(gene_of_interest=GOI, pathogenic=TRUE){
  message("Requesting known disease causing variants from ClinVar..")
  url<-paste0("https://www.ncbi.nlm.nih.gov/clinvar?term=", gene_of_interest,
              "[sym]AND%22clinsig%20pathogenic%22[Properties]")
  url<-paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=", gene_of_interest,"[gene]&retmax=500&retmode=json")
  res<- GET(url)
  # if (res$status_code!=200){
  #   message("Error communicating with ClinVar, please try again later. Error code: ",res$status_code)
  #   stop()
  # } else {
  #   message("Success..")
  # }
  stop_for_status(res)

  rawToChar(res$content)
  variantsGOI <- fromJSON(rawToChar(res$content))
  dfGOIVars<-as.data.frame(variantsGOI[4])
  ## TO DO


}








#' Process variants from HuVar output
#'
#' @param filename path to the HuVar output saved data
#' @param phenotypes_to_keep list object of characters detailing the disease phenotypes to be kept
#'
#' @return dataframe of cleaned and filtered HuVar single nucleotide variants
#' @export
#'
#' @examples
process_HuVarOutput<-function(filename=huvarfile, phenotypes_to_keep="none"){
  if(file.exists(filename)){
    message("Loading HuVar variant file..")
    huVarInfo <- read.csv(filename)
  } else {
    message("\n")
    stop('WARNING: HuVar Variant file does not exist, please check.')
  }
  # clean variant names
  huVarInfo<-subset(huVarInfo, huVarInfo$variant!="-")
  #keep specific phenotypes
  #if (phenotypes_to_keep!="none"){
  huVarInfo<-filter(huVarInfo, grepl(paste(phenotypes_to_keep, collapse='|'), Disease...Tissue..COSMIC.))
  #}
  huVarInfo$cDNA_pos<-gsub("c.", replacement = "", x=huVarInfo$variant)
  huVarInfo$cDNA_pos<-str_first_number(huVarInfo$cDNA_pos)
  huVarInfo$freq<-1
  huVarInfo$phenotype<-ifelse(huVarInfo$Disease...Tissue..COSMIC.=="Leber congenital amaurosis 8", 0, 1)
  huVarInfo$published<-1

  return(huVarInfo)
  }



#' Get variants from NIH Clinical Tables
#'
#' @param gene_name character vector of HGNC symbol gene name
#'
#' @return dataframe of cleaned and filtered disease causing variants for the gene of interest
#' @export
#'
#' @examples
get_gene_variants<-function(gene_name=GOI){
  # get variants for GOI from https://clinicaltables.nlm.nih.gov/apidoc/variants/v4/doc.html#params

  # https://clinicaltables.nlm.nih.gov/demo.html?db=variants
  message("Requesting known disease causing variants from ClinicalTables..")
  url<-paste0("https://clinicaltables.nlm.nih.gov/api/variants/v4/search?terms=", GOI,
              "&sf=GeneSymbol&maxList=500&df=AlternateAllele,AlleleID,AminoAcidChange,Chromosome,ChromosomeAccession,Cytogenetic,dbSNP,GeneID,GeneSymbol,GenomicLocation,hgnc_id,hgnc_id_num,HGVS_c,HGVS_exprs,HGVS_p,Name,NucleotideChange,phenotypes,phenotype,PhenotypeIDS,PhenotypeList,RefSeqID,ReferenceAllele,Start,Stop,Type,VariationID")
  res<- GET(url)

  rawToChar(res$content)
  variantsGOI <- fromJSON(rawToChar(res$content))
  dfGOIVars<-as.data.frame(variantsGOI[4])
  names(dfGOIVars)<-c("AlternateAllele","AlleleID","AminoAcidChange","Chromosome","ChromosomeAccession","Cytogenetic","dbSNP","GeneID","GeneSymbol","GenomicLocation","hgnc_id","hgnc_id_num","HGVS_c","HGVS_exprs","HGVS_p","Name","NucleotideChange","phenotypes","phenotype","PhenotypeIDS","PhenotypeList","RefSeqID","ReferenceAllele","Start","Stop","Type","VariationID")

  dfGOIVars<-subset(dfGOIVars, dfGOIVars$Type=="single nucleotide variant")
  #subset X17 based on containing c.
  dfGOIVars<-dfGOIVars[ grep("c", dfGOIVars$NucleotideChange),]
  # get cDNA_pos
  dfGOIVars$cDNA_pos<-gsub("c.", replacement = "", x=dfGOIVars$NucleotideChange)
  # remove any snps starting with symbol, ie neg, or *
  dfGOIVars$cDNA_pos[!grepl("^[[:digit:]]+",dfGOIVars$cDNA_pos)] <-NA
  dfGOIVars<-subset(dfGOIVars, !is.na(dfGOIVars$cDNA_pos))

  dfGOIVarsTidy<-data.frame(variant= dfGOIVars$NucleotideChange,
                            cDNA_pos=str_first_number(dfGOIVars$cDNA_pos),
                            aa_change=dfGOIVars$AminoAcidChange,
                            aa_pos=str_first_number(dfGOIVars$AminoAcidChange),
                            freq=1,
                            published=1,
                            phenotype=sample(c("0", "1"), nrow(dfGOIVars), replace=TRUE)
  )
  message("Indentified ", dim(dfGOIVarsTidy)[1], " disease causing single nucleotide variants..")

  return(dfGOIVarsTidy)
}







#' Rescale variants based upon new intronic sizes (fixed)
#'
#' @param variantInfo `dataframe` of variant information for gene of interest
#' @param transcript_grange `GRanges` object of the transcript model for the gene of interest
#' @param minIntronSize `numeric` vector to fix all `introns` to for representation
#'
#' @return
#' @export
#'
#' @examples
rescale_variants<-function(variantInfo, transcript_grange, minIntronSize){
  #exon totals
  exonTotals<-list()
  count<-0
  for (i in 1:length(transcript_grange)){
    exonTotals[i]<-width(transcript_grange)[i] + count
    count=count + width(transcript_grange)[i]
  }

  if (minIntronSize!="no"){

  for (j in seq(dim(variantInfo)[1])){
    SNPtoRecal<-variantInfo$cDNA_pos[j]
    #check SNP in which exon
    exonCount=1
    while (SNPtoRecal > exonTotals[exonCount]) {
      exonCount = exonCount+1
    }
    # print (exonCount)
    SNPtoRecal<-SNPtoRecal+((exonCount-1)*minIntronSize)
    variantInfo$cDNA_pos[j]<-SNPtoRecal
  }
  }

  variantGR <- GRanges(unique(seqnames(transcript_grange)), IRanges(variantInfo$cDNA_pos, width=1, score=variantInfo$freq,SNPsideID=ifelse(variantInfo$phenotype=="1", "top", "bottom"),names=variantInfo$variant, Type=variantInfo$Type))
  #variantGR$color <- ifelse(variantInfo$published=="1", "red", "black")
  #variantGR$label <- as.character(1:length(variantGR))
  #variantGR$alpha <- sample(100:255, nrow(variantInfo), replace = TRUE)/255
  #variantGR$label.col <- ifelse(variantGR$alpha>0.5, "white", "black")
  # sample.gr$SNPsideID <- sample(c("top", "bottom"),
  #                               length(sample.gr),
  #                               replace=TRUE)
  return(variantGR)
}









#' Rescale introns
#'
#' @param grange `GRanges` object of the transcript model for the gene of interest
#' @param minIntronSize `numeric` vector to fix all `introns` to for representation
#'
#' @return rescaled `GRanges` object with updated `intron` sizes
#' @export
#'
#' @examples
rescale_introns<-function(grange, minIntronSize){
  #granges_reduce_introns(grange, minIntronSize)
  message("Moving start to zero based positions..")
  CDSstart<-start(grange)[1]
  start(grange)<-start(grange)-CDSstart
  end(grange)<-end(grange)-CDSstart

  if (minIntronSize!="no"){
    message("Rescaling introns to size: ", minIntronSize)

    #exon totals
    exonTotals<-list()
    count<-0
    for (i in 1:length(grange)){
      exonTotals[i]<-width(grange)[i] + count
      count=count + width(grange)[i]
    }

    for (i in 2:length(grange)[1]){
      width<-width(grange)[i]
      start(grange)[i]<-end(grange)[i-1]+minIntronSize
      end(grange)[i]<-start(grange)[i]+width-1
    }
    }  else {
      message("Maintaining introns at original sizes.")
    }


  return(grange)

}


###
# create variant track
###
create_variant_track<-function(dfVariants, transcript_structure){
  # total exon length
  exonTotals<-list()
  count<-0
  for (i in 1:dim(transcript_structure)[1]){
    exonTotals[i]<-width(new)[i] + count
    count=count + width(new)[i]

  }





}



#' Get EnsemblID and Canonical transcript ID
#'
#' @param gene_of_interest `character` vector of the HGNC symbol of the gene of interest
#'
#' @return `transcript` and `exon` models
#' @export
#'
#' @examples
get_EnsemblIDs<-function(gene_of_interest=GOI){
  # support human for now
  message("Identifying EnsemblID for gene: ", gene_of_interest)
  cachedTxDataFile<-paste0(gene_of_interest, "_EnsemblTxData.gz")
  cachedExDataFile<-paste0(gene_of_interest, "_EnsemblExData.gz")

  if (file.exists(cachedTxDataFile)){
    message("Cached EnsemblTxData file found.. Loading..")
    idconvert<-read_csv(cachedTxDataFile, show_col_types = FALSE)
  } else {
    message("Contacting Ensembl, please be patient.")
    ensembl <- tryCatch({
      useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    },  error = function(e) {
      useEnsembl(biomart = "genes",
                 dataset = "hsapiens_gene_ensembl"
      )
    })
    #ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    message("Querying Ensembl database..")
      idconvert<-getBM(attributes = c('entrezgene_id', 'external_gene_name', "ensembl_gene_id","chromosome_name", "start_position","end_position","strand","hgnc_symbol","ucsc","band", "ensembl_transcript_id", "transcript_is_canonical", "ensembl_exon_id"),
                       filters = 'external_gene_name',
                       values = gene_of_interest,
                       mart = ensembl)
    write_csv(file = cachedTxDataFile, x = idconvert)
    message("Saved Ensembl Tx data as cache.")
  }

  #remove non canonical chromosome assemblies
  idconvert<-subset(idconvert, nchar(as.character(chromosome_name)) <= 3)
  # identify canonical transcript
  idconvert<-subset(idconvert, transcript_is_canonical=="1")
  idconvert<- idconvert %>%
    dplyr::select(-entrezgene_id) %>%
    unique()
  ensemblTXID<-unique(idconvert$ensembl_transcript_id)
  message("Canonical transcript: ", ensemblTXID)
  # get exon info for GOI
  if (file.exists(cachedExDataFile)){
    message("Cached EnsemblExData file found.. Loading..")
    exon_info<-read_csv(cachedExDataFile, show_col_types = FALSE)
  } else {
    message("Querying Ensembl Ex database..")
    exon_info<-getBM(attributes = c("external_gene_name",	"ensembl_exon_id","ensembl_transcript_id","transcript_count", "rank" ,"strand","chromosome_name", "exon_chrom_start", "exon_chrom_end","cds_start", "cds_end", "cds_length"),
                     filters = 'ensembl_transcript_id',
                     values = ensemblTXID,
                     mart = ensembl)
    write_csv(file = cachedExDataFile, x = exon_info)
    message("Saved Ensembl Ex data as cache.")
  }

  message("Number of exons: ", dim(exon_info)[1])
  #check for untranslated exons
  exon_info<-subset(exon_info, nchar(as.character(chromosome_name)) <= 3)
  if (sum(is.na(exon_info$cds_start)) > 0){
    message(sum(is.na(exon_info$cds_start)), " untranslated exon(s) identified and removed." )
    exon_info<-subset(exon_info, !is.na(cds_start))
  }
  # sort output by exon rank
  exon_info<- exon_info %>%
    arrange(rank)

  return(exon_info)
}



#' Title
#'
#' @param transcript_of_interest
#'
#' @return
#' @export
#'
#' @examples
get_Ensembl_peptide<-function(gene_of_interest){
  # support human for now
  message("Identifying Ensembl protein information for gene: ", gene_of_interest)
  cachedPrDataFile<-paste0(gene_of_interest, "_EnsemblPrData.gz")

  if (file.exists(cachedPrDataFile)){
    message("Cached EnsemblPrData file found.. Loading..")
    ens_protein<-read_csv(cachedPrDataFile, show_col_types = FALSE)
  } else {
    message("Contacting Ensembl, please be patient.")
    ensembl <- tryCatch({
      useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    },  error = function(e) {
      useEnsembl(biomart = "genes",
                 dataset = "hsapiens_gene_ensembl"
      )
    })
    message("Querying Ensembl database..")
    ens_gene<-getBM(attributes = c("external_gene_name", "ensembl_gene_id","chromosome_name", "hgnc_symbol", "ensembl_transcript_id", "transcript_is_canonical"),
                     filters = 'external_gene_name',
                     values = gene_of_interest,
                     mart = ensembl)
    # identify canonical transcript
    ens_gene<-subset(ens_gene, transcript_is_canonical=="1")
    ensemblTXID<-unique(ens_gene$ensembl_transcript_id)
    # get peptide info
    ens_protein<-getBM(attributes = c('external_gene_name', "ensembl_gene_id", "ensembl_transcript_id", "transcript_is_canonical","chromosome_name",
                                     "interpro", "interpro_short_description", "interpro_start", "interpro_end","ensembl_peptide_id"
                                     #, "smart", "smart_start", "smart_end"
                                   ),
                     filters = 'ensembl_transcript_id',
                     values = ensemblTXID,
                     mart = ensembl)
    ensemblPEPID<-unique(ens_protein$ensembl_peptide_id)
    message("Canonical protein: ", ensemblPEPID)
    ensemblPEPseq<-getSequence(id = ensemblPEPID, type = "ensembl_peptide_id", seqType = "peptide", mart = ensembl)
    ensemblPEPlength<-nchar(ensemblPEPseq)[1]
    ens_protein<-rbind(ens_protein[1,], ens_protein)
    ens_protein$interpro[1]<-"peptide"
    ens_protein$interpro_short_description[1]<-"peptide"
    ens_protein$interpro_start[1]<-"1"
    ens_protein$interpro_end[1]<-ensemblPEPlength
    write_csv(file = cachedPrDataFile, x = ens_protein)
    message("Saved Ensembl protein data as cache.")
  }

  # sort output by interpro feature start
  ens_protein<- ens_protein %>%
    arrange(interpro_start)

    return(ens_protein)
}








#' Read variant data
#'
#' @param varFile `path` to local variant file
#' @param do_count default of `FALSE`, function to count the frequency of observations of each variant
#'
#' @return variant `dataframe`
#' @export
#'
#' @examples
read_variant_data<-function(varFile, do_count=FALSE){
  if(file.exists(varFile)){
    message("Loading variant file..")
    dfVariant <- read.csv(varFile)
  } else {
    message("\n")
    stop('WARNING: Variant file does not exist, please check.')
  }
  # clean frequencies
  dfVariant$cDNA_pos_original<-dfVariant$cDNA_pos
  message("Checking variant file format..")
  stopifnot("`variant` must be a character." = is.character(dfVariant$variant))
  stopifnot("`cDNA_pos` must be a numeric" = is.numeric(dfVariant$cDNA_pos))
  stopifnot("`freq` must be a numeric" = is.numeric(dfVariant$freq))
  stopifnot("`phenotype` must be a binary" = all(dfVariant$phenotype %in% 0:1))
  stopifnot("`published` must be a binary" = all(dfVariant$published %in% 0:1))

  message("Scoring variant frequencies..")
  # get scoring of frequency
  if (do_count=="TRUE"){
    dfVariant<-do.call(rbind,
                       by(dfVariant,
                          dfVariant,
                          function(x){
                            data.frame(unique(x),"Score"=nrow(x))
                          }
                       )
    )
  } else {
    message(" - skipped")
  }
  return(dfVariant)
}





#' Title
#' https://stackoverflow.com/questions/4350440/split-data-frame-string-column-into-multiple-columns#4351918
#' @param column
#' @param pattern
#' @param into_prefix
#'
#' @return
#' @export
#'
#' @examples
split_into_multiple <- function(column, pattern = ", ", into_prefix){
  cols <- str_split_fixed(column, pattern, n = Inf)
  # Sub out the ""'s returned by filling the matrix to the right, with NAs which are useful
  cols[which(cols == "")] <- NA
  cols <- as.tibble(cols)
  # name the 'cols' tibble as 'into_prefix_1', 'into_prefix_2', ..., 'into_prefix_m'
  # where m = # columns of 'cols'
  m <- dim(cols)[2]

  names(cols) <- paste(into_prefix, 1:m, sep = "_")
  return(cols)
}



#
#
# after <- dfClinVar %>%
#   bind_cols(split_into_multiple(.$PhenotypeList, "\\|", "PhenotypeList"))# %>%
#   # selecting those that start with 'type_' will remove the original 'type' column
#   #dplyr::select(attr, starts_with("PhenotypeList_"))
#
# #remove notprovided
# is.na(after) <- after == "not provided"
