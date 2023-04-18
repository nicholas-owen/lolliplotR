#' Gene Variant Lolliplot
#' ~~~~~~~~~~~~~~~~~~~~~~
#'
#' Author: Nicholas Owen
#' Github: github.com/nicholas-owen
#' Email:
#' Version: 0.001a
#' Date: 2023-03-24
#'
#' Description:
#'
#'
#'

# load dependancies
library(bionickomics)
options(warn = -1)
load_pack(biomaRt)
load_pack(GenomicRanges)
load_pack(GenomicFeatures)
load_pack(trackViewer)
load_pack(RColorBrewer)
load_pack(httr)
load_pack(jsonlite)
load_pack(tibble)
load_pack(strex) # https://stackoverflow.com/questions/14003568/extract-first-number-from-string
load_pack(readr)
source("./bionickomics_utils.r") # to be included in bionickomics
cat("\014") # clear console



# Set gene of interest as HGNC symbol
GOI<-"CRB1"

# Set drawing parameters
minIntronSize<-500 # to rescale the intronic regions to a smaller representation, set to "no" for keeping introns original sizes -broken, needs intron sizes not removed
exonHeight<-0.075 # numeric representing the size of the boxes for exonic region
exonPalette<-"Spectral" # RColorBrewer palette from which all exons will be filled from
varPalette<-"Set3"
varPhenotype<-c("ret","leber")

# Examples

# 1. Using ClinVar to obtain all pathogenic variants reported for a gene of interest (GOI)
# 1.1
dfVars<-get_ClinVar_variants(gene_of_interest = GOI)

# 1.2.1
create_gene_lolliplot(gene_name = GOI, dfVars, variant_input_type="df", exonHeight, minIntronSize, varPalette=varPalette, color_by_data="type")

#1.2.2
create_gene_lolliplot(gene_name = GOI, dfVars, variant_input_type="df", exonHeight, minIntronSize, varPalette=varPalette, color_by_data="published")

#1.3.1
create_peptide_lolliplot(gene_name = GOI, dfVars, variant_input_type="df", exonHeight, varPalette=varPalette, color_by_data="type")

#1.3.2
create_peptide_lolliplot(gene_name = GOI, dfVars, variant_input_type="df", exonHeight, varPalette=varPalette, color_by_data="published")











# 2. Using local variants file for gene of interest (GOI)
# 2.1 HuVarDB variant info
varFile<-"D:/Box Sync/Scripts/R/General/CRB1_HuVarDB.csv"
create_gene_lolliplot(gene_name = GOI, varFile, variant_input_type="file", exonHeight, minIntronSize, varPalette=varPalette, color_by_data="type")

# 2.2.i Using local variants file for gene of interest (GOI), colored by variant novelty
varFile<-"./variant_info.csv"
create_gene_lolliplot(gene_name = GOI, varFile, variant_input_type="file", exonHeight, minIntronSize, varPalette=varPalette, color_by_data="published")
# 2.2.ii Using local variants file for gene of interest (GOI), colored by variant type
create_gene_lolliplot(gene_name = GOI, varFile, variant_input_type="file", exonHeight, minIntronSize, varPalette=varPalette, color_by_data="type")

# 3. Download variants from Clinical Tables (NIH)
variantInfoDL<-get_gene_variants(gene_name = GOI)
create_gene_lolliplot(gene_name = GOI, variantInfoDL, variant_input_type="df", exonHeight, minIntronSize, varPalette=varPalette, color_by_data="published")


# Accessory functions

# 1. Downloading data from HuVar Database






###################################

# varFile<-"D:/Box Sync/Scripts/R/General/variant_info.csv"
varFile<-"D:/Box Sync/Scripts/R/General/CRB1_HuVarDB.csv"
varFile<-"D:/Box Sync/Scripts/R/General/PAX6_HuVarDB.csv"
GOI<-"PAX6"
minIntronSize<-50
exonHeight<-0.075
exonPalette<-"Spectral"
#outFile<-"testSMN1.svg"
# huvarfile<-"./CRB1_var_info_HuVarbase.csv"
# phenotypesToKeep<-c("Retinitis","Leber","Retinal")
#
# HuVarDB<-process_HuVarOutput(filename = huvarfile, phenotypes_to_keep = phenotypesToKeep)
# write_csv(HuVarDB, file = "CRB1_HuVarDB.csv.gz")
#


# use provided GOI variant file
create_gene_lolliplot(gene_name = GOI, varFile, variant_input_type="file", exonHeight, minIntronSize)


# use downloaded GOI variants
variantInfoDL<-get_gene_variants(gene_name = GOI)
create_gene_lolliplot(gene_name = GOI, variantInfoDL, variant_input_type="df", exonHeight, minIntronSize)






























# TO DO
# - protein domains
#https://www.ensembl.info/2015/06/01/biomart-or-how-to-access-the-ensembl-data-from-r/
# domain_location_ENSG00000198763 <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','ensembl_peptide_id','interpro','interpro_start','interpro_end','pfam','pfam_start','pfam_end'), filters ='ensembl_gene_id', values ="ENSG00000198763", mart = ensembl)

# install_bitbucket("ibi_group/disgenet2r")
# load_pack(disgenet2r)
# get_disgenet_api_key(email = "n.owen@ucl.ac.uk", password = "****")
# results <- gene2disease( gene = c( "CRB1"), verbose = TRUE, api_key = "7a24d5091759a68c4196d25f02a2bbdc6a14964e")

