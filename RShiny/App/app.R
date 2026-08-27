# Seurat object with the final celltype annotation made in 2026 during Lode's internship under Clint's guidance.
# This object holds only the RNA layer, as that's the assay used for the downstream analysis. The object holds several annotations, of which the celltype_final_2026
# is the final one. The object holds several integrated layers (scANVI, scVI and several Harmony integrations)
# Of these integrations, the Harmony one of Samuel's integration was used for downstream processing. The UMAP of this integration is also present in the object.
# Furthermore, the RNA assay holds 3 layers: counts, data and X, which all hold the counts (So they are redundant?), no scaled data is present, as that would bloat
# the object.

# The important metadata in the object, is: orig.ident (holds experiment data), seurat_clusters (leiden clusters, or louvain?), experiment (holds the experiment perforned), treatment (Test or WT), tech (A or B, referring to 3' or 5' sequencing I think, only Toxo has a different sequencing method), leiden_1 and leiden_2 (referring to leiden resolution used), and celltype_final_2026 

# Additional information calculated in the downstream analysis were: 
# - Conserved markers (significant across all celltypes) (C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Downstream_processing/conserved_markers_summary.pdf)
# - Marker list, just all genes with a LogFC, pval and how many times significant (C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Downstream_processing/all_ranked_markers_summary.xlsx)
# - TF scores using the pseudobulk samples (C:\Users\irc\Desktop\Internship Bioinformatics 2025-2026\Lode\Downstream_processing\TF)
# - GSEA and ORA analysis for each celltype in each experiment (C:\Users\irc\Desktop\Internship Bioinformatics 2025-2026\Lode\Downstream_processing\GSEA and C:\Users\irc\Desktop\Internship Bioinformatics 2025-2026\Lode\Downstream_processing\ORA)
# - DESeq2 results from LNP vs WT, Toxo vs LNP and Toxo Test vs WT (C:\Users\irc\Desktop\Internship Bioinformatics 2025-2026\Lode\Downstream_processing\DESeq2)

#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
# The "New names: * `` -> `...1` output when running the tool is from automatic renaming of unnamed columns in the excel file, and is nothing to worry about
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#


# Import libraries
library("renv")
library("shiny")
library("shinythemes")
library("plotly")
library("bslib")
library("bootstrap")
library("shinyWidgets")
library("DT")
library("shinycssloaders")
library("readxl")
library("dplyr")
library("decoupleR")
library("UCell")

# prep_app_data.R
# Run this ONCE, wherever the full Seurat object lives, to extract only what
# the Shiny app actually needs. The output files are what get deployed to the server.

# SeuratObj <- readRDS("C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/RShiny/App/RShiny_tool_object.rds")
# 
# umap <- SeuratObj@reductions$RNA_harmony_umap_samuel@cell.embeddings
# meta <- SeuratObj@meta.data
# expr <- GetAssayData(SeuratObj, assay = "RNA", layer = "data")
# 
# # save each separately, so the app can load only what it needs at startup
# # (all three are needed, but keeping them separate is cleaner than one big list)
# saveRDS(umap, "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/RShiny/App/data/app_umap.rds")
# saveRDS(meta, "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/RShiny/App/data/app_meta.rds")
# saveRDS(expr, "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/RShiny/App/data/app_expr.rds")
# 
# message("Done. File sizes:")
# for (f in c("app_umap.rds", "app_meta.rds", "app_expr.rds")) {
#   path <- file.path("C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/RShiny/App/data", f)
#   message(f, ": ", round(file.info(path)$size / 1e6, 1), " MB")
# }
# 
# rm(SeuratObj)
# gc()



# IF YOU TRY TO RUN THE APP, AND AN ERROR RETURNS, then run the following line:
# source("C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/RShiny/App/app.R")
# this will return the real error, instead of the standard shiny error, which is not informative

################################################################################################################################################################
##################### Splitting up the data for improved speed #################################################################################################

umap <- readRDS("C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/RShiny/App/data/app_umap.rds")
meta <- readRDS("C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/RShiny/App/data/app_meta.rds")
expr <- readRDS("C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/RShiny/App/data/app_expr.rds")

# genes
genes <- rownames(expr)

condition <- c(
  paste("Orig.ident:", unique(meta$orig.ident)),
  paste("Treatment:", unique(meta$treatment)),
  paste("Experiment:", unique(meta$experiment))
)

metadata <- c("celltype",
              "cluster",
              "treatment",
              "experiment",
              "orig.ident")
################################################################################################################################################################
##################### Excel file loading #######################################################################################################################

######################### FOR MARKER GENES BLOCK ###############################################################################################################

# read excel sheet with the markers
sig_sheet_names <- readxl::excel_sheets("C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/Downstream_processing/all_ranked_markers_summary.xlsx")

# loop over all sheets in the excel sheet and check if essential columns are present?
sig_data_list <- setNames(
  lapply(sig_sheet_names, function(sn) {
    df <- readxl::read_excel("C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/Downstream_processing/all_ranked_markers_summary.xlsx", sheet = sn)
    
    # define the required columns
    required_cols <- c("n_experiments_significant", "max_padj", "min_logFC")
    missing_cols <- setdiff(required_cols, colnames(df))
    
    # check if columns are missing
    if (length(missing_cols) > 0) {
      stop("Sheet '", sn, "' is missing expected column(s): ",
           paste(missing_cols, collapse = ", "))
    }
    # return dataframe
    df
  }),
  sig_sheet_names
)

##################### DESeq2 Excel file loading ###############################################################################################################

# define the right files where to find the DESeq2 data
deseq2_files <- c(
  "LNP"           = "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/Downstream_processing/DESeq2/DEseq2_results_LNP.xlsx",
  "Toxo Test vs WT" = "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/Downstream_processing/DESeq2/DEseq2_results_Toxo_test_wt.xlsx",
  "Toxo vs Other" = "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/Downstream_processing/DESeq2/DEseq2_results_Toxo_vs_other.xlsx"
)

# check if files exist
for (f in deseq2_files) {
  if (!file.exists(f)) stop("DESeq2 file not found: ", f)
}


# define essential columns of the file, "gene" is unnamed in the excel file, but will be named here in the script
deseq2_required_cols <- c("gene", "baseMean", "log2FoldChange",
                          "lfcSE", "stat", "pvalue", "padj")

# find the sheet names in each excel file 
# lapply does it to each of the 3 files
# setNames ensures that the names are the same as the ones defined in the deseq2_files before
deseq2_sheet_names <- setNames(
  lapply(deseq2_files, readxl::excel_sheets), # readxl reads the excel file given
  names(deseq2_files)
)

# now actually load the data, by looping over the names in deseq2_files
deseq2_data_list <- setNames(
  lapply(names(deseq2_files), function(fname) {
    
    # this is a loop
    # get filepath
    fpath <- deseq2_files[[fname]]
    
    # get sheets for that file
    sheets <- deseq2_sheet_names[[fname]]
    
    # loop over every sheet
    setNames(
      lapply(sheets, function(sn) {
        
        # read one particular sheet
        df <- readxl::read_excel(fpath, sheet = sn)
        
        # rename first column to "gene"
        colnames(df)[1] <- "gene"
        
        # check if columns are missing
        missing_cols <- setdiff(deseq2_required_cols, colnames(df))
        if (length(missing_cols) > 0) {
          stop("File '", fname, "', sheet '", sn, "' missing column(s): ",
               paste(missing_cols, collapse = ", "))
        }
        
        # return the dataframe
        df
      }),
      sheets
    )
  }),
  names(deseq2_files)
)

##################### GSEA Excel file loading ###################################################################################################################

# define the path of the GSEA file (which holds the inferred pathways from the DESeq2 results of Toxo vs other immunogenic groups)
gsea_file <- "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/Downstream_processing/GSEA/gsea_results.xlsx"  

# check if file exists
if (!file.exists(gsea_file)) {
  stop("GSEA file not found: ", gsea_file)
}

# read file
gsea_sheet_names <- readxl::excel_sheets(gsea_file)

# define columns
gsea_required_cols <- c("Name", "Term", "ES", "NES", "NOM p-val",
                        "FDR q-val", "FWER p-val", "Tag %", "Gene %", "Lead_genes")
# loops through the sheets, reads them and checks if all columns are present
gsea_data_list <- setNames(
  lapply(gsea_sheet_names, function(sn) {
    df <- readxl::read_excel(gsea_file, sheet = sn)
    missing_cols <- setdiff(gsea_required_cols, colnames(df))
    if (length(missing_cols) > 0) {
      stop("GSEA sheet '", sn, "' missing column(s): ",
           paste(missing_cols, collapse = ", "))
    }
    # return dataframe
    df
  }),
  gsea_sheet_names
)

######################################### ORA results file loading #########################################################################################################
# ORA tests the following: of all genes called significant, is there a given pathway present in that gene set?
# GSEA takes a full ranked list and looks if a pathway's gene set skews towards one end
# the ORA results we stored, took the results from the marker excel file (with celltypes as sheets and experiments as columns)
# so this ORA analysis tried to find celltype specific modules across all conditions for a celltype
# This tab is maybe not useful

# can add multiple files here
ora_files <- c(
  "Log2FC = 0.5 and n_signif = 1" = "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/Downstream_processing/ORA/ORA_padj0.05_logFC0.5_nExp1.xlsx",
  "Log2FC = 0.5 and n_signif = 6" = "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/Downstream_processing/ORA/ORA_padj0.05_logFC0.5_nExp6.xlsx",
  "Log2FC = 0.5 and n_signif = 7" = "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/Downstream_processing/ORA/ORA_padj0.05_logFC0.5_nExp7.xlsx",
  "Log2FC = 1 and n_signif = 1" = "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/Downstream_processing/ORA/ORA_padj0.05_logFC1_nExp1.xlsx",
  "Log2FC = 1 and n_signif = 6" = "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/Downstream_processing/ORA/ORA_padj0.05_logFC1_nExp6.xlsx",
  "Log2FC = 1 and n_signif = 7" = "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/Downstream_processing/ORA/ORA_padj0.05_logFC1_nExp7.xlsx"
)

# check if file exists
for (f in ora_files) {
  if (!file.exists(f)) stop("ORA file not found: ", f)
}

# define required columns
ora_required_cols <- c("Gene_set", "Term", "Overlap", "P-value", "Adjusted P-value",
                       "Old P-value", "Old Adjusted P-value", "Odds Ratio",
                       "Combined Score", "Genes")

# return the names of the sheets in the files
ora_sheet_names <- setNames(
  lapply(ora_files, readxl::excel_sheets),
  names(ora_files)
)

# loop through the files, then through the sheets and check if all columns exist
ora_data_list <- setNames(
  lapply(names(ora_files), function(fname) {
    fpath <- ora_files[[fname]]
    sheets <- ora_sheet_names[[fname]]
    setNames(
      lapply(sheets, function(sn) {
        df <- readxl::read_excel(fpath, sheet = sn)
        missing_cols <- setdiff(ora_required_cols, colnames(df))
        if (length(missing_cols) > 0) {
          stop("ORA file '", fname, "', sheet '", sn, "' missing column(s): ",
               paste(missing_cols, collapse = ", "))
        }
        df
      }),
      sheets
    )
  }),
  names(ora_files)
)

##################### TF z-score Excel file loading ############################################################################################################

# file with z-scores of the inferred TFs
tf_file <- "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/Downstream_processing/TF/TF_activity_zscores_by_celltype_pdata.xlsx"  

# check if it exists
if (!file.exists(tf_file)) {
  stop("TF z-score file not found: ", tf_file)
}

# define sheet names
tf_sheet_names <- readxl::excel_sheets(tf_file)

tf_gene_col <- "...1"  # first column is unnamed, and readxl will name this "...1"

# Full set of possible experiment columns across all sheets combined.
# Not every sheet will contain every one of these -- that's expected.
tf_all_possible_experiments <- c("CITEseq_Final", "CITEseq_LNP_CpG_LNPs", "CITEseq_LNP_WT",
                                 "CITEseq_LNP_eLNPs", "CITEseq_LNP_pIC", "CITEseq_LNP_pIC_LNPs",
                                 "CITEseq_Notch", "CITEseq_Test", "CITEseq_Toxo_Test",
                                 "CITEseq_Toxo_WT")

# Load each sheet, only requiring the gene column to be present.
# For each sheet, also record which experiment columns it actually has.
tf_data_list <- setNames(
  lapply(tf_sheet_names, function(sn) {
    df <- readxl::read_excel(tf_file, sheet = sn)
    
    # check if columns are present
    if (!(tf_gene_col %in% colnames(df))) {
      stop("TF sheet '", sn, "' is missing the gene/TF column '", tf_gene_col, "'.")
    }
    
    # check which experiments exist for this celltype
    present_experiments <- intersect(tf_all_possible_experiments, colnames(df))
    
    if (length(present_experiments) == 0) {
      stop("TF sheet '", sn, "' has no recognized experiment columns. ",
           "Check column naming against tf_all_possible_experiments.")
    }
    
    list(
      data = df,
      experiments = present_experiments
    )
  }),
  tf_sheet_names
)

all_tfs <- unique(unlist(lapply(tf_data_list, function(x) x$data[[tf_gene_col]])))


# ################################################################################################################################################################
# ##################### CollecTRI network loading ################################################################################################################
# 
# # Downloaded the collectri network manually, as decoupleR and OmnipathR gave issues
# # 'https://zenodo.org/records/19773408' is the link
# collectri_raw_path  <- "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/RShiny/App/data/CollecTRI2.tsv.gz"
# 
# # save a cache of the network
# collectri_cache_path <- "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/RShiny/App/data/collectri_net_mouse.rds"
# 
# # check if file exists
# if (file.exists(collectri_cache_path)) {
#   
#   collectri_net <- readRDS(collectri_cache_path)
#   
# } else {
#   
#   if (!file.exists(collectri_raw_path)) {
#     stop("CollecTRI raw file not found: ", collectri_raw_path)
#   }
#   
# # PARSING THE FILE  
#   # Line 1 is a metadata/comment line; real header starts on line 2
#   raw <- readr::read_tsv(
#     collectri_raw_path,
#     skip = 1,
#     col_types = readr::cols(.default = "c")   # read everything as character; we only need 2 columns
#   )
#   
#   # columns in the file
#   tf_col     <- "Transcription Factor (Associated Gene Name)"
#   target_col <- "Target Gene (Associated Gene Name)"
#   
#   # check if columns are missing from the file
#   missing_cols <- setdiff(c(tf_col, target_col), colnames(raw))
#   if (length(missing_cols) > 0) {
#     stop("CollecTRI raw file missing expected column(s): ", paste(missing_cols, collapse = ", "))
#   }
#   
#   # Simple human -> mouse symbol casing conversion (first letter caps, rest lowercase). (not very robust)
#   # NOTE: this is an approximation, not a true ortholog map -- a small number of
#   # genes have mouse symbols that differ from a simple case conversion.
#   to_mouse_case <- function(x) {
#     paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
#   }
#   
#   # create the collectri network from the file
#   collectri_net <- raw %>%
#     select(source = !!tf_col, target = !!target_col) %>%
#     filter(!is.na(source), !is.na(target), source != "", target != "") %>%
#     distinct(source, target) %>%
#     mutate(
#       source = to_mouse_case(source),
#       target = to_mouse_case(target)
#     )
#   
#   # save a cache of the collectri network that was built, to not always redo this part
#   saveRDS(collectri_net, collectri_cache_path)
# }
# 
# # sort network based on sources
# collectri_tfs <- sort(unique(collectri_net$source))
# 
# # Sanity check: how much overlap exists between network targets and your dataset's genes.
# # If this is near-zero, the casing conversion likely isn't matching your gene symbols --
# # verify before trusting downstream results.
# .collectri_overlap_check <- length(intersect(collectri_net$target, genes))
# message("CollecTRI target genes overlapping with dataset: ", .collectri_overlap_check, # this results in 10884 out of 13200 genes in the network having a match
#         " out of ", length(unique(collectri_net$target)), " unique network targets.")
# 
# ################################################################################################################################################################
# ##################### UCell rankings pre-computation ############################################################################################################
# 
# # this ranking step takes a lot of time, so running it for each query would increase runtime, which is why we will cache these results after computing them once
# ucell_ranks_cache <- "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/RShiny/App/data/ucell_ranks.rds"
# 
# # if a precomputed ranking exists, use it
# if (file.exists(ucell_ranks_cache)) {
#   ucell_ranks <- readRDS(ucell_ranks_cache)
#   message("FOUND RANKED CACHE, WILL USE THIS")
#   # if such a ranking does not exist, calculate it (takes time)
# } else {
#   message("Calculating UCell's rankings. Takes time. Results are cached for future use")
#   ucell_ranks <- StoreRankings_UCell(expr)
#   saveRDS(ucell_ranks, ucell_ranks_cache)
# }

################################################################################################################################################################
##################### CollecTRI network loading for AddModuleScore() pipeline ################################################################################################################
# because the previous loading and processing did not include the contribution (stimularoty or inhibitory) interaction between TF and gene
# we need to separate positive and negative interactions, so they do not end up in the same signature

collectri_raw_path  <- "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/RShiny/App/data/CollecTRI2.tsv.gz"

# save a cache of the network
collectri_cache_path <- "C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/RShiny/App/data/collectri_net_mouse.rds"

# check if file exists
if (file.exists(collectri_cache_path)) {

  collectri_net <- readRDS(collectri_cache_path)

} else {

  if (!file.exists(collectri_raw_path)) {
    stop("CollecTRI collectri_raw_path file not found: ", collectri_raw_path)
  }
  
  # PARSING THE FILE
    # Line 1 is a metadata/comment line; real header starts on line 2
    raw <- readr::read_tsv(
      collectri_raw_path,
      skip = 1,
      col_types = readr::cols(.default = "c")   # read everything as character; we only need 2 columns
    )
  
  # columns holding sign/direction evidence, one per source database
  sign_cols <- c(
    "[ExTRI2] Sign", "[TFactS] Sign", "[GOA] Sign", "[SIGNOR] Sign", "[NTNU Curated] Sign"
  )
  
    # columns in the file
    tf_col     <- "Transcription Factor (Associated Gene Name)"
    target_col <- "Target Gene (Associated Gene Name)"

    # check if columns are missing from the file
    missing_cols <- setdiff(c(tf_col, target_col), colnames(raw))
    if (length(missing_cols) > 0) {
      stop("CollecTRI raw file missing expected column(s): ", paste(missing_cols, collapse = ", "))
    }

    # Simple human -> mouse symbol casing conversion (first letter caps, rest lowercase). (not very robust)
    # NOTE: this is an approximation, not a true ortholog map -- a small number of
    # genes have mouse symbols that differ from a simple case conversion.
    to_mouse_case <- function(x) {
      paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
    }

  missing_sign_cols <- setdiff(sign_cols, colnames(raw))
  if (length(missing_sign_cols) > 0) {
    stop("CollecTRI raw file missing expected sign column(s): ", paste(missing_sign_cols, collapse = ", "))
  }
  
  # normalize all evidence tokens (across databases and pipe-separated multi-values)
  # into one shared vocabulary, then classify each row as "activating", "repressing",
  # "conflicting" (both seen) or "unknown" (no usable evidence)
  classify_sign <- function(row_values) {
    tokens <- unlist(strsplit(row_values, "\\|"))
    tokens <- trimws(tokens)
    tokens <- tokens[tokens != "" & !is.na(tokens)]
    
    tokens <- dplyr::case_when(
      tokens %in% c("UP", "+")             ~ "UP",
      tokens %in% c("DOWN", "-")           ~ "DOWN",
      TRUE                                  ~ "UNKNOWN"   # covers UNKNOWN, Unknown, NA-as-text, anything else
    )
    
    has_up   <- any(tokens == "UP")
    has_down <- any(tokens == "DOWN")
    
    dplyr::case_when(
      has_up  && !has_down ~ "activating",
      has_down && !has_up  ~ "repressing",
      has_up  &&  has_down ~ "conflicting",
      TRUE                  ~ "unknown"
    )
  }
  
  # combine all sign columns row-wise into one string per row, then classify
  raw$combined_sign_evidence <- do.call(paste, c(raw[sign_cols], sep = "|"))
  raw$sign <- vapply(raw$combined_sign_evidence, classify_sign, character(1))
  
  message("Sign classification breakdown: ")
  print(table(raw$sign))
  
  # create the collectri network from the file, now keeping sign
  collectri_net <- raw %>%
    select(source = !!tf_col, target = !!target_col, sign) %>%
    filter(!is.na(source), !is.na(target), source != "", target != "") %>%
    distinct(source, target, .keep_all = TRUE) %>%
    mutate(
      source = to_mouse_case(source),
      target = to_mouse_case(target)
    )
  
  saveRDS(collectri_net, collectri_cache_path)
}

# sort network based on sources
collectri_tfs <- sort(unique(collectri_net$source))

# Sanity check: how much overlap exists between network targets and your dataset's genes.
# If this is near-zero, the casing conversion likely isn't matching your gene symbols --
# verify before trusting downstream results.
.collectri_overlap_check <- length(intersect(collectri_net$target, genes))
message("CollecTRI target genes overlapping with dataset: ", .collectri_overlap_check, # this results in 10884 out of 13200 genes in the network having a match
        " out of ", length(unique(collectri_net$target)), " unique network targets.")

################################################################################
##################### Manual AddModuleScore implementation #######################
# Reimplements Seurat::AddModuleScore() directly on the sparse `expr` matrix,
# since SeuratObj was removed from memory at startup. Same algorithm:
# bins genes by average expression, samples matched control genes per bin,
# and returns (mean signature expression - mean control expression) per cell.

# the Seurat function AddModuleScore() requires a SeuratObject, which we will not load to improve efficiency of this script
compute_module_score <- function(expr_mat, features, nbin = 24, ctrl = 100, seed = 1) {
  
  set.seed(seed)
  
  # average expression of every gene across all cells, used for binning
  data_avg <- Matrix::rowMeans(expr_mat)
  
  # order genes by average expression, then cut into nbin equal-sized bins
  # (tiny random jitter breaks exact ties, matching Seurat's approach)
  data_avg_sorted <- data_avg[order(data_avg)]
  jittered <- data_avg_sorted + stats::rnorm(length(data_avg_sorted)) / 1e30
  data_cut <- cut(
    x = jittered,
    breaks = stats::quantile(jittered, probs = seq(0, 1, length.out = nbin + 1)),
    labels = FALSE,
    include.lowest = TRUE
  )
  names(data_cut) <- names(data_avg_sorted)
  
  # for each signature gene, sample `ctrl` control genes from the same bin
  ctrl_use <- character(0)
  for (gene in features) {
    if (!(gene %in% names(data_cut))) next
    gene_bin <- data_cut[gene]
    bin_genes <- names(data_cut)[data_cut == gene_bin]
    n_sample <- min(ctrl, length(bin_genes))
    ctrl_use <- c(ctrl_use, sample(bin_genes, size = n_sample, replace = FALSE))
  }
  ctrl_use <- unique(ctrl_use)
  
  # mean expression of signature genes vs. matched control genes, per cell
  features_present <- intersect(features, rownames(expr_mat))
  feature_scores <- Matrix::colMeans(expr_mat[features_present, , drop = FALSE])
  ctrl_scores     <- Matrix::colMeans(expr_mat[ctrl_use, , drop = FALSE])
  
  score <- feature_scores - ctrl_scores
  
  list(
    score           = score,
    features_present = features_present,
    features_missing = setdiff(features, features_present),
    n_control_genes  = length(ctrl_use)
  )
}

################################################################################################################################################################
########################### Info Table | Home Page #############################################################################################################
# X-axis (columns)
# get column names
report_cols <- unique(as.character(meta$orig.ident))

# set to character and remove duplicates
experiment_map <- tapply(
  as.character(meta$experiment),
  as.character(meta$orig.ident),
  function(x) unique(x)[1]
)

# set to character and remove duplicates
WT_map <- tapply(
  as.character(meta$treatment),
  as.character(meta$orig.ident),
  function(x) unique(x)[1]
)

# Y-axis (rows)
# define the rows of the table, mostly lay out here
conditions <- c(as.character(
  tags$div(
    tags$img(src = "Orig_idents.png", height = "30px")," Orig.idents"
  )),
  as.character(
    tags$div(tags$img(src = "mouse_icon.png", height = "30px"),
             " Treatment")),
  as.character(
    tags$div(
      tags$img(src = "lab.png", height = "30px")," Experiment CITEseq:"
    ))
)
tbl <- data.frame(
  Condition = conditions,
  matrix(
    "",
    nrow = length(conditions),
    ncol = length(report_cols),
    dimnames = list(NULL, report_cols)
  ),
  check.names = FALSE)
# Fill rows
tbl[1, report_cols] <- report_cols                      # Orig.idents
tbl[2, report_cols] <- WT_map[report_cols]             # Treatment
tbl[3, report_cols] <- experiment_map[report_cols]     # Experiment

################################################################################################################################################################
############## Base of UI ######################################################################################################################################

# this is just HTML layout
ui <- page_navbar(
  theme = shinytheme("united"),
  title = img(src = "irc_logo_transparant.png", height = "30px"),
  window_title = "cDC1 cell atlas",
  bg = "#8EE5EE",
  inverse = TRUE,
  
  header = tags$head(
    tags$style(HTML("
      .Title {
        font-size: 32px;
        margin-bottom: 5px;
            }
            .body {
              font-family: sans-serif;
            }
            .card{
              position: relative;
              width: 100%;
              max-width: 800px;
              height: 700px;
              margin: 2em auto;
              padding: 1em;
              overflow: visible;
              display: flex;
              flex-direction: column;
         }
         
            .table{
              width: 1200px;
              margin-left:0px;
        }

          .card .plotly {
              flex: 1;
              min-height: 0;
        }

          .input{
            position:relative;
            width: 100%;
            padding: 2em;
            margin: 2em auto;
            
          }
          .metaplots{
          margin:50px;
          padding; 50px;
          }
            
  "))
 ),
                   
################## Home Page ###################################################################################################################################
nav_panel(title = "Home",
         tags$div(
           class = "Title",
           h1("Home")
         ),
         
         # top explanation here
         tags$div(

          p("This dataset holds approximately 137062 cells and 36048 genes, all different states of the cDC1 differentiation process.", tags$br(),
            "This atlas was integrated with Harmony, and holds 9 different experimental groups (WT, Toxo, LNP) from 27 different samples.", tags$br(),
            "Following preprocessing, the cell annotation was finalized using a combination of dotplots, Clint's original annotation of the", tags$br(),
            "individual datasets, and Samuel's original annotation. Next, further downstream analysis was performed on pseudobulk data, which", tags$br(),
            "was made from mouse x celltype pseudobulking. The downstream analyses included:", tags$br(),
            "- DESeq2 results from Toxo Test vs WT, Toxo vs LNP and LNP vs WT", tags$br(),
            "- Markers given certain thresholds", tags$br(),
            "- Inferred TFs based on pseudobulk data", tags$br(),
            "- GSEA and ORA analysis results", tags$br(),
            "- Module finding using UCell", tags$br(),
            ),
          
          p("Important metadata is:", tags$br(),
            "- Celltype: Early Immature, Late Immature, Early Mature, Late Mature, Pre_cDC1, Proliferating_cDC1, cDC1_engulfing_RBC, Other cDC1s and Low_quality_Immature_cDC1", tags$br(),
            "- Samples: JVE008, JVE010, SAM016, SAM05, SAM06, SAM2, SAM3, VBO004, VBO005, VBO006, VBO007, VBO008, VBO009, VBO010, VBO011, VBO012", tags$br(),
            "- Experimental groups: CITEseq_Toxo (has WT and Test group), CITEseq_Notch, CITEseq_Final, CITEseq_Test, CITEseq_LNP_WT, CITEseq_LNP_eLNPs, CITEseq_LNP_pIC_LNPs, CITEseq_LNP_CpG_LNPs, CITEseq_LNP_pIC",
          ),
         div(class ="table",
             tableOutput("report_table")
         )
    )
),
                   
################################################################################################################################################################
############################## Gene plots ######################################################################################################################


# lay out of the gene plots tab
nav_panel(title = "Gene Plots",
          # top explanation
          p("On this page, you can select a gene to check its expression across the dataset, as well as explore the associated violin plot.", tags$br(),
            "You can also view the different subsets of the data on the bottom right plot"),
            card( style = "margin-bottom: 5px;",
                  max_height = "150px",
                  selectizeInput(
                    inputId = "gene",
                    label = "Select gene for feature plot",
                    choices = NULL,
                    selected = NULL,
                    multiple = FALSE,
                    options = list(
                      placeholder = 'Type a gene...',
                      create = TRUE,   # allows typing custom values
                      dropdownParent = "body"
                    ))
          ),
          
          layout_columns(
            
            # feature plot
            card(  
              withSpinner(
                plotlyOutput("feature_plot", height = "700px",width = "700px")
              )),
            
            # violin plot
            card(
              card_header("Violin Plot"),
              plotOutput("ViolinPlot", height = "700px"),
              
              # dimplot
              selectizeInput(
                inputId = "Dimplot",
                choices = metadata,
                label = "Select MetaData for Dimplot",
                selected = NULL,
                multiple = FALSE,
                options = list(
                  placeholder = 'Type a gene...',
                  create = TRUE,   # allows typing custom values
                  dropdownParent = "body"
                )),
              
              plotOutput("DimPlotMeta", height = "500")
              
            ),col_widths = 12)
),

################################################################################################################################################################
################# Cell Metadata ################################################################################################################################

# cell metadata tab lay out. Here, you can plot a gene on a subset of the data of your choice
nav_panel(title = "Cell Metadata", 
          tags$div( class = "input",
            p("Here, you can select a gene to be showed on a subset of choice. Subsets can be samples, treatments or tech (WT or Test)",
              "Each selected subset is displayed as a separate UMAP, with the colour indicating the selected metadata variable.", tags$br(),
              "You can select multiple samples/experiments/treatments to be subsetted, resulting in multiple UMAPs for one particular gene.
              However, selecting more than 7-8 subsets, results in distortion of some UMAPs in the RShiny tool itself.", tags$br(),
              "This problem is fixed when opening the tool in the Firefox browser."),
                    layout_columns(
                      style = "max-width: 1000px; margin: 0 auto;",     
                      
                      # gene choice
                      selectizeInput(
                        inputId = "meta",
                        label = "Select Gene to analysis",
                        choices = NULL,
                        selected = NULL,
                        multiple = F,
                        options = list(
                          placeholder = 'Type a gene...',
                          create = TRUE   # allows typing custom values
                        )),
                      
                      # subset choice
                      selectizeInput(
                        inputId = "cond",
                        label = "Select or type the subset",
                        choices = condition,
                        selected = NULL,
                        multiple = TRUE,
                        options = list(
                          placeholder = 'Type a subset...',
                          create = TRUE,   # allows typing custom values
                          dropdownParent = "body"
                        ),)
                      
                    ),
                    # choose what to put at the bottom right UMAP
                    radioButtons(
                      inputId = "colour_by",
                      label = "Colour by",
                      choices = c(
                        "Gene expression" = "gene",
                        "Cell type" = "celltype",
                        "Cluster" = "cluster",
                        "Treatment" = "treatment",
                        "Experiment" = "experiment",
                        "Orig.ident" = "orig.ident"
                      ),
                      selected = "gene",
                      inline = TRUE
                    ),col_widths= c(4,4,4)),
          
          actionButton("start", "Generate plots"),
          
          uiOutput("metaplots_ui")
),
################################################################################################################################################################
######################################### Marker gene filter ###################################################################################################

# marker gene filter according to given thresholds by the user
sig_search_tab <- tabPanel(
  "Gene Significance Search",
  tags$div(
    style = "margin-bottom: 15px;",
    h4("Marker Gene Search"),
    p("Filter genes by significance thresholds per celltype, or view genes conserved across all experiments for the selected celltype.", tags$br(),
      "Conserved implies a maximal p-value of 0.05, minimal logFC of 1 and significant across all experiments for that celltype.", tags$br(),
      "The ranked list was made by performing a Wilcoxon rank sum test on the raw counts of the whole dataset.")
  ),
  sidebarLayout(
    sidebarPanel(
      
      # select celltype of interest
      selectInput(
        "sig_celltype", "Cell type:",
        choices = sig_sheet_names
      ),
      
      # filtering metric (pval)
      numericInput(
        "sig_maxpval", "Max adjusted p-value:",
        value = 0.05, min = 0, max = 1, step = 0.001
      ),
      
      # filtering metric (logFC)
      numericInput(
        "sig_minlogfc", "Min logFC:",
        value = 1, step = 0.05
      ),
      
      # filtering metric (amount of experiments in which the gene should be significant)
      numericInput(
        "sig_minexp", "Min. number of experiments significant:",
        value = 1, min = 1, step = 1
      ),
      
      # perform search
      actionButton("sig_search", "Search", class = "btn-primary"),
      br(), br(),
      downloadButton("sig_download", "Download results (.csv)"),
      br(), br(), hr(),
      actionButton("conserved_search", "Show Conserved Markers", class = "btn-success"),
      br(), br(),
      downloadButton("conserved_download", "Download conserved markers (.csv)")
    ),
    mainPanel(
      DTOutput("sig_table"),
      hr(),
      h4("Conserved Markers"),
      DTOutput("conserved_table")
    )
  )
),

################################################################################################################################################################
######################################### Forest plot LogFC  ###################################################################################################
# interpretation of the DESeq2 results, using a forest plot
# the user specifies a list of genes, and which DE the user wants

deseq2_forest_tab <- tabPanel(
  "DESeq2 Forest Plot",
  p("Select a DESeq2 result file and one or more comparisons (celltype/condition), then", tags$br(),
  "enter a list of genes to generate stacked forest plots of log2 fold changes."), 
  br(),
  p("The DESeq2 comparisons are the following: Toxoplasma treatment vs its WT group, Toxoplasma treatment vs one ", tags$br(),
  "other immunogenic (LNP) treatment(specify which one) and the LNP treatments individually vs their WT group."), br(),
  
  p("The file specifies which DESeq2 result (LNP holds the comparison of the individual LNP groups", tags$br(),
  "against their WT group, Toxo Test vs WT is the DESeq2 of the Toxoplasma treatment versus its WT, and", tags$br(),
  "the Toxo vs Other file holds the comparison of the Toxoplasma treatment versus the individual LNP groups.", tags$br(),
  "Then, you select the results from which DESeq2 analysis in that file you want plotted as forest plot.", tags$br(),
  "By default, all are selected. For example, in 'LNP', the CITEseq_LNP_CpG_LN_Early Mature holds the DE", tags$br(),
  "results between LNP treatment (LNP CpG LNPs) and the LNP WT, for the Early Mature celltype."),
  
  sidebarLayout(
    sidebarPanel(
      
      # select DESeq2 analysis of choice
      selectInput(
        "deseq2_file", "DESeq2 result file:",
        choices = names(deseq2_files)
      ),
      
      # select which specific comparison you want (like in the Toxo Test vs WT analysis, which celltype)
      selectizeInput(
        "deseq2_comparison", "Comparison(s) / celltype(s):",
        choices = NULL,   # populated server-side based on file
        multiple = TRUE,
        options = list(placeholder = "Select one or more comparisons...")
      ),
      
      # Which genes you want plotted
      textAreaInput(
        "deseq2_genes", "Gene list (one per line, or comma-separated):",
        rows = 6,
        placeholder = "Xcr1\nClec9a\nBatf3"
      ),
      actionButton("deseq2_plot_btn", "Generate Forest Plots", class = "btn-primary")
    ),
    mainPanel(
      withSpinner(uiOutput("deseq2_forest_plot_ui")),
      br(),
      h4("Combined results across selected comparisons"),
      DTOutput("deseq2_forest_table")
    )
  )
),

################################################################################################################################################################
######################################### GSEA #################################################################################################################
# lay out of GSEA plot
gsea_tab <- tabPanel(
  "GSEA Results",
  p("Select a GSEA comparison (celltype/experiment) to view its enriched pathways, ",
    "ranked by normalized enrichment score (NES) or significance. These GSEA results are derived from the DESeq2 results of the Toxo vs other immunogenic conditions.
    This information is derived from the DESeq2 between the Toxoplasma treatment group and the other individual LNP treatment groups. So the CITEseq_LNP_CpG_LNPs_Late Matu comparison,
    holds the inferred pathways (from GSEA), of the DE results between the Toxoplasma treatment group and the LNP CpG LNPs treatment, and this for the early mature celltype.
    A positive NES refers to a pathway enriched in that celltype, for the Toxoplasma treatment group. A negative NES, means that that pathway is enriched in the respective LNP treatment
    group, for that specific celltype."),
  tags$br(),
  p("FDR or false discovery rate relates to the percentage of false positives in the results, and is set by default to 0,25 (reasoning in the original GSEA paper)"),
  br(),
  p("Naming is cutoff, because excel sheets have a character limit of 31. So for example for Late Mature cells, the Toxoplasma vs LNP_CpG_LNPs will be called CITEseq_LNP_CpG_LNPs_Late MatuR"),
  sidebarLayout(
    sidebarPanel(
      
      # select which comparison (Toxo vs which LNP) you want
      selectInput(
        "gsea_sheet", "Comparison (celltype / experiment):",
        choices = gsea_sheet_names
      ),
      
      # sort by what
      selectInput(
        "gsea_sort_by", "Sort pathways by:",
        choices = c(
          "NES (high to low)"        = "NES_desc",
          "NES (low to high)"        = "NES_asc",
          "FDR q-val (most sig.)"    = "FDR_asc"
        )
      ),
      
      # cutoff value of FDR
      numericInput(
        "gsea_fdr_cutoff", "FDR q-val threshold:",
        value = 1, min = 0, max = 1, step = 0.01
      ),
      helpText("GSEA conventionally uses FDR < 0.25 (looser than typical DE",
               "thresholds like 0.05), because it tests many overlapping",
               "gene sets. Adjust as needed for your use case.")
    ),
    mainPanel(
      DTOutput("gsea_table"),
      br(),
      downloadButton("gsea_download", "Download this comparison's results (.csv)"),
      br(), br(),
      h4("Summary: Term, NES, FDR q-val"),
      DTOutput("gsea_summary_table")
    )
  )
),


################################################################################################################################################################
######################################### ORA results ##########################################################################################################

ora_tab <- tabPanel(
  "ORA Results",
  p("Select an ORA threshold set and celltype to view enriched pathways ", tags$br(),
    "based on the gene set over-representation analysis. Overlap dictates how many genes in the dataset overlap with the database genes for that pathway.", tags$br(),
    "P-values indicate statistical significance (FIX: what are the 'old p-values'), and odds ratio measures how strongly the selected gene set is", tags$br(),
    "associated with a particular pathway, compared to what you might expect by chance."),
  br(),
  p("These results were calculated from the marker genes list. THis list was made using a Wilxocon rank test for all genes in each experiment separately,", tags$br(),
  "for each celltype individually. 6 ORA analyses were performed in total, 3 for LogFC 0.5 and 3 for LogFC 1, with 1, 6 or 7 significant",tags$br(),
  "experiments (meaning the genes involved needed to be significant across n amount of experimental groups, before passing this filtering"),
  br(),
  p("So in short, this ORA tried to find enriched pathways in a specific celltype across all experimental conditions. The ORA used genes which", tags$br(),
    "followed these thresholds: LogFC >= 1, pval < 0.05 and significant in at least 5 experimental groups."),
  
  sidebarLayout(
    sidebarPanel(
      
      # select input ORA file (6 of them can be chosen, see before)
      selectInput(
        "ora_file", "Specify file:",
        choices = names(ora_files)
      ),
      
      # select the celltype you want enriched pathways for
      selectInput(
        "ora_sheet", "Celltype:",
        choices = NULL  # populated server-side based on file
      ),
      
      # sort by what
      selectInput(
        "ora_sort_by", "Sort pathways by:",
        choices = c(
          "Adjusted P-value (most sig.)" = "adjp_asc",
          "Odds Ratio (high to low)"     = "or_desc",
          "Combined Score (high to low)" = "cs_desc"
        )
      ),
      checkboxInput(
        "ora_sig_only", "Filter by adjusted p-value threshold", value = FALSE
      ),
      numericInput(
        "ora_padj_cutoff", "Adjusted p-value threshold:",
        value = 0.05, min = 0, max = 1, step = 0.01
      )
    ),
    mainPanel(
      DTOutput("ora_table"),
      br(),
      downloadButton("ora_download", "Download this comparison's results (.csv)"),
      br(), br(),
      h4("Summary: Term, Odds Ratio, Adjusted P-value"),
      DTOutput("ora_summary_table")
    )
  )
),

################################################################################################################################################################
######################################### TF inference #########################################################################################################
# TFs inferred from the pseudobulk data, for each experimental group individually for each celltype
# this resulted in an excel file which has celltypes as sheets, and then the experimental conditions as columns and TFs as rows
# the data in the file are the z-scores of each TF (= number of standard deviations the activity score of that TF differs from the mean activity score of that TF
# across all experimental groups), this can be dominated by certain groups, so keep this in mind
tf_tab <- tabPanel(
  "TF Inference",
  p("Select a celltype and one or more transcription factors to view their inferred activity (z-score) across experimental conditions, based on ", tags$br(),
    "pseudobulk data."),
  br(),
  p("Z-score is how many standard deviations the activity score, of that inferred TF, differs from the mean activity score", tags$br(),
    "of that TF, across all experimental groups. Negative means lower than the mean, positive means higher than the mean."),
  
    sidebarLayout(
      sidebarPanel(
        selectInput(
          "tf_celltype", "Celltype:",
          choices = tf_sheet_names
        ),
        radioButtons(
          "tf_mode", "Mode:",
          choices = c("Top TFs by experiment" = "top_n",
                      "Manual TF selection"    = "manual"),
          selected = "top_n"
        ),
        conditionalPanel(
          condition = "input.tf_mode == 'top_n'",
          selectInput(
            "tf_reference_experiment", "Rank TFs by z-score in:",
            choices = NULL  # populated server-side, based on celltype
          ),
          numericInput(
            "tf_top_n", "Number of top TFs:",
            value = 50, min = 1, max = 200, step = 1
          )
        ),
        conditionalPanel(
          condition = "input.tf_mode == 'manual'",
          selectizeInput(
            "tf_selected", "Transcription factor(s):",
            choices = NULL,
            multiple = TRUE,
            options = list(placeholder = "Type a TF name...", dropdownParent = "body")
          )
        ),
        helpText("z-scores are calculated per TF across all experimental groups; ",
                 "values reflect relative activity, not absolute expression.")
      ),
      mainPanel(
        withSpinner(plotOutput("tf_heatmap", height = "700px")),
        br(),
        DTOutput("tf_table"),
        br(),
        downloadButton("tf_download", "Download TF data (.csv)")
      )
    )
  ),

# ################################################################################################################################################################
# ############################## Module finding via TF (UCell) ###################################################################################################

# ucell_tab <- tabPanel(
#   "UCell TF Score",
#   p("Select one or more transcription factors. Their target genes are retrieved from the ", tags$br(),
#     "CollecTRI regulatory network, intersected with genes present in this dataset, and used ", tags$br(),
#     "to compute a UCell module score per cell, shown on the UMAP. Note, UCell needs some ", tags$br(),
#     "time to calculate module scores, so be patient."),
#   p("UCell evaluates gene signatures in single-cell databases. THis signature score is based on the Mann-Whitney U statistic.", tags$br(),
#     "It first ranks all genes for each cell in the dataset (this is precomputed for this tool), and then it looks at the ranks of the subset the user provides,", tags$br(),
#     "and calculates the statistic for that gene set against the rest of the genes in the cell. Scores are then normalized", tags$br(),
#     "Each cell is scored on its own. The test statistic measures how consistent the signature genes are ranked at the very top of that cell's expression list,", tags$br(),
#     "compared to all the other genes. A score of 1 means that all genes in the signature are the highest-expressing of that cell. 0.5 is random distribution,", tags$br(),
#     "while 0 means that all genes in the signature are the lowest expressed genes"
#   ),
# 
#   h1("UNDER CONSTRUCTION"),
#   sidebarLayout(
#     sidebarPanel(
#       selectizeInput(
#         "ucell_tfs", "Transcription factor(s):",
#         choices = collectri_tfs,
#         multiple = TRUE,
#         options = list(placeholder = "Type a TF name...", dropdownParent = "body")
#       ),
#       radioButtons(
#         "ucell_mode", "Scoring mode:",
#         choices = c("Combined score (all TFs pooled into one signature)" = "combined",
#                     "Separate score per TF" = "separate"),
#         selected = "combined"
#       ),
#       actionButton("ucell_run", "Calculate UCell Score", class = "btn-primary"),
#       helpText("A signature needs at least 2 target genes present in the dataset ",
#                "to be scored.")
#     ),
#     mainPanel(
#       withSpinner(uiOutput("ucell_umap_ui")),
#       br(),
#       h4("Target gene summary"),
#       DTOutput("ucell_gene_table")
#     )
#   )
# ),


# ################################################################################################################################################################
# ############################## Module finding (AddModuleScore()) ###################################################################################################

# note: We had major problems with the runtime while performing UCell on the whole dataset, as well as taking up all of the CPU and the RAM, so we
# will use the AddModuleScore() of Seurat instead
# this function takes the average expression of the signature genes in each cell, and subtracts the average expression of a randomly sampled control gene set
# the result is a score centered near 0, positive means the signature genes are expressed higher than expression matched background
# negative means lower

# the workflow goes as follows:
# Seurat takes every gene in the dataset and sorts them into expression bins (default 24), based on average expression across all cells
# so you end up with low expression bins, middle, higher...
# then for the signature genes, it finds which bins they fall into. Then for each signature gene, Seurat randomly samples genes from the same expression bin
# These control genes are not biologically related, they just have similar expression to the signature gene
# then for each cell, a per cell score is calculated [(avg expression of signature genes in that cell) - (avg expression of matched control genes in that cell)']
# so a positive score will mean that your gene set's expression is meaningfully higher than what you'd expect from a random set of genes at that same overall expression level

# # Seurat explanation: Calculate the average expression levels of each program (cluster) on single cell level, 
# # subtracted by the aggregated expression of control feature sets. All analyzed features are binned based on averaged expression, 
# # and the control features are randomly selected from each bin.

module_score_tab <- tabPanel(
  "Gene Set Score",
  p("Enter a set of genes to compute a per-cell module score, shown on the UMAP. ", tags$br(),
    "The score is the average expression of your genes minus the average expression of ", tags$br(),
    "a matched set of control genes (sampled from genes with similar overall expression levels), ", tags$br(),
    "following the same method as Seurat's AddModuleScore."),
  sidebarLayout(
    sidebarPanel(
      textAreaInput(
        "module_genes", "Gene list (one per line, or comma-separated):",
        rows = 6,
        placeholder = "Xcr1\nClec9a\nBatf3"
      ),
      numericInput(
        "module_ctrl", "Control genes per signature gene:",
        value = 100, min = 10, max = 500, step = 10
      ),
      actionButton("module_run", "Calculate Module Score", class = "btn-primary"),
      helpText("At least 2 genes from your list need to be present in the dataset.")
    ),
    mainPanel(
      withSpinner(plotOutput("module_umap_plot", height = "600px")),
      br(),
      h4("Gene summary"),
      DTOutput("module_gene_table")
    )
  )
),

# ################################################################################################################################################################
# ############################## Module finding via TF #####################################################################################

module_tf_tab <- tabPanel(
  "TF Module Score (Seurat)",
  p("Select one or more transcription factors. Their activating target genes are retrieved ", tags$br(),
    "from the CollecTRI network (repressed/inhibitory and conflicting/unknown-direction targets ", tags$br(),
    "are excluded), intersected with genes present in this dataset, and scored using a ", tags$br(),
    "Seurat-style module score (average signature expression minus matched control genes).", tags$br(),
    "A positive score means the signature gene set has a higher expression than that of a sampled gene set from the data"),
  sidebarLayout(
    sidebarPanel(
      selectizeInput(
        "module_tf_tfs", "Transcription factor(s):",
        choices = NULL,
        multiple = TRUE,
        options = list(placeholder = "Type a TF name...", dropdownParent = "body")
      ),
      radioButtons(
        "module_tf_mode", "Scoring mode:",
        choices = c("Combined score (all TFs pooled into one signature)" = "combined",
                    "Separate score per TF" = "separate"),
        selected = "combined"
      ),
      numericInput(
        "module_tf_ctrl", "Control genes per signature gene:",
        value = 100, min = 10, max = 500, step = 10
      ),
      actionButton("module_tf_run", "Calculate Module Score", class = "btn-primary"),
      helpText("A signature needs at least 2 activating target genes present in the dataset ",
               "to be scored.")
    ),
    mainPanel(
      withSpinner(uiOutput("module_tf_umap_ui")),
      br(),
      h4("Target gene summary"),
      DTOutput("module_tf_gene_table")
    )
  )
),

################################################################################################################################################################
############################## Contact #########################################################################################################################
  nav_panel(title = "Contact", 
           h3("Creator Cell Atlas:"),
           p("Lode Van de Vreken (intern): lodevandevreken@gmail.com"),
           h3("Mentor Cell Atlas:"),
           p("Clint De Nolf"),
           
           tags$div( class = "input"),
           
           tags$div()
  ),
  
################################################################################################################################################################
################################ Links #########################################################################################################################
  nav_spacer(),
  nav_menu(
   title = "Links",
   align = "right",
   nav_item(tags$a("Posit", href = "https://posit.co")),
   nav_item(tags$a("Shiny", href = "https://shiny.posit.co"))
  )
)
################################################################################################################################################################
############################      SERVER     ###################################################################################################################
################################################################################################################################################################
# Define server logic
server <- function(input, output,session) {
  
  ####################################### Home Table############################################################################################################
  output$report_table <- renderTable(
    tbl,
    striped = TRUE,
    bordered = TRUE,
    spacing = "s",
    # Makes it so that text within the table will also be treated as a html elements
    sanitize.text.function = function(x) x
  )
  
  ######################################### Gene plots##########################################################################################################
  output$selected_gene <- renderPrint({input$gene})
  
  # take feature plot and put it into output
  output$feature_plot <- renderPlotly({ 
    
    req(input$gene)
    
    # UMAP
    df <- data.frame(
      UMAP_1 = umap[, "RNAharmonyumapsamuel_1"],
      UMAP_2 = umap[, "RNAharmonyumapsamuel_2"]
    )
    
    # set expression and celltype in the dataframe
    df$expr <- as.numeric(expr[input$gene, ])
    df$celltype <- meta$celltype_final_2026
    
    # To make Hover Text of cells
    df$hover <- paste(
      rownames(df),
      "<br>Expression:",
      round(df$expr, 2),
      "<br>Celltype:",df$celltype
    )
    # Making a Plotly plot that resembles a featureplot
    plot_ly(
      df,
      x = ~UMAP_1,
      y = ~UMAP_2,
      color = ~expr,
      type = "scattergl",
      mode = "markers",
      text = ~hover,
      #Removing Coordinates
      hoverinfo= "text",
      marker = list(size = 3)
    )
  })
  
  # create the output of the violin plot, showing expression of the gene of interest for each celltype
  output$ViolinPlot <- renderPlot({
    
    req(input$gene)
    
    df_violin <- data.frame(
      expr= as.numeric(expr[input$gene,]),
      celltype = meta$celltype_final_2026
    )
    
    # make the violin plot
    ggplot(df_violin, aes(x = celltype, y= expr,fill = celltype)) +
      geom_violin()+
      scale_fill_brewer(palette = "Set3") +
      theme_bw()+  
      theme(
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)
      )
  })

  # make the dimplot, which is the UMAP with a metadata group of interest
  output$DimPlotMeta <- renderPlot({
    
    req(input$Dimplot)
    
    df <- data.frame(
      UMAP_1 = umap[, "RNAharmonyumapsamuel_1"],
      UMAP_2 = umap[, "RNAharmonyumapsamuel_2"]
    )
    
    # define this metadata in the dataframe
    df$celltype <- meta$celltype_final_2026
    df$cluster <- meta$leiden_2
    df$treatment <- meta$treatment
    df$experiment <- meta$experiment
    df$orig.ident <- meta$orig.ident
    
    df$MetaDataColumn <- switch(
      input$Dimplot,
      "celltype" = df$celltype,
      "cluster" = df$cluster,
      "treatment" = df$treatment,
      "experiment" = df$experiment,
      "orig.ident" = df$orig.ident
    )
    
    # actually make the dimplot
    ggplot(df,aes(x = UMAP_1, y = UMAP_2, color=MetaDataColumn ))+
      geom_point(size= 0.5)+
      theme_classic()+
      labs(x= "UMAP_1", y= "UMAP_2")
  })
  
  
######################################### Cell metadata ########################################################################################################
  ##############################################################################################################################################################
  
  # this tab allows you to plot a gene's expression on different subsets of the data
  
  # Recalculates when the "generate plots" button is clocked
  PlotNameList <- eventReactive(input$start, {
    req(input$cond)
    input$cond
  })
  
  # build page layout, ony plotlyOutput box per selected subset, since the number of subsets is not known in advance
  output$metaplots_ui <- renderUI({
    
    plots <- PlotNameList()
    req(plots)
    tagList(
      
      # loop over plots and generate a div + plotlyOutput for each, naming them metaplot_1, metaplot_2 etc
      lapply(seq_along(plots), function(i) {
        
        div(
          style = "
          width: 1200px;
          height: 650px;
          margin: 20px auto;
        ",
          h4(plots[i]),
          
          withSpinner(
            plotlyOutput(
              outputId = paste0("metaplot_", i),
              width = "1200px",
              height = "600px"
            )
          )
        )
      })
    )
  })
  
  # This is where the actual plot logic occurs
  # it builds one shared dataframe (UMAP coords and metadata), then loops over each subset and creates its own renderPlotly
  observeEvent(PlotNameList(), {
    
    plots <- PlotNameList()
    
    req(input$meta)
    req(input$colour_by)
    
    # Create dataframe
    df <- data.frame(
      
      # UMAp coordinates
      UMAP_1 = umap[, "RNAharmonyumapsamuel_1"],
      UMAP_2 = umap[, "RNAharmonyumapsamuel_2"]
    )
    
    # define the metadata
    df$expr <- as.numeric(expr[input$meta, ])
    df$celltype <- meta$celltype_final_2026
    df$cluster <- meta$leiden_2
    df$treatment <- meta$treatment
    df$experiment <- meta$experiment
    df$orig.ident <- meta$orig.ident
    
    # Select colouring variable
    df$groupby <- switch(
      input$colour_by,
      "gene" = df$expr,
      "celltype" = df$celltype,
      "cluster" = df$cluster,
      "treatment" = df$treatment,
      "experiment" = df$experiment,
      "orig.ident" = df$orig.ident
    )
    
    
    # Create each plot
    lapply(seq_along(plots), function(i) {
      
      local({ # wraps each loop iteration (if not used, all closures created in lapply would share the same ii value, so every plot would show the last subset instead of its own)
        
        ii <- i
        subset_name <- plots[ii]
        
        output[[paste0("metaplot_", ii)]] <- renderPlotly({
          
          # parse selection
          parts <- strsplit(subset_name, ":", fixed = TRUE)[[1]]
          
          validate(
            need(
              length(parts) == 2,
              paste("Invalid subset:", subset_name)
            )
          )
          
          column <- tolower(trimws(parts[1]))
          word <- trimws(parts[2])
          
          validate(
            need(
              column %in% colnames(df),
              paste("Unknown metadata column:", column)
            )
          )
          
          # filtering step
          df_Filter <- df[df[[column]] == word, , drop = FALSE]
          
          validate(
            need(
              nrow(df_Filter) > 0,
              paste("No cells found for:", subset_name)
            )
          )
          
          # make actual plot
          p <- plot_ly(
            data = df_Filter,
            x = ~UMAP_1,
            y = ~UMAP_2,
            color = ~groupby,
            type = "scattergl",
            mode = "markers",
            marker = list(size = 3)
          )
          
          # more lay out
          p %>%
            layout(
              title = list(
                text = subset_name,
                x = 0.02
              ),
              
              xaxis = list(
                title = "UMAP 1",
                scaleanchor = "y",
                scaleratio = 1
              ),
              
              yaxis = list(
                title = "UMAP 2"
              ),
              
              margin = list(
                l = 60,
                r = 180,
                t = 70,
                b = 60
              )
            )
        })
      })
    })
  })

######################################### Marker gene filter ###################################################################################################
  
  # runs when "Search" is clicked and pulls the selected excel sheet from sig_data_list and then filters by the three thresholds
  sig_results <- eventReactive(input$sig_search, {
    req(input$sig_celltype)
    
    df <- sig_data_list[[input$sig_celltype]]
    validate(need(!is.null(df), "No data found for this celltype."))
    
    filtered <- df %>%
      filter(
        max_padj                  <= input$sig_maxpval,
        min_logFC                 >= input$sig_minlogfc,
        n_experiments_significant >= input$sig_minexp
      ) %>%
      arrange(max_padj)
    
    filtered
  }, ignoreNULL = FALSE)
  
  gene_col <- "names"
  
  # derives the gene-name column from sig_results() for a clean table display
  sig_genes_only <- reactive({
    filtered <- sig_results()
    validate(need(gene_col %in% colnames(filtered),
                  paste0("Expected gene column '", gene_col, "' not found in this sheet.")))
    data.frame(gene = filtered[[gene_col]])
  })
  
  # table shows names
  output$sig_table <- renderDT({
    datatable(
      sig_genes_only(),
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  # download results to a csv
  output$sig_download <- downloadHandler(
    filename = function() {
      paste0("significant_genes_", input$sig_celltype, "_",
             Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv2(sig_results(), file, row.names = FALSE)
    }
  )
  
  ################## Conserved marker filter ####################
  
  # same as before, just that there are hard cutoffs here, no influence of the user here
  conserved_results <- eventReactive(input$conserved_search, {
    req(input$sig_celltype)
    
    df <- sig_data_list[[input$sig_celltype]]
    validate(need(!is.null(df), "No data found for this celltype."))
    
    # Per-experiment p-value columns are named "padj_EXPERIMENTNAME"
    experiment_pval_cols <- grep("^padj_", colnames(df), value = TRUE)
    n_experiments_present <- length(experiment_pval_cols)
    
    validate(need(n_experiments_present > 0,
                  "Could not detect per-experiment p-value columns - check column naming."))
    
    df %>%
      filter(
        max_padj                  < 0.05,
        min_logFC                 > 1,
        n_experiments_significant == n_experiments_present
      ) %>%
      arrange(max_padj)
  })
  
  conserved_genes_only <- reactive({
    filtered <- conserved_results()
    validate(need(gene_col %in% colnames(filtered),
                  paste0("Expected gene column '", gene_col, "' not found in this sheet.")))
    data.frame(gene = filtered[[gene_col]])
  })
  
  output$conserved_table <- renderDT({
    datatable(
      conserved_genes_only(),
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  output$conserved_download <- downloadHandler(
    filename = function() {
      paste0("conserved_markers_", input$sig_celltype, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv2(conserved_results(), file, row.names = FALSE)
    }
  )
  
######################################### Forest plot LogFC DESeq2  ############################################################################################
  
  # Update comparison choices when file changes; default to all comparisons selected
  # when the user picks a DESeq2 file, then the comparison dropdown is updated with that files actual sheet name
  observeEvent(input$deseq2_file, {
    req(input$deseq2_file)
    
    # sheet name choices for the chosen file
    choices <- deseq2_sheet_names[[input$deseq2_file]]
    updateSelectizeInput(
      session, "deseq2_comparison",
      choices = choices,
      selected = choices  # default: all comparisons in the file
    )
  })
  
  # on button click, this parses the gene list (given as input), and loops ober ebery selected comparison sheet and filters to the requested genes
  deseq2_plot_data <- eventReactive(input$deseq2_plot_btn, {
    req(input$deseq2_file, input$deseq2_comparison, input$deseq2_genes)
    
    gene_list <- trimws(unlist(strsplit(input$deseq2_genes, "[,\n]+")))
    gene_list <- gene_list[gene_list != ""]
    validate(need(length(gene_list) > 0, "Please enter at least one gene."))
    
    all_missing <- character(0)
    
    combined <- lapply(input$deseq2_comparison, function(comp) {
      df <- deseq2_data_list[[input$deseq2_file]][[comp]]
      if (is.null(df)) return(NULL)
      
      filtered <- df %>% filter(.data[["gene"]] %in% gene_list)
      if (nrow(filtered) == 0) return(NULL)
      
      missing_here <- setdiff(gene_list, filtered[["gene"]])
      if (length(missing_here) > 0) {
        all_missing <<- c(all_missing, paste0(missing_here, " (", comp, ")"))
      }
      
      # here it calculates the 95% CI for every gene's logFC and logFCSE
      filtered %>%
        mutate(
          comparison = comp,
          ci_lower = log2FoldChange - 1.96 * lfcSE,
          ci_upper = log2FoldChange + 1.96 * lfcSE,
          significant = ifelse(!is.na(padj) & padj < 0.05, "Significant", "Not significant")
        )
    }) %>%
      # everything is stacked into one data frame here
      bind_rows()
    
    validate(need(nrow(combined) > 0,
                  "None of the entered genes were found in any selected comparison."))
    
    if (length(all_missing) > 0) {
      showNotification(
        paste("Genes not found in some comparisons:", paste(all_missing, collapse = ", ")),
        type = "warning", duration = 10
      )
    }
    
    combined %>%
      mutate(comparison = factor(comparison, levels = input$deseq2_comparison)) %>%
      arrange(comparison, log2FoldChange)
  })
  
  # Dynamic plot height: scale with number of comparisons and genes so facets stay readable
  deseq2_plot_height <- reactive({
    df <- deseq2_plot_data()
    n_comparisons <- length(unique(df$comparison))
    n_genes <- length(unique(df$gene))
    max(400, n_comparisons * (100 + n_genes * 25))
  })
  
  output$deseq2_forest_plot_ui <- renderUI({
    plotOutput("deseq2_forest_plot", height = paste0(deseq2_plot_height(), "px"))
  })
  
  # draws the actual forest plot, a point plus its CI per gene
  output$deseq2_forest_plot <- renderPlot({
    plot_df <- deseq2_plot_data()
    
    ggplot(plot_df, aes(x = log2FoldChange,
                        y = reorder(gene, log2FoldChange),
                        color = significant)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper), size = 0.6) +
      scale_color_manual(values = c("Significant" = "firebrick",
                                    "Not significant" = "grey60")) +
      facet_wrap(~ comparison, ncol = 1, scales = "free_y") +
      labs(
        x = "Log2 Fold Change (95% CI)",
        y = NULL,
        color = "padj < 0.05",
        title = paste0(input$deseq2_file, " - selected comparisons")
      ) +
      theme_bw(base_size = 14) +
      theme(strip.text = element_text(face = "bold"))
  })
  
  # create the table, showing the numeric results
  output$deseq2_forest_table <- renderDT({
    datatable(
      deseq2_plot_data() %>%
        select(comparison, gene, baseMean, log2FoldChange, lfcSE, pvalue, padj),
      options = list(pageLength = 15, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      formatRound(columns = c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"), digits = 4)
  })
  
######################################### GSEA results #########################################################################################################
  
  # pulls the selected sheet from the file
  gsea_selected <- reactive({
    req(input$gsea_sheet)
    
    # pulls the sheet from the file
    df <- gsea_data_list[[input$gsea_sheet]]
    validate(need(!is.null(df), "No data found for this comparison."))
    df
  })
  
  # filters the user's FDR cutoff and then sorts by which option was picked below
  gsea_filtered_sorted <- reactive({
    df <- gsea_selected() %>%
      
      # actual filtering
      filter(.data[["FDR q-val"]] < input$gsea_fdr_cutoff)
    
    # sort by which option was chosen
    df <- switch(
      input$gsea_sort_by,
      "NES_desc" = df %>% arrange(desc(NES)),
      "NES_asc"  = df %>% arrange(NES),
      "FDR_asc"  = df %>% arrange(.data[["FDR q-val"]])
    )
    
    df
  })
  
  # full table with the GSEA columns, plus a condensed term/NES/FDR summary
  output$gsea_table <- renderDT({
    datatable(
      gsea_filtered_sorted() %>%
        select(Term, NES, ES, `NOM p-val`, `FDR q-val`, `FWER p-val`,
               `Tag%`, `Gene%`, Lead_genes),
      options = list(pageLength = 20, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      formatRound(columns = c("NES", "ES", "NOM p-val", "FDR q-val", "FWER p-val"),
                  digits = 4)
  })
  
  output$gsea_summary_table <- renderDT({
    datatable(
      gsea_filtered_sorted() %>%
        select(Term, NES, `FDR q-val`),
      options = list(pageLength = 20, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      formatRound(columns = c("NES", "FDR q-val"), digits = 4)
  })
  
  # possible download of the table
  output$gsea_download <- downloadHandler(
    filename = function() {
      paste0("GSEA_", input$gsea_sheet, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv2(gsea_filtered_sorted(), file, row.names = FALSE)
    }
  )
  
######################################### ORA results ##########################################################################################################

  # repopulates the celltype dropdown whenever the user switches ORA threshold file
  observeEvent(input$ora_file, {
    req(input$ora_file)
    updateSelectInput(
      session, "ora_sheet",
      choices = ora_sheet_names[[input$ora_file]]
    )
  })
  
  # pulls the chosen file and sheet combination from ora_data_list
  ora_selected <- reactive({
    req(input$ora_file, input$ora_sheet)
    df <- ora_data_list[[input$ora_file]][[input$ora_sheet]]
    validate(need(!is.null(df), "No data found for this file/celltype."))
    df
  })
  
  # also reactive, optionally filters by adjusted p-val and then sorts by p-val, odds ratio or combined score
  ora_filtered_sorted <- reactive({
    df <- ora_selected()
    
    if (isTRUE(input$ora_sig_only)) {
      df <- df %>% filter(.data[["Adjusted P-value"]] < input$ora_padj_cutoff)
    }
    
    df <- switch(
      input$ora_sort_by,
      "adjp_asc" = df %>% arrange(.data[["Adjusted P-value"]]),
      "or_desc"  = df %>% arrange(desc(`Odds Ratio`)),
      "cs_desc"  = df %>% arrange(desc(`Combined Score`))
    )
    
    df
  })
  
  # same table summary as before (GSEA)
  output$ora_table <- renderDT({
    datatable(
      ora_filtered_sorted() %>%
        select(Term, Overlap, `P-value`, `Adjusted P-value`,
               `Old P-value`, `Old Adjusted P-value`,
               `Odds Ratio`, `Combined Score`, Genes),
      options = list(pageLength = 20, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      formatRound(columns = c("P-value", "Adjusted P-value", "Old P-value",
                              "Old Adjusted P-value", "Odds Ratio", "Combined Score"),
                  digits = 4)
  })
  
  output$ora_summary_table <- renderDT({
    datatable(
      ora_filtered_sorted() %>%
        select(Term, `Odds Ratio`, `Adjusted P-value`),
      options = list(pageLength = 20, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      formatRound(columns = c("Odds Ratio", "Adjusted P-value"), digits = 4)
  })
  
  output$ora_download <- downloadHandler(
    filename = function() {
      paste0("ORA_", input$ora_file, "_", input$ora_sheet, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv2(ora_filtered_sorted(), file, row.names = FALSE)
    }
  )
######################################### TF inference #########################################################################################################

  # Update experiment choices (for top-N mode) and TF choices (for manual mode)
  # whenever celltype changes, update 2 dropdowns based only on what's actually present for that celltype (reference experiment choice and manual TF choices)
  observeEvent(input$tf_celltype, {
    req(input$tf_celltype)
    entry <- tf_data_list[[input$tf_celltype]]
    
    updateSelectInput(
      session, "tf_reference_experiment",
      choices = entry$experiments
    )
    
    updateSelectizeInput(
      session, "tf_selected",
      choices = entry$data[[tf_gene_col]],
      server = TRUE
    )
  })
  
  # Determine which TFs to display, depending on mode
  tf_chosen <- reactive({
    entry <- tf_data_list[[input$tf_celltype]]
    validate(need(!is.null(entry), "No data found for this celltype."))
    
    # if top TFs are selected, show these
    if (input$tf_mode == "top_n") {
      req(input$tf_reference_experiment, input$tf_top_n)
      validate(need(input$tf_reference_experiment %in% entry$experiments,
                    "Selected experiment is not present for this celltype."))
      
      entry$data %>%
        mutate(abs_z = abs(.data[[input$tf_reference_experiment]])) %>%
        arrange(desc(abs_z)) %>% # take absolute z score for the top listing, not just positive or negative
        slice_head(n = input$tf_top_n) %>%
        pull(.data[[tf_gene_col]])
      
      # else, manual input is given, and these need to be searched
    } else {
      req(input$tf_selected)
      validate(need(length(input$tf_selected) > 0, "Select at least one transcription factor."))
      input$tf_selected
    }
  })
  
  # filter the celltype's data down to the chosen TFs, keeping only the experiment columns actually present for that celltype
  tf_selected_data <- reactive({
    entry <- tf_data_list[[input$tf_celltype]]
    chosen_tfs <- tf_chosen()
    
    entry$data %>%
      filter(.data[[tf_gene_col]] %in% chosen_tfs) %>%
      rename(TF = !!tf_gene_col) %>%
      select(TF, all_of(entry$experiments))
  })
  
  # shows which experiment columns exist for the current celltype
  tf_present_experiments <- reactive({
    req(input$tf_celltype)
    tf_data_list[[input$tf_celltype]]$experiments
  })
  
  # this make the format long (Experiment and z-score columns)
  tf_long <- reactive({
    tf_selected_data() %>%
      tidyr::pivot_longer(
        cols = all_of(tf_present_experiments()),
        names_to = "Experiment",
        values_to = "z_score"
      )
  })
  
  # Order TFs on the heatmap by their z-score in the reference experiment (top-N mode only),
  # so the "top hit" is visually at one end rather than alphabetical
  # the list is sorted alphabetically, while the heatmap is sorted top to bottom from scores
  tf_heatmap_order <- reactive({
    if (input$tf_mode == "top_n") {
      entry <- tf_data_list[[input$tf_celltype]]
      entry$data %>%
        filter(.data[[tf_gene_col]] %in% tf_chosen()) %>%
        arrange(desc(abs(.data[[input$tf_reference_experiment]]))) %>%
        pull(.data[[tf_gene_col]])
    } else {
      tf_chosen()
    }
  })
  
  # creates the heatmap
  output$tf_heatmap <- renderPlot({
    plot_df <- tf_long()
    plot_df$TF <- factor(plot_df$TF, levels = rev(tf_heatmap_order()))
    
    # find max z-value
    z_max <- max(abs(plot_df$z_score), na.rm = TRUE)  # symmetric range around 0
    
    subtitle_txt <- if (input$tf_mode == "top_n") {
      paste0("Top ", input$tf_top_n, " TFs by |z-score| in ", input$tf_reference_experiment)
    } else {
      "Manually selected TFs"
    }
    
    # actual creation
    ggplot(plot_df, aes(x = Experiment, y = TF, fill = z_score)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(
        low = "blue", mid = "white", high = "red", midpoint = 0,
        # define scale
        limits = c(-z_max, z_max)
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8)
      ) +
      labs(x = NULL, y = NULL, fill = "z-score",
           title = paste0("TF activity - ", input$tf_celltype),
           subtitle = subtitle_txt)
  })
  
  # create a numeric table
  output$tf_table <- renderDT({
    datatable(
      tf_selected_data(),
      options = list(pageLength = 15, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      formatRound(columns = tf_present_experiments(), digits = 3)
  })
  
  # allow csv output
  output$tf_download <- downloadHandler(
    filename = function() {
      mode_tag <- if (input$tf_mode == "top_n") {
        paste0("top", input$tf_top_n, "_", input$tf_reference_experiment)
      } else {
        "manual_selection"
      }
      paste0("TF_zscores_", input$tf_celltype, "_", mode_tag, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv2(tf_selected_data(), file, row.names = FALSE)
    }
  )
  
# ############################## Module finding via TF (UCell) ###################################################################################################
# 
#   # on a button click, it subsets the collectri_net to the selected TFs' regulatory edges, meaning the target nodes
#   ucell_result <- eventReactive(input$ucell_run, {
#     req(input$ucell_tfs)
#     t_start <- Sys.time()
#     message("UCell run started at: ", t_start)
#     
#     net_sub <- collectri_net %>% filter(source %in% input$ucell_tfs)
#     validate(need(nrow(net_sub) > 0,
#                   "No target genes found for the selected TF(s) in CollecTRI."))
#     message("net_sub built: ", nrow(net_sub), " rows, at +",
#             round(difftime(Sys.time(), t_start, units = "secs"), 1), "s")
#     
#     # if you want a score for all TF targets combined, pool all target genes into 1 signature
#     if (input$ucell_mode == "combined") {
#       
#       target_genes <- unique(net_sub$target)
#       present <- intersect(target_genes, genes)
#       missing_n <- length(setdiff(target_genes, genes))
#       
#       # check if there are enough genes for the signature
#       validate(need(length(present) >= 2,
#                     paste0("Only ", length(present), " target gene(s) present in dataset -- ",
#                            "need at least 2 for a UCell score. Try selecting more/different TFs.")))
#       
#       signature_list <- list("Combined_Score" = present)
#       
#       gene_summary <- data.frame(
#         Signature          = "Combined_Score",
#         TFs_included        = paste(input$ucell_tfs, collapse = ", "),
#         n_targets_total     = length(target_genes),
#         n_targets_present   = length(present),
#         n_targets_missing   = missing_n
#       )
#     
#       # this is the separate mode, which calculates a signature for each individual TFs' target nodes, instead of pooling
#     } else {
#       message("Performing UCell IN SEPARATE MODE")
#       
#       sig_by_tf <- lapply(input$ucell_tfs, function(tf) {
#         targets <- unique(net_sub$target[net_sub$source == tf])
#         intersect(targets, genes)
#       })
#       names(sig_by_tf) <- input$ucell_tfs
#       
#       n_present_per_tf <- sapply(sig_by_tf, length)
#       keep <- n_present_per_tf >= 2
#       
#       if (any(!keep)) {
#         showNotification(
#           paste0("Skipped (fewer than 2 target genes present): ",
#                  paste(names(sig_by_tf)[!keep], collapse = ", ")),
#           type = "warning", duration = 8
#         )
#       }
#       
#       sig_by_tf <- sig_by_tf[keep]
#       validate(need(length(sig_by_tf) > 0,
#                     "None of the selected TFs had at least 2 target genes present in the dataset."))
#       
#       signature_list <- sig_by_tf
#       
#       gene_summary <- data.frame(
#         Signature        = names(sig_by_tf),
#         TFs_included      = names(sig_by_tf),
#         n_targets_total   = sapply(names(sig_by_tf), function(tf) {
#           length(unique(net_sub$target[net_sub$source == tf]))
#         }),
#         n_targets_present = sapply(sig_by_tf, length)
#       )
#     }
#     
#     # use pre-computed ranking, instead of calculating it again
#     message("Signature built, entering ScoreSignatures_UCell at +",
#             round(difftime(Sys.time(), t_start, units = "secs"), 1), "s")
#     
#     scores <- ScoreSignatures_UCell(features = signature_list,
#                                     precalc.ranks = ucell_ranks)
# 
#     message("ScoreSignatures_UCell finished at +",
#             round(difftime(Sys.time(), t_start, units = "secs"), 1), "s")
#     
#     list(scores = scores, gene_summary = gene_summary)
#   })
#   
#   # attach UMAP coordinates to the resulting per-cell UCell scores
#   message("PLOTTING UCell RESULTS")
#   ucell_plot_df <- reactive({
#     res <- ucell_result()
#     score_df <- as.data.frame(as.matrix(res$scores))
#     colnames(score_df) <- gsub("_UCell$", "", colnames(score_df))
#     
#     df <- data.frame(
#       UMAP_1 = umap[, "RNAharmonyumapsamuel_1"],
#       UMAP_2 = umap[, "RNAharmonyumapsamuel_2"]
#     )
#     cbind(df, score_df)
#   })
#   
#   # size the plot area based on how many signature panels there are (one for combined, multiple for separate)
#   output$ucell_umap_ui <- renderUI({
#     req(ucell_result())
#     n_panels <- ncol(as.data.frame(ucell_result()$scores))
#     height <- if (input$ucell_mode == "combined") 600 else max(500, ceiling(n_panels / 2) * 400)
#     plotOutput("ucell_umap_plot", height = paste0(height, "px"))
#   })
#   
#   # reshapes to long format and draws one UMAP panel per signature, coloured by UCell score
#   output$ucell_umap_plot <- renderPlot({
#     df <- ucell_plot_df()
#     score_cols <- setdiff(colnames(df), c("UMAP_1", "UMAP_2"))
#     
#     df_long <- df %>%
#       tidyr::pivot_longer(cols = all_of(score_cols), names_to = "Signature", values_to = "Score")
#     
#     ggplot(df_long, aes(x = UMAP_1, y = UMAP_2, color = Score)) +
#       geom_point(size = 0.5) +
#       scale_color_gradientn(colors = c("lightgrey", "blue", "darkblue")) +
#       facet_wrap(~ Signature) +
#       theme_classic(base_size = 13) +
#       labs(title = "UCell module score", x = "UMAP_1", y = "UMAP_2")
#   })
#   
#   # output a per-signature summary table
#   output$ucell_gene_table <- renderDT({
#     req(ucell_result())
#     datatable(
#       ucell_result()$gene_summary,
#       options = list(pageLength = 10, scrollX = TRUE),
#       rownames = FALSE
#     )
#   })
  
############################## Module finding via TF (AddModuleScore()) ###################################################################################################

  module_result <- eventReactive(input$module_run, {
    req(input$module_genes)
    
    gene_list <- trimws(unlist(strsplit(input$module_genes, "[,\n]+")))
    gene_list <- gene_list[gene_list != ""]
    validate(need(length(gene_list) > 0, "Please enter at least one gene."))
    
    result <- compute_module_score(
      expr_mat = expr,
      features = gene_list,
      ctrl = input$module_ctrl
    )
    
    validate(need(length(result$features_present) >= 2,
                  paste0("Only ", length(result$features_present), " gene(s) from your list ",
                         "are present in the dataset -- need at least 2 to compute a score.")))
    
    if (length(result$features_missing) > 0) {
      showNotification(
        paste("Genes not found in dataset:", paste(result$features_missing, collapse = ", ")),
        type = "warning", duration = 8
      )
    }
    
    result
  })
  
  module_plot_df <- reactive({
    res <- module_result()
    
    df <- data.frame(
      UMAP_1 = umap[, "RNAharmonyumapsamuel_1"],
      UMAP_2 = umap[, "RNAharmonyumapsamuel_2"],
      Score  = res$score
    )
    df
  })
  
  output$module_umap_plot <- renderPlot({
    df <- module_plot_df()
    
    ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = Score)) +
      geom_point(size = 0.5) +
      scale_color_gradientn(colors = c("lightgrey", "blue", "darkblue")) +
      coord_fixed() +
      theme_classic(base_size = 13) +
      labs(title = "Gene set module score", x = "UMAP_1", y = "UMAP_2")
  })
  
  output$module_gene_table <- renderDT({
    res <- module_result()
    
    summary_df <- data.frame(
      Metric = c("Genes present", "Genes missing", "Control genes used"),
      Value  = c(
        paste(res$features_present, collapse = ", "),
        if (length(res$features_missing) > 0) paste(res$features_missing, collapse = ", ") else "None",
        res$n_control_genes
      )
    )
    
    datatable(
      summary_df,
      options = list(pageLength = 5, dom = 't'),
      rownames = FALSE
    )
  }) 
  
######################################### Module Score (Seurat AddModuleScore-style) #########################################################################

  
  observe({
    updateSelectizeInput(session, "module_tf_tfs", choices = collectri_tfs, server = TRUE)
  })
  
  module_tf_result <- eventReactive(input$module_tf_run, {
    req(input$module_tf_tfs)

    # subset network to selected TFs, keeping only activating relationships
    net_sub <- collectri_net %>%
      filter(source %in% input$module_tf_tfs, sign == "activating")

    validate(need(nrow(net_sub) > 0,
                  "No activating target genes found for the selected TF(s) in CollecTRI."))

    if (input$module_tf_mode == "combined") {

      target_genes <- unique(net_sub$target)
      present <- intersect(target_genes, genes)
      missing_n <- length(setdiff(target_genes, genes))

      validate(need(length(present) >= 2,
                    paste0("Only ", length(present), " activating target gene(s) present in dataset -- ",
                           "need at least 2 for a module score. Try selecting more/different TFs.")))

      score_result <- compute_module_score(
        expr_mat = expr,
        features = present,
        ctrl = input$module_tf_ctrl
      )

      score_matrix <- data.frame(Combined_Score = score_result$score)

      gene_summary <- data.frame(
        Signature          = "Combined_Score",
        TFs_included        = paste(input$module_tf_tfs, collapse = ", "),
        n_targets_activating = length(target_genes),
        n_targets_present    = length(present),
        n_targets_missing    = missing_n
      )

    } else {

      sig_by_tf <- lapply(input$module_tf_tfs, function(tf) {
        targets <- unique(net_sub$target[net_sub$source == tf])
        intersect(targets, genes)
      })
      names(sig_by_tf) <- input$module_tf_tfs

      n_present_per_tf <- sapply(sig_by_tf, length)
      keep <- n_present_per_tf >= 2

      if (any(!keep)) {
        showNotification(
          paste0("Skipped (fewer than 2 activating target genes present): ",
                 paste(names(sig_by_tf)[!keep], collapse = ", ")),
          type = "warning", duration = 8
        )
      }

      sig_by_tf <- sig_by_tf[keep]
      validate(need(length(sig_by_tf) > 0,
                    "None of the selected TFs had at least 2 activating target genes present in the dataset."))

      score_list <- lapply(sig_by_tf, function(target_set) {
        compute_module_score(expr_mat = expr, features = target_set, ctrl = input$module_tf_ctrl)$score
      })
      score_matrix <- as.data.frame(score_list)
      colnames(score_matrix) <- names(sig_by_tf)

      gene_summary <- data.frame(
        Signature            = names(sig_by_tf),
        TFs_included          = names(sig_by_tf),
        n_targets_activating  = sapply(names(sig_by_tf), function(tf) {
          length(unique(net_sub$target[net_sub$source == tf]))
        }),
        n_targets_present     = sapply(sig_by_tf, length)
      )
    }

    list(scores = score_matrix, gene_summary = gene_summary)
  })

  module_tf_plot_df <- reactive({
    res <- module_tf_result()

    df <- data.frame(
      UMAP_1 = umap[, "RNAharmonyumapsamuel_1"],
      UMAP_2 = umap[, "RNAharmonyumapsamuel_2"]
    )
    cbind(df, res$scores)
  })

  output$module_tf_umap_ui <- renderUI({
    req(module_tf_result())
    n_panels <- ncol(module_tf_result()$scores)
    height <- if (input$module_tf_mode == "combined") 600 else max(500, ceiling(n_panels / 2) * 400)
    plotOutput("module_tf_umap_plot", height = paste0(height, "px"))
  })

  output$module_tf_umap_plot <- renderPlot({
    df <- module_tf_plot_df()
    score_cols <- setdiff(colnames(df), c("UMAP_1", "UMAP_2"))

    df_long <- df %>%
      tidyr::pivot_longer(cols = all_of(score_cols), names_to = "Signature", values_to = "Score")

    ggplot(df_long, aes(x = UMAP_1, y = UMAP_2, color = Score)) +
      geom_point(size = 0.5) +
      scale_color_gradientn(colors = c("lightgrey", "blue", "darkblue")) +
      coord_fixed() +
      facet_wrap(~ Signature) +
      theme_classic(base_size = 13) +
      labs(title = "TF module score (activating targets only)", x = "UMAP_1", y = "UMAP_2")
  })

  output$module_tf_gene_table <- renderDT({
    req(module_tf_result())
    datatable(
      module_tf_result()$gene_summary,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)