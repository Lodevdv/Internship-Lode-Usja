library(Matrix)
library(UCell)
library(Seurat)

options(warn = 2)  # turns warnings into errors so we get a traceback immediately

SeuratObj <- readRDS("C:/Users/irc/Desktop/Internship Bioinformatics 2025-2026/Lode/Internship-Lode-Usja/RShiny/App/RShiny_tool_object.rds")
data_layer <- GetAssayData(SeuratObj, assay = "RNA", layer = "data")
data_test <- data_layer[, sample(ncol(data_layer), 130000)]

t0 <- Sys.time()
scores <- ScoreSignatures_UCell(
  data_test,
  features = list(test = c("Xcr1", "Cxcl9", "Cxcl10", "Fscn1")),
  maxRank = 1500,
  BPPARAM = BiocParallel::SnowParam(workers = 3)
)
message("Took: ", round(difftime(Sys.time(), t0, units = "secs"), 1), " sec")