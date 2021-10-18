library(matrixStats)
library(MatrixGenerics)
library(BiocGenerics)
library(stats4)
library(BiocGenerics)
library(parallel)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(SingleCellExperiment)
library(scRNAseq)

# 宸ヤ綔鐩綍
work_dir <- "C:/Users/Administrator/Desktop/bkb" 

# tcga瀵瑰簲鑲跨槫鏌ヨ
project <- "TCGA-KIRC"
data_category <- "Transcriptome Profiling"
data_type <- "Gene Expression Quantification"
workflow_type <- "HTSeq - Counts"
legacy <- FALSE


# 璁剧疆涓哄綋鍓嶅伐浣滅洰褰?
setwd(work_dir)
getwd()

# 鏁版嵁涓嬭浇鏌ヨ
DataDirectory <- paste0(work_dir,"/GDC/",gsub("-","_",project))
FileNameData <- paste0(DataDirectory, "_","RNAseq_HTSeq_Counts",".rda")

# 鏁版嵁鎯呭喌涓嬭浇
query <- GDCquery(project = project,
                  data.category = data_category,
                  data.type = data_type, 
                  workflow.type = workflow_type,
                  legacy = legacy)

# 鎬绘暟鏌ヨ
samplesDown <- getResults(query,cols=c("cases"))
cat("Total sample to download:", length(samplesDown))

# 鑲跨槫鏌ヨ
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown, typesample = "TP")
cat("Total TP samples to down:", length(dataSmTP))

# 姝ｅ父鏌ヨ
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,typesample = "NT")
cat("Total NT samples to down:", length(dataSmNT))


# 涓村簥鏁版嵁涓嬭浇
GDCdownload(query = query,
            directory = DataDirectory,files.per.chunk=6, method='client')

# data 赋值
data <- GDCprepare(query = query, 
                   save = TRUE, 
                   directory =  DataDirectory,
                   save.filename = FileNameData)


data_expr <- assay(data)
dim(data_expr)
expr_file <- paste0(DataDirectory, "_","All_HTSeq_Counts",".txt")
write.table(data_expr, file = expr_file, sep="\t", row.names =T, quote = F)