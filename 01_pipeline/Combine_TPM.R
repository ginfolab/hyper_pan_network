################################################################################
#####################################################整合所有的TPM 并进行初步预处理
################################################################################
#############这部分要修改
path <- "/mnt/gpfsA/home/zb-jiang/ncbi/public/sra/Ecoli/"
setwd(path)
folder_paths <- list.dirs(path, recursive = FALSE)
folder_names <- basename(folder_paths)
library(readxl)
library(openxlsx)
wb <- createWorkbook()
#############这部分要修改
cluster_info <- read.table("/mnt/gpfsA/home/zb-jiang/genome/Ecoli_pangenome/info.txt",header = T,sep = "\t")
for (variable in 1:length(folder_names)) {
  current_path <- paste(folder_paths[variable],"/TPM.txt",sep = "")
  
  if (!file.exists(current_path)) {
    next  # 文件不存在，跳到下一个循环
  }
  TPM <- read.table(current_path,header = T,sep = "\t")[,-1]
  rownames(TPM) <- cluster_info$clustid
  filtered_TPM <- TPM[rowSums(TPM >= 4) > 0, ]
  filtered_TPM <- filtered_TPM[,colSums(filtered_TPM >= 4) > 1000 ]
  if (ncol(filtered_TPM) > 30) {
    filtered_TPM <- filtered_TPM[, 1:30]
  }
  if (ncol(filtered_TPM) < 15 ) {
    next
  }
  addWorksheet(wb, sheetName = folder_names[variable])
  writeData(wb, sheet = folder_names[variable], filtered_TPM, startCol = 1, startRow = 1, colNames = TRUE, rowNames = T)
}
#############这部分要修改
existing_excel_path <- "Escherichia_coli_TPM_new.xlsx"
saveWorkbook(wb, existing_excel_path)