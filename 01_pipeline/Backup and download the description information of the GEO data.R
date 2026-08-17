# 安装和加载 GEOquery 包
library(GEOquery)
setwd("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli")
# 假设你的 GSE 列表是一个向量，例如：
gse_list <- read.table("./GSE_name.txt",header = T,sep = "\t")  # 这里替换为你的 GSE ID 列表
gse_list <- gse_list$GSE_ID
# 创建一个空的 data frame 来存储 GSE 信息
gse_summary <- data.frame(GSE_ID = character(),
                          Summary = character(),
                          Strain = character(),
                          stringsAsFactors = FALSE)

# 遍历 GSE 列表并获取每个 GSE 的信息
for (gse in gse_list) {
  gse_data <- getGEO(gse, GSEMatrix = TRUE)  # 获取 GSE 数据
  
  # 获取 GSE 的元数据
  gse_info <- pData(phenoData(gse_data[[1]]))  # 获取第一个数据集的元数据
  
  # 获取摘要信息（如果有的话）
  summary_info <- ifelse(length(gse_data[[1]]@experimentData@other[["summary"]]) > 0, gse_data[[1]]@experimentData@other[["summary"]], "No summary available")  # 检查摘要是否存在
  # 获取品系名称（一般可以在 `characteristics_ch1` 列中找到）
  strain_info <- ifelse(length(gse_info$characteristics_ch1) > 0, paste(unique(gse_info$`strain:ch1`), collapse = ", "), "No strain info available")  # 检查品系信息是否存在
  
  # 将数据添加到 data frame
  gse_summary <- rbind(gse_summary, data.frame(GSE_ID = gse,
                                               Summary = summary_info,
                                               Strain = strain_info,
                                               stringsAsFactors = FALSE))
}

saveRDS(gse_summary,"gse_summary.RDS")


library(GEOquery)

setwd("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli")

# 读取 GSE ID 列表
gse_list <- read.table("./GSE_name.txt", header = TRUE, sep = "\t")$GSE_ID

# 初始化保存表格
gse_summary <- data.frame(
  GSE_ID = character(),
  Summary = character(),
  Strain = character(),
  Keyword = character(),
  Genotype_info = character(),
  Growth_info = character(),
  stringsAsFactors = FALSE
)

# 保存失败的GSE
failed_gse <- c()

# 设置较长的超时时间，避免连接中断
options(timeout = 300)

# 遍历 GSE 列表
for (gse in gse_list) {
  message("Processing: ", gse)
  
  tryCatch({
    gse_data <- getGEO(gse, GSEMatrix = TRUE, destdir = tempdir())  # 强制下载到临时文件夹，避免缓存出错
    
    gse_info <- pData(phenoData(gse_data[[1]]))
    expData <- gse_data[[1]]@experimentData
    
    # 提取 summary
    summary_info <- if (!is.null(expData@other[["summary"]])) {
      expData@other[["summary"]]
    } else {
      "No summary available"
    }
    
    # 提取 keyword
    keyword_info <- if (!is.null(expData@other[["overall_design"]])) {
      paste(expData@other[["overall_design"]], collapse = ", ")
    } else {
      "No overall_design available"
    }
    
    # 提取 strain（试多个字段）
    strain_info <- "No strain info available"
    if ("strain:ch1" %in% colnames(gse_info)) {
      strain_info <- paste(unique(gse_info$`strain:ch1`), collapse = "; ")
    } else if ("characteristics_ch1" %in% colnames(gse_info)) {
      strain_info <- paste(unique(gse_info$characteristics_ch1), collapse = "; ")
    }
    
    # 提取 genotype 信息
    genotype_info <- "No genotype info"
    if ("genotype:ch1" %in% colnames(gse_info)) {
      genotype_info <- paste(unique(gse_info$`genotype:ch1`), collapse = "; ")
    } else if ("characteristics_ch1" %in% colnames(gse_info)) {
      geno_related <- grep("genotype", gse_info$characteristics_ch1, value = TRUE, ignore.case = TRUE)
      if (length(geno_related) > 0) {
        genotype_info <- paste(unique(geno_related), collapse = "; ")
      }
    }
    
    # 提取 growth 信息
    growth_info <- "No growth info"
    if ("growth:ch1" %in% colnames(gse_info)) {
      growth_info <- paste(unique(gse_info$`growth:ch1`), collapse = "; ")
    } else if ("characteristics_ch1" %in% colnames(gse_info)) {
      grow_related <- grep("growth", gse_info$characteristics_ch1, value = TRUE, ignore.case = TRUE)
      if (length(grow_related) > 0) {
        growth_info <- paste(unique(grow_related), collapse = "; ")
      }
    }
    
    # 追加信息到表格
    gse_summary <- rbind(
      gse_summary,
      data.frame(
        GSE_ID = gse,
        Summary = summary_info,
        Strain = strain_info,
        Keyword = keyword_info,
        Genotype_info = genotype_info,
        Growth_info = growth_info,
        stringsAsFactors = FALSE
      )
    )
    
  }, error = function(e) {
    warning("❌ Failed to process ", gse, ": ", conditionMessage(e))
    failed_gse <<- c(failed_gse, gse)  # 收集失败的 GSE ID
  })
}

# 保存结果
saveRDS(gse_summary, "gse_summary.RDS")
write.csv(gse_summary, "gse_summary.csv", row.names = FALSE)

# 保存失败列表
if (length(failed_gse) > 0) {
  writeLines(failed_gse, "gse_failed_list.txt")
  message("⚠️ 有 ", length(failed_gse), " 个 GSE 下载失败，已记录到 gse_failed_list.txt")
} else {
  message("✅ 所有 GSE 都成功处理")
}
