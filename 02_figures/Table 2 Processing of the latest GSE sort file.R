# 加载必要的包
library(readxl)
library(writexl)
library(dplyr)
library(purrr)
library(GEOquery) # 用于从 NCBI 抓取数据

# =========================================================================
# 第一部分：极速提取精确的样本信息 (NAR 强制要求的 Used_GSM_IDs & Sample_Count)
# =========================================================================
tpm_file <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/New_pangenome_Escherichia_coli_TPM5.xlsx"
dataset_names <- excel_sheets(tpm_file)

message("开始从 TPM 矩阵中提取实际使用的样本列表...")
sample_info_list <- map_dfr(dataset_names, function(gse) {
  # 【核心技巧】只读取第一行（列名），不加载庞大的矩阵，极大地节省内存和时间
  col_names <- colnames(read_excel(tpm_file, sheet = gse, n_max = 1))
  
  # 剔除第一列的 "Gene_ID"（或类似名称），剩下的全是我们真正用来建网的样本号 (SRR/GSM)
  actual_samples <- col_names[-1] 
  
  data.frame(
    Dataset = gse,
    Sample_Count = length(actual_samples),
    Used_GSM_IDs = paste(actual_samples, collapse = ", "),
    stringsAsFactors = FALSE
  )
})


# =========================================================================
# 第二部分：复原丢失的 GEO 数据抓取代码 (获取 Sequencing_Platform 等信息)
# =========================================================================
message("开始从 NCBI 抓取 GEO 元数据 (由于有106个数据集，请稍候几分钟)...")
geo_fetch_list <- map_dfr(dataset_names, function(gse) {
  message("  Fetching metadata for: ", gse)
  
  # 【核心技巧】GSEMatrix = FALSE 确保只极速下载文本信息，不会下载动辄几个 G 的表达矩阵
  geo_meta <- tryCatch(getGEO(gse, GSEMatrix = FALSE), error = function(e) NULL)
  
  if (is.null(geo_meta)) {
    return(data.frame(Dataset = gse, Sequencing_Platform = NA, stringsAsFactors = FALSE))
  }
  
  header <- Meta(geo_meta)
  
  # 提取测序平台 (如 GPL11154, Illumina HiSeq 等)。部分数据集存在多个平台，用逗号连接
  platforms <- ifelse(is.null(header$platform_id), NA, paste(header$platform_id, collapse = ", "))
  
  # 备注：如果你需要重新抓取摘要或标题，也可以在这里加上
  # summary_text <- paste(header$summary, collapse = " ")
  # pmid <- paste(header$pubmed_id, collapse = ", ")
  
  data.frame(
    Dataset = gse,
    Sequencing_Platform = platforms, # NAR 强制要求的新增列
    stringsAsFactors = FALSE
  )
})


# =========================================================================
# 第三部分：与原有的 NLP 和 Clade 分组数据进行终极大合并
# =========================================================================
# 1. 读取你现有的、极其珍贵的 NLP 和已整理信息表
original_gse_data <- read_excel("/Users/jiangzhenbo/Desktop/Rwork/test/GSE_info.xlsx")

# 2. 从热图中提取 Clade 分组 (保留你原有的逻辑)
data_order <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Heatmap_AutoCluster_Qc_FDR0.05_TOP0.05.RDS")
pdf(NULL) 
ht_drawn <- draw(data_order) 
dev.off() 

r_order_list <- row_order(ht_drawn)
original_rownames <- rownames(data_order@matrix)

dataset_clade_info <- map_dfr(seq_along(r_order_list), ~ data.frame(
  Dataset = original_rownames[r_order_list[[.x]]],
  Clade = paste0("Clade_", .x),
  stringsAsFactors = FALSE
))

# 3. 终极左连接：确保你的 Table 2 既有结构，又有新老数据
final_table2 <- dataset_clade_info %>%
  left_join(sample_info_list, by = "Dataset") %>%    # 加入精确样本和数量
  left_join(geo_fetch_list, by = "Dataset") %>%      # 加入新抓取的测序平台
  left_join(original_gse_data, by = "Dataset") %>%   # 完美保留你的所有 NLP 等信息
  # 调整列的顺序，把 NAR 最关心的列放在最前面，让审稿人一目了然
  select(
    Dataset, Clade, Sample_Count, Sequencing_Platform, Used_GSM_IDs,
    everything() # 把你原来表格里的所有列（摘要、NLP词条等）自动排在后面
  )

# 4. 输出最终成品
output_path <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/paper图片/Table2_GSE_info_sorted_NAR_Ready.xlsx"
write_xlsx(final_table2, path = output_path)

message(">>> 恭喜！NAR 格式的完美 Table 2 已成功生成！路径: ", output_path)