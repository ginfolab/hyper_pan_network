#clustid2genename2Bname2describe
library(dplyr)
library(stringr)
table1 <- read.table("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Ec_gene_Name2.txt", header=T, sep="\t", quote="")

table2 <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Module_Gene_Region_Annotation.RDS")
MBGD_anno <- read.table("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/MBGD_anno.txt", header=T, sep="\t", quote="") %>%
  mutate(Descri = paste0(category,": ",Descri))

final_table <- table2 %>%
  left_join(table1 , by = c("Gene" = "clustid"))   %>%
  mutate(Module = str_replace(Module, "S1_M", "Module_")) 


anno_map <- MBGD_anno %>%
  select(category, Descri) %>%
  # 确保没有重复项，防止 join 后数据行数膨胀
  distinct(category, .keep_all = TRUE)

# 2. 定义需要替换的四列
target_cols <- c("KEGG", "COG", "TIGR", "MBGD")

# 3. 循环进行替换
final_table_translated <- final_table

for (col_name in target_cols) {
  final_table_translated <- final_table_translated %>%
    # 左连接翻译表
    left_join(anno_map, by = setNames("category", col_name)) %>%
    # 用查到的 Descri 替换原有的列，如果没查到则保留原值
    mutate(!!sym(col_name) := ifelse(!is.na(Descri), Descri, !!sym(col_name))) %>%
    # 删掉临时生成的 Descri 列，为下一轮循环做准备
    select(-Descri)
}

# 第一次运行需要先安装包（装过就不用了）
# install.packages("writexl")
library(writexl)

# 保持 .xlsx 后缀名不变
output_filename <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/paper图片/Table5 Genomic_Network_Positions_and_Annotations.xlsx"
output_path <- file.path(output_filename)

# 使用专用的 write_xlsx 函数，完美导出真正的 Excel 格式！
write_xlsx(final_table_translated, path = output_path)

message(paste(">>> 真正的 Excel 文件已成功导出至:", output_path))

