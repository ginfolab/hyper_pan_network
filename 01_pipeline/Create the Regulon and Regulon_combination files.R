################1. 首先创建Regulon & Operon 的网络 ################
setwd("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/")
library(dplyr)
library(tidyverse)
gene_name <- read.table("Gene_transform.txt",header = T,sep = "\t")

#Operon
Operon_file <- readRDS("Operon_family_Cluster2name.RDS") %>%
  rename(Operon = Operon_family,
         Operon_clustid = sp_EC_id)  %>%
  group_by(Operon) %>%
  summarise(Operon_clustid = paste(Operon_clustid, collapse = ", ")) %>%
  ungroup() %>%
  mutate(Operon_gene_count = str_count(Operon_clustid, ",") + 1)

#Regulon(gene numbers >=2)
Regulon_file <- read.table("Regulon_family.txt",header = T,sep = "\t")  %>%
  inner_join(gene_name, by = "gene",relationship = "many-to-many") %>%
  mutate(clustid = paste0("Ec_", clustid)) %>%
  select(-gene) %>%
  rename(Regulon = regulatorName,
         Regulon_clustid = clustid) %>%
  group_by(Regulon) %>%
  summarise(Regulon_clustid = paste(Regulon_clustid, collapse = ", ")) %>%
  ungroup() %>%
  mutate(Regulon_gene_count = str_count(Regulon_clustid, ",") + 1)%>%
  filter(Regulon_gene_count >=2)

#Regulon_combination
RC_file <- read.table("Regulon_family_combination_final.txt",header = T,sep = "\t")  %>%
  select(-gene_pairs) %>%
  mutate(clustid = str_split(clustid, ",")) %>%   # 拆成 list
  mutate(clustid = lapply(clustid, function(x) str_c("Ec_", str_trim(x)))) %>%  # 每个 gene 加 Ec_
  mutate(clustid = sapply(clustid, function(x) str_c(x, collapse = ", "))) %>% # 合并回去
  mutate(term = if_else(str_detect(term, "&"), term, paste0(term,"(unique)")))

names(RC_file) <- c("RC","RC_clustid","RC_gene_count")

# 计算 Regulon 和 Operon 之间的基因交集
#1. Operon2RC and Operon2uniqueRegulon
Operon2RC <- expand_grid( Operon_file, RC_file) %>%
  mutate(
    overlap_genes = map2(
      str_split(Operon_clustid, ", "), 
      str_split(RC_clustid, ", "), 
      ~ intersect(.x, .y) %>% paste(collapse = ", ")
    ) %>% 
      map_chr(~ ifelse(length(.x) > 0, .x, NA))  # 用NA替代没有重叠基因的结果
  ) %>%
  filter(!is.na(overlap_genes) & overlap_genes != "") %>%  # 删除overlap_genes为空的行
  mutate(
    overlap_genes_count = str_count(overlap_genes, ",") + 1
  ) %>%
  mutate(
    Edge_type = if_else(str_detect(RC, "&"), "Operon2RC", "Operon2uniqueR"),
    Jaccard_similarity = overlap_genes_count / (Operon_gene_count + RC_gene_count - overlap_genes_count),
    Precision_mapping = overlap_genes_count / Operon_gene_count
  )%>%
  rename(Term1 = Operon,
         Term2 = RC)  %>%
  select(Term1, Term2,Edge_type,overlap_genes_count,Jaccard_similarity,Precision_mapping) 

#2. RC- Regulon
RC2Regulon <- expand_grid(RC_file, Regulon_file) %>%
  mutate(
    overlap_genes = map2(
      str_split(RC_clustid, ", "), 
      str_split(Regulon_clustid, ", "), 
      ~ intersect(.x, .y) %>% paste(collapse = ", ")
    ) %>% 
      map_chr(~ ifelse(length(.x) > 0, .x, NA))  # 用NA替代没有重叠基因的结果
  ) %>%
  filter(!is.na(overlap_genes) & overlap_genes != "") %>%  # 删除overlap_genes为空的行
  mutate(
    overlap_genes_count = str_count(overlap_genes, ",") + 1
  ) %>%
  mutate(
    Edge_type = if_else(str_detect(RC, "&"), "RC2Regulon", "Runique2Regulon"),
    Jaccard_similarity = overlap_genes_count / (Regulon_gene_count + RC_gene_count - overlap_genes_count),
    Precision_mapping = overlap_genes_count / RC_gene_count
  )%>%
  rename(Term1 = RC,
         Term2 = Regulon)  %>%
  select(Term1, Term2,Edge_type,overlap_genes_count,Jaccard_similarity,Precision_mapping ) 
result_final <- rbind(Operon2RC,RC2Regulon)
write.table(result_final,"../Gephi/Cyto/Regulon&Operon_Edge.txt",sep = "\t",row.names = F,quote = F)
saveRDS(result_final,"/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Regulon&Operon_Edge.RDS")

#制作node文件
names(Regulon_file) <- c("Node","members","Size")
Regulon_file$Type <- "Regulon"
names(Operon_file) <- c("Node","members","Size")
Operon_file$Type <- "Operon"
names(RC_file) <- c("Node","members","Size")
RC_file <- RC_file %>% mutate(Type =  if_else(str_detect(Node, "&"), "RC", "Runique"))

Node_file <- rbind(Regulon_file,Operon_file,RC_file)
write.table(Node_file,"../Gephi/Cyto/Regulon&Operon_Node.txt",sep = "\t",row.names = F,quote = F)
################1. 首先创建Regulon & Operon 的网络 ################


################2. 其次创建最新版本Regulon table################
###############制作RC_multi——tabl e####################
library(dplyr)
library(purrr)
library(tidyr)
library(combinat)
library(tidyverse)
library(stringr)
Function_term_table <- read.table("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Regulon&Operon_Node.txt",header = T,sep = "\t")
names(Function_term_table)  <- c("Function_term", "Term_clustid",  "Term_size" , "Type" )
All_RC_table <- Function_term_table[Function_term_table$Type == "RC",]  %>%
  mutate(Function_term = str_replace_all(Function_term, "&", " & "),
         Type = "RC_only",
         Function_term_TF_num = str_count(Function_term,"&")+1) 
RC_number <- dim(All_RC_table)[1]
All_Regulon_table <- Function_term_table[Function_term_table$Type == "Regulon",]

# 初始化结果表格
All_combination_results <- tibble()  
for (i in 1:RC_number) {
  # 拆分组合名称
  Current_regulon_candidate <- strsplit(All_RC_table$Function_term[i], " & ")[[1]]
  # 获取对应的 regulon 表格
  Current_regulon_table <- All_Regulon_table %>%
    filter(Function_term %in% Current_regulon_candidate)
  # 生成基因向量
  Regulon_list <- Current_regulon_table %>%
    mutate(Gene_vector = strsplit(Term_clustid, ", ") %>% map(trimws))
  n <- nrow(Regulon_list)
  if (n > 2) {  # 至少需要3个 regulon 才能组合
    # 每组组合内部的计算
    Regulon_combination_df <- map_dfr(2:(n-1), function(k) {
      combn(1:n, k, simplify = FALSE) %>%
        map_dfr(function(idx) {
          gene_sets <- map(idx, ~ Regulon_list$Gene_vector[[.x]])
          common_genes <- reduce(gene_sets, intersect)
          if (length(common_genes) > 1) {
            tibble(
              Function_term = paste(Regulon_list$Function_term[idx], collapse = " & "),
              Term_clustid = paste(common_genes, collapse = ","),
              Term_size = length(common_genes),
              Type = "Regulon_combination",
              Source = All_RC_table$Function_term[i]  # 可选：记录原始组合名
            )
          } else {
            NULL
          }
        })
    })
    # 合并进总表格
    All_combination_results <- bind_rows(All_combination_results, Regulon_combination_df)
  }
}

All_combination_results_final <- All_combination_results %>%
  left_join(All_RC_table %>% select(Function_term, Term_size) %>% rename(Source_term_size =Term_size), by = c("Source" = "Function_term")) %>%
  filter(Term_size != Source_term_size) %>%
  distinct() 
# 保存结果
saveRDS(All_combination_results_final,"/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/All_combination_results.RDS")
# 如果需要读取之前保存的结果，可以使用以下代码
#All_combination_results_final <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/All_combination_results.RDS")
RC_muti_table <- All_combination_results_final %>%
  mutate(
    Function_term_count = str_count(Function_term, "&") + 1,
    Source_count = str_count(Source, "&") + 1
  ) %>%
  filter(Function_term_count != Source_count) %>%
  group_by(Term_clustid) %>%
  mutate(row_count = n()) %>%
  summarise(
    Term_size = first(Term_size),
    Source = first(Source),
    Function_term = if (first(row_count) == 1) {
      paste0("{", first(Function_term), "} RC_muti")
    } else {
      # 多行时，拆分所有 Function_term 后去重再合并
      all_terms <- unique(unlist(str_split(Function_term, " & ")))
      all_terms_combined <- paste(all_terms, collapse = " & ")
      paste0("{", all_terms_combined, "} RC_muti")
    },
    .groups = "drop"
  ) %>%
  mutate(
    Type = "RC_muti",
    Term_clustid = str_replace_all(Term_clustid, ",", ", ")
  ) %>%
  select(Function_term, Term_size, Type, Term_clustid)


############### 再制作RC_only ####################
Function_term_table <- read.table("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Regulon&Operon_Node.txt",header = T,sep = "\t")
names(Function_term_table)  <- c("Function_term", "Term_clustid",  "Term_size" , "Type" )
Runique_table <- Function_term_table[Function_term_table$Type == "Runique",] %>%
  mutate(
    Type = "RC_only",
    Function_term = str_replace(Function_term, "\\(unique\\)", ""))

All_RC_table <- Function_term_table[Function_term_table$Type == "RC",]  %>%
  mutate(Function_term = str_replace_all(Function_term, "&", " & "),
         Type = "RC_only") 

RC_only_table <- rbind(Runique_table,All_RC_table) %>%
  mutate(Function_term = paste0("{",Function_term,"} RC_only")) 


############### 最后制作Regulon Operon 再合并 去重####################
Regulon_Operon_candidate <- c("Regulon", "Operon")
Regulon_Operon_table <- Function_term_table[Function_term_table$Type %in% Regulon_Operon_candidate,]
Function_term_table <- rbind(RC_muti_table,RC_only_table,Regulon_Operon_table)


priority <- c("Operon", "RC_only", "RC_muti", "Regulon")
test1 <- Function_term_table %>% mutate(lable = str_replace_all(Term_clustid,"Ec_","")) %>%  # 对 Term_clustid 进行处理，去掉 Ec_ 前缀
  mutate(
    label_sorted = str_split(lable, ",\\s*") %>%      # 拆成字符向量
      map(~ sort(as.numeric(.))) %>%                  # 转成数字后排序
      map_chr(~ paste0("Ec_", .) %>% paste(collapse = ", ")) # 加 Ec_ 再拼回字符串
  )%>%
  mutate(Type_priority = match(Type, priority)) %>%   # 转换成优先级数值，数值越小优先级越高
  #  group_by(label_sorted) %>%
  #slice_min(Type_priority, n = 1, with_ties = FALSE) %>%  # 每组保留优先级最高的那一行  with_ties = FALSE：如果有多个行的最小值相同，只保留其中一行（默认是第一行），不保留“并列”的其它行。
  #ungroup() %>%
  select(-Type_priority,-Term_clustid,-lable)  %>%
  rename(Term_clustid =   label_sorted)

saveRDS(test, "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/最新的RegulonDB基因注释(dupli).RDS")
################2. 其次创建最新版本RC table################



###########################加入closed Itemset list ############################
library(data.table)
library(stringr)
library(dplyr)
# 读取并处理 Function_term_table
Function_term_table <-  as.data.table(readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/最新的RegulonDB基因注释(dupli).RDS")) %>% select(Function_term, Term_clustid,Term_size, Type)
Function_term_table[, Term_clustid := str_split(as.character(Term_clustid), ",\\s*")]

# 展开 Function_term_table：每行一个基因
Function_gene_dt <- Function_term_table[, .(gene = unlist(Term_clustid)), by = .(Function_term)]

# 处理 closed_itemset
closed_itemset <- as.data.table(readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Escherichia_coli_closed_U15.RDS"))
closed_itemset[, itemset_name := paste0("{", Itemsets, "}")]
closed_itemset[, Itemset_clustid := str_split(as.character(Itemsets), ",\\s*")]

# 展开 closed_itemset：每行一个基因
closed_gene_dt <- closed_itemset[, .(gene = unlist(Itemset_clustid)), by = .(itemset_name)]

# 用 gene 进行连接：找到所有功能项和 closed itemset 的交集
overlap_dt <- merge(Function_gene_dt, closed_gene_dt, by = "gene", allow.cartesian = TRUE)

# 统计每个 Function_term 和 itemset_name 的重叠基因数
overlap_count <- overlap_dt[, .N, by = .(Function_term, itemset_name)]

# 筛选重叠基因数 >= 2 的组合
overlap_filtered <- overlap_count[N >= 2]

# 汇总为 list
overlap_summary <- overlap_filtered[, .(
  overlap_itemsets = list(itemset_name),
  overlap_count = .N
), by = Function_term]

# 合并回 Function_term_table
Function_term_table_final <- merge(Function_term_table, overlap_summary, by = "Function_term", all.x = TRUE)

# 替换 NA 为空 list / 0
Function_term_table_final[is.na(overlap_itemsets), overlap_itemsets := list(list())]
Function_term_table_final[is.na(overlap_count), overlap_count := 0]

Function_term_table_final <- Function_term_table_final %>%
  rename(
    Closed_itemsets = overlap_itemsets, 
    Closed_count = overlap_count)

# 从 closed_itemset 中提取需要的列（itemset_name 和 Universality）
universality_info <- closed_itemset[, .(itemset_name, Universality)]

# 合并到 overlap_filtered 中，添加 Universality 信息
overlap_with_uni <- merge(overlap_filtered, universality_info, by = "itemset_name", all.x = TRUE)

# 聚合：每个 Function_term 对应一组 Universality 列表
universality_summary <- overlap_with_uni[, .(
  Closed_Universality = list(Universality)
), by = Function_term]

# 合并回最终表
Function_term_table_final <- merge(
  Function_term_table_final,
  universality_summary,
  by = "Function_term",
  all.x = TRUE
)

# 替换 NA 为 list()
Function_term_table_final[is.na(Closed_Universality), Closed_Universality := list(list())]
Function_term_table_final_closed <- Function_term_table_final


# 展开 Function_term_table：每行一个基因
#Function_gene_dt <- Function_term_table[, .(gene = unlist(Term_clustid)), by = .(Function_term)]

# 处理 closed_itemset
maximal_itemset <- as.data.table(readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Escherichia_coli_maximal_U15.RDS"))
maximal_itemset[, itemset_name := paste0("{", Itemsets, "}")]
maximal_itemset[, Itemset_clustid := str_split(as.character(Itemsets), ",\\s*")]

# 展开 closed_itemset：每行一个基因
maximal_gene_dt <- maximal_itemset[, .(gene = unlist(Itemset_clustid)), by = .(itemset_name)]

# 用 gene 进行连接：找到所有功能项和 closed itemset 的交集
overlap_dt <- merge(Function_gene_dt, maximal_gene_dt, by = "gene", allow.cartesian = TRUE)

# 统计每个 Function_term 和 itemset_name 的重叠基因数
overlap_count <- overlap_dt[, .N, by = .(Function_term, itemset_name)]

# 筛选重叠基因数 >= 2 的组合
overlap_filtered <- overlap_count[N >= 2]

# 汇总为 list
overlap_summary <- overlap_filtered[, .(
  overlap_itemsets = list(itemset_name),
  overlap_count = .N
), by = Function_term]

# 合并回 Function_term_table
Function_term_table_final <- merge(Function_term_table, overlap_summary, by = "Function_term", all.x = TRUE)

# 替换 NA 为空 list / 0
Function_term_table_final[is.na(overlap_itemsets), overlap_itemsets := list(list())]
Function_term_table_final[is.na(overlap_count), overlap_count := 0]

Function_term_table_final <- Function_term_table_final %>%
  rename(
    Maximal_itemsets = overlap_itemsets, 
    Maximal_count = overlap_count)
universality_info <- maximal_itemset[, .(itemset_name, Universality)]

# 合并到 overlap_filtered 中，添加 Universality 信息
overlap_with_uni <- merge(overlap_filtered, universality_info, by = "itemset_name", all.x = TRUE)

# 聚合：每个 Function_term 对应一组 Universality 列表
universality_summary <- overlap_with_uni[, .(
  Maximal_Universality = list(Universality)
), by = Function_term]

# 合并回最终表
Function_term_table_final <- merge(
  Function_term_table_final,
  universality_summary,
  by = "Function_term",
  all.x = TRUE
)

# 替换 NA 为 list()
Function_term_table_final[is.na(Maximal_Universality), Maximal_Universality := list(list())]

Function_term_table_final_maximal <- Function_term_table_final %>% select(Function_term, Maximal_itemsets, Maximal_count, Maximal_Universality)

Function_term_table_final <- Function_term_table_final_closed %>%
  left_join(Function_term_table_final_maximal, by = "Function_term") 

Function_term_table_final <- Function_term_table_final %>%
  mutate(Term_clustid = sapply(Term_clustid, function(x) paste(x, collapse = ", ")))

saveRDS(Function_term_table_final, "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/最新的RegulonDB基因注释(dupli).RDS")



