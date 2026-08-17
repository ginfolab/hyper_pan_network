####4. map regulonDB到 Module################
#########预先处理步骤############

# 设置 Type 的优先级顺序
#创建 Module2Function
type_priority <- c("Operon", "RC_only", "RC_muti", "Regulon")
Function_priority <- c("Best_mapping", "High_F_score", "Best_Precision", "High_Precision","Best_Recall","High_Recall")
library(dplyr)
library(purrr)
library(tidyr)
library(combinat)
library(tidyverse)
Module_relation_node1 <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Module_relation_node1.RDS")  # %>% #已经删除了没用的S3 模块了  但是有 Regulon 的模块，注意要删除 
 filter(Module_type %in% c("S1", "S3"))  #删除Regulon 模块
Module_num   <- nrow(Module_relation_node1)
Function_term_table <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/最新的RegulonDB基因注释.RDS")
#动态赋值每个Function_term_num 
Function_term_num_df <- as.data.frame(table(Function_term_table$Type)) %>%
  rename(Function_term = Var1,
         Function_term_num = Freq)
for (i in seq_len(nrow(Function_term_num_df))) {
  var_name <- paste0(Function_term_num_df$Function_term[i],"_num")
  var_value <- Function_term_num_df$Function_term_num[i]
  assign(var_name, var_value)
}

#统计所有的mapping
Module2Function  <-  expand_grid(Module_relation_node1,Function_term_table) %>%
  mutate(
    overlap_genes = map2(
      str_split(Clustid, ","),
      str_split(Term_clustid, ", "),
      ~ intersect(.x, .y) %>% paste(collapse = ",")
    )
  )%>%
  filter(overlap_genes != "") %>%  # 删除没有重叠基因的行
  mutate(overlap_genes_count = str_count(overlap_genes, ",") + 1) %>%
  mutate(
    Precision = overlap_genes_count / Module_size,
    Recall = overlap_genes_count / Term_size,
    F_score = if_else(Precision + Recall == 0, 0, 2 * Precision * Recall / (Precision + Recall))
  ) 
saveRDS(Module2Function, "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Module2Function.RDS")
#########当我还没有确定Low_cutoff   我需要画一个动态的low cut 的 coverage的折线图 ############
Bar_statistics_df <- data.frame()
Final_statistics_df <- data.frame()
Pie_statistics_df <- data.frame()
Total_module_num <- nrow(Module_relation_node1)
High_cutoff = 0.8
for (Low_cutoff in seq(0.1, 1, by = 0.02)) {
  Module2Function_candidate <- Module2Function %>%
    filter(F_score >= Low_cutoff ) %>%
    rowwise() %>%
    mutate(
      Function_mapping = case_when(
        F_score == 1 ~ "Best_mapping",
        F_score >= High_cutoff & F_score < 1 ~ "High_F_score",
        F_score >= Low_cutoff  & F_score < High_cutoff & Precision == 1 ~ "Best_Precision",
        F_score >= Low_cutoff  & F_score < High_cutoff & Precision >= High_cutoff ~ "High_Precision",
        F_score >= Low_cutoff  & F_score < High_cutoff & Recall == 1 ~ "Best_Recall",
        F_score >= Low_cutoff  & F_score < High_cutoff & Recall >= High_cutoff ~ "High_Recall"
      ),
      Weight = max(c(F_score, Recall, Precision), na.rm = TRUE)
    ) %>%
    ungroup() %>% 
    rename(
      Function_term_clustid = Term_clustid,
      Module_clustid = Clustid,
      Module = Id,
      Function_term_type = Type
    )
  
  #做第1个折线统计表
  current_F_table <- Module2Function_candidate %>%
    select(Function_term, Function_term_type,Module ) %>%
    distinct()
  
  Current_statistics_df <- as.data.frame(table(current_F_table$Function_term_type)) %>%
    rename(Function_term = Var1,
           Mapping_num = Freq) %>%
    left_join(Function_term_num_df, by = "Function_term") %>%
    mutate(Mapping_Ratio = Mapping_num / Function_term_num,
           Fscore_cut = Low_cutoff) #%>%
  #select(-Function_term_num)
  Final_statistics_df <- rbind(Final_statistics_df,Current_statistics_df)
  
  #做第2个Pie统计表
  current_Pie_table <- Module2Function_candidate %>%
    group_by(Module) %>%
    arrange(desc(F_score)) %>%
    slice(1) %>%
    ungroup() %>%
    select(Module, Function_term_type ) %>%
    distinct()
  
  Current_pie_df <- as.data.frame(table(current_Pie_table$Function_term_type)) %>%
    rename(Function_term = Var1,
           Mapping_num = Freq) %>%
    right_join(Function_term_num_df, by = "Function_term") %>%
    mutate(Mapping_num = if_else(is.na(Mapping_num), 0, Mapping_num)
    ) %>%
    bind_rows(
      tibble(
        Function_term = "No_mapping",
        Mapping_num = Module_num - sum(.$Mapping_num)
      )
    )%>%
    mutate(Fscore_cut = Low_cutoff) %>%
    select(-Function_term_num)
  Pie_statistics_df <- rbind(Pie_statistics_df,Current_pie_df)
  
  #做第3个barplot统计表
  current_bar_table <- Module2Function_candidate %>%
    count(Module) %>%
    mutate(Category = ifelse(n > 1, "Repeated", "Single")) %>%
    count(Category) %>%
    rename(Count = n) %>%
    mutate(Proportion = round(Count / Total_module_num, 3)) %>%
    { 
      bind_rows(
        .,
        tibble(
          Category = "No_mapping",
          Count = Total_module_num - sum(.$Count),
          Proportion = round(1 - sum(.$Proportion), 3)
        )
      )
    } %>%
    mutate(Fscore_cut = Low_cutoff)
 
  
  Bar_statistics_df <- rbind(Bar_statistics_df, current_bar_table)
  
}
saveRDS(Final_statistics_df, "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Final_statistics_df4Fscore.RDS")
Pie_statistics_df_final <- Pie_statistics_df %>%
  filter(sapply(Fscore_cut, function(x) any(near(x, seq(0.1, 1.0, by = 0.1)))))
saveRDS(Pie_statistics_df_final, "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Pie_statistics_df.RDS")
saveRDS(Bar_statistics_df, "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Bar_statistics_df.RDS")



################开始绘制图片 ################
#制作 模块 和 Function term 的size 山脊图
# 读取数据
Final_statistics_df <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Final_statistics_df4Fscore.RDS")
Pie_statistics_df_final <- readRDS( "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Pie_statistics_df.RDS")
Bar_statistics_df <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Bar_statistics_df.RDS")
library(ggplot2)
library(tidyr)
library(patchwork)
library(dplyr)
library(ggrepel)
library(ggridges) 
# 设置颜色映射
fill_colors <- c(
  "Operon" = "#9467bd",
  "RC_muti" = "darkred",
  "RC_only" = "#ff7f0e",
  "Regulon" = "#2ca02c"
)

  Guides_line <- Bar_statistics_df %>%
  filter(Category %in% c("Single")) %>%
  group_by(Fscore_cut) %>%
  summarise(Count_sum = sum(Count), .groups = "drop") %>%
  arrange(desc(Count_sum)) %>%
  slice(1) %>%
  pull(Fscore_cut)


Module_data <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Module_relation_node1.RDS")   %>% 
  filter(Module_type == "S1") %>% #删除S3 模块
  select(Id,  Module_size) %>%
  mutate(Type = "Module") 

RegulonDB_data <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/最新的RegulonDB基因注释.RDS") %>%
  select("Function_term","Term_size", "Type")%>%
  rename(Id = Function_term, Module_size = Term_size) 
Module_RegulonDB_data <- rbind(Module_data,RegulonDB_data) 


#第三个 画山脊图
# 计算每个 Type 的数量
label_df <- Module_RegulonDB_data %>%
  count(Type) %>%
  mutate(label = paste0("n = ", n))

p3 <- ggplot(Module_RegulonDB_data, aes(x = Module_size, y = Type, fill = Type)) +
  geom_density_ridges(scale = 1, alpha = 0.6) +
  scale_x_log10() +
  geom_text(data = label_df, aes(x = 1, y = Type, label = label), 
            inherit.aes = FALSE, hjust = 1, vjust = -0.5, size =5) +
  scale_fill_manual(values = c(
    fill_colors,
    "Module" = "gray40"
  )) +
  theme_minimal() +
  labs(
    title = "Log-scaled Module & Function term Size Distribution",
    x = "log10(Module/Term Size)",
    y = "Type"
  ) +
  theme(legend.position = "none")



#第二个 画饼图
# 假设 Pie_statistics_df 已存在并包含 Fscore_cut, Function_term, Mapping_num
# # 绘图数据准备过滤 Fscore_cut 为 0.X 的情况（确保是数值并在 0~1 范围内）
Pie_statistics_df_final  <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Pie_statistics_df.RDS")
plot_data <- Pie_statistics_df_final %>%
  mutate(Fscore_cut = as.numeric(Fscore_cut)) %>%
  arrange(Fscore_cut) %>%
  slice(which(Fscore_cut %in% sort(unique(Fscore_cut))[1:10])) %>%
  group_by(Fscore_cut) %>%
  mutate(
    percent = Mapping_num / sum(Mapping_num),
    label_pos = cumsum(Mapping_num) - Mapping_num / 2,
    percent_label = paste0(round(percent * 100), "%")
  )%>%
  mutate(percent_label2 = ifelse(percent >= 0.01, percent_label, NA)) %>%
 filter(!is.na(percent_label2)) %>%
  ungroup() 

library(patchwork)

# 第一个画不同function term 的mapping 折线图
Final_statistics_df <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Final_statistics_df4Fscore.RDS")
p1 <- ggplot(Final_statistics_df, aes(x = Fscore_cut, y = Mapping_Ratio, color = Function_term)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  scale_color_manual(
    values = c(
      fill_colors
    )
  ) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(
    title = "Distribution of Functional Terms Mapped to Modules(ManytoMany) Across F-score Cutoffs",
    x = "Fscore Cutoff",
    y = "Mapping Ratio",
    color = "Function Term"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  ) +
  geom_vline(xintercept = Guides_line, linetype = "dotted", color = "gray40", size = 0.9)

# 饼图
p2 <- ggplot(plot_data, aes(x = "", y = Mapping_num, fill = Function_term)) +
  geom_col(width = 1, color = "grey90") +
  coord_polar(theta = "y") +
  facet_wrap(~ Fscore_cut, ncol = 5) +
  scale_fill_manual(values = c(fill_colors,
                               "No_mapping" = "gray40")) +
  theme_void() +
  theme(
    strip.text = element_text(size = 10),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(title = "Function Term")) +
  geom_label_repel(
    aes(label = percent_label2, y = label_pos),
    segment.color = "white",
    color = "black",
    size = 3,
    show.legend = FALSE
  )+
  labs(title = "Functional Term Mapping in Module(1to1) Across F-score Cutoffs")


#堆积 折线 图
# 确保 Category 是有序 factor，顺序与 fill 手动指定一致
Bar_statistics_df <- Bar_statistics_df %>%
  mutate(Category = factor(Category, levels = c("No_mapping", "Repeated", "Single")))

# 提取 Fscore_cut = Guides_line 的标注数据，并计算堆积位置
label_data1 <- Bar_statistics_df %>%
  filter(Fscore_cut == Guides_line) %>%
 arrange(desc(Category))%>%
  mutate(
    y_pos = cumsum(Proportion) - Proportion / 2,
    label = paste0(Category, ": ", Count)
  )

label_data2 <- Bar_statistics_df %>%
  filter(Fscore_cut == 0.8) %>%
  arrange(desc(Category))%>%
  mutate(
    y_pos = cumsum(Proportion) - Proportion / 2,
    label = paste0(Category, ": ", Count)
  )


p4 <- ggplot(Bar_statistics_df, aes(x = Fscore_cut, y = Proportion, fill = Category)) +
  geom_area(position = "stack") +
  geom_vline(xintercept = Guides_line, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "black") +
  geom_text(
    data = label_data1,
    aes(x =Guides_line, y = y_pos, label = label),
    inherit.aes = FALSE,
    hjust = -0.1,
    size = 3.5
  )   +
  annotate("text", 
           x = Guides_line, 
           y = 1.02, 
           label = paste0("Cutoff = ", round(Guides_line, 2)), 
           angle = 90, 
           vjust = -0.5, 
           hjust = 2.5,
           size = 3.5, 
           color = "black")  +
  
  geom_text(
    data = label_data2,
    aes(x =0.8, y = y_pos, label = label),
    inherit.aes = FALSE,
    hjust = -0.1,
    size = 3.5
  )   +
  annotate("text", 
           x = 0.8, 
           y = 1.02, 
           label = paste0("Cutoff = ", round(0.8, 2)), 
           angle = 90, 
           vjust = -0.5, 
           hjust = 2.5,
           size = 3.5, 
           color = "black")  +
  scale_fill_manual(values = c("Single" = "#1f77b4", "Repeated" = "#ff7f0e", "No_mapping" = "gray40")) +
  labs(
    title = "Proportion of Modules mapping with Single / Repeated / None Functional Terms",
    x = "Fscore Cutoff",
    y = "Proportion",
    fill = "Category"
  ) +
  xlim(min(Bar_statistics_df$Fscore_cut), max(Bar_statistics_df$Fscore_cut) + 0.1) +  # 给文字留点右边距
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )


# 组合，p1 在上，p2 在下
combined_plot <- (p1 | p4) / (p2 | p3) + plot_layout(heights = c(1, 1))
print(combined_plot)
################开始绘制图片 ################






############### 桑基图制作 ######################
###当我确定了Low_cutoff
library(dplyr)

Module2Function <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Module2Function.RDS")  %>%
      filter(Module_type %in% c("S1", "S3")) #删除Regulon 模块
# 读取 Module2Function 数据#读取function term 之后map到module#读取function term 之后 跟新的RC进行整合
High_cutoff = 0.8        # 设置高 F-score 截断值  
Low_cutoff = 0.48         # 设置低 F-score 截断值    这个截断值由之前的那个动态折线图决定
# 设置 Type 的优先级顺序
Module2Function_candidate <- Module2Function %>%
  filter(F_score >= High_cutoff | (Precision >= High_cutoff & F_score >= Low_cutoff)| (Recall >= High_cutoff & F_score >= Low_cutoff) ) %>%
  rowwise() %>%
  mutate(
    Function_mapping = case_when(
      F_score == 1 ~ "Best_mapping",
      F_score >= High_cutoff & F_score < 1 ~ "High_F_score",
      #F_score >= Low_cutoff  & F_score < High_cutoff & Precision == 1 ~ "Best_Precision",
      F_score >= Low_cutoff  & F_score < High_cutoff & Precision >= High_cutoff ~ "High_Precision",
      # F_score >= Low_cutoff  & F_score < High_cutoff & Recall == 1 ~ "Best_Recall",
      F_score >= Low_cutoff  & F_score < High_cutoff & Recall >= High_cutoff ~ "High_Recall",
    ),
    Weight = max(c(F_score, Recall, Precision), na.rm = TRUE)
  ) %>%
  ungroup() %>% 
  rename(
    Function_term_clustid = Term_clustid,
    Module_clustid = Clustid,
    Module = Id,
    Function_term_type = Type
  )

Function_priority <- c("Best_mapping", "High_F_score", "High_Precision","High_Recall")

#留作以后注释用
saveRDS(Module2Function_candidate,paste0("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/Module2Function_candidate4RegulonDB(F-score = H",High_cutoff,"_L",Low_cutoff,").RDS"))

#############制作RegulonDB 中的 Regulon > RCmuti > RC_only > Operon 的关系框架############
#########以Function term 为点绘制泛网络
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
RegulonDB_file <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/最新的RegulonDB基因注释.RDS")
Operon_file <- RegulonDB_file %>%   filter(Type == "Operon") %>% select(-Type) %>%  rename( Operon = Function_term,  Operon_size = Term_size,   Operon_clustid = Term_clustid)
RCmuti_file <- RegulonDB_file %>% filter(Type == "RC_muti")%>%  select(-Type) %>%  rename( RCmuti =  Function_term,  RCmuti_size = Term_size,   RCmuti_clustid = Term_clustid)
Regulon_file <- RegulonDB_file %>%   filter(Type == "Regulon")%>% select(-Type) %>% rename(Regulon = Function_term,  Regulon_size = Term_size,  Regulon_clustid = Term_clustid)
RConly_file <- RegulonDB_file %>%  filter(Type == "RC_only") %>% select(-Type) %>% rename( RConly =  Function_term,  RConly_size = Term_size,   RConly_clustid = Term_clustid)

Node_file <- data.frame()
###########链接文件 成为EDGE ###############
Opero2RConly <-   expand_grid(Operon_file, RConly_file) %>%
  mutate(
    overlap_genes = map2(
      str_split(Operon_clustid, ", "),
      str_split(RConly_clustid, ", "),
      ~ intersect(.x, .y) %>% paste(collapse = ", ")
    ) %>%
      map_chr(~ ifelse(length(.x) > 0, .x, NA))  # 用NA替代没有重叠基因的结果
  ) %>%
  filter(!is.na(overlap_genes) & overlap_genes != "") %>%  # 删除overlap_genes为空的行
  mutate(
    overlap_genes_count = str_count(overlap_genes, ",") + 1
  ) %>%
  mutate(Jaccard = overlap_genes_count / (RConly_size + Operon_size - overlap_genes_count),
         Edge_type = "Opero2RConly",
         Recall = overlap_genes_count / RConly_size) %>%
  select(-Operon_clustid, -RConly_clustid,-overlap_genes,Operon_size,RConly_size)  %>%
  mutate(O2RCo = paste0(Operon , " - ", RConly))
Opero2RConly_modi <- Opero2RConly %>% select(TAG = O2RCo, Jaccard,Edge_type )

RConly2RCmuti <-   expand_grid(RConly_file,RCmuti_file ) %>%
  mutate(
    overlap_genes = map2(
      str_split(RCmuti_clustid, ", "),
      str_split(RConly_clustid, ", "),
      ~ intersect(.x, .y) %>% paste(collapse = ", ")
    ) %>%
      map_chr(~ ifelse(length(.x) > 0, .x, NA))  # 用NA替代没有重叠基因的结果
  ) %>%
 filter(!is.na(overlap_genes) & overlap_genes != "") %>%  # 删除overlap_genes为空的行
  mutate(
    overlap_genes_count = str_count(overlap_genes, ",") + 1
  ) %>%
  mutate(Jaccard = overlap_genes_count / (RConly_size + RCmuti_size - overlap_genes_count),
         Edge_type = "RConly2RC_muti",
         Recall = overlap_genes_count / RCmuti_size,
         RConly_TF = str_remove_all(RConly, "\\} RC_only|\\{"),# 删除 RConly 中的无用字符串
         RConly_TF = str_split(RConly_TF, " & "),
         RCmuti_TF = str_remove_all(RCmuti, "\\} RC_muti|\\{"), # 删除 RCmuti 中的无用字符串
         RCmuti_TF = str_split(RCmuti_TF, " & ")
  ) %>%
  filter(
    map2_lgl(RConly_TF, RCmuti_TF, ~ all(.y %in% .x))        #保留 RConly_TF 中完全包含 RCmuti_TF 的行
  ) %>%
  select(-RCmuti_clustid, -RConly_clustid,-overlap_genes,-RConly_TF,-RCmuti_TF) %>%
  mutate(RCo2RCm = paste0(RConly, " - ", RCmuti))
RConly2RCmuti_modi <- RConly2RCmuti %>% select(TAG = RCo2RCm, Jaccard,Edge_type )

RCmuti2Regulon <- expand_grid(RCmuti_file,Regulon_file ) %>%
  mutate(
    overlap_genes = map2(
      str_split(Regulon_clustid, ", "),
      str_split(RCmuti_clustid, ", "),
      ~ intersect(.x, .y) %>% paste(collapse = ", ")
    ) %>%
      map_chr(~ ifelse(length(.x) > 0, .x, NA))  # 用NA替代没有重叠基因的结果
  ) %>%
  filter(!is.na(overlap_genes) & overlap_genes != "") %>%  # 删除overlap_genes为空的行
  mutate(
    overlap_genes_count = str_count(overlap_genes, ",") + 1
  ) %>%
  mutate(Jaccard = overlap_genes_count / (RCmuti_size + Regulon_size - overlap_genes_count),
         Edge_type = "RCmuti2Regulon",
         RCmuti_TF = str_remove_all(RCmuti, "\\} RC_muti|\\{"), # 删除 RCmuti 中的无用字符串
         RCmuti_TF = str_split(RCmuti_TF, " & "),
         Recall = overlap_genes_count / Regulon_size
  ) %>%
  
  filter(
    map2_lgl(RCmuti_TF, Regulon, ~ all(.y %in% .x))        #保留 RConly_TF 中完全包含 RCmuti_TF 的行
  ) %>%
  select(-Regulon_clustid, -RCmuti_clustid,-overlap_genes, -RCmuti_TF)   %>%  
  mutate(RCm2R = paste0(RCmuti, " - ", Regulon))
RCmuti2Regulon_modi <- RCmuti2Regulon %>% select(TAG = RCm2R, Jaccard,Edge_type )

# 将所有的 pairwies edge 文件合并并记录Jaccard 和Edge_type
allmodi <- bind_rows(Opero2RConly_modi, RConly2RCmuti_modi, RCmuti2Regulon_modi)
saveRDS(allmodi, "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/RegulonDB_edge_modi.RDS")


############将所有的edge文件合并 （忽略）############
Opero2RConly_modi <- Opero2RConly %>%
  select(Operon, RConly, O2RCo)  
RConly2RCmuti_modi <- RConly2RCmuti %>%
  select(RConly, RCmuti, RCo2RCm) 
RCmuti2Regulon_modi <- RCmuti2Regulon %>%
  select(RCmuti, Regulon, RCm2R) 
library(dplyr)
ALL_join_file <- Opero2RConly_modi %>%
  full_join(RConly2RCmuti_modi, by = "RConly", relationship = "many-to-many") %>%
  full_join(RCmuti2Regulon_modi, by = "RCmuti", relationship = "many-to-many")  %>%
  mutate(overlap_genes_count = NA_real_,   Jaccard = NA_real_,   Edge_type = NA_character_,   Recall = NA_real_) %>% # 初始化新列为 NA 先占个地方 
  select("Operon","RConly","RCmuti", "Regulon" , "O2RCo","RCo2RCm","RCm2R",
         "overlap_genes_count", "Jaccard", "Edge_type", "Recall")

# 第一步：准备 key 到信息的映射表
O2RCo_map <- Opero2RConly %>%
  select(O2RCo, overlap_genes_count, Jaccard, Edge_type, Recall)
RCo2RCm_map <- RConly2RCmuti %>%
  select(RCo2RCm, overlap_genes_count, Jaccard, Edge_type, Recall)
RCm2R_map <- RCmuti2Regulon %>%
  select(RCm2R, overlap_genes_count, Jaccard, Edge_type, Recall)

# 第2步：O2RCo 不为空的情况 
ALL_join_file1 <- ALL_join_file %>%
  left_join(O2RCo_map, by = c("O2RCo")) %>%
  mutate(
    overlap_genes_count = coalesce(overlap_genes_count.x, overlap_genes_count.y),         #coalesce() 是一个非常有用的函数，用于返回多个列中第一个非缺失值（非 NA）
    Jaccard = coalesce(Jaccard.x, Jaccard.y),   
    Edge_type = coalesce(Edge_type.x, Edge_type.y),
    Recall = coalesce(Recall.x, Recall.y)
  ) %>%
  select(-ends_with(".x"), -ends_with(".y"))

# 第3步：O2RCo 是空值且 RCo2RCm 非空
ALL_join_file2 <- ALL_join_file1 %>%
  left_join(
    RCo2RCm_map %>% rename_with(~ paste0(., "_rco2rcm"), -RCo2RCm),    #将 RCo2RCm_map 中除了 RCo2RCm 这一列以外的所有列都加上后缀 _rcm2r
    by = c("RCo2RCm")
  ) %>%
  mutate(
    overlap_genes_count = if_else(is.na(overlap_genes_count), overlap_genes_count_rco2rcm, overlap_genes_count),
    Jaccard = if_else(is.na(Jaccard), Jaccard_rco2rcm, Jaccard),
    Edge_type = if_else(is.na(Edge_type), Edge_type_rco2rcm, Edge_type),
    Recall = if_else(is.na(Recall), Recall_rco2rcm, Recall)
  ) %>%
  select(-ends_with("_rco2rcm"))

# 第4步：前两步为空但 RCm2R 非空
ALL_join_file3 <- ALL_join_file2 %>%
  left_join(
    RCm2R_map %>% rename_with(~ paste0(., "_rcm2r"), -RCm2R),  #将 RCm2R_map 中除了 RCm2R 这一列以外的所有列都加上后缀 _rcm2r
    by = c("RCm2R")
  ) %>%
  mutate(
    overlap_genes_count = if_else(is.na(overlap_genes_count), overlap_genes_count_rcm2r, overlap_genes_count),
    Jaccard = if_else(is.na(Jaccard), Jaccard_rcm2r, Jaccard),
    Edge_type = if_else(is.na(Edge_type), Edge_type_rcm2r, Edge_type),
    Recall = if_else(is.na(Recall), Recall_rcm2r, Recall)
  ) %>%
  select(-ends_with("_rcm2r"))

saveRDS(ALL_join_file3, "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/Regulondataframe4layer.RDS")
############将所有的edge文件合并（忽略）############  
#############制作RegulonDB 中的 Regulon > RCmuti > RC_only > Operon 的关系框架############



###############1. 选择并填充候选线路###############
library(dplyr)
High_cutoff = 0.8        # 设置高 F-score 截断值  
Low_cutoff = 0.48         # 设置低 F-score 截断值    这个截断值由之前的那个动态折线图决定
Function_priority <- c("Best_mapping", "High_F_score", "High_Precision","High_Recall")
Module2Function_candidate <- readRDS(paste0("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/Module2Function_candidate4RegulonDB(F-score = H",High_cutoff,"_L",Low_cutoff,").RDS")) %>%
 # select(Module, Function_term) #%>%
 mutate(TAG = paste0(Function_term,"_(",Module,")")) # 将 Function_term 列拆分为列表

#保留候选Terms 线路
candidate_terms <-  unique(Module2Function_candidate$Function_term)

candidate_map  <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/Regulondataframe4layer.RDS") %>%      
rowwise() %>%
  mutate(
    match_vec = list(c_across(c(Operon, RConly, RCmuti, Regulon)) %in% candidate_terms),
    match_sum = sum(unlist(match_vec)),
    has_adjacent_match = any(diff(which(unlist(match_vec))) == 1)
  ) %>%
  ungroup() %>%
  #filter(match_sum >= 1) #%>% # 保留至少一个匹配的行
filter(match_sum >= 2 & has_adjacent_match) #%>%      保留至少2个匹配的行 而且还要连续
  #select(-match_vec, -match_sum, -has_adjacent_match)

##读取整合candidate cascade文件，因为包含了很多空值 所以都要用随机数代替空值
# 自定义函数：替换 NA 为唯一的 pseudo ID
replace_na_with_pseudo <- function(vec, layer) {
  na_idx <- which(is.na(vec))
  vec[na_idx] <- paste0("pseudo_", layer, "_", seq_along(na_idx))
  return(vec)
}     # 自定义函数：替换NA为唯一的 pseudo ID
candidate_map <- candidate_map %>%
  mutate(
    Operon = replace_na_with_pseudo(Operon, "Operon"),
    RConly = replace_na_with_pseudo(RConly, "RConly"),
    RCmuti = replace_na_with_pseudo(RCmuti, "RCmuti"),
    Regulon = replace_na_with_pseudo(Regulon, "Regulon")
  )
###############选择并填充候选线路###############




#####################2. 创建Node文件#########################
######读取MOdule的注释RegulonDB信息  . 制作map#################
library(igraph)
library(tidygraph)
library(stringr)
library(dplyr)
library(tidyr)
library(purrr)  

# 选择 Function_term 和 F_score 列，并按 F_score 降序排列，保留每个 Function_term 的最高 F_score （保留每个Function term的唯一值）
mapping_F <- Module2Function_candidate %>%
  select(Function_term, F_score,Function_mapping,Module,Module_type)  %>%
  group_by(Function_term) %>%
  arrange(desc(F_score)) %>%
  slice(1) 

#做点文件 修改 nodes 表，添加是否为 pseudo 节点的标识
level_order <- c("all","Operon", "RConly", "RCmuti", "Regulon")
RegulonDB_file <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/最新的RegulonDB基因注释.RDS")
nodes <- candidate_map %>%
  pivot_longer(cols = c(Operon, RConly, RCmuti, Regulon), 
               names_to = "Type", 
               values_to = "Function_term") %>%
  unique() %>%
 # mutate(pseudo = if_else(Function_term %in% candidate_terms, FALSE , TRUE )) %>%
  mutate(pseudo = if_else(str_detect(Function_term,"pseudo_"), TRUE , FALSE )) %>%
  left_join(RegulonDB_file %>% select(Function_term, Term_size), by = "Function_term")  %>%
  mutate(annotation = if_else(Function_term %in% Module2Function_candidate$Function_term, 
                              "Yes", 
                              "No"),
         Type = factor(Type, levels = level_order),
         node.short_name = if_else(
           str_length(Function_term) > 5,
           paste0(str_sub(Function_term, 1, 5), "..."),
           Function_term)) %>%   # 先在这里创建 node.short_name
  select(node.name = Function_term, 
         node.branch = Type, 
         node.size = Term_size,
         node.short_name,
         color = annotation, 
         node.level = Type,
         pseudo)%>%
  add_row(node.name = "all",node.branch = "all", node.size = 1000,
          node.short_name = "RegulonDB", color = "No", node.level = "all", pseudo = TRUE)  %>%
  rowwise() %>%
  mutate(node.size = if_else(
    is.na(node.size) & node.branch == "Operon", runif(1, 1, 2),
    if_else(
      is.na(node.size) & node.branch == "RConly", runif(1, 2, 3),
      if_else(
        is.na(node.size) & node.branch == "RCmuti", runif(1, 3, 4),
        if_else(
          is.na(node.size) & node.branch == "Regulon", runif(1, 4, 5),
          node.size
        )
      )
    )
  )) %>%
  ungroup()  %>%
  mutate(node.anno = if_else(color == "Yes",1, 0.1))  %>%
  left_join(mapping_F, by = c("node.name" = "Function_term"))  %>%
  mutate(F_score = if_else(is.na(F_score), 0, F_score),
         Function_mapping = if_else(is.na(Function_mapping), "No_mapping", Function_mapping),
         Module = if_else(is.na(Module), "No_mapping", Module),
         Module_type = if_else(is.na(Module_type), "No_mapping", Module_type)) %>%
  unique()



# 更新 node.level 为数字
nodes$node.level <- match(nodes$node.level , level_order) 
nodes$color <- factor(nodes$color, levels = unique(nodes$color))


#做EDGE文件  重新构建 edges：同时去除 pseudo 节点边连接（后面画图时会 filter ） 
Edge_file <- bind_rows(
  candidate_map %>%  select(from = Operon   , to = RConly) %>%  mutate(Edge_type = "Operon2RConly"),
  candidate_map %>%  select(from =RConly   , to = RCmuti) %>%  mutate(Edge_type = "RConly2RCmuti"),
  candidate_map %>%  select(from =RCmuti , to =  Regulon ) %>%  mutate(Edge_type = "RCmuti2Regulon")
) %>%
  unique()

Edge_file2 <- Edge_file %>%
 select(from, to) %>%
  mutate(from = as.character(from),
         to = as.character(to))

# 最内层连中心点    # 创建一个新的 tibble，表示从 "all" 到 Operons 的连接
Operons <- nodes %>%
  filter(node.level == match("Operon", level_order)) %>%
  pull(node.name)

to_center <- tibble(
  from =  "all",
  to =  Operons
)

#pairwies edge 文件，记录Jaccard 和Edge_type
Edge_modi <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/RegulonDB_edge_modi.RDS")
edges <- bind_rows(Edge_file2, to_center) %>%
  #filter(!(from == "RegulonDB" & to == "RegulonDB")) %>%
  drop_na() %>%
  unique()  %>%
  mutate(TAG = paste(from, to, sep = " - ") ) %>%
  left_join(Edge_modi, by = "TAG")  %>%
  mutate(Jaccard = if_else(is.na(Jaccard), 0, Jaccard),
         Edge_type = if_else(is.na(Edge_type), "None", Edge_type),
         Edge_transpancy  = case_when(
           from %in% candidate_terms & to %in% candidate_terms ~ 1,
           from %in% candidate_terms | to %in% candidate_terms ~ 0.1,
           TRUE ~ 0.0001
         )
  )


graph <- tbl_graph(nodes, edges)%>%
  activate(edges) %>%
  mutate(Jaccard = edges$Jaccard,
         Edge_type = Edge_type,
         Edge_transpancy = Edge_transpancy
        )


#fill_colors <- c(
#  "Operon" = "#9467bd",
#  "RConly" = "#ff7f0e",   # 注意这里用的是你之前代码里的 RConly 名称
#  "RCmuti" = "darkred",
#  "Regulon" = "#2ca02c"
#)



fill_colors <- c(
  "Best_mapping" = "darkred",
  "High_F_score" = "brown1",   # 注意这里用的是你之前代码里的 RConly 名称
  #"Best_Recall" = "darkgreen",
  "High_Recall" = "greenyellow",
  #"Best_Precision" = "darkblue",
  "High_Precision" = "steelblue1"
)



library(ggraph)

nodes$Function_mapping <- factor(
  nodes$Function_mapping,
  levels = Function_priority
)

# 创建层级标签数据
layer_labels <- data.frame(
  x = c(0.9, 0.6,  0.35, 0),  # 半径估计
  y = 0,
  label = c("Regulon",  "RC_muti", "RC_only", "Operon")
)

library(ggforce)
ggraph(graph, layout = 'tree', circular = TRUE) +                                     #partition tree 
  geom_edge_diagonal(aes(color = Edge_type,
                         alpha = Edge_transpancy ,
                         edge_width = Jaccard,
                         filter = !node1.pseudo & !node2.pseudo & node1.node.name != "all" & node2.node.name != "all",
                         )
  ) + 
  geom_node_point(aes(size = node.size, 
                      fill = Function_mapping,   
                      alpha = F_score,
                      shape = Module_type, 
                      filter = !pseudo
                      ),
                        
                  color = "transparent") + 
  # 加上标签：只显示 color 为 Yes 的节点
  geom_node_text(
    aes(label = Module, filter = color == "Yes"),                   #node.short_name  Module
    size = 1,
    repel = TRUE,
    check_overlap = TRUE,
    max.overlaps = 10 # 不限制重叠标签数量
  ) + 

  scale_fill_manual(values = fill_colors) +  
  scale_shape_manual(values = c("S1" = 21, "S2" = 22, "S3" = 23, "No_mapping" = 24))+
  scale_size(range = c(1,10)) + 
  scale_edge_width(range = c(0.01, 0.5)) +  
  scale_edge_color_manual(values = c(
    "Opero2RConly" =  "#9467bd",
    "RConly2RC_muti" = "palegoldenrod",
    "RCmuti2Regulon" = "darkcyan",
    "None" = "grey100"
  )) + 
  theme_void() +  
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA)
  )+
# ⭐ 添加圆弧标注：对每一层（node.branch）做环状标记
  geom_text(data = layer_labels, 
            aes(x = x, y = y, label = label), 
            inherit.aes = FALSE,
            angle = 0, hjust = 0.5, 
            size = 5,
            color = "darkred",
            alpha = 0.4)  +
  
  # ✅ 设置图例，强制显示 size / fill / alpha
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 5)),
    size = guide_legend(override.aes = list(shape = 21, fill = "grey")),
    alpha = guide_legend(override.aes = list(shape = 21, fill = "grey")),
    shape = guide_legend()  # shape 也保留图例
  )


#说几点阐述一下 我们的通路只关注有联系的点。如果是孤立的点。我们不关注 


#最后创建一个function Term mapping 统计表
library(dplyr)

High_cutoff = 0.8        # 设置高 F-score 截断值  
Low_cutoff = 0.48         # 设置低 F-score 截断值    这个截断值由之前的那个动态折线图决定
Function_priority <- c("Best_mapping", "High_F_score", "High_Precision","High_Recall")
Module2Function_candidate <- readRDS(paste0("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/Module2Function_candidate4RegulonDB(F-score = H",High_cutoff,"_L",Low_cutoff,").RDS")) %>%
  # select(Module, Function_term) #%>%
  mutate(TAG = paste0(Function_term,"_(",Module,")")) # # 将 Function_term 列拆分为列表
 

Candidate_relation <- Module2Function_candidate %>%
  mutate(Module_info = paste0(Module, " (", round(F_score,2), ")")) %>%
  select(Function_term,Module_info) %>%
  group_by(Function_term) %>%
  summarise(Module_list = paste0(Module_info, collapse = ","), .groups = "drop")


candidate_map  <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/Regulondataframe4layer.RDS") %>%      
  rowwise() %>%
  mutate(
    match_vec = list(c_across(c(Operon, RConly, RCmuti, Regulon)) %in% candidate_terms),
    match_sum = sum(unlist(match_vec)),
    has_adjacent_match = any(diff(which(unlist(match_vec))) == 1)
  ) %>%
  ungroup() %>%
  #filter(match_sum >= 1) #%>% # 保留至少一个匹配的行
  filter(match_sum >= 2 & has_adjacent_match) %>%
  select(Operon,RConly,RCmuti,Regulon) %>%
  left_join(Candidate_relation, by = c("Operon" = "Function_term")) %>%
  rename(Operon_mapping_Module = Module_list) %>%
  left_join(Candidate_relation, by = c("RConly" = "Function_term")) %>%
  rename(RConly_mapping_Module = Module_list) %>%
  left_join(Candidate_relation, by = c("RCmuti" = "Function_term")) %>%
  rename(RCmuti_mapping_Module = Module_list) %>%
  left_join(Candidate_relation, by = c("Regulon" = "Function_term")) %>%
  rename(Regulon_mapping_Module = Module_list) %>%
  mutate(
    Operon_mapping_Module = ifelse(is.na(Operon_mapping_Module), "No mapping", Operon_mapping_Module),
    RConly_mapping_Module = ifelse(is.na(RConly_mapping_Module), "No mapping", RConly_mapping_Module),
    RCmuti_mapping_Module = ifelse(is.na(RCmuti_mapping_Module), "No mapping", RCmuti_mapping_Module),
    Regulon_mapping_Module = ifelse(is.na(Regulon_mapping_Module), "No mapping", Regulon_mapping_Module)
  ) 
  
saveRDS(candidate_map,paste0("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/candidate_module_map.RDS"))


