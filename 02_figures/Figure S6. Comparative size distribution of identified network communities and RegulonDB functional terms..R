#制作 模块 和 Function term 的size 山脊图
# 读取数据
Module_data <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Module_relation_node1.RDS")   %>% 
  filter(Module_type == "S1") %>% #删除S3 模块
  select(Id,  Module_size) %>%
  mutate(Type = "Module") 

RegulonDB_data <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/RegulonDB/最新的RegulonDB基因注释.RDS") %>%
  select("Function_term","Term_size", "Type")%>%
  rename(Id = Function_term, Module_size = Term_size) 
Module_RegulonDB_data <- rbind(Module_data,RegulonDB_data) 

#画 山脊图
library(dplyr)
library(ggplot2)
library(ggridges)

# 计算每个 Type 的数量
label_df <- Module_RegulonDB_data %>%
  count(Type) %>%
  mutate(label = paste0("n = ", n))

# 绘图
p3 <- ggplot(Module_RegulonDB_data, aes(x = Module_size, y = Type, fill = Type)) +
  geom_density_ridges(scale = 1, alpha = 0.6) +
  scale_x_log10() +
  geom_text(data = label_df, aes(x = 1, y = Type, label = label), 
            inherit.aes = FALSE, hjust = 1, vjust = -0.5, size =5) +
  scale_fill_manual(values = c(
    "Operon" = "#9467bd",
    "RC_muti" = "darkred",
    "RC_only" = "#ff7f0e",
    "Regulon" = "#2ca02c",
    "Module" = "gray40"
  )) +
  theme_minimal() +
  labs(
    title = "Log-scaled Module & Function term Size Distribution",
    x = "log10(Module Size)",
    y = "Type"
  ) +
  theme(legend.position = "none")

