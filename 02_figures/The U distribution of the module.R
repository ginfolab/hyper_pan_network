library(dplyr)
library(ggplot2)
library(dplyr)
library(ggrepel) # 引入 ggrepel 包来防止文字重叠
Gene2U_original <- read.table("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Escherichia_coli_single_network_node_U15Gephi.txt",header=T,sep="\t")
Gene2U <- Gene2U_original %>% select(Gene = id, Universality =  Weight, Module = Step_1)

Module2Size_original <- read.table('/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Module(Cyto)_relation_node(最新).txt', header=T, sep="\t")
Module2Size <- Module2Size_original %>% select(Module = Node, Module_size = Size)

Module_region <- read.table("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Module_region.txt", header=T, sep="\t")



Module2gene_table <- Module2Size %>% left_join(Gene2U, by = "Module")  %>% left_join(Module_region, by = "Module") %>%
  select(Gene, Universality, Module, Module_size, Region)

base_region_colors <- c("Region_1"="#CC6195", 
                        "Region_2"="#9EA433", 
                        "Region_3"="#6D5AC0", 
                        "Region_4"="#35A595")

library(ggplot2)
library(dplyr)

# 定义计算众数的函数
get_mode <- function(x) {
  # 去除 NA 值
  x <- na.omit(x)
  # 找出所有不重复的值
  ux <- unique(x)
  # 找出出现次数最多的那个值
  ux[which.max(tabulate(match(x, ux)))]
}

# 计算每个Module的平均Universality（散点图一个点代表一个Module）
module_summary <- Module2gene_table %>%
  group_by(Module, Region, Module_size) %>%
  summarize(Mean_Universality = mean(Universality, na.rm = TRUE),
            Mode_Universality = get_mode(Universality)) %>%
  ungroup()  %>%
  mutate(Module = gsub("^S1_M", "Module_", Module)) 



# 画散点图
p <- ggplot(module_summary, aes(x = Module_size, y = Mean_Universality, color = Region)) +
  geom_point(alpha = 0.8, size = 3) +
  
  # 关键修改：只对 Module_size > 10 的数据添加 Module 标签
  geom_text_repel(data = subset(module_summary, Module_size > 10), 
                  aes(label = Module), 
                  size = 3,            # 标签字体大小
                  max.overlaps = 20,   # 允许的最大重叠次数（可调）
                  show.legend = FALSE) + # 不把标签加到图例里
  
  scale_color_manual(values = base_region_colors) +
  scale_x_log10() + 
  theme_bw() +
  labs(title = "Module Size vs. Mean Universality",
       x = "Module Size (log scale)",
       y = "Mean Universality") +
  theme(panel.grid.minor = element_blank())

print(p)