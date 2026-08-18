################3. 最后创建泛网络模块和Item set的关系网络 ################
#############################
#1.创建模块和模块的包含关系#
library(dplyr)
library(tidyverse)
Module_info <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/clustid2Module2Step.RDS")
Module_info_final <- data.frame()

for (variable in 2:length(names(Module_info))) {
  Final_module_info  <- Module_info[,c(1,variable)]
  names(Final_module_info) <- c("Clustid","Module")
  Final_module_info_merged <- Final_module_info %>%
    group_by(Module) %>%
    summarise(
      Clustid = paste(Clustid, collapse = ","),
      Module_size = n()  # 计算每个 Module 的 Clustid 数量
    ) %>%
    ungroup()%>%
    mutate(Module_type = paste0("S",variable - 1))%>%
    mutate(Interval = paste0("<[",variable - 1,",",variable,"]>"))
  Module_info_final <- rbind(Module_info_final,Final_module_info_merged)
}

test_module <- Module_info_final %>% select(-Interval)

Module_relation <- expand_grid(
  test_module %>% rename(Module1 = Module, Clustid1 = Clustid, Size1 = Module_size, Type1 = Module_type),
  test_module %>% rename(Module2 = Module, Clustid2 = Clustid, Size2 = Module_size, Type2 = Module_type)
) %>%
  filter(Module1 != Module2) %>%  # 过滤掉相同模块的组合
  mutate(
    overlap_genes = map2(
      str_split(Clustid1, ","),
      str_split(Clustid2, ","),
      ~ intersect(.x, .y) %>% paste(collapse = ",")
    )
  ) %>%
  filter(overlap_genes != "") %>%  # 删除没有重叠基因的行
  select(Module1, Module2, overlap_genes)  %>%
  mutate(Module1_trimmed  = str_remove(Module1, "_.*"  )) %>%
  mutate(Module2_trimmed  = str_remove(Module2, "_.*"  ) )%>%
  mutate(trimmed_final  = paste0(Module1_trimmed,"_",Module2_trimmed))
#给module设置分组
S3_group <- data.frame(
  Module =unique(Module_info$Step_3),
  Group = unique(Module_info$Step_3)
)
Group_info <- Module_relation %>%      
  filter(trimmed_final %in% c("S1_S3", "S2_S3")) %>%    
  select(Module1,Module2)%>%  
  rename(Module = Module1,
         Group = Module2 ) %>%  
  bind_rows(S3_group)

Module_info_final <- Module_info_final %>% 
  left_join(Group_info,by = "Module") %>% 
  rename(Id = Module) %>% 
  mutate(Group = if_else(is.na(Group) | Group == "", Id, Group)) 
#saveRDS(Module_info_final,"/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Module_info_final.RDS")

Module_relation_final <- Module_relation %>%        
  #filter(!trimmed_final %in% c("S1_S3", "S3_S1", "S2_S1", "S3_S2"))  %>%  
  filter(trimmed_final %in% c("S1_S3"))  %>%  
  mutate(overlap_genes_count = str_count(overlap_genes, ",") + 1,
         Interval = if_else(trimmed_final == "S1_S2", "<[1,2]>", "<[2,3]>"),
  )%>%  
  select(-Module1_trimmed,-Module2_trimmed,-overlap_genes,trimmed_final) %>%  
  rename(Source = Module1,
         Target = Module2)
#删除无变化Module
S3_frequency <- as.data.frame(table(Module_relation_final$Target))
Module_delelation <- Module_delelation <- as.character(S3_frequency[S3_frequency$Freq == 1,1])
Module_relation_node <- Module_info_final %>%
  filter(!Id %in% Module_delelation,
         !Module_type %in% "S2") %>%
  select(-Interval)
#write.table(Module_relation_node,"/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Module(Cyto)_relation_node.txt",sep = "\t",quote = F,row.names = F)
saveRDS(Module_relation_node,"/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Pan_network_Module_Info(删除无变化module).RDS")
Modulesize <- test_module %>%
  filter(Module_type != "S2")%>%
  select("Module","Module_size")

Module_relation_edge <- Module_relation_final %>%
  mutate(Target = if_else(Target %in% Module_delelation, Source, Target))  %>%
  mutate(Edge_type = if_else(Target == Source,"Self_loop","Normal"))%>%
  left_join(Modulesize,by = c("Source"= "Module")) %>%
  rename(Source_size = Module_size)%>%
  left_join(Modulesize,by = c("Target"= "Module")) %>%
  rename(Target_size = Module_size)%>%
  mutate(Source2Target_ratio = Source_size/Target_size)

#write.table(Module_relation_edge,"/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Module(Cyto)_relation_edge.txt",sep = "\t",quote = F,row.names = F)
saveRDS(Module_relation_edge,"/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Module_relation_edge(hierarchy).RDS")