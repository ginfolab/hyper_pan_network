# ==============================================================================
# 终极全景映射评估图：Module视角 + Operon视角 + Regulon视角 双环饼图拼图
# ==============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggnewscale)
  library(ggrepel)
  library(patchwork) # 拼图神器
})

# ==============================================================================
# 1. 全局参数与数据读取
# ==============================================================================
base_dir <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli"
Perfect_F_Cut  <- 0.8  
minOverlapSize <- 1   

message(">>> 正在读取基础数据...")
# 1. 模块与功能的映射表
Module2Function <- readRDS(file.path(base_dir, "Gephi/Module2Function(Remove_noK12).RDS"))

# 2. 模块基本信息 (用于 S7)
Module_info <- readRDS(file.path(base_dir, "Gephi/Module_relation_node1.RDS")) %>%
  filter(Module_type == "S1") %>% select(Id)

# 3. 所有的 Operon 和 Regulon 背景库
All_Annotations <- readRDS(file.path(base_dir, "RegulonDB/最新的RegulonDB基因注释.RDS"))

# ==============================================================================
# Part 1: Module 视角图 (原 Figure S7)
# ==============================================================================
message(">>> 正在处理 Module 视角数据 (Figure S7)...")
Priority_mod <- c("Operon", "RC_closed", "RC_open", "Regulon", "No_mapping")
Mapping_order_mod <- c("Perfect_Match", "Module_is_Subset", "Module_Contains_Term", "Partial Overlap", "No_mapping")

Module2Term <- Module2Function %>%
  select(Id, Module_size, Function_term, Type, Term_size, overlap_genes_count, Precision, Recall, F_score) %>%
  mutate(Type = case_when(Type == "RC_only" ~ "RC_closed", Type == "RC_muti" ~ "RC_open", TRUE ~ Type)) %>%
  mutate(Mapping_state = case_when(
    F_score >= Perfect_F_Cut ~ "Perfect_Match",
    Precision >= Perfect_F_Cut ~ "Module_is_Subset",
    Recall >= Perfect_F_Cut ~ "Module_Contains_Term",
    overlap_genes_count >= minOverlapSize ~ "Partial Overlap",
    TRUE ~ "No_mapping"
  )) %>%
  filter(Mapping_state != "No_mapping") %>%
  mutate(Type = factor(Type, levels = Priority_mod), Mapping_state = factor(Mapping_state, levels = Mapping_order_mod)) %>%
  group_by(Id) %>% arrange(Mapping_state, desc(F_score), Type) %>% slice(1) %>% ungroup() %>%
  select(Id, Type, Mapping_state, F_score) %>% right_join(Module_info, by = "Id") %>%
  mutate(Mapping_state = replace_na(as.character(Mapping_state), "No_mapping"),
         Type = replace_na(as.character(Type), "No_mapping")) %>%
  mutate(Type = if_else(Mapping_state == "No_mapping", "No_mapping", Type),
         Type = factor(Type, levels = Priority_mod),
         Mapping_state = factor(Mapping_state, levels = Mapping_order_mod))

inner_mod <- Module2Term %>% count(Type) %>% mutate(fraction = n/sum(n), percentage = round(fraction*100,1), ymax = cumsum(fraction), ymin = lag(ymax, default=0), label_pos = (ymax+ymin)/2)
outer_mod <- Module2Term %>% count(Type, Mapping_state) %>% mutate(Mapping_state = factor(Mapping_state, levels = Mapping_order_mod)) %>% group_by(Type) %>% mutate(fraction = n/sum(n), percentage = round(fraction*100,1), ymax_rel = cumsum(fraction), ymin_rel = lag(ymax_rel, default=0)) %>% left_join(inner_mod %>% select(Type, ymin_type = ymin, ymax_type = ymax), by = "Type") %>% mutate(ymin = ymin_type + (ymax_type - ymin_type)*ymin_rel, ymax = ymin_type + (ymax_type - ymin_type)*ymax_rel, label_pos = (ymax+ymin)/2) %>% ungroup()

Color_Type_mod <- c("Operon"="lightsalmon3", "RC_closed"="lightblue", "RC_open"="cornflowerblue", "Regulon"="dodgerblue4", "No_mapping"="grey50")
Color_Map_mod <- c("Perfect_Match"="firebrick4", "Module_is_Subset"="darkkhaki", "Module_Contains_Term"="cadetblue", "Partial Overlap"="cornsilk", "No_mapping"="grey50")

p_module <- ggplot() +
  geom_rect(data=outer_mod, aes(ymin=ymin, ymax=ymax, xmin=3.2, xmax=3.5, fill=Mapping_state), color="white", linewidth=0.3) + scale_fill_manual(values=Color_Map_mod, guide="none") +
  new_scale_fill() + geom_rect(data=inner_mod, aes(ymin=ymin, ymax=ymax, xmin=1, xmax=3, fill=Type), color="white", linewidth=0.4) + scale_fill_manual(values=Color_Type_mod, guide="none") +
  geom_text(data=inner_mod, aes(x=2.7, y=label_pos, label=paste0(Type, "\n(n=", n, ", ", percentage, "%)")), size=3.5, family="Arial", fontface="bold") +
  geom_text_repel(data=outer_mod %>% filter(n>=25), aes(x=3.3, y=label_pos, label=paste0(Mapping_state, "\n(n=", n, ")")), size=3.5, family="Arial", fontface="bold", min.segment.length=0, force=5, direction="y", segment.size=0.4, segment.color="grey50", nudge_x=0.8, max.overlaps=Inf) +
  coord_polar(theta="y") + xlim(1, 4.8) + theme_void() +
  labs(title = "A. Module Mapping State", subtitle = "What annotations do modules capture?") +
  theme(plot.title = element_text(hjust=0.5, face="bold", size=16), plot.subtitle = element_text(hjust=0.5, size=12, color="grey30", margin=margin(b=10)))

# ==============================================================================
# Part 2: Operon 视角图 
# ==============================================================================
message(">>> 正在处理 Operon 视角数据...")
Mapping_order_op <- c("Perfect_Match", "Operon_is_Subset", "Operon_Contains_Module", "Partial Overlap", "No_mapping")

Operon_all <- All_Annotations %>% filter(Type == "Operon") %>% select(Function_term) %>% distinct()
Operon_Final <- Module2Function %>% filter(Type == "Operon") %>%
  select(Id, Function_term, Term_size, overlap_genes_count, Precision, Recall, F_score) %>%
  mutate(Mapping_state = case_when(F_score >= Perfect_F_Cut ~ "Perfect_Match", Recall >= Perfect_F_Cut ~ "Operon_is_Subset", Precision >= Perfect_F_Cut ~ "Operon_Contains_Module", overlap_genes_count >= minOverlapSize ~ "Partial Overlap", TRUE ~ "No_mapping")) %>%
  filter(Mapping_state != "No_mapping") %>% mutate(Mapping_state = factor(Mapping_state, levels = Mapping_order_op)) %>%
  group_by(Function_term) %>% arrange(Mapping_state, desc(F_score)) %>% slice(1) %>% ungroup() %>%
  right_join(Operon_all, by = "Function_term") %>%
  mutate(Type = "Operon", Mapping_state = replace_na(as.character(Mapping_state), "No_mapping"), Mapping_state = factor(Mapping_state, levels = Mapping_order_op))

inner_op <- Operon_Final %>% count(Type) %>% mutate(fraction = n/sum(n), percentage = round(fraction*100,1), ymax = cumsum(fraction), ymin = lag(ymax, default=0), label_pos = (ymax+ymin)/2)
outer_op <- Operon_Final %>% count(Type, Mapping_state) %>% mutate(Mapping_state = factor(Mapping_state, levels = Mapping_order_op)) %>% group_by(Type) %>% mutate(fraction = n/sum(n), percentage = round(fraction*100,1), ymax_rel = cumsum(fraction), ymin_rel = lag(ymax_rel, default=0)) %>% left_join(inner_op %>% select(Type, ymin_type = ymin, ymax_type = ymax), by = "Type") %>% mutate(ymin = ymin_type + (ymax_type - ymin_type)*ymin_rel, ymax = ymin_type + (ymax_type - ymin_type)*ymax_rel, label_pos = (ymax+ymin)/2) %>% ungroup()

Color_Map_op <- c("Perfect_Match"="firebrick4", "Operon_is_Subset"="darkkhaki", "Operon_Contains_Module"="cadetblue", "Partial Overlap"="cornsilk", "No_mapping"="grey50")

p_operon <- ggplot() +
  geom_rect(data=outer_op, aes(ymin=ymin, ymax=ymax, xmin=3.2, xmax=3.5, fill=Mapping_state), color="white", linewidth=0.3) + scale_fill_manual(values=Color_Map_op, guide="none") +
  new_scale_fill() + geom_rect(data=inner_op, aes(ymin=ymin, ymax=ymax, xmin=1, xmax=3, fill=Type), color="white", linewidth=0.4) + scale_fill_manual(values=c("Operon"="lightsalmon3"), guide="none") +
  geom_text(data=inner_op, aes(x=2, y=label_pos, label=paste0("Total ", Type, "s\n(n=", n, ")")), size=5, family="Arial", fontface="bold") +
  geom_text_repel(data=outer_op %>% filter(n>0), aes(x=3.3, y=label_pos, label=paste0(Mapping_state, "\n(n=", n, ", ", percentage, "%)")), size=3.5, family="Arial", fontface="bold", min.segment.length=0, force=5, direction="y", segment.size=0.4, segment.color="grey50", nudge_x=0.8, max.overlaps=Inf) +
  coord_polar(theta="y") + xlim(1, 4.8) + theme_void() +
  labs(title = "B. Operon-centric Mapping State", subtitle = "How completely are Operons captured?") +
  theme(plot.title = element_text(hjust=0.5, face="bold", size=16), plot.subtitle = element_text(hjust=0.5, size=12, color="grey30", margin=margin(b=10)))

# ==============================================================================
# Part 3: Regulon 视角图 
# ==============================================================================
message(">>> 正在处理 Regulon 视角数据...")
Mapping_order_reg <- c("Perfect_Match", "Regulon_is_Subset", "Regulon_Contains_Module", "Partial Overlap", "No_mapping")

Regulon_all <- All_Annotations %>% filter(Type == "Regulon") %>% select(Function_term) %>% distinct()
Regulon_Final <- Module2Function %>% filter(Type == "Regulon") %>%
  select(Id, Function_term, Term_size, overlap_genes_count, Precision, Recall, F_score) %>%
  mutate(Mapping_state = case_when(F_score >= Perfect_F_Cut ~ "Perfect_Match", Recall >= Perfect_F_Cut ~ "Regulon_is_Subset", Precision >= Perfect_F_Cut ~ "Regulon_Contains_Module", overlap_genes_count >= minOverlapSize ~ "Partial Overlap", TRUE ~ "No_mapping")) %>%
  filter(Mapping_state != "No_mapping") %>% mutate(Mapping_state = factor(Mapping_state, levels = Mapping_order_reg)) %>%
  group_by(Function_term) %>% arrange(Mapping_state, desc(F_score)) %>% slice(1) %>% ungroup() %>%
  right_join(Regulon_all, by = "Function_term") %>%
  mutate(Type = "Regulon", Mapping_state = replace_na(as.character(Mapping_state), "No_mapping"), Mapping_state = factor(Mapping_state, levels = Mapping_order_reg))

inner_reg <- Regulon_Final %>% count(Type) %>% mutate(fraction = n/sum(n), percentage = round(fraction*100,1), ymax = cumsum(fraction), ymin = lag(ymax, default=0), label_pos = (ymax+ymin)/2)
outer_reg <- Regulon_Final %>% count(Type, Mapping_state) %>% mutate(Mapping_state = factor(Mapping_state, levels = Mapping_order_reg)) %>% group_by(Type) %>% mutate(fraction = n/sum(n), percentage = round(fraction*100,1), ymax_rel = cumsum(fraction), ymin_rel = lag(ymax_rel, default=0)) %>% left_join(inner_reg %>% select(Type, ymin_type = ymin, ymax_type = ymax), by = "Type") %>% mutate(ymin = ymin_type + (ymax_type - ymin_type)*ymin_rel, ymax = ymin_type + (ymax_type - ymin_type)*ymax_rel, label_pos = (ymax+ymin)/2) %>% ungroup()

Color_Map_reg <- c("Perfect_Match"="firebrick4", "Regulon_is_Subset"="darkkhaki", "Regulon_Contains_Module"="cadetblue", "Partial Overlap"="cornsilk", "No_mapping"="grey50")

p_regulon <- ggplot() +
  geom_rect(data=outer_reg, aes(ymin=ymin, ymax=ymax, xmin=3.2, xmax=3.5, fill=Mapping_state), color="white", linewidth=0.3) + scale_fill_manual(values=Color_Map_reg, guide="none") +
  new_scale_fill() + geom_rect(data=inner_reg, aes(ymin=ymin, ymax=ymax, xmin=1, xmax=3, fill=Type), color="white", linewidth=0.4) + scale_fill_manual(values=c("Regulon"="dodgerblue4"), guide="none") +
  geom_text(data=inner_reg, aes(x=2, y=label_pos, label=paste0("Total ", Type, "s\n(n=", n, ")")), size=5, family="Arial", fontface="bold") +
  geom_text_repel(data=outer_reg %>% filter(n>0), aes(x=3.3, y=label_pos, label=paste0(Mapping_state, "\n(n=", n, ", ", percentage, "%)")), size=3.5, family="Arial", fontface="bold", min.segment.length=0, force=5, direction="y", segment.size=0.4, segment.color="grey50", nudge_x=0.8, max.overlaps=Inf) +
  coord_polar(theta="y") + xlim(1, 4.8) + theme_void() +
  labs(title = "C. Regulon-centric Mapping State", subtitle = "How completely are Regulons captured?") +
  theme(plot.title = element_text(hjust=0.5, face="bold", size=16), plot.subtitle = element_text(hjust=0.5, size=12, color="grey30", margin=margin(b=10)))

# ==============================================================================
# Part 4: 终极拼图与导出 (借助 patchwork)
# ==============================================================================
message(">>> 正在使用 patchwork 组合三张图表...")

# 将三张图横向排列，并且在最上方加上总标题
combined_plot <- p_module + p_operon + p_regulon + 
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Comprehensive Mapping Evaluation of Core-Network Modules",
    theme = theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 20)))
  )

# 保存为超宽的高清 PDF
output_pdf <- file.path(base_dir, "paper图片/FigureS7_Mapping_Triple_Rings.pdf")
ggsave(output_pdf, plot = combined_plot, width = 24, height = 8, device = cairo_pdf)

message(">>> 完美！三环拼图已成功保存至: ", output_pdf)