##################################
# Representative gene expression heatmap for malignant metaprograms
# Author: Kevin Johnson
# Date Updated: 2026.04.10
##################################

library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(EnvStats)
library(scales)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
out_data_dir <- file.path(proj_dir, "processed_data/rna/")
script_dir  <- file.path(proj_dir, "scripts")

setwd(proj_dir)
heatCols <- readRDS(file.path(script_dir, "utils", "heatCols.RDS"))
color.scheme <- colorRampPalette(c(heatCols))(n=333)

# Read in the cell state classification based on p-value assignment.
care_state_md <- read.delim(file = paste0(proj_dir, "/results/scoring/caremut_malignant_cell_state_assignment.txt"), sep="\t", header = TRUE)

# Relabel some of the features.
care_state_md <- care_state_md %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "Oligo.",
                                    `IDHmut-noncodel` = "Astro."),
         patient_id = sapply(strsplit(care_id, "T"), "[[", 1),
         timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2)), 
         State = recode(State, `MP_AC1_MUT` = "AC-like",
                             `MP_OPC_MUT` = "OPC-like",
                             `MP_NPC_MUT` = "NPC-like",
                             `MP_MES_MUT` = "MES-like",
                             `MP_AC2_MUT` = "AC-like",
                             "Undifferentiated" = "Undifferentiated")) 


# Load in the CAREmut UMI produced at the beginning of the project.
umi_data_all <- readRDS(paste0(proj_dir, "/data/snrna/care_mut_umi_data_all_20230729.RDS"))
names(umi_data_all)

# malignant metaprograms
mut_mp <- read.table(paste0(proj_dir, "/results/metaprograms/care_mut_selected_malignant_metaprograms.csv"), sep = ",", header = T, row.names = 1)
mut_mp_long <- mut_mp %>% 
  pivot_longer(cols=c(MP_AC1:MP_CC), names_to = "mp", values_to = "gene") %>% 
  arrange(mp) %>% 
  filter(mp%in%c("MP_AC1", "MP_OPC", "MP_MES", "MP_NPC"))


# Define the MP column that maps to each state
state_mp_map <- c(
  "AC-like"         = "MP_AC1_MUT",
  "MES-like"        = "MP_MES_MUT",
  "OPC-like"        = "MP_OPC_MUT",
  "NPC-like"        = "MP_NPC_MUT"
)

# Extract representative cells with clear metaprogram expression profiles
set.seed(42)
md_malignant_scored <- care_state_md %>%
  filter(State %in% names(state_mp_map)) %>%
  group_by(idh_codel_subtype, State) %>%
  group_modify(~ {
    mp_col <- state_mp_map[.y$State]
    top20 <- .x %>%
      arrange(desc(.data[[mp_col]])) %>%
      slice_head(n = ceiling(0.2 * nrow(.x)))
    top20 %>% sample_n(min(200, nrow(top20)), replace = FALSE)
  }) %>%
  ungroup()


md_malignant_undiff <- care_state_md %>%
  filter(State == "Undifferentiated") %>%
  group_by(idh_codel_subtype) %>%
  group_modify(~ {
    bottom20 <- .x %>%
      mutate(mp_variance = apply(select(., MP_AC1_MUT, MP_MES_MUT, MP_OPC_MUT, MP_NPC_MUT), 1, var)) %>%
      arrange(mp_variance) %>%
      slice_head(n = ceiling(0.2 * nrow(.x)))
    bottom20 %>% sample_n(min(200, nrow(bottom20)), replace = FALSE)
  }) %>%
  ungroup()

md_malignant <- bind_rows(md_malignant_scored, md_malignant_undiff)

# Examine the distribution. There are 200 for each MP assigned and each subtype (n = 400 per MP)
table(md_malignant$idh_codel_subtype, md_malignant$State)
# Distributed across samples
table(md_malignant$SampleID)

# Define a function that subsets columns based on CellID
subset_columns <- function(mat, cell_ids) {
  col_idx <- which(colnames(mat) %in% cell_ids)
  mat[, col_idx, drop = FALSE]
}

# Apply the function to each element of the list
umi_data_all_malignant <- map(umi_data_all, subset_columns, cell_ids = md_malignant$CellID)

# Convert the list of dgCMatrix objects to a single matrix
cm <- do.call(cbind, umi_data_all_malignant)
cm <- as.matrix(cm)

umi2upm <- function(m) {
  count_sum <- colSums(m)
  upm_data <- (t(t(m)/count_sum)) * 1e+06
  upm_data
}

cm <- umi2upm(cm)
dim(cm)


# log-transform and center    
mp_genes <- unique(mut_mp_long$gene)
mp_gene_names <- mp_genes[mp_genes %in% rownames(cm)]
m <- cm[mp_gene_names, ]
m <- log2(m / 10 + 1)
m <- apply(m, 1, function(x) x - mean(x))

M_new2        <- m
M_new2        <- apply(M_new2, 2, rev)
M_meltII      <-  reshape2::melt(t(M_new2)) 
M_meltII$Var2 <- factor(M_meltII$Var2)

# Create a new variable MP_score
care_state_md <- care_state_md %>%
  mutate(MP_score = State)

M_meltII_order <- M_meltII %>% 
  inner_join(care_state_md, by=c("Var2"="CellID")) %>% 
  left_join(mut_mp_long, by=c("Var1"="gene")) %>% 
  dplyr::select(-mp) %>% 
  distinct() %>% 
  arrange(State, MP_score) %>%
  mutate(Var2 = factor(Var2, levels = unique(Var2))) %>% 
  dplyr::select(Var1:value, State, isCC, MP_score) %>% 
  left_join(mut_mp_long, by=c("Var1"="gene"))

cell_state_order <- c("NPC-like",
                      "OPC-like",
                      "Undifferentiated",
                      "AC-like",
                      "MES-like")

mp_order <- c("MP_NPC",
              "MP_OPC" ,
              "MP_AC1",
              "MP_MES")

M_meltII_order$Var1 <- factor(M_meltII_order$Var1, levels = unique(mut_mp_long$gene))
M_meltII_order$State <- factor(M_meltII_order$State, levels = cell_state_order)
M_meltII_order$mp <- factor(M_meltII_order$mp, levels = mp_order)

# Selected metaprogram genes to visualize
y_labels_to_display <- c("OPCML", "DLL3", "OLIG1", # OPC-like
                         "CD44","NAMPT", "GAP43", # MES-like
                         "ID4", "AQP4", "GFAP", # AC-like
                         "SOX11", "SOX4", "TOX3") # NPC-like

y_labels_to_display <- c("CD44","NAMPT", "GAP43")


hm_plot_panel <- ggplot(data = M_meltII_order %>% 
                          mutate(State = recode(State, `Undifferentiated` = "Undiff."),
                                 mp = recode(mp, `MP_AC1` = "MP_AC")) %>% 
                          filter(State%in%c("NPC-like",
                                            "OPC-like",
                                            "AC-like",
                                            "MES-like",
                                            "Undiff."), mp%in%c("MP_NPC",
                                                                "MP_OPC" ,
                                                                "MP_AC",
                                                                "MP_MES")), aes(x=Var2, y=Var1, fill=value, color=value)) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(-4,4), low=color.scheme[1:111],  mid =color.scheme[112:222], high = color.scheme[223:333], midpoint = 0, oob=squish, name=NULL) +                                
  scale_fill_gradient2(limits=c(-4,4), low=color.scheme[1:111],  mid =color.scheme[112:222], high = color.scheme[223:333], midpoint = 0, oob=squish, name=NULL)  +
  theme(panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.title = element_text(size=8), legend.text = element_text(size = 8), legend.text.align = 0.5, legend.justification = "bottom" ) + 
  theme(axis.text.y=element_text(size = 12)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_grid(mp~State, scales="free", space = "free") +
  theme(panel.spacing = unit(0, "lines")) +
  labs(x="Cells", y="", color="Log2\nexpression") +
  scale_y_discrete(breaks = y_labels_to_display)

ggsave(paste0(fig_dir, "state_metaprogram_heatmap.pdf"), hm_plot_panel, width = 5, height = 4, bg="transparent")

ggsave(paste0(fig_dir, "state_metaprogram_heatmap.png"), hm_plot_panel, width = 5, height = 4, bg="transparent")


# END