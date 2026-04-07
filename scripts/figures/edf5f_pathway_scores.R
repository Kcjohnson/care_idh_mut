##################################
# Create a visualization for pathway-based metaprogram scores across malignant cell states.
# Author: Kevin Johnson
# Date Updated: 2026.04.06
##################################

library(tidyverse)

caremut_malignant <- caremut_md %>% 
  filter(CellType=="Malignant")

cell_ids <- rownames(scores)

caremut_malignant_new <- caremut_malignant %>% 
  anti_join(scores, by = )

load("/vast/palmer/pi/verhaak/kcj28/care_mut/data/pmp/CARE_IDHmut_10x_Tumor_scores_PMP_240903_10_10_10.RData")
pmp_scores <- scores
colnames(pmp_scores)[2:12] <- c("PMP1_Mitochrondria", "PMP2_Chromatin", "PMP3_AC", "PMP4_NEU", "PMP5_Morphogenesis", "PMP6_CellCycle", "PMP7_CalciumSignaling", "PMP8_AC_Cilia", "PMP9_Immune", "PMP10_GlycolysisStress", "PMP11_Angiogenesis")




care_state_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/scoring/caremut_malignant_cell_state_assignment.txt", sep="\t", header = TRUE)

care_state_md_pmp <- care_state_md %>% 
  inner_join(pmp_scores, by="CellID") %>% 
  mutate(cell_state = recode(State, `MP_AC1_MUT` = "AC-like",
                      `MP_OPC_MUT` = "OPC-like",
                      `MP_NPC_MUT` = "NPC-like",
                      `MP_MES_MUT` = "MES-like",
                      `MP_AC2_MUT` = "AC-like",
                      "Undifferentiated" = "Undifferentiated")) 

mp_public_scores_state_avg <- care_state_md_pmp %>% 
  dplyr::select(CellID, care_id, cell_state, PMP1_Mitochrondria:PMP11_Angiogenesis) %>% 
  pivot_longer(cols = c(PMP1_Mitochrondria:PMP11_Angiogenesis), names_to = "signatures", values_to = "scores") %>% 
  group_by(care_id, signatures, cell_state) %>% 
  summarise(avg_score = median(scores),
            cell_counts = n()) %>% 
  mutate(patient_id = sapply(strsplit(care_id, "T"), "[[", 1),
         timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2))) %>% 
  filter(timepoint!="T3") %>% 
  ungroup() %>% 
  arrange(care_id, signatures, cell_state) %>% 
  filter(cell_counts > 25)


pdf(file = "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/figures/pathway_score_by_state.pdf", height = 5, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(mp_public_scores_state_avg %>% 
         filter(signatures%in%c("PMP10_GlycolysisStress", "PMP1_Mitochrondria")) %>% 
                  mutate(signatures = recode(signatures, `PMP10_GlycolysisStress` = "PMP Glycolysis/Stress",
                                            `PMP1_Mitochrondria` = "PMP Mitochondria energy production")), aes(x=cell_state, y=avg_score, fill=cell_state)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  facet_wrap(.~signatures, scales="free_y") +
  plot_theme +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        strip.text = element_text(size=8)) +
  #stat_n_text() +
  stat_compare_means(method="kruskal", label="p.format") +
  labs(y="Median pathway metaprogram score\nper tumor (min. 25 cells)", x="Malignant state") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  guides(fill=FALSE)
dev.off()
