#################################
## Create a final metaprogram table based on initial metaprogram derivcation and iterative NMF for testing in the Undifferentiated compartment
## Author: Kevin Johnson
## Updated:2026.04.02
#################################

library(dplyr)

# Gene expression metaprograms that were defined across all malignant cells.
care_mut_mps <- read.table("results/metaprograms/malignant_downsampled/meta_programs_generated_for_all_mutants_n74_min_groupsize5_2026-04-01.csv", sep = ",", row.names = 1, header = T)
colnames(care_mut_mps) <- c("MP1_OPC", "MP2_AC", "MP3_RP", "MP4_CC", "MP5_MES", "MP6_AC2", "MP7_LQ", "MP8_MES2", "MP9_AC3", "MP10_Stress", "MP11_Mix", "MP12_Mix")

# Metaprograms defined for undifferentiated malignant cells
undiff_mps <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/metaprograms/undifferentiated_downsampled/meta_programs_generated_for_undifferentiated_n71_min_groupsize5_2026-04-02.csv", sep = ",", row.names = 1, header = T)
colnames(undiff_mps) <- paste0(colnames(undiff_mps) , "_Undiff.")
undiff_mps_novel <- undiff_mps %>% 
  dplyr::select(MP_NPC = MP_4_Undiff.)

# Select the states used for classification
final_metaprograms <- care_mut_mps %>% 
  bind_cols(undiff_mps_novel) %>% 
  dplyr::select(MP_AC1 = MP2_AC, MP_AC2 = MP6_AC2, MP_MES = MP5_MES, MP_OPC = MP1_OPC, MP_NPC, MP_CC = MP4_CC)


write.csv(final_metaprograms, file = "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/metaprograms/care_mut_selected_malignant_metaprograms.csv")


# Check overlap with CARE IDH-wildtype metaprograms:
care_wt_mps <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/misc/carewt_malignant_metaprograms.txt", sep = "\t", header = T)
colnames(care_wt_mps) <- gsub("MP_", "CARE_IDHwt_MP", colnames(care_wt_mps))

sum(final_metaprograms$MP_AC1%in%care_wt_mps$CARE_IDHwt_MP4_AC)
sum(final_metaprograms$MP_AC2%in%care_wt_mps$CARE_IDHwt_MP4_AC)
sum(final_metaprograms$MP_MES%in%care_wt_mps$CARE_IDHwt_MP6_MES)
sum(final_metaprograms$MP_NPC%in%care_wt_mps$CARE_IDHwt_MP7_NPC)
sum(final_metaprograms$MP_CC%in%care_wt_mps$CARE_IDHwt_MP3_CC)

### END ###