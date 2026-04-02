##############################
### Examine chromatin potential and the impact of T1 ATAC on T2 RNA
### Author: Kevin Johnson
### Updated: 2025.10.18
##############################

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# Chromatin TF motif activity scores for malignant cells that also have RNA
tf_df <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/archr_care_celltype_tf_motif_activity_zscore_20240505.txt", header = TRUE)

# ATAC data scored for the RNA-based metaprograms
atac_df <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/atac_malignant_module_scores.txt", header = TRUE)
colnames(atac_df) <- gsub("Module.MP_", "ATAC_module_", colnames(atac_df))
colnames(atac_df) <- gsub("_MUT", "", colnames(atac_df))
# Reformatting so that it's easier to structure the downstream analysis 
atac_df <- atac_df %>% 
  dplyr::select(CellID, ATAC_module_AC1:ATAC_module_NPC)

# Read in the cell state classification based on p-value assignment.
rna_df <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/classification/caremut_all_select_state_assignment_20240416.txt", sep="\t", header = TRUE)

# Relabel some of the features.
rna_df <- rna_df %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "IDH-O",
                                    `IDHmut-noncodel` = "IDH-A"),
         cell_state = recode(State, `MP_AC1_MUT` = "AC-like",
                             `MP_OPC_MUT` = "OPC-like",
                             `MP_NPC_MUT` = "NPC-like",
                             `MP_MES_MUT` = "MES-like",
                             `MP_AC2_MUT` = "AC-like",
                             "Undifferentiated" = "Undifferentiated")) %>% 
  dplyr::select(-MP_CC_MUT)
colnames(rna_df) <- gsub("MP_", "RNA_module_", colnames(rna_df))
colnames(rna_df) <- gsub("_MUT", "", colnames(rna_df))

# Combined RNA module and ATAC module scores for the metaprograms.
merged_df <- rna_df %>%
  inner_join(atac_df, by = "CellID") 

# Filtering to only those cases that longitudinal T1 and T2 data to keep it simple and a minimum of 25 cells at both time points
filtered_df <- merged_df %>%
  filter(timepoint %in% c("T1", "T2")) %>%
  group_by(patient_id, timepoint) %>%
  filter(n() >= 25) %>%
  ungroup()

# Keep only patients with both timepoints
eligible <- filtered_df %>%
  group_by(patient_id) %>%
  summarise(n_timepoints = n_distinct(timepoint)) %>%
  filter(n_timepoints == 2) %>%
  pull(patient_id)

filtered_df <- filtered_df %>%
  filter(patient_id %in% eligible)

# For each sample and cell state, take the average of the different module scores.
agg_df <- filtered_df %>%
  group_by(patient_id, timepoint, cell_state) %>%
  summarise(
    across(starts_with("RNA_module_"), mean, na.rm = TRUE),
    across(starts_with("ATAC_module_"), mean, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

agg_wide <- agg_df %>%
  pivot_wider(
    id_cols = c(patient_id, cell_state),
    names_from = timepoint,
    values_from = c(starts_with("RNA_module_"), starts_with("ATAC_module_")),
    names_sep = "_"
  )

programs <- c("AC1", "AC2", "MES", "OPC", "NPC")

# For each lineage (L) exam
results <- map_dfr(programs, function(L) {
  atac_col <- paste0("ATAC_module_", L, "_T1")
  rna_col <- paste0("RNA_module_", L, "_T2")
  rna_col_T1 <- paste0("RNA_module_", L, "_T1")
  
  df <- agg_wide %>%
    filter(!is.na(.data[[atac_col]]), !is.na(.data[[rna_col]]))
  
  if (nrow(df) > 3) {
    cor_res <- suppressWarnings(cor.test(df[[atac_col]], df[[rna_col]], method = "spearman"))
    lm_res <- lm(df[[rna_col]] ~ df[[atac_col]] + df[[rna_col_T1]])
    tibble(
      lineage = L,
      rho = cor_res$estimate,
      pval = cor_res$p.value,
      lm_pval = summary(lm_res)$coefficients[2, 4],
      n_pairs = nrow(df)
    )
  } else {
    tibble(lineage = L, rho = NA, pval = NA, lm_pval = NA, n_pairs = nrow(df))
  }
})

ggplot(agg_wide, aes(x = ATAC_module_AC1_T1, y = RNA_module_AC1_T2)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  facet_wrap(~cell_state) +
  labs(x = "T1 ATAC module (AC1)", y = "T2 RNA module (AC1)",
       title = "Chromatin accessibility at T1 vs RNA expression at T2") +
  theme_minimal(base_size = 12)

ggplot(agg_wide, aes(x = ATAC_module_AC2_T1, y = RNA_module_AC2_T2)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  facet_wrap(~cell_state) +
  labs(x = "T2 ATAC module (AC1)", y = "T2 RNA module (AC1)",
       title = "Chromatin accessibility at T2 vs RNA expression at T2") +
  theme_minimal(base_size = 12)

# There appears to be a weak effect for AC2 T1 ATAC and AC2 T2 RNA.
ggplot(results, aes(x = lineage, y = rho, fill = lineage)) +
  geom_col() +
  geom_text(aes(label = sprintf("n=%d\np=%.2g", n_pairs, pval)), vjust = -0.5) +
  labs(
    title = "T1 ATAC Accessibility vs T2 RNA Expression",
    y = "Spearman ρ (T1 ATAC → T2 RNA)",
    x = "Lineage module"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## Test whether the global profiles (cell state independent) can explain a link between T1 ATAC and T2 RNA
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
agg_global <- filtered_df %>%
  group_by(patient_id, timepoint) %>%
  summarise(
    across(starts_with("RNA_module_"), mean, na.rm = TRUE),
    across(starts_with("ATAC_module_"), mean, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

agg_global_wide <- agg_global %>%
  pivot_wider(
    id_cols = patient_id,
    names_from = timepoint,
    values_from = c(starts_with("RNA_module_"), starts_with("ATAC_module_")),
    names_sep = "_"
  )

results_global <- purrr::map_dfr(programs, function(L) {
  atac_col <- paste0("ATAC_module_", L, "_T1")
  rna_col  <- paste0("RNA_module_", L, "_T2")
  rna_col_T1 <- paste0("RNA_module_", L, "_T1")
  
  df <- agg_global_wide %>%
    filter(!is.na(.data[[atac_col]]), !is.na(.data[[rna_col]]))
  
  if (nrow(df) > 3) {
    cor_res <- suppressWarnings(cor.test(df[[atac_col]], df[[rna_col]], method = "spearman"))
    lm_res  <- lm(df[[rna_col]] ~ df[[atac_col]] + df[[rna_col_T1]])
    
    tibble(
      lineage = L,
      rho = cor_res$estimate,
      pval = cor_res$p.value,
      lm_pval = summary(lm_res)$coefficients[2, 4],
      n_pairs = nrow(df)
    )
  } else {
    tibble(lineage = L, rho = NA, pval = NA, lm_pval = NA, n_pairs = nrow(df))
  }
})

# It seems that with this analysis, AC2 ATAC profile is associated with RNA expression at T2.
results_global


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## Test whether cellular proportions at T2 can be predicted based on mean TF activity for ALL malignant cells
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# ---- compute composition per patient × state × timepoint ----
prop_df <- filtered_df %>%
  filter(timepoint %in% c("T1", "T2")) %>%
  count(patient_id, timepoint, cell_state) %>%
  group_by(patient_id, timepoint) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  complete(
    patient_id,
    timepoint,
    cell_state,
    fill = list(counts = 0, prop = 0)) %>% 
  select(patient_id, cell_state, timepoint, prop)

prop_df$prop <- prop_df$prop*100 

# reshape to wide format (T1 and T2 proportions)
prop_wide <- prop_df %>%
  pivot_wider(
    names_from = timepoint,
    values_from = prop,
    names_prefix = "Prop_"
  )

comp_df <- agg_global %>%
  filter(timepoint == "T1") %>%
  select(patient_id, starts_with("ATAC_module_"), starts_with("RNA_module_")) %>%
  left_join(prop_wide, by = c("patient_id"))

programs <- c("AC1", "AC2", "MES", "OPC", "NPC")

comp_df$cell_state <- gsub("-like", "", comp_df$cell_state)

comp_results <- map_dfr(programs, function(L) {
  atac_col <- paste0("ATAC_module_", L)
  rna_col  <- paste0("RNA_module_", L)
  
  # Skip if columns are missing
  if (!all(c(atac_col, rna_col, "Prop_T2") %in% names(comp_df))) {
    message("Skipping ", L, ": columns missing")
    return(tibble(lineage = L, coef = NA, pval = NA, n = NA))
  }
  
  # Filter and prepare data
  df <- comp_df %>%
    filter(
      !is.na(.data[[atac_col]]),
      !is.na(.data[[rna_col]]),
      !is.na(Prop_T2),
      cell_state == L,
    ) %>%
    mutate(
      ATAC = .data[[atac_col]],
      RNA  = .data[[rna_col]]
    )
  
  # Skip if too few samples
  if (nrow(df) < 5) {
    return(tibble(lineage = L, coef = NA, pval = NA, n = nrow(df)))
  }
  
  # Fit linear model
  fit <- lm(Prop_T2 ~ ATAC + RNA, data = df)
  s <- summary(fit)
  
  # Extract the ATAC coefficient and its p-value
  if (!"ATAC" %in% rownames(s$coefficients)) {
    return(tibble(lineage = L, coef = NA, pval = NA, n = nrow(df)))
  }
  
  atac_coef <- s$coefficients["ATAC", "Estimate"]
  atac_pval <- s$coefficients["ATAC", "Pr(>|t|)"]
  
  tibble(
    lineage = L,
    coef = atac_coef,
    pval = atac_pval,
    n = nrow(df)
  )
})

comp_results

# Test AC1/2 manually due to issues with AC1/AC2 and AC-like collapsing
comp_df_ac <- comp_df %>% 
  filter(cell_state=="AC")
fit_ac1 <- lm(Prop_T2 ~ ATAC_module_AC1 + RNA_module_AC1, data = comp_df_ac)
summary(fit_ac1)
# Again, AC2 seems to have some association, but not statistically significant
fit_ac2 <- lm(Prop_T2 ~ ATAC_module_AC2 + RNA_module_AC2, data = comp_df_ac)
summary(fit_ac2)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### TF activity
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

merged_tf <- rna_df %>%
inner_join(tf_df, by = "CellID")

filtered_tf <- merged_tf %>%
  filter(timepoint %in% c("T1", "T2")) %>%
  group_by(patient_id, timepoint) %>%
  filter(n() >= 25) %>%
  ungroup()

agg_tf <- filtered_tf %>%
  filter(patient_id%in%eligible) %>% 
  group_by(patient_id, timepoint) %>%
  summarise(
    across(starts_with("RNA_module_"), mean, na.rm = TRUE),
    across(matches("^[A-Z0-9_]+_[0-9]+$"), mean, na.rm = TRUE),  # captures motif cols like TFAP2B_1
    n_cells = n(),
    .groups = "drop"
  )

agg_tf_wide <- agg_tf %>%
  pivot_wider(
    id_cols = c(patient_id),
    names_from = timepoint,
    values_from = c(starts_with("RNA_module_"), matches("^[A-Z0-9_]+_[0-9]+$")),
    names_sep = "_"
  )
names(agg_tf_wide)[1:10]

motif_cols <- grep("_[0-9]+_T1$", names(agg_tf_wide), value = TRUE)

tf_results <- map_dfr(motif_cols, function(motif_col_T1) {
  motif_name <- sub("_T1$", "", motif_col_T1)
  
  map_dfr(programs, function(L) {
    rna_T2 <- paste0("RNA_module_", L, "_T2")
    rna_T1 <- paste0("RNA_module_", L, "_T1")
    
    # skip if missing columns
    if (!all(c(rna_T1, rna_T2, motif_col_T1) %in% names(agg_tf_wide))) {
      return(NULL)
    }
    
    df <- agg_tf_wide %>%
      filter(!is.na(.data[[motif_col_T1]]),
             !is.na(.data[[rna_T1]]),
             !is.na(.data[[rna_T2]])) %>%
      mutate(
        TF_T1  = .data[[motif_col_T1]],
        RNA_T1 = .data[[rna_T1]],
        RNA_T2 = .data[[rna_T2]]
      )
    
    if (nrow(df) < 5) return(NULL)
    
    # correlation and linear model
    cor_res <- suppressWarnings(cor.test(df$TF_T1, df$RNA_T2, method = "spearman"))
    fit     <- lm(RNA_T2 ~ TF_T1 + RNA_T1, data = df)
    s       <- summary(fit)
    
    tibble(
      motif     = motif_name,
      lineage   = L,
      rho       = unname(cor_res$estimate),
      pval_rho  = cor_res$p.value,
      coef_TF   = s$coefficients["TF_T1", "Estimate"],
      pval_TF   = s$coefficients["TF_T1", "Pr(>|t|)"],
      n         = nrow(df),
      r2        = s$r.squared
    )
  })
})

tf_results <- tf_results %>%
  group_by(lineage) %>%
  mutate(padj_TF = p.adjust(pval_TF, method = "fdr")) %>%
  ungroup()

tf_results_clean <- tf_results %>%
  # Rename AC1 → AC and remove AC2
  mutate(lineage = ifelse(lineage == "AC1", "AC", lineage)) %>%
  filter(lineage != "AC2") %>%
  # Rank by adjusted p-value within each lineage
  group_by(lineage) %>%
  mutate(rank = rank(padj_TF, ties.method = "first")) %>%
  ungroup()

# Identify top 3 per lineage
top_tfs <- tf_results_clean %>%
  group_by(lineage) %>%
  slice_min(padj_TF, n = 3, with_ties = FALSE)

library(ggrepel)
source("/vast/palmer/pi/verhaak/kcj28/care_mut/scripts/misc/plot_theme.R")

# Volcano plot
pdf(paste0("/vast/palmer/pi/verhaak/kcj28/care_mut/results/figures/archr/s1_tf_activity_predicts_s1_expression.pdf"), width=4.5, height=4.5, useDingbats = FALSE, bg = "transparent")
ggplot(tf_results_clean, aes(x = coef_TF, y = -log10(padj_TF))) +
  geom_point(aes(alpha = 0.7)) +
  #geom_text_repel(
  #  data = top_tfs,
  #  aes(label = motif),
  #  size = 3,
  #  color = "black",
  #  box.padding = 0.5,
  #  max.overlaps = Inf
  #) +
  guides(alpha=FALSE) +
  facet_wrap(~lineage, scales = "free_y") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme_bw(base_size = 8) +
  labs(
    #title = "TF activity (Sample 1) assoc. RNA module (Sample 2)",
    title = "S1 TF activity vs. S2 RNA expression\nModel: lm(RNA_S2 ~ TF_S2 + RNA_S1)",
    x = "Effect size (beta)",
    y = "-log10(adjusted p-value)"
  )
dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Examine TF activity at T1 on T2 proportions while controlling for T1 proportions
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

tf_prop_df <- merged_tf %>%
  filter(patient_id%in%eligible) %>% 
  filter(timepoint %in% c("T1", "T2")) %>%
  count(patient_id, timepoint, cell_state) %>%
  group_by(patient_id, timepoint) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  complete(
    patient_id,
    timepoint,
    cell_state,
    fill = list(counts = 0, prop = 0)) %>% 
  select(patient_id, cell_state, timepoint, prop)


prop_wide <- tf_prop_df %>%
  pivot_wider(names_from = timepoint, values_from = prop, names_prefix = "Prop_")

agg_tf_T1 <- agg_tf %>%
  filter(timepoint == "T1") %>%
  select(patient_id, starts_with("RNA_module_"), matches("^[A-Z0-9_]+_[0-9]+$"))

comp_tf <- agg_tf_T1 %>%
  left_join(prop_wide, by = c("patient_id"))


motif_cols <- grep("^[A-Z0-9_]+_[0-9]+$", names(comp_tf), value = TRUE)

cell_states <- unique(comp_tf$cell_state)


tf_prop_results <- map_dfr(motif_cols, function(motif_col_T1) {
  motif_name <- sub("_T1$", "", motif_col_T1)
  
  map_dfr(cell_states, function(cs) {
    # subset to this cell state
    df <- comp_tf %>%
      filter(cell_state == cs) %>%
      filter(!is.na(.data[[motif_col_T1]]),
             !is.na(Prop_T1),
             !is.na(Prop_T2)) %>%
      mutate(
        TF_T1 = .data[[motif_col_T1]]
      )
    
    # skip underpowered subsets
    if (nrow(df) < 5) return(NULL)
    
    # correlation
    cor_res <- suppressWarnings(cor.test(df$TF_T1, df$Prop_T2, method = "spearman"))
    
    # linear model: Prop_T2 ~ TF_T1 + Prop_T1
    fit <- lm(Prop_T2 ~ TF_T1 + Prop_T1, data = df)
    s <- summary(fit)
    
    tibble(
      motif = motif_name,
      cell_state = cs,
      rho = unname(cor_res$estimate),
      pval_rho = cor_res$p.value,
      coef_TF = s$coefficients["TF_T1", "Estimate"],
      pval_TF = s$coefficients["TF_T1", "Pr(>|t|)"],
      coef_Prop_T1 = s$coefficients["Prop_T1", "Estimate"],
      pval_Prop_T1 = s$coefficients["Prop_T1", "Pr(>|t|)"],
      n = nrow(df),
      r2 = s$r.squared
    )
  })
})

tf_prop_results %>% 
  filter(pval_rho < 0.05, pval_TF < 0.05)

tf_prop_results <- tf_prop_results %>%
  group_by(cell_state) %>%
  mutate(padj_TF = p.adjust(pval_TF, method = "BH")) %>%
  ungroup()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Examine TF activity at T1 on T2 RNA module while controlling for T1 RNA module
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
merged_tf <- rna_df %>%
  inner_join(tf_df, by = "CellID")

# Keep only T1/T2 malignant cells and patients with >=25 cells per time point
filtered_tf <- merged_tf %>%
  filter(patient_id%in%eligible) %>% 
  filter(timepoint %in% c("T1", "T2")) %>%
  group_by(patient_id, timepoint) %>%
  filter(n() >= 25) %>%
  ungroup()

# Aggregate *globally per patient per time point* (drop cell_state)
agg_tf <- filtered_tf %>%
  group_by(patient_id, timepoint) %>%
  summarise(
    across(starts_with("RNA_module_"), mean, na.rm = TRUE),
    across(matches("^[A-Z0-9_]+_[0-9]+$"), mean, na.rm = TRUE),  # motif columns
    n_cells = n(),
    .groups = "drop"
  )

# Pivot to wide format (patient-level)
agg_tf_wide <- agg_tf %>%
  pivot_wider(
    id_cols = patient_id,
    names_from = timepoint,
    values_from = c(starts_with("RNA_module_"), matches("^[A-Z0-9_]+_[0-9]+$")),
    names_sep = "_"
  )
names(agg_tf_wide)[1:10]


# Identify relevant columns
motif_cols <- grep("_[0-9]+_T1$", names(agg_tf_wide), value = TRUE)
rna_T2_cols <- grep("^RNA_module_.*_T2$", names(agg_tf_wide), value = TRUE)

tf_results_all <- map_dfr(rna_T2_cols, function(rna_T2) {
  rna_T1 <- sub("_T2$", "_T1", rna_T2)
  lineage_name <- sub("^RNA_module_|_T2$", "", rna_T2)
  
  # Skip if missing T1 column
  if (!rna_T1 %in% names(agg_tf_wide)) return(NULL)
  
  map_dfr(motif_cols, function(motif_col_T1) {
    motif_name <- sub("_T1$", "", motif_col_T1)
    
    df <- agg_tf_wide %>%
      filter(!is.na(.data[[motif_col_T1]]),
             !is.na(.data[[rna_T2]]),
             !is.na(.data[[rna_T1]]))
    
    if (nrow(df) < 5) return(NULL)
    
    cor_res <- suppressWarnings(cor.test(df[[motif_col_T1]], df[[rna_T2]], method = "spearman"))
    lm_res  <- lm(df[[rna_T2]] ~ df[[motif_col_T1]] + df[[rna_T1]])
    
    tibble(
      lineage = lineage_name,
      motif   = motif_name,
      rho     = cor_res$estimate,
      pval_cor = cor_res$p.value,
      beta_motif = summary(lm_res)$coefficients[2, 1],
      pval_motif = summary(lm_res)$coefficients[2, 4],
      beta_T1    = summary(lm_res)$coefficients[3, 1],
      pval_T1    = summary(lm_res)$coefficients[3, 4],
      n = nrow(df)
    )
  })
})

top_motifs <- tf_results_all %>% 
  filter(pval_cor < 0.05, pval_motif < 0.05)


# Get a better idea of what's expected by looking at T1 ATAC vs T1 RNA
motif_cols_T1 <- grep("^[A-Z0-9_]+_[0-9]+_T1$", names(agg_tf_wide), value = TRUE)
rna_cols_T1   <- grep("^RNA_module_.*_T1$",      names(agg_tf_wide), value = TRUE)

# All motif × RNA pairs at T1 (Spearman; add linear if you like)
t1_pairwise <- purrr::map_dfr(rna_cols_T1, function(rna_col) {
  lineage <- sub("^RNA_module_|_T1$", "", rna_col)
  purrr::map_dfr(motif_cols_T1, function(motif_col) {
    df <- agg_tf_wide %>%
      filter(!is.na(.data[[motif_col]]), !is.na(.data[[rna_col]]))
    if (nrow(df) < 5) return(NULL)
    
    ct <- suppressWarnings(cor.test(df[[motif_col]], df[[rna_col]], method = "spearman"))
    
    tibble(
      lineage    = lineage,
      motif_full = motif_col,
      motif      = sub("_T1$", "", sub("_[0-9]+_T1$", "", motif_col)),  # clean name
      rho        = unname(ct$estimate),
      pval       = ct$p.value,
      n          = nrow(df)
    )
  })
}) %>%
  group_by(lineage) %>%
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  ungroup()

top_t1 <- t1_pairwise %>%
  filter(!is.na(rho), rho>0) %>%
  group_by(lineage) %>%
  slice_max(order_by = abs(rho), n = 15, with_ties = FALSE) %>%
  ungroup()

### delta RNA2-RNA1 ###

tf_results_delta <- map_dfr(rna_T2_cols, function(rna_T2) {
  # define paired T1 column and lineage name
  rna_T1 <- sub("_T2$", "_T1", rna_T2)
  lineage_name <- sub("^RNA_module_|_T2$", "", rna_T2)
  
  # skip if missing T1 column
  if (!rna_T1 %in% names(agg_tf_wide)) return(NULL)
  
  map_dfr(motif_cols, function(motif_col_T1) {
    motif_name <- sub("_T1$", "", motif_col_T1)
    
    df <- agg_tf_wide %>%
      filter(!is.na(.data[[motif_col_T1]]),
             !is.na(.data[[rna_T2]]),
             !is.na(.data[[rna_T1]])) %>%
      mutate(
        dRNA = .data[[rna_T2]] - .data[[rna_T1]],   # the RNA change
        TF_T1 = .data[[motif_col_T1]]
      )
    
    if (nrow(df) < 5) return(NULL)
    
    # correlation: T1 TF vs ΔRNA
    cor_res <- suppressWarnings(cor.test(df$TF_T1, df$dRNA, method = "spearman"))
    
    # linear model: ΔRNA ~ TF_T1
    lm_res <- lm(dRNA ~ TF_T1, data = df)
    s <- summary(lm_res)
    
    tibble(
      lineage = lineage_name,
      motif   = motif_name,
      rho     = unname(cor_res$estimate),
      pval_cor = cor_res$p.value,
      beta_motif = s$coefficients["TF_T1", "Estimate"],
      pval_motif = s$coefficients["TF_T1", "Pr(>|t|)"],
      n = nrow(df),
      r2 = s$r.squared
    )
  })
})


top_motifs_delta <- tf_results_delta %>% 
  filter(pval_cor < 0.05, pval_motif < 0.05)

top_motifs_delta <- tf_results_delta %>%
  group_by(lineage) %>%
  mutate(padj_motif = p.adjust(pval_motif, method = "BH")) %>%
  ungroup() %>% 
  filter(padj_motif<0.05)

### END ###