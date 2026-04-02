##############################
### Examine the overlap in upregulated and downregulated genes amongst malignant cells
### Author: Kevin Johnson
### Updated: 2025.10.18
##############################

## Investigate the overlap in increased and decreased gene accessibility in a cell state-independent manner

library(tidyverse)
library(readr)
library(purrr)

path <- "/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/state_controlled_dag/malignant"

files <- list.files(path, pattern = "longitudinal_dag\\.txt$", full.names = TRUE)

read_and_annotate <- function(f) {
  fname <- basename(f)
  parts <- strsplit(fname, "_x_")[[1]]
  first_part <- strsplit(parts[1], "_")[[1]]
  
  patient_id <- first_part[1]
  cell_state <- first_part[2]
  sample_tested <- parts[2]
  sample_background <- str_remove(parts[3], "_longitudinal_dag\\.txt")
  
  df <- read_tsv(
    f,
    col_types = cols(
      seqnames = col_character(),
      start = col_double(),
      end = col_double(),
      strand = col_character(),
      name = col_character(),
      idx = col_double(),
      Log2FC = col_double(),
      FDR = col_double(),
      MeanDiff = col_double()
    )
  )
  
  df %>%
    mutate(
      patient_id = patient_id,
      cell_state = cell_state,
      sample_tested = sample_tested,
      sample_background = sample_background,
      accessibility_change = case_when(
        Log2FC > 0 ~ "Increased",
        Log2FC < 0 ~ "Decreased",
        TRUE ~ "No change"
      )
    )
}

# Safely wrap
safe_read <- safely(read_and_annotate)

# Apply to all files
results_list <- map(files, safe_read)

# Separate successes and failures
successes <- results_list |> keep(~ is.null(.x$error)) |> map("result")
failures  <- results_list |> keep(~ !is.null(.x$error))

if (length(failures) > 0) {
  cat("The following files failed to load:\n")
  print(map_chr(failures, ~ .x$error))
}

# Combine successfully read files
all_results <- bind_rows(successes)

# Example summary
summary_table <- all_results %>%
  group_by(patient_id, cell_state, accessibility_change) %>%
  summarise(n = n(), .groups = "drop")

print(summary_table)



##
# Optionally filter for significant changes only
filtered_results <- all_results %>%
  filter(FDR < 0.05)  # adjust threshold if desired

# Summarize within each sample
sample_summary <- filtered_results %>%
  group_by(patient_id, cell_state, sample_tested) %>%
  summarise(
    n_total = n(),
    n_increased = sum(accessibility_change == "Increased"),
    n_decreased = sum(accessibility_change == "Decreased"),
    pct_increased = 100 * n_increased / n_total,
    pct_decreased = 100 * n_decreased / n_total,
    .groups = "drop"
  )

# Mean percent across the cohort
cohort_summary <- sample_summary %>%
  summarise(
    mean_pct_increased = mean(pct_increased, na.rm = TRUE),
    mean_pct_decreased = mean(pct_decreased, na.rm = TRUE)
  )

print(cohort_summary)


###

# Filter for significant increases
increased_genes <- all_results %>%
  filter(
    FDR < 0.05,
    Log2FC > 0,
    !is.na(name)
  )

# Count how many patients showed increase per gene per cell_state
gene_patient_counts <- increased_genes %>%
  distinct(patient_id, cell_state, name) %>%
  group_by(cell_state, name) %>%
  summarise(
    n_patients_increased = n_distinct(patient_id),
    .groups = "drop"
  )

#  Count total patients evaluated per cell_state
patients_per_state <- all_results %>%
  distinct(patient_id, cell_state) %>%
  group_by(cell_state) %>%
  summarise(
    n_patients_total = n_distinct(patient_id),
    .groups = "drop"
  )

# Merge the two and compute proportion
recurrent_genes <- gene_patient_counts %>%
  left_join(patients_per_state, by = "cell_state") %>%
  mutate(
    prop_patients = n_patients_increased / n_patients_total
  ) %>%
  filter(n_patients_increased >= 3) %>%   # threshold for recurrence
  arrange(desc(prop_patients), desc(n_patients_increased))

print(recurrent_genes)

recurrent_genes_top_2 <- recurrent_genes %>% 
  filter(prop_patients > 0.2)
table(recurrent_genes_top_25$name)

recurrent_genes %>% 
  group_by(name) %>% 
  dplyr::summarise(counts = n()) %>% 
  arrange(desc(counts))

## Decreased accessibility

# Filter for significant decreases
decreased_genes <- all_results %>%
  filter(
    FDR < 0.05,
    Log2FC < 0,
    !is.na(name)
  )

# Count how many patients show decrease per gene per cell_state
gene_patient_counts_dec <- decreased_genes %>%
  distinct(patient_id, cell_state, name) %>%
  group_by(cell_state, name) %>%
  summarise(
    n_patients_decreased = n_distinct(patient_id),
    .groups = "drop"
  )

# Count total patients evaluated per cell_state
patients_per_state <- all_results %>%
  distinct(patient_id, cell_state) %>%
  group_by(cell_state) %>%
  summarise(
    n_patients_total = n_distinct(patient_id),
    .groups = "drop"
  )

# Merge and compute proportion decreased
recurrent_genes_dec <- gene_patient_counts_dec %>%
  left_join(patients_per_state, by = "cell_state") %>%
  mutate(
    prop_patients = n_patients_decreased / n_patients_total
  ) %>%
  filter(n_patients_decreased >= 3) %>%   # threshold for recurrence
  arrange(desc(prop_patients), desc(n_patients_decreased))

print(recurrent_genes_dec)

tmp <- recurrent_genes_dec %>% 
  group_by(name) %>% 
  dplyr::summarise(counts = n()) %>% 
  arrange(desc(counts))

# Combine both if needed first (optional)
recurrent_genes_all <- bind_rows(
  recurrent_genes %>% mutate(direction = "Increased"),
  recurrent_genes_dec %>% mutate(direction = "Decreased")
)

# Count in how many cell states each gene appears per direction
shared_genes <- recurrent_genes_all %>%
  group_by(direction, name) %>%
  summarise(
    n_cell_states = n_distinct(cell_state),
    cell_states = paste(sort(unique(cell_state)), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_cell_states > 1) %>%
  arrange(desc(n_cell_states))

print(shared_genes)


###


# Filter to significant genes only
sig_results <- all_results %>%
  filter(FDR < 0.05, !is.na(name))
table(sig_results$patient_id, sig_results$cell_state)

# Count cell states per patient
patients_with_multi_states <- sig_results %>%
  distinct(patient_id, cell_state) %>%
  count(patient_id) %>%
  filter(n > 2)

# Restrict to those patients
multi_state_results <- sig_results %>%
  semi_join(patients_with_multi_states, by = "patient_id")

# For each patient × gene, determine direction consistency
gene_direction_summary <- multi_state_results %>%
  group_by(patient_id, name) %>%
  summarise(
    n_cell_states = n_distinct(cell_state),
    directions = paste0(sign(Log2FC), collapse = ","),
    all_positive = all(Log2FC > 0),
    all_negative = all(Log2FC < 0),
    consistent = all_positive | all_negative,
    .groups = "drop"
  ) %>%
  filter(n_cell_states >= 2)   # only genes tested in ≥2 states

# Compute per-patient proportion of consistent genes
patient_consistency <- gene_direction_summary %>%
  group_by(patient_id) %>%
  summarise(
    n_genes_multi_state = n(),
    n_consistent = sum(consistent),
    prop_consistent = n_consistent / n_genes_multi_state,
    .groups = "drop"
  ) %>%
  arrange(desc(prop_consistent))

print(patient_consistency)


###
# 1️⃣ Filter for significant genes
sig_results <- all_results %>%
  filter(FDR < 0.05, !is.na(name))

table(sig_results$patient_id, sig_results$cell_state)

# 2️⃣ Function to compute pairwise overlap metrics within one patient
pairwise_overlap <- function(df_patient) {
  states <- unique(df_patient$cell_state)
  if (length(states) < 2) return(NULL)
  
  combn(states, 2, simplify = FALSE, FUN = function(pair) {
    s1 <- pair[1]
    s2 <- pair[2]
    
    df1 <- df_patient %>% filter(cell_state == s1) %>% select(name, Log2FC)
    df2 <- df_patient %>% filter(cell_state == s2) %>% select(name, Log2FC)
    
    # shared genes
    shared <- inner_join(df1, df2, by = "name", suffix = c("_s1", "_s2"))
    n_shared <- nrow(shared)
    n_s1 <- nrow(df1)
    n_s2 <- nrow(df2)
    
    # direction consistency
    n_same_dir <- sum(sign(shared$Log2FC_s1) == sign(shared$Log2FC_s2), na.rm = TRUE)
    prop_same_dir <- ifelse(n_shared > 0, n_same_dir / n_shared, NA)
    
    tibble(
      patient_id = df_patient$patient_id[1],
      state1 = s1,
      state2 = s2,
      n_s1 = n_s1,
      n_s2 = n_s2,
      n_shared = n_shared,
      prop_shared_s1 = n_shared / n_s1,
      prop_shared_s2 = n_shared / n_s2,
      prop_shared_avg = n_shared / mean(c(n_s1, n_s2)),
      n_same_dir = n_same_dir,
      prop_same_dir = prop_same_dir
    )
  }) %>% bind_rows()
}

# 3️⃣ Apply to each patient
patient_overlap <- sig_results %>%
  group_by(patient_id) %>%
  group_split() %>%
  map_dfr(pairwise_overlap)

# 4️⃣ Optionally summarize across pairs per patient
patient_summary <- patient_overlap %>%
  group_by(patient_id) %>%
  summarise(
    mean_prop_shared = mean(prop_shared_avg, na.rm = TRUE),
    mean_prop_same_dir = mean(prop_same_dir, na.rm = TRUE),
    .groups = "drop"
  )

print(patient_summary)



####

# 1️⃣ Filter to significant differential accessibility genes
sig_results <- all_results %>%
  filter(FDR < 0.05, !is.na(name))

# 2️⃣ Function to compute direction-specific overlaps for one patient
pairwise_overlap_directional <- function(df_patient) {
  states <- unique(df_patient$cell_state)
  if (length(states) < 2) return(NULL)
  
  combn(states, 2, simplify = FALSE, FUN = function(pair) {
    s1 <- pair[1]
    s2 <- pair[2]
    
    df1 <- df_patient %>% filter(cell_state == s1) %>% select(name, Log2FC)
    df2 <- df_patient %>% filter(cell_state == s2) %>% select(name, Log2FC)
    
    # shared genes
    shared <- inner_join(df1, df2, by = "name", suffix = c("_s1", "_s2"))
    n_shared <- nrow(shared)
    n_s1 <- nrow(df1)
    n_s2 <- nrow(df2)
    
    # Direction breakdown
    n_joint_increased <- sum(shared$Log2FC_s1 > 0 & shared$Log2FC_s2 > 0)
    n_joint_decreased <- sum(shared$Log2FC_s1 < 0 & shared$Log2FC_s2 < 0)
    n_discordant <- n_shared - n_joint_increased - n_joint_decreased
    
    # Proportions
    prop_joint_increased <- ifelse(n_shared > 0, n_joint_increased / n_shared, NA)
    prop_joint_decreased <- ifelse(n_shared > 0, n_joint_decreased / n_shared, NA)
    prop_discordant <- ifelse(n_shared > 0, n_discordant / n_shared, NA)
    
    tibble(
      patient_id = df_patient$patient_id[1],
      state1 = s1,
      state2 = s2,
      n_s1 = n_s1,
      n_s2 = n_s2,
      n_shared = n_shared,
      jaccard = n_shared / (n_s1 + n_s2 - n_shared),
      prop_shared_avg = n_shared / mean(c(n_s1, n_s2)),
      n_joint_increased = n_joint_increased,
      n_joint_decreased = n_joint_decreased,
      n_discordant = n_discordant,
      prop_joint_increased = prop_joint_increased,
      prop_joint_decreased = prop_joint_decreased,
      prop_discordant = prop_discordant
    )
  }) %>% bind_rows()
}

# 3️⃣ Apply across all patients
patient_overlap_dir <- sig_results %>%
  group_by(patient_id) %>%
  group_split() %>%
  map_dfr(pairwise_overlap_directional)

# 4️⃣ Optional: summarize per patient (average across pairs)
patient_summary_dir <- patient_overlap_dir %>%
  group_by(patient_id) %>%
  summarise(
    mean_prop_shared = mean(prop_shared_avg, na.rm = TRUE),
    mean_prop_joint_increased = mean(prop_joint_increased, na.rm = TRUE),
    mean_prop_joint_decreased = mean(prop_joint_decreased, na.rm = TRUE),
    mean_prop_discordant = mean(prop_discordant, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_prop_shared))

print(patient_summary_dir)


###

library(ggplot2)
library(tidyr)

ggplot(patient_overlap_dir, aes(x = state1, y = state2, size = prop_shared_avg, fill = prop_joint_increased - prop_joint_decreased)) +
  geom_point(shape = 21, color = "black") +
  facet_wrap(~patient_id) +
  scale_fill_gradient2(low = "#377eb8", mid = "white", high = "#e41a1c", midpoint = 0,
                       name = "Directional bias\n(+ = increased)") +
  scale_size_continuous(name = "Proportion shared") +
  theme_minimal(base_size = 12) +
  labs(title = "Pairwise cross-state similarity of accessibility shifts")



# Suppose you have a filtered set of genes shared across multiple states in one patient
subset_data <- all_results %>%
  filter(patient_id == "SJ07", FDR < 0.05) %>%
  select(name, cell_state, Log2FC)

ggplot(subset_data, aes(x = cell_state, y = reorder(name, Log2FC), fill = Log2FC)) +
  geom_tile() +
  scale_fill_gradient2(low = "#377eb8", mid = "white", high = "#e41a1c", midpoint = 0) +
  theme_minimal(base_size = 10) +
  labs(title = "Example patient: directionality of shared accessibility across cell states",
       x = "Cell state", y = "Gene (sorted by mean Log2FC)")


ggplot(patient_summary_dir, aes(x = mean_prop_shared, y = mean_prop_joint_increased - mean_prop_joint_decreased, color = patient_id)) +
  geom_point(size = 4) +
  geom_vline(xintercept = mean(patient_summary_dir$mean_prop_shared), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "Mean proportion of shared accessible genes",
    y = "Directional bias (Increased − Decreased)",
    title = "Global vs. state-specific accessibility trends"
  ) +
  theme_minimal(base_size = 13)


### Jaccard ####
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

# Define Jaccard function
jaccard_index <- function(a, b) {
  inter <- length(intersect(a, b))
  union <- length(union(a, b))
  if (union == 0) return(NA_real_)
  inter / union
}

# For each patient, compute pairwise Jaccard for same-direction genes
compute_jaccard_patient <- function(df_patient) {
  states <- unique(df_patient$cell_state)
  if (length(states) < 2) return(NULL)
  
  combn(states, 2, simplify = FALSE, FUN = function(pair) {
    s1 <- pair[1]; s2 <- pair[2]
    
    df1 <- df_patient %>% filter(cell_state == s1, FDR < 0.05)
    df2 <- df_patient %>% filter(cell_state == s2, FDR < 0.05)
    
    # Split by direction
    inc1 <- df1 %>% filter(Log2FC > 0) %>% pull(name)
    inc2 <- df2 %>% filter(Log2FC > 0) %>% pull(name)
    dec1 <- df1 %>% filter(Log2FC < 0) %>% pull(name)
    dec2 <- df2 %>% filter(Log2FC < 0) %>% pull(name)
    
    tibble(
      patient_id = df_patient$patient_id[1],
      state1 = s1,
      state2 = s2,
      jaccard_increased = jaccard_index(inc1, inc2),
      jaccard_decreased = jaccard_index(dec1, dec2),
      n_inc_s1 = length(inc1),
      n_inc_s2 = length(inc2),
      n_dec_s1 = length(dec1),
      n_dec_s2 = length(dec2)
    )
  }) %>% bind_rows()
}

# Apply to all patients
jaccard_results <- all_results %>%
  group_by(patient_id) %>%
  group_split() %>%
  map_dfr(compute_jaccard_patient)


patient_summary <- jaccard_results %>%
  group_by(patient_id) %>%
  summarise(
    mean_jaccard_inc = mean(jaccard_increased, na.rm = TRUE),
    mean_jaccard_dec = mean(jaccard_decreased, na.rm = TRUE),
    n_states = n_distinct(c(state1, state2)),
    total_dags = sum(c(n_inc_s1, n_dec_s1, n_inc_s2, n_dec_s2), na.rm = TRUE) / 2, # rough average
    .groups = "drop"
  )


ggplot(patient_summary, aes(x = total_dags, y = mean_jaccard_inc, size = n_states)) +
  geom_point(color = "#e41a1c", alpha = 0.8) +
  geom_point(aes(y = mean_jaccard_dec), color = "#377eb8", alpha = 0.8) +
  labs(
    x = "Number of differentially accessible genes (total)",
    y = "Mean Jaccard index (same-direction overlap)",
    title = "Shared accessibility shifts across cell states per patient",
    subtitle = "Red = increased, Blue = decreased"
  ) +
  theme_minimal(base_size = 13)


ggplot(jaccard_results, aes(x = state1, y = state2, fill = jaccard_increased)) +
  geom_tile(color = "white") +
  facet_wrap(~patient_id, scales = "free") +
  scale_fill_gradient(low = "white", high = "#e41a1c", name = "Jaccard\nIncreased") +
  theme_minimal(base_size = 10) +
  labs(title = "Cross-state overlap of increased accessibility (Jaccard index)")



jaccard_long <- jaccard_results %>%
  select(patient_id, state1, state2,
         jaccard_increased, jaccard_decreased) %>%
  pivot_longer(cols = starts_with("jaccard_"),
               names_to = "direction", values_to = "jaccard") %>%
  mutate(direction = recode(direction,
                            jaccard_increased = "Increased accessibility",
                            jaccard_decreased = "Decreased accessibility"))

# Continuous blue-red palette, capped between 0 and 1
scale_fill_jaccard <- scale_fill_gradientn(
  colours = c("#f7fbff", "#6baed6", "#08306b"),
  limits = c(0, 1),
  name = "Jaccard index"
)

ggplot(jaccard_long, aes(x = state1, y = state2, fill = jaccard)) +
  geom_tile(color = "grey85") +
  facet_grid(patient_id ~ direction) +
  scale_fill_gradientn(colours = c("#f7fbff", "#6baed6", "#08306b"),
                       limits = c(0, 1),
                       name = "Jaccard index") +
  coord_equal() +
  labs(
    x = "Cell state 1",
    y = "Cell state 2",
    title = "Cross-state reproducibility of accessibility shifts",
    subtitle = "Each tile shows same-direction Jaccard index between cell states per patient"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

