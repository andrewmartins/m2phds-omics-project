library(dplyr)
library(ggplot2)
library(tidydr)

# loading data
participants = read.csv(here::here("data","participants.csv"))
participants

vaccine_hist = read.csv(here::here("data", "participant_vaccineHistory.csv"))

flow_outcome = read.csv(here::here("data","publicData_ex_vivo_flow.csv"))
# loading rnaSeq data for the first time
rnaSeq_2024 = read.csv(here::here("data","publicData_rnaseq_2024_UGA.csv"))
rnaSeq_2024 = rnaSeq_2024[rnaSeq_2024$timepoint == 0,]
rnaSeq_2024

# filtering on patients who have a 7 day plasmablast signal
flow_outcome = flow_outcome%>%
  filter(
    name == "Antibody-secreting_cells_(ASC)",
    timepoint %in% c(0, 7)
  ) %>% 
  distinct(subject, timepoint, .keep_all = TRUE) %>%  # drop the duplicate baselines
  select(participant_id, timepoint, value) %>%
  pivot_wider(
    id_cols     = participant_id,
    names_from  = timepoint,
    values_from = value,
    names_glue  = "{.value}_d{timepoint}"
  ) %>%
  filter(!is.na(value_d7)) %>%
  mutate(
    across(
      starts_with("value_"), log, .names = "log_{.col}"
    )
  ) %>%
  rename(value_d7 = plasmablast_frequency)
  

flow_outcome

# identifying the same patients with an exposure signal in RNAseq
common_ids <- intersect(
  flow_outcome$participant_id,
  rnaSeq_2024$participant_id
)
length(common_ids)

# filtering the RNAseq data to only include participants with a 7 day plasmablast signal and an exposure signal
rnaSeq_2024 = rnaSeq_2024[rnaSeq_2024$participant_id %in% common_ids,]

# filtering on participants who have a 7 day plasmablast signal and an exposure signal in RNAseq
participants = participants[participants$participant_id %in% common_ids,]

# filtering on participants' vaccine history who have a plasmablast signal at 7 days
vaccine_hist = vaccine_hist[vaccine_hist$subject %in% participants$subject,]

# joining the flow_outcome data to the participants data to get the subject column
participants = participants %>%
  inner_join(flow_outcome, by = "participant_id")

# aggregating the vaccine history data to get the number of prior vaccines for each participant
vaccine_hist_agg = vaccine_hist %>%
  group_by(subject) %>%
  summarise(
    n_prior_vax = sum(received %in% c("Yes_selfReported","Yes_Study")) - 1
  )

# calculating median number of prior vaccinations to impute vaccine history
med_n_vax = median(vaccine_hist_agg$n_prior_vax, na.rm = TRUE)

# calculating median day 0 plasmablast frequency to impute missing data
med_d0_pf = median(participants$log_value_d0, na.rm = TRUE)

# joining the vaccine history data to the participants data
participants = left_join(participants, vaccine_hist_agg, by = "subject")

# median imputation
participants_cleaned = participants %>%
  mutate(n_prior_vax = ifelse(is.na(n_prior_vax), med_n_vax, n_prior_vax))

participants_cleaned = participants_cleaned %>%
  mutate(log_value_d0 = ifelse(is.na(log_value_d0), med_d0_pf, log_value_d0))


saveRDS(participants_cleaned, here::here("data", "participants_cleaned.rds"))
saveRDS(rnaSeq_2024, here::here("data", "rnaseq.rds"))
saveRDS(vaccine_hist, here::here("data", "vaccine_hist.rds"))