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



flow_outcome

# filtering on patients who have a 7 day plasmablast signal
flow_outcome = flow_outcome%>%
  filter(
    name == "Antibody-secreting_cells_(ASC)",
    timepoint == 7
  )

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
  inner_join(flow_outcome, by = "participant_id") %>%
  select(-subject.x) %>%
  rename(subject = subject.y)

# aggregating the vaccine history data to get the number of prior vaccines for each participant
vaccine_hist_agg = vaccine_hist %>%
  group_by(subject) %>%
  summarise(
    n_prior_vax = sum(received %in% c("Yes_selfReported","Yes_Study"))
  )

# joining the vaccine history data to the participants data
participants = left_join(participants, vaccine_hist_agg, by = "subject")
participants_cleaned = participants %>%
  mutate(n_prior_vax = coalesce(n_prior_vax,0))

saveRDS(participants_cleaned, here::here("data", "participants_cleaned.rds"))
saveRDS(rnaSeq_2024, here::here("data", "rnaseq.rds"))
saveRDS(vaccine_hist, here::here("data", "vaccine_hist.rds"))