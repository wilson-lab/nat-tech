##############################################
### Work out edit credits for BANC neurons ###
##############################################

# Load libraries
library(bancr)
library(tidyverse)

# Read data
data <- readxl::read_xlsx("data/hartman_et_al_2024/bm_an.xlsx")

# Get edits
b.changes <- bancr::banc_change_log(data$root_id)

# Map to known users
b.changes.users <- b.changes %>%
  dplyr::mutate(user_affiliation = dplyr::case_when(
    grepl("Laia|Ruai|Borja",user_name) ~ "Aelysia",
    TRUE ~ user_affiliation
  )) %>%
  dplyr::group_by(user_name) %>%
  dplyr::mutate(edits = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(user_name, user_affiliation, edits) %>%
  dplyr::arrange(dplyr::desc(edits))

# Reorder by total edits
b.changes.users <- b.changes.users %>%
  mutate(user_name = reorder(user_name, -edits))

# Create the bar plot
ggplot(b.changes.users, aes(x = user_name, y = edits, fill = user_affiliation)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Accumulated Edits by User",
    x = "User Name",
    y = "Total Edits",
    fill = "User Affiliation"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotates x-axis labels for better readability

# Save
readr::write_csv(b.changes.users, file = "data/hartman_et_al_2024/bm_an_credit.csv")
