
message("Saving Image")
save.image("combine_sceptre_pwr_sim.rda")
#stop()


# opening log file to collect all messages, warnings and errors
message("Opening log file")
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


suppressPackageStartupMessages({
  library(tidyverse)
})


# Combine each file into one tsv and save it
message("Combining power analysis outputs...")

# Initialize an empty list to store data frames
data_list <- list()

# Loop over input files
for (file_path in snakemake@input) {
  # Read the current file
  current_data <- read_tsv(file_path, col_types = cols())
  
  # Append the data frame to the list
  data_list[[length(data_list) + 1]] <- current_data
}

# Combine all data frames into one
message("Combining all data frames into one")
combined_data <- do.call(rbind, data_list)

# Remove duplicate rows based on "response_id", "grna_target", and "rep" columns
message("Removing duplicate rows")
filt_combined_data <- combined_data %>%
  distinct(response_id, grna_target, rep, .keep_all = TRUE)
message(paste0((nrow(combined_data) - nrow(filt_combined_data))/20, " duplicate rows removed"))

# Write the combined data frame to the output file
write_tsv(filt_combined_data, snakemake@output[[1]])

message("Power analysis outputs combined successfully.")



# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)
