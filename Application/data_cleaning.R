#!/usr/bin/env Rscript
# data_cleaning.R
# This script reads raw symmetric bluetooth proximity data (bt_symmetric.csv),
# filters out noise, defines continuous co-location events from discrete scans,
# and outputs clean events and actor mappings for modeling.

rm(list=ls())
library(data.table)

# 1. Load raw proximity data
# bt_symmetric.csv contains Bluetooth signal strength (RSSI) scans over time between actors
proximity <- fread("Application/Data/bt_symmetric.csv")
colnames(proximity) = c("time", "sender", "receiver", "signal_strength")

# 2. Filter data
# Delete empty scans (indicated by receiver == -1 and signal_strength == 0)
proximity = proximity[!(receiver == -1 & signal_strength == 0)]
# Delete scans with receivers outside of the study/experiment (receiver == -2)
proximity = proximity[!(receiver == -2)]

# 3. Ensure undirected dyad property
# Since proximity events are undirected, we enforce that sender < receiver (i.e. from < to)
# to keep dyads uniquely identified and consistent.
tmp_ind = proximity$sender > proximity$receiver
tmp = proximity$sender
proximity$sender[tmp_ind] = proximity$receiver[tmp_ind]
proximity$receiver[tmp_ind] = tmp[tmp_ind]

# Construct a unique key for each dyad
proximity$from_to = paste(proximity$sender, proximity$receiver)

# Sort chronologically and by dyad to prepare for consecutive signal tracking
proximity = proximity[order(time)]
proximity = proximity[order(from_to)]

# Filter out very weak/noisy signals (keep signal strength stronger than -90 dBm)
proximity = proximity[signal_strength > -90]

# Extract unique actors present in the proximity data
unique_actors = unique(c(proximity$sender, proximity$receiver))


# 4. Define Co-location Event Windows (tmp_fun)
# This function identifies when a series of discrete RSSI scans represents
# a continuous co-location event, and outputs its start and end times.
tmp_fun = function(x, min) {
  # Calculate time differences between consecutive scans for this dyad
  x$diff = c(0, diff(x$time))
  
  # A gap larger than the threshold (min = 10 minutes) indicates the end of an event
  x$is_change = x$diff > min
  
  # The first scan in the sequence always marks a change/new event start
  x$is_change[1] = TRUE
  
  # Determine change states before and after each scan
  x$lag_change = c(TRUE, x$is_change[-length(x$is_change)])
  x$future_change = c(x$is_change[-1], TRUE)
  
  # Start of co-location: marked by a change where the next scan does not immediately break the event
  x$is_start = x$is_change & !x$future_change
  
  # End of co-location: no change here, but the next scan marks a new event (gap > threshold)
  x$is_end = !x$is_change & x$future_change
  
  # Keep only the start and end boundary scans of each event
  x = x[is_start == TRUE | is_end == TRUE]
  return(x)
}

# Apply the event-extraction function to each dyad (.SD) with a 10-minute threshold (10*60 seconds)
proximity <- proximity[, tmp_fun(.SD, min = 10*60), by = from_to]

# 5. Format and standardise event variables
# Type 1 = event start (incidence), Type 0 = event end (dissolution)
proximity$type = as.numeric(proximity$is_start == 1)
proximity_events = proximity[,.(time, from = sender, to = receiver, type)]

# Map actor IDs to sequential integers starting at 1
proximity_events$from = match(proximity_events$from, unique_actors)
proximity_events$to = match(proximity_events$to, unique_actors)
proximity_events = proximity_events[!is.na(from) & !is.na(to)]

# Sort the final dataset chronologically
proximity_events = proximity_events[order(time)]

# Set start time of the experiment as t = 0
proximity_events$time_scaled = as.numeric(proximity_events$time - min(proximity_events$time))

# Re-enforce undirected dyad sorting (from < to) after ID mapping
tmp_ind = proximity_events$from > proximity_events$to
tmp = proximity_events$from
proximity_events$from[tmp_ind] = proximity_events$to[tmp_ind]
proximity_events$to[tmp_ind] = tmp[tmp_ind]

# Prepare node ID mapping table
actor_data <- data.frame(1:length(unique_actors), unique_actors)

# 6. Save clean outputs
write.table(proximity_events, file = "Application/Data/events.csv", sep = ",", row.names = FALSE)
write.table(actor_data, file = "Application/Data/actor_data.csv", sep = ",", row.names = FALSE)

