library(here)
library(tidyverse)

df = readr::read_csv(here::here("config/samples_old.csv"))

# rename columns

df_new = df %>% 
  dplyr::rename(sample = run,
                of_start = start_of,
                of_end = end_of,
                no_start = start_no,
                no_end = end_no,
                fps = frame_rate) %>% 
  tidyr::separate(.,
                  col = sample,
                  into = c("date", "time", NA, NA, NA),
                  sep = "_",
                  remove = F) %>% 
  dplyr::mutate(adj_top = NA,
                adj_bottom = NA,
                adj_left = NA,
                adj_right = NA,
                of_video_length = of_end - of_start,
                no_video_length = no_end - no_start,
                intensity_floor = 0,
                intensity_ceiling = 155,
                area_floor = 80,
                area_ceiling = 200,
                cab_coords_of_q1 = NA,
                cab_coords_of_q2 = NA,
                cab_coords_of_q3 = NA,
                cab_coords_of_q4 = NA,
                cab_coords_no_q1 = NA,
                cab_coords_no_q2 = NA,
                cab_coords_no_q3 = NA,
                cab_coords_no_q4 = NA
                ) %>% 
  dplyr::select(-ends_with("dims")) %>% 
  dplyr::select(sample, date, time, tank_side, of_start, of_end, no_start, no_end, fps,
                everything()) %>% 
  readr::write_csv(here::here("config/samples.csv"))
