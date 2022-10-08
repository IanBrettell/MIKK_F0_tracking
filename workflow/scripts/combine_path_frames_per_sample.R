# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = as.list(list.files("/hps/software/users/birney/ian/repos/MIKK_F0_tracking/book/figs/paths_frame_grabs/0.05/dist_angle/15/open_field",
                        full.names = T))


## True
IN = snakemake@input
OUT_PNG = snakemake@output[["png"]]
OUT_PDF = snakemake@output[["pdf"]]
  

# Put into data frame

df = tibble::tibble(PATH = unlist(IN)) %>% 
  dplyr::mutate(BASENAME = basename(PATH)) %>% 
  tidyr::separate(col = "BASENAME",into = c(NA, NA, "SAMPLE", NA, "RUN", NA),sep = "_")



fig = cowplot::ggdraw() +
  cowplot::draw_image(df %>% 
                        dplyr::filter(SAMPLE == "22-1" & RUN == "A") %>% 
                        dplyr::pull(PATH),
                      x = 0, y = 0.66, width = 0.5, height = 0.33) +
  cowplot::draw_image(df %>% 
                        dplyr::filter(SAMPLE == "22-1" & RUN == "B") %>% 
                        dplyr::pull(PATH),x = 0.5, y = 0.66, width = 0.5, height = 0.33) +
  cowplot::draw_image(df %>% 
                        dplyr::filter(SAMPLE == "18-2" & RUN == "A") %>% 
                        dplyr::pull(PATH),x = 0, y = 0.33, width = 0.5, height = 0.33) +
  cowplot::draw_image(df %>% 
                        dplyr::filter(SAMPLE == "18-2" & RUN == "B") %>% 
                        dplyr::pull(PATH),x = 0.5, y = 0.33, width = 0.5, height = 0.33) +
  cowplot::draw_image(df %>% 
                        dplyr::filter(SAMPLE == "10-1" & RUN == "A") %>% 
                        dplyr::pull(PATH),x = 0, y = 0, width = 0.5, height = 0.33) +
  cowplot::draw_image(df %>% 
                        dplyr::filter(SAMPLE == "10-1" & RUN == "B") %>% 
                        dplyr::pull(PATH),x = 0.5, y = 0, width = 0.5, height = 0.33) +
  cowplot::draw_plot_label(c("A", "B", "C"),x = c(0,0,0), y = c(1, 0.66, 0.33)) +
  cowplot::draw_plot_label(c("22-1 (David)", "18-2 (Elsa)", "10-1 (Janeway)"),x = c(0.34,0.35,0.29), y = c(1, 0.67, 0.35),size = 13,colour = c("#FB737A", "#FF66A6", "#F8766D"))


ggsave(OUT_PNG,
       fig,
       device = "png",
       width = 5.5,
       height = 8,
       units = "in")

ggsave(OUT_PDF,
       fig,
       device = "pdf",
       width = 5.5,
       height = 8,
       units = "in")
