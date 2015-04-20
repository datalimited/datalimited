library("dplyr")

sci_names <- ram_ts[,c("stockid", "scientificname")]
sci_names <- sci_names[!duplicated(sci_names), ]

ramts <- readRDS("data/ramts.rds") %>%
  inner_join(sci_names)

# spp_categories is built into the datalimited package:
ramts <- dplyr::inner_join(ramts, spp_categories, by = "scientificname")

ramts <- ramts %>% filter(!is.na(c_touse))
# stocks <- sort(unique(ramts$stocklong))

stocks <- sort(unique(ramts$stocklong))
stocks <- stocks[-1] # let's not start on a broken stock
