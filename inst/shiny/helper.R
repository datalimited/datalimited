ramts <- readRDS("data/ramts.rds")
ramts <- ramts %>% filter(!is.na(c_touse)) # TODO won't this create gaps in years?
# stocks <- sort(unique(ramts$stocklong))

stocks <- sort(unique(ramts$stocklong))
stocks <- stocks[-1]
# stocks <- c(
  # "Black Grouper Gulf of Mexico",
  # "Black oreo West end of Chatham Rise",
#  "Dover sole Gulf of Alaska")
