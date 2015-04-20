library(dplyr)
library(ggplot2)
orig <- readRDS("../ensembles/generated-data/ram-orig-fits.rds")
orig <- dplyr::filter(orig, method == "Costello") %>%
  rename(stockid = stock)

load_all()
new <- ram_prm_dat
new$pred <- predict_prm(ram_prm_dat)

q <- dplyr::inner_join(dplyr::select(orig, stockid, year, b2bmsy),
  dplyr::select(new, stockid, year, pred)) %>%
  arrange(stockid, year)

# orig2 <- dplyr::inner_join(orig, select(new, stockid, bbmsy, year))

ggplot(q[1:300, ], aes(b2bmsy, pred)) + geom_point() +
  facet_wrap(~stockid, scales = "free")
# ggplot(new[1:300, ], aes(bbmsy, pred)) + geom_point() +
#   facet_wrap(~stockid, scales = "free")
# ggplot(orig2[1:300, ], aes(bbmsy, b2bmsy)) + geom_point() +
#   facet_wrap(~stockid, scales = "free")
