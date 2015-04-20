# Fit Costello et al.-style panel regression model to a cached version of the
# RAM database

library(dplyr)
devtools::load_all()
d <- dplyr::inner_join(ramts, spp_categories, by = "scientificname")

ram_prm_dat <- plyr::ddply(d, "stockid", function(x) {
  format_prm(year = x$year, catch = x$catch, bbmsy = x$bbmsy_ram,
    species_cat = x$spp_category[1L])
}) %>%
  dplyr::inner_join(dplyr::select(d, -bbmsy_ram, -catch))

devtools::use_data(ram_prm_dat, overwrite = TRUE)

ram_prm_model <- fit_prm(ram_prm_dat)
devtools::use_data(ram_prm_model, overwrite = TRUE)

# TODO reduce the size of this lm() object
