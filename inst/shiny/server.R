library("shiny")
library("datalimited")
library("ggplot2")

shinyServer(
  function(input, output) {

    output$plot <- renderPlot({
      x <- comsir(ct = blue_gren$ct, yr = blue_gren$yr, k = 800, r = 0.6,
        nsim = 5e4, a = 0.8, x = 0.5, n_posterior = 5e3, start_r = c(0.2, 1.0),
        logistic_model = TRUE, obs = FALSE, norm_k = FALSE, norm_r = FALSE,
        logk = TRUE, norm_x = FALSE, norm_a = FALSE, normal_like = FALSE)
      bbmsy <- reshape2::dcast(x$quantities, sample_id ~ yr,
        value.var = "bbmsy")[,-1] # convert long to wide format
      bbmsy_out <- summarize_bbmsy(bbmsy, log = TRUE)
      bbmsy_out$year <- seq_len(nrow(bbmsy_out))

      ggplot(bbmsy_out, aes(year, bbmsy_q50)) + geom_line() +
        geom_ribbon(aes(ymin = bbmsy_q25, ymax = bbmsy_q75), alpha = 0.2) +
        geom_hline(yintercept = 1, lty = 2)
    })
  }
)
