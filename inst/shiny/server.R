library("shiny")
library("datalimited")
library("ggplot2")
library("dplyr")




shinyServer(
  function(input, output) {


    output$plot_comsir <- renderPlot({
      dat <- filter(ramts, stocklong == input$stock)
      x <- comsir(ct = dat$c_touse, yr = dat$year, k = 800, r = 0.6,
        nsim = 3e4, a = 0.8, x = 0.5, n_posterior = 1e3, start_r = input$start_r,
        logistic_model = TRUE, obs = FALSE, norm_k = FALSE, norm_r = FALSE,
        logk = TRUE, norm_x = FALSE, norm_a = FALSE, normal_like = FALSE)
      bbmsy <- reshape2::dcast(x$quantities, sample_id ~ yr,
        value.var = "bbmsy")[,-1] # convert long to wide format
      bbmsy_out <- summarize_bbmsy(bbmsy, log = TRUE)
      bbmsy_out$year <- dat$year

      ggplot(bbmsy_out, aes(year, bbmsy_q50)) + geom_line() +
        geom_ribbon(aes(ymin = bbmsy_q25, ymax = bbmsy_q75), alpha = 0.2) +
        geom_hline(yintercept = 1, lty = 2)
    })

    output$plot_cmsy <- renderPlot({
      dat <- filter(ramts, stocklong == input$stock)
      cmsy_out <- cmsy(
        yr             = dat$year,
        ct             = dat$c_touse,
        prior_log_mean = dat$log_mean[1],
        prior_log_sd   = dat$log_sd[1],
        start_r        = input$start_r,
        sig_r          = 0.05,
        reps           = 1e4)
      # if(!is.null(cmsy_out)) {
      bbmsy <- cmsy_out$biomass[, -1] / cmsy_out$quantities$bmsy
      bbmsy[is.infinite(bbmsy)] <- NA  # TODO investigate the -Inf values
      # TODO switch log to TRUE:
      bbmsy_out <- summarize_bbmsy(bbmsy, na.rm = TRUE, log = FALSE) # TODO investigate the negative values
      bbmsy_out$year <- dat$year
      ggplot(bbmsy_out, aes(year, bbmsy_q50)) + geom_line() +
        geom_ribbon(aes(ymin = bbmsy_q25, ymax = bbmsy_q75), alpha = 0.2) +
        geom_hline(yintercept = 1, lty = 2)
    })
  }
)
