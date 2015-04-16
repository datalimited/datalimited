library("shiny")
library("datalimited")
library("ggplot2")
library("dplyr")

shinyServer(
  function(input, output) {

    #     output$plot_comsir <- renderPlot({
    #       dat <- filter(ramts, stocklong == input$stock)
    #       x <- comsir(ct = dat$c_touse, yr = dat$year, k = 800, r = 0.6,
    #         nsim = 3e4, a = 0.8, x = 0.5, n_posterior = 1e3, start_r = resilience(NA),
    #         logistic_model = TRUE, obs = FALSE, norm_k = FALSE, norm_r = FALSE,
    #         logk = TRUE, norm_x = FALSE, norm_a = FALSE, normal_like = FALSE)
    #       bbmsy <- reshape2::dcast(x$quantities, sample_id ~ yr,
    #         value.var = "bbmsy")[,-1] # convert long to wide format
    #       bbmsy_out <- summarize_bbmsy(bbmsy, log = TRUE)
    #       bbmsy_out$year <- dat$year
    #
    #       ggplot(bbmsy_out, aes(year, bbmsy_q50)) + geom_line() +
    #         geom_ribbon(aes(ymin = bbmsy_q25, ymax = bbmsy_q75), alpha = 0.2) +
    #         geom_ribbon(aes(ymin = bbmsy_q2.5, ymax = bbmsy_q97.5), alpha = 0.1) +
    #         geom_hline(yintercept = 1, lty = 2) +
    #         geom_line(data = dat, aes(year, b_bmsy_touse), colour = "blue")
    #     })

    output$plot_cmsy <- renderPlot({
      dat <- filter(ramts, stocklong == input$stock)
      cmsy_out <- cmsy(
        yr             = dat$year,
        ct             = dat$c_touse,
        #prior_log_mean = dat$log_mean[1],
        #prior_log_sd   = dat$log_sd[1],
        prior_log_mean = log(input$prior_mean),
        prior_log_sd   = log(input$prior_sd),
        start_r        = resilience(input$resilience),
        sig_r          = input$sig_r,
        reps           = input$reps,
        interbio       = input$interbio,
        revise_bounds  = input$revise_bounds,
        interyr_index  = input$interyr_index)

      bbmsy <- reactive({
        validate(
          need(!is.null(cmsy_out), "Insufficient samples drawn")
        )
        cmsy_out$biomass[, -1] / cmsy_out$bmsy
      })

      bbmsy_out <- summarize_bbmsy(bbmsy(), log = TRUE)
      bbmsy_out$year <- dat$year
      ggplot(bbmsy_out, aes(year, bbmsy_q50)) + geom_line(lwd = 1) +
        geom_ribbon(aes(ymin = bbmsy_q25, ymax = bbmsy_q75), alpha = 0.2) +
        geom_ribbon(aes(ymin = bbmsy_q2.5, ymax = bbmsy_q97.5), alpha = 0.1) +
        geom_hline(yintercept = 1, lty = 2) + theme_bw() + xlab("Year") +
        ylab(expression(B/B[MSY])) +
        geom_line(data = dat, aes(year, b_bmsy_touse), colour = "red", lwd = 1) +
        theme( plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
    })
  }
)
