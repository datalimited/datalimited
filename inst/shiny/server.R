library("shiny")
library("datalimited")
library("ggplot2")
library("dplyr")

shinyServer(
  function(input, output) {

    # a generic function we'll use to plot B/Bmsy time series
    plot_bbmsy <- function(est_dat, orig_dat) {

      this_theme <-  theme_bw() + theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

      p1 <- ggplot(est_dat, aes(year, bbmsy_q50)) + geom_line(lwd = 1) +
        geom_ribbon(aes(ymin = bbmsy_q25, ymax = bbmsy_q75), alpha = 0.2) +
        geom_ribbon(aes(ymin = bbmsy_q2.5, ymax = bbmsy_q97.5), alpha = 0.1) +
        geom_hline(yintercept = 1, lty = 2) + xlab("Year") +
        ylab(expression(B/B[MSY])) + ylim(0, 3) +
        geom_line(data = orig_dat, aes(year, b_bmsy_touse), colour = "red", lwd = 1) +
        this_theme + ylab("Catch")

      p2 <- ggplot(orig_dat, aes(year, c_touse), colour = "black", lwd = 1) +
        geom_line() + this_theme
      gridExtra::grid.arrange(p1, p2, ncol = 1)
    }


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
        reps           = input$cmsy_reps,
        interbio       = input$interbio,
        revise_bounds  = input$revise_bounds,
        interyr_index  = input$interyr_index)

      bbmsy_cmsy <- reactive({
        validate(
          need(!is.null(cmsy_out), "Insufficient samples drawn")
        )
        cmsy_out$biomass[, -1] / cmsy_out$bmsy
      })

      bbmsy_out <- summarize_bbmsy(bbmsy_cmsy(), log = TRUE)
      bbmsy_out$year <- dat$year
      plot_bbmsy(bbmsy_out, dat)
    })

    output$plot_comsir <- renderPlot({
      print(input$comsir_a) # doesn't re-render without these?
      print(input$comsir_x)
      print(input$comsir_cv)
      dat <- filter(ramts, stocklong == input$stock)
      comsir_out <- comsir(ct = dat$c_touse,
        yr = dat$year,
        k = input$comsir_k,
        r = input$comsir_r,
        nsim = input$comsir_reps,
        a = input$comsir_a,
        x = input$comsir_x,
        cv = input$comsir_cv,
        n_posterior = input$comsir_n_posterior,
        start_r = resilience(input$resilience),
        logistic_model = input$comsir_logistic,
        obs = input$comsir_obs,
        norm_k = input$comsir_norm_k,
        norm_r = input$comsir_norm_r,
        logk = input$comsir_logk,
        norm_x = input$comsir_norm_x,
        norm_a = input$comsir_norm_a,
        normal_like = input$comsir_normal_like)

      bbmsy_comsir <- reactive({
        validate(
          need(!is.null(comsir_out), "Insufficient samples drawn")
        )
        bbmsy <- reshape2::dcast(comsir_out$quantities, sample_id ~ yr,
          value.var = "bbmsy")[,-1] # convert long to wide format
        bbmsy
      })

      bbmsy_out <- summarize_bbmsy(bbmsy_comsir(), log = TRUE)
      bbmsy_out$year <- dat$year

      plot_bbmsy(bbmsy_out, dat)
    })


    output$plot_prm <- renderPlot({
      dat <- filter(ramts, stocklong == input$stock)

      dat_formatted <- format_prm(year = dat$year, catch = dat$c_touse,
        bbmsy = dat$b_bmsy_touse, species_cat = dat$spp_category[1L])

      out <- predict_prm(dat_formatted, ci = TRUE, level = 0.95)
      out$year <- dat_formatted$year
      out <- rename(out, bbmsy_q50 = fit,
        bbmsy_q2.5 = lower,
        bbmsy_q97.5 = upper)

      out2 <- predict_prm(dat_formatted, ci = TRUE, level = 0.50)
      out2$year <- dat_formatted$year
      out2 <- rename(out2,
        bbmsy_q25 = lower,
        bbmsy_q75 = upper) %>%
        select(-fit)
      out <- inner_join(out, out2)

      plot_bbmsy(out, dat)
    })


  }
)
