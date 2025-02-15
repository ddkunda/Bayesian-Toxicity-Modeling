library(shiny)
library(rstan)
library(ggplot2)
library(dplyr)
library(boot)  
library(parallel)

# For faster Stan compilation (caching) and parallel sampling:
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---------------------------------------------------------------------------
# STAN MODEL: Binomial likelihood with logistic link
#
#   p_tox(dose[k]) = inv_logit(alpha + beta * dose[k])
#   y[k] ~ binomial(n[k], p_tox)
#
#   Priors:
#     alpha ~ Normal(prior_alpha_mu, prior_alpha_sd)
#     beta  ~ Normal(prior_beta_mu,  prior_beta_sd)
# ---------------------------------------------------------------------------
stan_code <- "
data {
  int<lower=1> K;               // number of dose groups
  array[K] real dose;           // dose values
  array[K] int y;               // number of toxicities at dose group k
  array[K] int n;               // total mice at dose group k
  
  real prior_alpha_mu;
  real<lower=0> prior_alpha_sd;
  real prior_beta_mu;
  real<lower=0> prior_beta_sd;
}
parameters {
  real alpha;
  real beta;
}
model {
  // Priors
  alpha ~ normal(prior_alpha_mu, prior_alpha_sd);
  beta  ~ normal(prior_beta_mu,  prior_beta_sd);
  
  // Binomial likelihood
  for (k in 1:K){
    real p = inv_logit(alpha + beta * dose[k]);
    y[k] ~ binomial(n[k], p);
  }
}
"

# Compile the Stan model once globally
stan_model_obj <- stan_model(model_code = stan_code)

# ---------------------------------------------------------------------------
# SHINY UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Bayesian Toxicity Modeling"),
  
  sidebarLayout(
    sidebarPanel(
      h4(HTML("<u>1. Input Data</u>")),
      helpText("Enter vectors of equal length, comma-separated."),
      textAreaInput("dose_levels", 
                    "Dose Levels (e.g., mg/kg)", 
                    "0, 5, 10, 20, 40", 
                    rows = 2),
      textAreaInput("tox_events",  
                    "Number of Toxicities at Each Dose", 
                    "0, 1, 3, 4, 9", 
                    rows = 2),
      textAreaInput("group_size",  
                    "Total Mice in Each Dose Group",     
                    "5, 5, 5, 5, 5", 
                    rows = 2),
      
      h4(HTML("<u>2. Priors</u>")),
      fluidRow(
        column(6, numericInput("alpha_mu", "Alpha prior mean", 0, step = 1)),
        column(6, numericInput("alpha_sd", "Alpha prior SD",   5, step = 1))
      ),
      fluidRow(
        column(6, numericInput("beta_mu", "Beta prior mean", 0, step = 0.5)),
        column(6, numericInput("beta_sd", "Beta prior SD",   2, step = 0.5))
      ),
      
      h4(HTML("<u>3. Maximum Tolerated Dose (MTD)</u>")),
      numericInput("tox_threshold", 
                   "Toxicity Probability Threshold (e.g., 0.3 = 30%)", 
                   value = 0.3, 
                   step = 0.05),
      helpText("We'll estimate the dose at which Probability(Toxicity) <= threshold."),
      
      actionButton("runAnalysis", "Run Bayesian Analysis")
    ),
    
    mainPanel(
      h4(HTML("<u>Posterior Summaries</u>")),
      verbatimTextOutput("posteriorSummary"),
      
      h4(HTML("<u>Posterior Dose-Toxicity Curve</u>")),
      plotOutput("toxCurvePlot", height = "400px"),
      
      h4(HTML("<u>Estimated MTD</u>")),
      
      verbatimTextOutput("mtdText"),
      
      h5(
        span("All questions can be sent to Christian Dide-Agossou, PhD:", style = "color:black; display: inline-block;"),
        span(htmlOutput("uicmt14"), style = "display: inline-block;")
      )
    )
  )
)

# ---------------------------------------------------------------------------
# SHINY SERVER
# ---------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # Run Bayesian model after user clicks "Run Bayesian Analysis"
  runModel <- eventReactive(input$runAnalysis, {
    # 1. Parse user data
    dose_vec <- as.numeric(strsplit(input$dose_levels, "[,\\s]+")[[1]])
    y_vec    <- as.integer(strsplit(input$tox_events,  "[,\\s]+")[[1]])
    n_vec    <- as.integer(strsplit(input$group_size,  "[,\\s]+")[[1]])
    
    # Basic validation
    validate(
      need(length(dose_vec) == length(y_vec) &&
             length(y_vec)  == length(n_vec),
           "Dose, #Toxicities, and Group Size vectors must all have the same length!"),
      need(all(y_vec <= n_vec),
           "Number of toxicities cannot exceed the group size in any dose!"),
      need(all(y_vec >= 0),
           "Toxicities must be non-negative."),
      need(all(n_vec > 0),
           "Group size must be > 0 for each dose.")
    )
    
    K <- length(dose_vec)
    
    # 2. Create data list for Stan
    stan_data_list3 <- list(
      K = K,
      dose = dose_vec,
      y    = y_vec,
      n    = n_vec,
      prior_alpha_mu = input$alpha_mu,
      prior_alpha_sd = input$alpha_sd,
      prior_beta_mu  = input$beta_mu,
      prior_beta_sd  = input$beta_sd
    )
    
    # 3. Fit with Stan
    fit <- sampling(
      object = stan_model_obj,
      data   = stan_data_list3,
      iter   = 3000,
      warmup = 1000,
      chains = 4,
      seed   = 42,
      control = list(adapt_delta = 0.95, max_treedepth = 15)
    )
    
    list(
      fit      = fit, 
      dose_vec = dose_vec,
      y_vec    = y_vec,
      n_vec    = n_vec
    )
  })
  
  # Posterior summary
  output$posteriorSummary <- renderPrint({
    req(runModel())
    fit_obj <- runModel()$fit
    print(summary(fit_obj, pars = c("alpha", "beta"),
                  probs = c(0.025, 0.5, 0.975))$summary)
  })
  
  # Posterior Dose-Toxicity Curve
  output$toxCurvePlot <- renderPlot({
    req(runModel())
    fit_obj  <- runModel()$fit
    dose_vec <- runModel()$dose_vec
    y_vec    <- runModel()$y_vec
    n_vec    <- runModel()$n_vec
    
    draws <- rstan::extract(fit_obj)
    
    # Create a dose grid from 0 up to ~ max(dose_vec)*1.2
    dose_grid <- seq(0, max(dose_vec)*1.2, length.out = 100)
    
    # p(dose) = inv_logit(alpha + beta*dose)
    p_mat <- sapply(dose_grid, function(d){
      alpha_draws <- draws$alpha
      beta_draws  <- draws$beta
      inv.logit(alpha_draws + beta_draws * d)
    })
    
    p_mean  <- apply(p_mat, 2, mean)
    p_lower <- apply(p_mat, 2, quantile, 0.025)
    p_upper <- apply(p_mat, 2, quantile, 0.975)
    
    df_plot <- data.frame(
      dose  = dose_grid,
      p     = p_mean,
      lower = p_lower,
      upper = p_upper
    )
    
    # Observed proportion of toxicity at each dose
    prop_vec <- y_vec / n_vec
    obs_dat  <- data.frame(
      dose = dose_vec,
      prop = prop_vec
    )
    
    ggplot() +
      geom_ribbon(
        data = df_plot,
        aes(x = dose, ymin = lower, ymax = upper),
        fill = "blue", alpha = 0.2
      ) +
      geom_line(
        data = df_plot,
        aes(x = dose, y = p),
        color = "blue", linewidth = 1
      ) +
      geom_point(
        data = obs_dat,
        aes(x = dose, y = prop),
        color = "red", size = 3, alpha = 0.6
      ) +
      labs(
        x = "Dose (e.g., mg/kg)",
        y = "Observed/Posterior Probability of Toxicity",
        title = "Posterior Toxicity Curve (Aggregated Binomial)"
      ) +
      theme_minimal(base_size = 14)
  })
  
  # Calculate MTD: dose where P(Toxicity) <= threshold
  output$mtdText <- renderPrint({
    req(runModel())
    fit_obj   <- runModel()$fit
    draws     <- rstan::extract(fit_obj)
    threshold <- input$tox_threshold
    
    alpha_draws <- draws$alpha
    beta_draws  <- draws$beta
    
    # Solve for dose: logit(threshold) = alpha + beta * dose
    dose_draws <- (qlogis(threshold) - alpha_draws) / beta_draws
    
    # Summarize
    mtd_median <- median(dose_draws)
    ci_lower   <- quantile(dose_draws, 0.025)
    ci_upper   <- quantile(dose_draws, 0.975)
    
    cat("Estimated MTD (median of posterior) where P(Toxicity) =",
        threshold, ":\n")
    cat(sprintf("  MTD: %.2f (e.g., mg/kg) [95%% CI: %.2f, %.2f]\n",
                mtd_median, ci_lower, ci_upper))
  })
  
  output$uicmt14 <- renderUI({
    email <- "christian.dideagossou@gmail.com"
    link <- paste0("mailto:", email)
    tags$a(href = link, email)  # Use tags$a to create the hyperlink
  })
}

shinyApp(ui = ui, server = server)
