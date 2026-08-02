suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
  library(ggthemes)
  library(DT)
})

find_input <- function() {
  candidates <- c(
    file.path("data", "candidate_scorecard_input_reference.csv"),
    file.path("..", "..", "doc", "data", "candidates",
      "candidate_scorecard_input_reference.csv"),
    file.path("doc", "data", "candidates",
      "candidate_scorecard_input_reference.csv")
  )
  found <- candidates[file.exists(candidates)]
  if (!length(found)) {
    stop("Could not locate candidate_scorecard_input_reference.csv")
  }
  found[[1L]]
}

scorecard <- read.csv(find_input(), check.names = FALSE,
  stringsAsFactors = FALSE)
scorecard$include <- as.logical(scorecard$include)
scorecard$raw_value <- as.numeric(scorecard$raw_value)

default_metrics <- unique(scorecard$metric[scorecard$include])
all_metrics <- unique(scorecard$metric)
all_cmps <- unique(scorecard$mp)

cmp_colors <- c(
  "HS+20 (MP29)" = "#C87A8A",
  "HS-20 (MP43)" = "#B6875B",
  "HSsym (MP45)" = "#909646",
  "HS-30 (MP47)" = "#55A067",
  "PR+20 (MP32)" = "#00A396",
  "PR-20 (MP44)" = "#409BBB",
  "PRsym (MP46)" = "#9189C7",
  "PR-30 (MP48)" = "#BE7AB4"
)

stock_metrics <- c(
  "SB / SB[MSY]", "F / F[MSY]", "VB / VB[2025]", "VB / VB[MSY]",
  "P(Kobe red)"
)
fishing_metrics <- c("Catch", "IACC", "Mean catch reduction")

normalize_preference <- function(rows) {
  rows$preference <- NA_real_
  for (metric_name in unique(rows$metric)) {
    use <- rows$metric == metric_name
    values <- rows$raw_value[use]
    direction <- unique(rows$preferred_direction[use])
    if (length(direction) != 1L) {
      stop("Inconsistent direction for metric: ", metric_name)
    }
    span <- max(values) - min(values)
    if (span == 0) {
      rows$preference[use] <- 1
    } else if (direction == "higher is better") {
      rows$preference[use] <- (values - min(values)) / span
    } else {
      rows$preference[use] <- (max(values) - values) / span
    }
  }
  rows
}

ui <- page_sidebar(
  title = "JM MSE candidate scorecard explorer",
  theme = bs_theme(version = 5, bootswatch = "cosmo"),
  sidebar = sidebar(
    width = 340,
    selectInput(
      "scheme", "Weighting scheme",
      choices = c(
        "Equal weights" = "equal",
        "Balanced fishing and stock condition" = "balanced",
        "Dispersion weighted (square-root CV)" = "dispersion",
        "Specify weights manually" = "custom"
      ),
      selected = "dispersion"
    ),
    checkboxGroupInput("cmps", "CMPs", choices = all_cmps,
      selected = all_cmps),
    checkboxGroupInput("metrics", "Performance metrics",
      choices = all_metrics, selected = default_metrics),
    conditionalPanel(
      condition = "input.scheme == 'custom'",
      helpText("Enter non-negative relative weights. They are normalized to sum to 100%."),
      uiOutput("custom_weights")
    ),
    hr(),
    helpText(paste(
      "Scores are relative to the selected CMP set. They describe trade-offs,",
      "not absolute acceptability or management preference."
    ))
  ),
  navset_card_tab(
    nav_panel(
      "Scores",
      plotOutput("score_plot", height = "480px"),
      DTOutput("score_table")
    ),
    nav_panel(
      "Performance quilt",
      plotOutput("quilt_plot", height = "560px")
    ),
    nav_panel(
      "Weights",
      DTOutput("weight_table")
    )
  )
)

server <- function(input, output, session) {
  selected_rows <- reactive({
    req(length(input$cmps) >= 2L, length(input$metrics) >= 1L)
    rows <- scorecard[
      scorecard$mp %in% input$cmps & scorecard$metric %in% input$metrics,
      , drop = FALSE
    ]
    validate(need(
      nrow(rows) == length(input$cmps) * length(input$metrics),
      "Every selected CMP must have every selected metric."
    ))
    normalize_preference(rows)
  })

  output$custom_weights <- renderUI({
    req(input$metrics)
    tagList(lapply(seq_along(input$metrics), function(i) {
      metric_name <- input$metrics[[i]]
      numericInput(
        paste0("custom_weight_", i), metric_name,
        value = 1, min = 0, step = 0.1
      )
    }))
  })

  metric_weights <- reactive({
    rows <- selected_rows()
    metrics <- input$metrics
    if (input$scheme == "equal") {
      weights <- rep(1, length(metrics))
    } else if (input$scheme == "balanced") {
      weights <- rep(0, length(metrics))
      stock <- metrics %in% stock_metrics
      fishing <- metrics %in% fishing_metrics
      if (any(stock)) weights[stock] <- 0.5 / sum(stock)
      if (any(fishing)) weights[fishing] <- 0.5 / sum(fishing)
      other <- !stock & !fishing
      if (any(other)) weights[other] <- 1 / sum(other)
    } else if (input$scheme == "dispersion") {
      cv <- vapply(metrics, function(metric_name) {
        values <- rows$raw_value[rows$metric == metric_name]
        stats::sd(values) / abs(mean(values))
      }, numeric(1))
      cv[!is.finite(cv)] <- 0
      weights <- sqrt(cv)
    } else {
      weights <- vapply(seq_along(metrics), function(i) {
        value <- input[[paste0("custom_weight_", i)]]
        if (is.null(value) || !is.finite(value) || value < 0) 0 else value
      }, numeric(1))
    }
    validate(need(sum(weights) > 0, "At least one metric weight must be positive."))
    data.frame(
      metric = metrics,
      scheme_weight = weights / sum(weights),
      stringsAsFactors = FALSE
    )
  })

  scores <- reactive({
    rows <- merge(selected_rows(), metric_weights(), by = "metric",
      all.x = TRUE, sort = FALSE)
    rows$contribution <- 100 * rows$preference * rows$scheme_weight
    totals <- aggregate(contribution ~ mp, rows, sum)
    names(totals)[2L] <- "score"
    totals <- totals[order(-totals$score, totals$mp), , drop = FALSE]
    totals$rank <- seq_len(nrow(totals))
    totals$mp <- factor(totals$mp, levels = totals$mp)
    totals
  })

  output$score_plot <- renderPlot({
    totals <- scores()
    ggplot(totals, aes(mp, score, fill = mp)) +
      geom_col(width = 0.72, show.legend = FALSE) +
      geom_text(aes(label = sprintf("%.1f", score)), vjust = -0.35,
        size = 4.4) +
      scale_fill_manual(values = cmp_colors, drop = FALSE) +
      scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0))) +
      labs(
        x = "CMP (ordered from highest score)",
        y = "Relative trade-off score (0-100)",
        subtitle = paste("Weighting:", input$scheme)
      ) +
      ggthemes::theme_few(base_size = 13) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
  })

  output$score_table <- renderDT({
    totals <- scores()
    out <- data.frame(
      Rank = totals$rank,
      CMP = as.character(totals$mp),
      Score = round(totals$score, 2),
      check.names = FALSE
    )
    datatable(out, rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  output$weight_table <- renderDT({
    weights <- metric_weights()
    weights$`Weight (%)` <- round(100 * weights$scheme_weight, 2)
    weights$scheme_weight <- NULL
    names(weights)[1L] <- "Metric"
    datatable(weights, rownames = FALSE,
      options = list(dom = "t", pageLength = 20))
  })

  output$quilt_plot <- renderPlot({
    rows <- selected_rows()
    rows$mp <- factor(rows$mp, levels = rev(input$cmps))
    rows$metric <- factor(rows$metric, levels = input$metrics)
    labels <- ifelse(
      rows$statistic == "C", sprintf("%.0f", rows$raw_value),
      ifelse(rows$statistic %in% c("SB0red", "PC270"),
        sprintf("%.1f%%", 100 * rows$raw_value),
        sprintf("%.2f", rows$raw_value))
    )
    ggplot(rows, aes(metric, mp, fill = preference)) +
      geom_tile(colour = "grey80", linewidth = 0.4) +
      geom_text(aes(label = labels), size = 4) +
      scale_fill_gradient(low = "#7952A8", high = "#F4EFF8",
        limits = c(0, 1), name = "Relative\npreference") +
      labs(x = NULL, y = NULL) +
      ggthemes::theme_few(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid = element_blank()
      )
  })
}

shinyApp(ui, server)
