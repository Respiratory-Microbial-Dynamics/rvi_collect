library(shiny)
library(DT)
library(readr)
library(dplyr)

options(shiny.maxRequestSize = 100 * 1024^2)

read_input <- function(path, delimiter, required) {
  data <- read_delim(path, delim = delimiter, col_types = cols(.default = col_character()),
                     trim_ws = TRUE, show_col_types = FALSE)
  missing <- setdiff(required, names(data))
  validate(need(!length(missing), paste("Missing columns:", paste(missing, collapse = ", "))))
  data
}

copy_script <- function(paths, root) {
  # Expand the home directory on the machine running Bash.
  root_assignment <- if (startsWith(root, "~/")) {
    paste0('root="$HOME"/', shQuote(substring(root, 3), type = "sh"))
  } else {
    paste0("root=", shQuote(if (nzchar(root)) root else "/", type = "sh"))
  }
  c("#!/usr/bin/env bash",
    'dest="${1:-.}"',
    root_assignment,
    'mkdir -p -- "$dest" || exit 1',
    "failed=0",
    "copy_fasta() {",
    '  local source="$root/$1" target="$dest/${1##*/}" status',
    '  if [[ ! -f "$source" ]]; then',
    '    printf "Missing FASTA file: %s\\n" "$source" >&2',
    '    failed=1',
    '  elif [[ -e "$target" || -L "$target" ]]; then',
    '    cmp -s -- "$source" "$target"',
    '    status=$?',
    '    if [[ "$status" -ne 0 ]]; then',
    '      printf "Destination file already exists or cannot be read: %s\\n" "$target" >&2',
    '      failed=1',
    '    fi',
    '  elif ! cp -- "$source" "$target"; then',
    '    failed=1',
    '  fi',
    "}",
    paste("copy_fasta", shQuote(paths, type = "sh")),
    'exit "$failed"')
}

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { background: #f2f5f5; color: #153433; }
    h2 { margin: 25px 0 20px; font-weight: 700; }
    .well { background: #fff; border-color: #d9e3e1; }
    .results { background: #fff; border: 1px solid #d9e3e1; border-radius: 8px; padding: 16px; }
    #results th, #results td { white-space: nowrap; }
    .downloads { display: flex; flex-wrap: wrap; gap: 8px; margin: 12px 0 18px; }
    .downloads .btn { color: #fff; background: #176258; border-color: #176258; }
  "))),
  titlePanel("RVI | Metadata and FASTA"),
  sidebarLayout(
    sidebarPanel(
      fileInput("metadata", "CUH metadata (.tsv)", accept = ".tsv"),
      radioButtons("view", "View",
                   c("CUH metadata only" = "metadata",
                     "Metadata + viral classification" = "combined"),
                   selected = "metadata"),
      fileInput("classification", "Viral classification (.csv)", accept = ".csv"),
      checkboxInput("deduplicate", "Keep the second entry for duplicate Addenbrookes AnonSampleIDs", TRUE),
      selectizeInput("origin", "Sample origin",
                     choices = c("Addenbrookes", "Control", "Other"),
                     selected = character(), multiple = TRUE),
      selectizeInput("samples", "Samples (seq_id)", choices = NULL, multiple = TRUE),
      selectizeInput("columns", "Visible columns", choices = NULL, multiple = TRUE,
                     options = list(placeholder = "All columns")),
      conditionalPanel(
        "input.view === 'combined'",
        selectizeInput("virus", "Virus", choices = NULL, multiple = TRUE),
        sliderInput("coverage", "Genome coverage (%)", min = 0, max = 100,
                    value = c(0, 100), step = 0.1),
        textInput("root", "FASTA source directory on the cluster",
                  "~/scratch/00_consensus/viral_lens_1.3.4")
      ),
      width = 3
    ),
    mainPanel(
      div(class = "results",
        p("Table search uses the visible columns."),
        checkboxInput("manual", "Export selected rows only", FALSE),
        textOutput("summary"),
        div(class = "downloads",
            downloadButton("metadata_download", "Metadata (.tsv)"),
            conditionalPanel("input.view === 'combined' && output.has_fasta",
                             downloadButton("script_download", "Collect FASTA (.sh)"))),
        DTOutput("results")
      ),
      width = 9
    )
  )
)

server <- function(input, output, session) {
  metadata <- reactive({
    validate(need(input$metadata, "Upload CUH_metadata.tsv."))
    read_input(input$metadata$datapath, "\t", c("seq_id", "Sample origin", "AnonSampleID"))
  })

  classification <- reactive({
    validate(need(input$classification, "Upload the viral classification CSV."))
    data <- read_input(input$classification$datapath, ",",
                       c("Sample_ID", "Virus_Taxon_ID", "Virus", "Reference_Taxon_ID",
                         "Percentage_of_Genome_Covered"))
    data |>
      filter(is.na(Virus_Taxon_ID) | Virus_Taxon_ID != "Virus_Taxon_ID") |>
      mutate(
        Percentage_of_Genome_Covered = suppressWarnings(as.numeric(Percentage_of_Genome_Covered)),
        .has_fasta = grepl("^[0-9]+_[0-9]+_[0-9]+$", Sample_ID) &
          grepl("^[0-9]+$", Reference_Taxon_ID),
        fasta_relative_path = if_else(
          .has_fasta,
          paste0(sub("_[^_]+$", "", Sample_ID), "/", Sample_ID, "/", Reference_Taxon_ID,
                 "/", Sample_ID, ".", Reference_Taxon_ID, ".consensus.fa"),
          NA_character_
        )
      ) |>
      select(-.has_fasta)
  })

  metadata_selection <- reactive({
    meta <- metadata()
    excluded_ids <- character()
    if (isTRUE(input$deduplicate)) {
      drop <- meta |>
        mutate(.metadata_row = row_number()) |>
        filter(`Sample origin` == "Addenbrookes", !is.na(AnonSampleID), nzchar(AnonSampleID)) |>
        group_by(AnonSampleID) |>
        filter(n() > 1L, row_number() != 2L) |>
        pull(.metadata_row)
      keep <- setdiff(seq_len(nrow(meta)), drop)
      excluded_ids <- setdiff(meta$seq_id[drop], meta$seq_id[keep])
      meta <- meta[keep, ]
    }
    list(data = meta, excluded_ids = excluded_ids)
  })

  combined <- reactive({
    meta <- metadata_selection()
    calls <- filter(classification(), !Sample_ID %in% meta$excluded_ids)
    full_join(meta$data, calls, by = c("seq_id" = "Sample_ID"), na_matches = "never",
              relationship = "many-to-many") |>
      relocate(any_of(c("seq_id", "AnonSampleID", "Virus", "Reference_Subtype",
                        "Sample_Subtype", "Flu_Segment", "Percentage_of_Genome_Covered",
                        "Received", "Sample origin", "AgeGroup")))
  })

  view_data <- reactive({
    if (identical(input$view, "combined")) combined() else metadata_selection()$data
  })

  observe({
    data <- view_data()
    columns <- names(data)
    samples <- sort(unique(data$seq_id))
    updateSelectizeInput(session, "columns", choices = columns,
                         selected = intersect(isolate(input$columns), columns), server = TRUE)
    updateSelectizeInput(session, "samples", choices = samples,
                         selected = intersect(isolate(input$samples), samples), server = TRUE)
  })

  observe({
    req(input$classification)
    viruses <- sort(unique(classification()$Virus))
    updateSelectizeInput(session, "virus", choices = viruses,
                         selected = intersect(isolate(input$virus), viruses), server = TRUE)
  })

  filtered <- reactive({
    data <- view_data()
    if (length(input$origin)) {
      origin <- data$`Sample origin`
      known <- c("Addenbrookes", "Control")
      origin[!is.na(origin) & nzchar(origin) & !origin %in% known] <- "Other"
      data <- data[origin %in% input$origin, , drop = FALSE]
    }
    if (length(input$samples)) data <- filter(data, seq_id %in% input$samples)
    if (identical(input$view, "combined")) {
      if (length(input$virus)) data <- filter(data, Virus %in% input$virus)
      coverage <- input$coverage
      if (length(coverage) == 2L && any(coverage != c(0, 100))) {
        data <- filter(data, between(Percentage_of_Genome_Covered, coverage[1], coverage[2]))
      }
    }
    data
  })

  output$results <- renderDT({
    data <- filtered()
    columns <- intersect(input$columns, names(data))
    if (length(columns)) data <- select(data, all_of(columns))
    datatable(data, rownames = FALSE, selection = "multiple",
              options = list(scrollX = TRUE, pageLength = 15))
  }, server = TRUE)

  selected <- reactive({
    data <- filtered()
    # Wait for DataTables to report which rows match its search across all pages.
    req(!is.null(input$results_rows_all))
    rows <- input$results_rows_all
    if (isTRUE(input$manual)) rows <- intersect(rows, input$results_rows_selected)
    data[intersect(rows, seq_len(nrow(data))), , drop = FALSE]
  })

  fasta_paths <- reactive({
    if (!identical(input$view, "combined")) return(character())
    unique(na.omit(selected()$fasta_relative_path))
  })

  output$has_fasta <- reactive(length(fasta_paths()) > 0L)
  outputOptions(output, "has_fasta", suspendWhenHidden = FALSE)

  output$summary <- renderText({
    data <- selected()
    summary <- paste(n_distinct(data$seq_id, na.rm = TRUE), "samples ·", nrow(data), "rows")
    if (identical(input$view, "combined")) {
      summary <- paste(summary, "·", length(fasta_paths()), "FASTA paths")
    }
    summary
  })

  output$metadata_download <- downloadHandler(
    filename = function() "selected_metadata.tsv",
    content = function(file) write_tsv(selected(), file, na = "")
  )

  output$script_download <- downloadHandler(
    filename = function() "collect_fasta.sh",
    content = function(file) {
      root <- trimws(input$root)
      paths <- fasta_paths()
      req(length(paths), nzchar(root), grepl("^(/|~/)", root))
      writeLines(copy_script(paths, root), file, useBytes = TRUE)
    }
  )
}

shinyApp(ui, server)
