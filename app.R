library(bslib)
library(ggplot2)
library(dplyr)
library(nanoparquet)
library(plotly)
library(tidyr)
library(shiny)

# Read in your data
purge_total_peptides <- nanoparquet::read_parquet(
  "data/2025_03_31_Purge_Combined_Total_Peptides.parquet"
)
sarco_total_peptides <- nanoparquet::read_parquet(
  "data/2025_03_31_Sarco_Combined_Total_Peptides.parquet"
)
all_total_peptides <- rbind(sarco_total_peptides, purge_total_peptides)
fasta <- nanoparquet::read_parquet(
  "data/2025-03-31_Sus_Scrofa_Total_fasta.parquet"
)

names <- all_total_peptides |>
  dplyr::select(c(Protein.Description, Protein.ID, Gene)) |>
  dplyr::distinct(Protein.ID, .keep_all = TRUE) |>
  dplyr::arrange(Gene, Protein.Description)

# Functions to filter data and get FASTA
get_protein_data <- function(input, day, fraction, protein) {
  selected <- input |>
    dplyr::filter(
      stringr::str_detect(Protein.ID, protein),
      Fraction == fraction,
      stringr::str_detect(Time.Point, day)
    )

  if (nrow(selected) < 1) {
    # Return a dummy data.frame if nothing is found
    selected <- data.frame(
      Protein.Description = "",
      Protein.ID = protein,
      Peptide.Sequence = "",
      Protein.Start = 0,
      Protein.End = 1,
      Is.Unique = FALSE,
      Total.AA = as.numeric(NA),
      Peptide.Position = as.numeric(NA),
      Fraction = fraction,
      Time.Point = paste0("Day", day)
    )
    return(selected)
  } else {
    selected <- selected |>
      dplyr::arrange(Protein.Description, Protein.Start, Protein.End) |>
      dplyr::left_join(y = fasta, by = "Protein.ID") |>
      dplyr::rowwise() |>
      dplyr::reframe(
        Protein.Description = Protein.Description,
        Protein.ID = Protein.ID,
        Peptide.Sequence = Peptide.Sequence,
        Protein.Start = Protein.Start,
        Protein.End = Protein.End,
        Is.Unique = Is.Unique,
        Total.AA = Total.AA,
        Peptide.Position = seq(Protein.Start, Protein.End),
        Fraction = fraction
      ) |>
      dplyr::mutate(Time.Point = paste0("Day", day))
    return(selected)
  }
}

get_protein_fasta <- function(input, protein) {
  input |>
    dplyr::filter(stringr::str_detect(Protein.ID, protein)) |>
    dplyr::select(-Total.AA) |>
    dplyr::mutate(Sequence = stringr::str_split(Sequence, "")) |>
    tidyr::unnest(Sequence) |>
    dplyr::mutate(Protein.Position = dplyr::row_number())
}

# Shared sidebar remains unchanged
share_sidebar <- sidebar(
  width = 300,
  div(
    htmltools::HTML(
      "Start by selecting the protein you are interested in visualizing. You can scroll or type to search for a:
<ul>
<li>Protein's Name (Desmin)</li>
<li>Protein's Gene Product (DES)</li>
<li>Protein's UniProt ID (P02540)</li>
</ul>"
    ),
    style = "display: inline;"
  ),
  shiny::selectInput(
    inputId = "sel_protein",
    label = "Protein Name, Gene, or UniProtID",
    choices = NULL,
    selected = NULL
  ),
  div(
    htmltools::HTML(
      'The <b><span style="color: red;">Red Dots</span></b> indicate the start of a Unique tryptic peptide.'
    )
  ),
  div(
    htmltools::HTML(
      'The <b><span style="color: black;">Black Dots</span></b> indicate the start of a Non-Unique tryptic peptide.'
    )
  ),
  div(
    htmltools::HTML(
      "<h4>Note</h4> Not all proteins are present; if your protein is missing, it was not identified in this experiment."
    )
  )
)

ui <- bslib::page_fillable(
  titlePanel("PepViewR: Visualization of Proteomic Data"),
  theme = bslib::bs_theme(
    version = 5,
    bg = "#ffffff",
    fg = "#101010",
    primary = "#e30000"
  ),
  fillable_mobile = TRUE,
  tags$style(HTML(
    '.box {
        display: flex;
        align-items: center;
        justify-content: center;
        align-content: center;
      }'
  )),
  navset_card_tab(
    # Welcome Tab
    nav_panel(
      class = ".box",
      h5("Welcome"),
      card("Introduction and links can go here.")
    ),
    # Protein Fraction Tab
    nav_panel(
      h5("Protein Fraction"),
      bslib::layout_sidebar(
        sidebar = share_sidebar,
        # IMPORTANT: Add an id to capture the selected sub-tab
        navset_card_underline(
          id = "fraction_tab",
          nav_panel(
            "Purge",
            card(
              uiOutput("text"),
              plotly::plotlyOutput("frac_purge"),
              max_height = 800
            )
          ),
          nav_panel(
            "Sarcoplasmic",
            card(
              uiOutput("text"),
              plotly::plotlyOutput("frac_sarco"),
              max_height = 800
            )
          )
        )
      )
    ),
    # Time Point Tab (placeholder for now)
    nav_panel(
      h5("Time Point"),
      bslib::layout_sidebar(
        sidebar = share_sidebar,
        navset_card_underline(
          # You can later add an id here if you want to branch on time point as well
          nav_panel("Day 1", card("Filler", max_height = 800)),
          nav_panel("Day 7", card("Filler", max_height = 800)),
          nav_panel("Day 14", card("Filler", max_height = 800))
        )
      )
    ),
    nav_spacer(),
    # About Tab
    nav_panel(
      class = ".box",
      h5("About"),
      card("About the authors and the Lonergan Lab.")
    )
  )
)

server <- function(input, output, session) {
  # Update the protein selection input choices
  shiny::updateSelectizeInput(
    session,
    inputId = "sel_protein",
    choices = paste0(
      names$Protein.Description,
      " (",
      names$Gene,
      "): ",
      names$Protein.ID
    ),
    server = TRUE,
    selected = "",
    options = list(
      maxItems = 1,
      placeholder = "Start here....",
      closeAfterSelect = TRUE
    )
  )

  # Reactive for FASTA sequence data
  selected_fasta <- reactive({
    req(input$sel_protein)
    selected_protein <- stringr::str_extract(input$sel_protein, "(?<=: )\\w+$")
    get_protein_fasta(fasta, protein = selected_protein)
  })

  output$text <- renderUI({
    req(input$sel_protein)
    if (stringr::str_detect(input$sel_protein, "020931560")) {
      return(NULL)
    } else {
      uni_id <- stringr::str_extract(input$sel_protein, "(?<=: )\\w+$")
      tags$div(
        style = "display: inline-block",
        tags$a(
          href = paste0("https://www.uniprot.org/uniprotkb/", uni_id, "/entry"),
          paste0("Visit UniProt for ", uni_id),
          target = "_blank"
        )
      )
    }
  })

  output$frac_purge <- plotly::renderPlotly({
    req(input$sel_protein)
    selected_protein <- stringr::str_extract(input$sel_protein, "(?<=: )\\w+$")

    day01 <- get_protein_data(
      all_total_peptides,
      "01",
      "Purge",
      selected_protein
    )
    day07 <- get_protein_data(
      all_total_peptides,
      "07",
      "Purge",
      selected_protein
    )
    day14 <- get_protein_data(
      all_total_peptides,
      "14",
      "Purge",
      selected_protein
    )
    data <- dplyr::bind_rows(day01, day07, day14)

    p <- ggplot(data, aes(x = Peptide.Position, y = 0)) +
      geom_point(aes(
        text = paste0("<b>Identified Peptide</b>: ", Peptide.Sequence)
      )) +
      scale_x_continuous(limits = c(0, max(data$Total.AA, na.rm = TRUE))) +
      scale_y_continuous(limits = c(-0.5, 1)) +
      geom_segment(
        data = data |>
          dplyr::distinct(
            Peptide.Sequence,
            Protein.Start,
            Time.Point,
            .keep_all = TRUE
          ),
        aes(x = Peptide.Position, y = 0, xend = Peptide.Position, yend = 0.1)
      ) +
      geom_point(
        data = data |>
          dplyr::distinct(
            Peptide.Sequence,
            Protein.Start,
            Time.Point,
            .keep_all = TRUE
          ),
        aes(
          x = Peptide.Position,
          y = 0.1,
          color = if_else(Is.Unique == "TRUE", "#ff0000", "#000000")
        )
      ) +
      scale_color_identity() +
      facet_wrap(
        ~ factor(
          Time.Point,
          levels = c("Day01", "Day07", "Day14"),
          labels = c(
            "Day 1 Postmortem",
            "Day 7 Postmortem",
            "Day 14 Postmortem"
          )
        ),
        ncol = 1
      ) +
      theme(
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour = "#000000", linewidth = 1)
      ) +
      ggtitle(label = input$sel_protein) +
      labs(x = "Amino Acid Residue Position (N- to C-Term)")

    # Optionally add FASTA annotation if available and short enough
    if (
      nrow(selected_fasta()) > 0 &&
        length(selected_fasta()$Protein.Position) < 2500
    ) {
      p <- p +
        geom_text(
          data = selected_fasta(),
          aes(
            x = Protein.Position,
            y = -0.5,
            label = Protein.Position,
            text = paste0("<b>Protein Position</b>: ", Protein.Position)
          ),
          family = "Courier",
          size = 3
        ) +
        geom_text(
          data = selected_fasta(),
          aes(
            x = Protein.Position,
            y = -0.3,
            label = Sequence,
            text = paste0("<b>Residue</b>: ", Sequence)
          ),
          family = "Courier",
          size = 3
        )
    }

    plotly::ggplotly(p, tooltip = c("text")) |>
      plotly::layout(
        hovermode = "x unified",
        xaxis = list(
          showspikes = TRUE,
          spikemode = "across",
          spikesnap = "cursor",
          spikethickness = 2,
          spikedash = "solid"
        )
      )
  })

  output$frac_sarco <- plotly::renderPlotly({
    req(input$sel_protein)
    selected_protein <- stringr::str_extract(input$sel_protein, "(?<=: )\\w+$")

    # Build data for Sarcoplasmic fraction
    day01 <- get_protein_data(
      all_total_peptides,
      "01",
      "Sarco",
      selected_protein
    )
    day07 <- get_protein_data(
      all_total_peptides,
      "07",
      "Sarco",
      selected_protein
    )
    day14 <- get_protein_data(
      all_total_peptides,
      "14",
      "Sarco",
      selected_protein
    )
    data <- dplyr::bind_rows(day01, day07, day14)

    p <- ggplot(data, aes(x = Peptide.Position, y = 0)) +
      geom_point(aes(
        text = paste0("<b>Identified Peptide</b>: ", Peptide.Sequence)
      )) +
      scale_x_continuous(limits = c(0, max(data$Total.AA, na.rm = TRUE))) +
      scale_y_continuous(limits = c(-0.5, 1)) +
      geom_segment(
        data = data |>
          dplyr::distinct(
            Peptide.Sequence,
            Protein.Start,
            Time.Point,
            .keep_all = TRUE
          ),
        aes(x = Peptide.Position, y = 0, xend = Peptide.Position, yend = 0.1)
      ) +
      geom_point(
        data = data |>
          dplyr::distinct(
            Peptide.Sequence,
            Protein.Start,
            Time.Point,
            .keep_all = TRUE
          ),
        aes(
          x = Peptide.Position,
          y = 0.1,
          color = if_else(Is.Unique == "TRUE", "#ff0000", "#000000")
        )
      ) +
      scale_color_identity() +
      facet_wrap(
        ~ factor(
          Time.Point,
          levels = c("Day01", "Day07", "Day14"),
          labels = c(
            "Day 1 Postmortem",
            "Day 7 Postmortem",
            "Day 14 Postmortem"
          )
        ),
        ncol = 1
      ) +
      theme(
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_rect(colour = "#000000", linewidth = 1)
      ) +
      ggtitle(label = input$sel_protein) +
      labs(x = "Amino Acid Residue Position (N- to C-Term)")

    # Optionally add FASTA annotation if available and short enough
    if (
      nrow(selected_fasta()) > 0 &&
        length(selected_fasta()$Protein.Position) < 2500
    ) {
      p <- p +
        geom_text(
          data = selected_fasta(),
          aes(
            x = Protein.Position,
            y = -0.5,
            label = Protein.Position,
            text = paste0("<b>Protein Position</b>: ", Protein.Position)
          ),
          family = "Courier",
          size = 3
        ) +
        geom_text(
          data = selected_fasta(),
          aes(
            x = Protein.Position,
            y = -0.3,
            label = Sequence,
            text = paste0("<b>Residue</b>: ", Sequence)
          ),
          family = "Courier",
          size = 3
        )
    }

    plotly::ggplotly(p, tooltip = c("text")) |>
      plotly::layout(
        hovermode = "x unified",
        xaxis = list(
          showspikes = TRUE,
          spikemode = "across",
          spikesnap = "cursor",
          spikethickness = 2,
          spikedash = "solid"
        )
      )
  })
}

shiny::shinyApp(ui = ui, server = server)
