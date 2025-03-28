library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(readr)
library(plotly)
library(tidyr)
library(markdown)
library(stringr)

purge_total_peptides <- readr::read_csv(
  file = "data/2025_03_28_Combined_Total_Peptides.csv"
)

purge_day01_total_peptides <- purge_total_peptides |>
  filter(Time.Point == "Day01")

purge_day07_total_peptides <- purge_total_peptides |>
  filter(Time.Point == "Day07")

purge_day14_total_peptides <- purge_total_peptides |>
  filter(Time.Point == "Day14")

fasta <- readr::read_csv(
  file = "data/2025-03-10_sus_scrofa_total_AA.csv"
)

get_protein_data <- function(input, day, protein) {
  selected <- input |>
    filter(str_detect(Protein.ID, protein))

  if (nrow(selected) == 0) {
    selected <- tibble(
      Protein.Description = "",
      Protein.ID = selected_protein,
      Peptide.Sequence = "",
      Protein.Start = 0,
      Protein.End = 1,
      Is.Unique = is.logical("FALSE"),
      Total.AA = 0,
      AA.Position = 0,
      Time.Point = paste0("Day", day)
    )
    return(selected)
  } else {
    selected <- selected |>
      arrange(Protein.Description, Protein.Start, Protein.End) |>
      left_join(y = fasta, by = "Protein.ID") |>
      rowwise() |>
      reframe(
        Protein.Description = Protein.Description,
        Protein.ID = Protein.ID,
        Peptide.Sequence = Peptide.Sequence,
        Protein.Start = Protein.Start,
        Protein.End = Protein.End,
        Is.Unique = Is.Unique,
        Total.AA = Total.AA,
        Peptide.Position = seq(Protein.Start, Protein.End)
      ) |>
      mutate(Time.Point = paste0("Day", day))

    return(selected)
  }
}

get_protein_fasta <- function(input, protein) {
  protein_fasta <- input |>
    filter(Protein.ID == protein) |>
    select(-Total.AA) |>
    mutate(Sequence = str_split(Sequence, "")) |>
    unnest(Sequence) |>
    mutate(Protein.Position = row_number())
}


ui <- bslib::page_fillable(
  h1("Pork Proteomic Results"),
  theme = bslib::bs_theme(bg = "#FFF", fg = "#101010"),
  fillable_mobile = TRUE,

  layout_sidebar(
    border = TRUE,
    border_color = "#000000",
    bg = "#e1e1e1",
    fg = "#101010",
    div(
      HTML(
        "These data illustrate the identification of specific tryptic peptides in the <b>Purge</b> or <b>Muscle Exudate</b> of fresh pork loins aged for either 1, 7, or 14 days postmortem."
      ),
      style = "display: inline;"
    ),
    sidebar = sidebar(
      width = 350,
      tags$h2("Welcome!"),
      div(
        HTML(
          "Start by selecting the protein you are interested in visualizing. You can scroll or type to search for a:
        <ul>
  <li>Protein's Name (Desmin)</li>
  <li>Protein's Gene Product (DES)</li> or
  <li>Protein's UniProt ID (P02540)</li>
</ul>"
        ),
        style = "display: inline;"
      ),
      shiny::selectInput(
        inputId = "sel_protein",
        label = "Protein Name, Gene, or ID",
        choices = NULL,
        selected = NULL
      ),
      div(
        HTML(
          'The <b><span style="color: red;">Red Dots</span></b> <i>above</i> the peptide indicate that is a Unique peptide.'
        )
      ),
      div(
        HTML(
          'The <b><span style="color: black;">Black Dots</span></b> <i> above</i> the peptide indicate that is Not a Unique peptide.'
        )
      ),
      div(
        HTML(
          "<h4>Note</h4> Not all possible proteins are present, thus if your protein of interest is not in the list, it was not identified in this experiment."
        )
      )
    ),
    card(
      plotly::plotlyOutput("main_plot"),
      max_height = 800
    )
  ),
)


server <- function(input, output, session) {
  names <- purge_total_peptides |>
    select(c(Protein.Description, Protein.ID, Gene)) |>
    distinct(Protein.ID, .keep_all = TRUE) |>
    arrange(Protein.Description)

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

  combine_day <- reactive({
    req(input$sel_protein)

    selected_protein <- input$sel_protein |>
      stringr::str_extract(pattern = "(?<=: )\\w+$")

    selected_fasta <- get_protein_fasta(
      input = fasta,
      protein = selected_protein
    )

    day01 <- get_protein_data(
      input = purge_day01_total_peptides,
      day = "01",
      protein = selected_protein
    )
    day07 <- get_protein_data(
      input = purge_day07_total_peptides,
      day = "07",
      protein = selected_protein
    )
    day14 <- get_protein_data(
      input = purge_day14_total_peptides,
      day = "14",
      protein = selected_protein
    )

    bind_rows(day01, day07, day14)
  })

  selected_fasta <- reactive({
    req(input$sel_protein)

    selected_protein <- input$sel_protein |>
      stringr::str_extract(pattern = "(?<=: )\\w+$")

    get_protein_fasta(fasta, protein = selected_protein)
  })

  output$main_plot <- plotly::renderPlotly({
    req(input$sel_protein)
    p <- ggplot(data = combine_day(), aes(x = Peptide.Position, y = 0)) +

      geom_point() +

      scale_x_continuous(limits = c(0, max(combine_day()$Total.AA))) +
      scale_y_continuous(limits = c(-0.5, 1)) +

      geom_segment(
        data = combine_day() |>
          distinct(
            Peptide.Sequence,
            Protein.Start,
            Time.Point,
            .keep_all = TRUE
          ),
        aes(x = Peptide.Position, y = 0, xend = Peptide.Position, yend = 0.1)
      ) +

      geom_point(
        data = combine_day() |>
          distinct(
            Peptide.Sequence,
            Protein.Start,
            Time.Point,
            .keep_all = TRUE
          ),
        aes(
          x = Peptide.Position,
          y = 0.1,
          text = paste0("<b>Identified Peptide</b>: ", Peptide.Sequence),
          color = if_else(Is.Unique == "TRUE", "#ff0000", "#000000")
        )
      ) +

      scale_color_identity() +

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
      ) +
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

      theme(
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "#aaaaaa"),
        strip.placement = "outside"
      ) +
      facet_wrap(~Time.Point, ncol = 1, axes = "all", axis.labels = "all_x") +
      ggtitle(
        label = input$sel_protein
      ) +

      labs(x = "Amino Acid Residue Position (C- to N-Term)")

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
