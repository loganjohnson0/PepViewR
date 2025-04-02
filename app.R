library(bslib)
library(ggplot2)
library(dplyr)
library(nanoparquet)
library(plotly)
library(tidyr)
library(shiny)


purge_total_peptides <- nanoparquet::read_parquet(
  file = "data/2025_03_31_Purge_Combined_Total_Peptides.parquet"
)
sarco_total_peptides <- nanoparquet::read_parquet(
  file = "data/2025_03_31_Sarco_Combined_Total_Peptides.parquet"
)

all_total_peptides <- rbind(sarco_total_peptides, purge_total_peptides)

fasta <- nanoparquet::read_parquet(
  file = "data/2025-03-31_Sus_Scrofa_Total_fasta.parquet"
)

names <- all_total_peptides |>
  dplyr::select(c(Protein.Description, Protein.ID, Gene)) |>
  dplyr::distinct(Protein.ID, .keep_all = TRUE) |>
  dplyr::arrange(Gene, Protein.Description)

get_protein_data <- function(input, day, fraction, protein) {
  selected <- input |>
    dplyr::filter(
      stringr::str_detect(Protein.ID, protein),
      Fraction == fraction,
      stringr::str_detect(Time.Point, day)
    )

  if (length(selected$Protein.ID) < 1) {
    selected <- data.frame(
      Protein.Description = "",
      Protein.ID = protein,
      Peptide.Sequence = "",
      Protein.Start = 0,
      Protein.End = 1,
      Is.Unique = is.logical("FALSE"),
      Total.AA = as.numeric(""),
      Peptide.Position = as.numeric(""),
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
  protein_fasta <- input |>
    dplyr::filter(stringr::str_detect(Protein.ID, protein)) |>
    dplyr::select(-Total.AA) |>
    dplyr::mutate(Sequence = stringr::str_split(Sequence, "")) |>
    tidyr::unnest(Sequence) |>
    dplyr::mutate(Protein.Position = dplyr::row_number())
}

share_sidebar <- sidebar(
  width = 300,
  div(
    htmltools::HTML(
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
    label = "Protein Name, Gene, or UniProtID",
    choices = NULL,
    selected = NULL
  ),
  div(
    htmltools::HTML(
      'The <b><span style="color: red;">Red Dots</span></b> <i>above</i> the peptide indicate the start of a Unique tryptic peptide.'
    )
  ),
  div(
    htmltools::HTML(
      'The <b><span style="color: black;">Black Dots</span></b> <i> above</i> the peptide indicate the start of a Non-Unique tryptic peptide.'
    )
  ),
  div(
    htmltools::HTML(
      "<h4>Note</h4> Not all possible proteins are present, thus if your protein of interest is not in the list, it was not identified in this experiment."
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

  # fmt: skip
  tags$style(HTML(
    '.box {
  display: flex;
  align-items: center;
  justify-content: center;
  align-content: center;
  }')),

  navset_card_tab(
    # Start of Welcome
    nav_panel(
      class = ".box",
      h5("Welcome"),
      card(
        "My idea was to have some sort of an introduction here.
      
      Maybe this would include some links to review papers.
      
      I was also thinking of including some images of either sarco vs. purge in how the experiment was conducted."
      )
    ),
    #####
    #
    # Start of Protein Fraction
    #
    #####
    nav_panel(
      h5("Protein Fraction"),
      class = ".box",
      bslib::layout_sidebar(
        sidebar = share_sidebar,
        navset_card_underline(
          # Start of Protein Fraction - Purge
          nav_panel(
            "Purge",
            card(
              uiOutput("text"),
              plotly::plotlyOutput("frac_purge"),
              max_height = 800
            )
          ),

          # Start of Protein Fraction - Sarcoplasmic
          nav_panel(
            "Sarcoplasmic",
            card("Space filler text for now until later.", max_height = 800)
          )
        )
      )
    ),
    #####
    #
    # Start of Time Point
    #
    #####
    nav_panel(
      h5("Time Point"),
      class = ".box",
      bslib::layout_sidebar(
        sidebar = share_sidebar,
        navset_card_underline(
          # Start of Time Point - Day 1
          nav_panel(
            "Day 1",
            card("Filler", max_height = 800)
          ),

          # Start of Time Point - Day 7
          nav_panel(
            "Day 7",
            card("Space filler text for now until later.", max_height = 800)
          ),

          # Start of Time Point - Day 14
          nav_panel(
            "Day 14",
            card("Space filler text for now until later.", max_height = 800)
          )
        )
      )
    ),
    nav_spacer(),

    #####
    #
    # Start of About
    #
    #####
    nav_panel(
      class = ".box",
      h5("About"),
      card(
        "Learn More Here.
      
      Was thinking this could be an 'about the authors' type of page and some info about the Lonergan Lab, etc."
      )
    )
  )
)


server <- function(input, output, session) {
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

  selected_fasta <- shiny::reactive({
    req(input$sel_protein)

    selected_protein <- stringr::str_extract(
      input$sel_protein,
      pattern = "(?<=: )\\w+$"
    )
    get_protein_fasta(fasta, protein = selected_protein)
  })

  combine_day <- shiny::reactive({
    req(input$sel_protein)

    selected_protein <- stringr::str_extract(
      input$sel_protein,
      pattern = "(?<=: )\\w+$"
    )

    day01 <- get_protein_data(
      input = all_total_peptides,
      day = "01",
      fraction = "Purge",
      protein = selected_protein
    )
    day07 <- get_protein_data(
      input = all_total_peptides,
      day = "07",
      fraction = "Purge",
      protein = selected_protein
    )
    day14 <- get_protein_data(
      input = all_total_peptides,
      day = "14",
      fraction = "Purge",
      protein = selected_protein
    )

    bind_rows(day01, day07, day14)
  })

  output$text <- renderUI({
    req(input$sel_protein)

    if (stringr::str_detect(input$sel_protein, "020931560")) {
      return(NULL)
    } else {
      uni_id <- stringr::str_extract(
        input$sel_protein,
        pattern = "(?<=: )\\w+$"
      )

      tags$div(
        sytle = "display: inline-block",
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

    p <- ggplot(data = combine_day(), aes(x = Peptide.Position, y = 0)) +

      geom_point(aes(
        text = paste0("<b>Identified Peptide</b>: ", Peptide.Sequence)
      )) +

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
        aes(
          x = Peptide.Position,
          y = 0,
          xend = Peptide.Position,
          yend = 0.1
        )
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
        ncol = 1,
        axes = "all",
        axis.labels = "all_x"
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

    if (length(selected_fasta()$Protein.Position) < 2500) {
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

      plotly::ggplotly(p, tooltip = c("text")) |>
        plotly::style(hoverinfo = "none", traces = c(4:9)) |>
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
    } else {
      p <- p +
        geom_text(
          data = combine_day(),
          aes(
            x = Peptide.Position,
            y = -0.5,
            label = Peptide.Position,
            text = paste0("<b>Protein Position</b>: ", Peptide.Position)
          ),
          family = "Courier",
          size = 3
        )

      plotly::ggplotly(p, tooltip = c("text")) |>
        plotly::style(hoverinfo = "none", traces = c(4:12)) |>
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
    }
  })
}

shiny::shinyApp(ui = ui, server = server)
