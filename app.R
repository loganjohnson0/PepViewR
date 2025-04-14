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
all_total_peptides <- bind_rows(sarco_total_peptides, purge_total_peptides)
fasta <- nanoparquet::read_parquet(
  "data/2025-03-31_Sus_Scrofa_Total_fasta.parquet"
)

get_protein_data <- function(input, day, fraction, protein) {
  selected <- input |>
    dplyr::filter(
      stringr::str_detect(Protein.ID, protein),
      Fraction == fraction,
      stringr::str_detect(Time.Point, day)
    )

  if (nrow(selected) < 1) {
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


fraction_sidebar <- sidebar(
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
    inputId = "fraction_protein",
    label = "Protein Name, Gene, or UniProtID",
    choices = NULL,
    selected = NULL
  ),
  uiOutput("fraction_text"),
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


time_sidebar <- sidebar(
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
    inputId = "time_protein",
    label = "Protein Name, Gene, or UniProtID",
    choices = NULL,
    selected = NULL
  ),
  uiOutput("time_text"),
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

      card(
        card_header("Motivation"),
        "Bioinformatic tools help generate meaningful outputs of large datasets from genomic, transcriptomic, proteomic, metabolomic, and other types of omic-based methods.
        The meat science research community has also employed these tools with samples of early and extended postmortem samples.
        One could argue results from these tools are helpful for interpreting these results in meat science; however, many of these tools provide results and outputs that are uninformative and out of context, given the samples are in a postmortem context.
        This work aimed to illustrate an alternative approach to visualizing and contextualizing results in meat science in a way that can be informative and useful for researchers and a broader audience."
      ),
      card(
        card_header("Experimental Background"),
        "These proteomic data were generated from samples selected based on 1 day postmortem pH from a commercial pork plant (Jess et al., In Press).
        The workflow diagram broadly illustrates the protein extraction methods used for this experiment.",
        imageOutput("illustration", height = "auto"),
        max_height = 400
      ),
    ),
    # Protein Fraction Tab
    nav_panel(
      h5("Protein Fraction"),
      bslib::layout_sidebar(
        sidebar = fraction_sidebar,

        navset_card_underline(
          id = "fraction_tab",
          nav_panel(
            "Purge",
            card(
              plotly::plotlyOutput("frac_purge"),
              max_height = 800
            )
          ),
          nav_panel(
            "Sarcoplasmic",
            card(
              plotly::plotlyOutput("frac_sarco"),
              max_height = 800
            )
          )
        )
      )
    ),
    # Time Point Tab
    nav_panel(
      h5("Time Point"),
      bslib::layout_sidebar(
        sidebar = time_sidebar,
        navset_card_underline(
          nav_panel(
            "Day 1",
            card(
              plotly::plotlyOutput("time_01"),
              max_height = 800
            )
          ),
          nav_panel(
            "Day 7",
            card(
              plotly::plotlyOutput("time_07"),
              max_height = 800
            )
          ),
          nav_panel(
            "Day 14",
            card(
              plotly::plotlyOutput("time_14"),
              max_height = 800
            )
          )
        )
      )
    ),
    nav_spacer(),
    # About Tab
    nav_panel(
      class = ".box",
      h5("About"),

      card(
        card_header("Logan Johnson"),
        layout_columns(
          col_widths = c(3, 9),
          imageOutput("lgj", height = "auto"),
          "Logan Johnson"
        )
      ),
      layout_columns(
        card(
          card_header("Elisabeth Huff-Lonergan"),
          layout_columns(
            col_widths = c(3, 9),
            imageOutput("ehl", height = "auto"),
            "Elisabeth Huff-Lonergan"
          )
        ),
        card(
          card_header("Steven Lonergan"),
          layout_columns(
            col_widths = c(3, 9),
            imageOutput("sml", height = "auto"),
            "Steven Lonergan"
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  output$illustration <- renderImage(
    {
      list(
        src = "www/graphical_workflow.png",
        contentType = "image/png",
        height = 200,
        alt = "Protein sample preparation graphic"
      )
    },
    deleteFile = FALSE
  )

  output$lgj <- renderImage(
    {
      list(
        src = "www/logan_johnson.jpg",
        contentType = "image/jpg",
        height = 200,
        alt = "Logan Johnson photo"
      )
    },
    deleteFile = FALSE
  )
  output$ehl <- renderImage(
    {
      list(
        src = "www/elisabeth_lonergan.jpeg",
        contentType = "image/jpeg",
        height = 200,
        alt = "Elisabeth Huff-Lonergan photo"
      )
    },
    deleteFile = FALSE
  )
  output$sml <- renderImage(
    {
      list(
        src = "www/steven_lonergan.png",
        contentType = "image/png",
        height = 200,
        alt = "Steven Lonergan photo"
      )
    },
    deleteFile = FALSE
  )

  names <- all_total_peptides |>
    dplyr::select(c(Protein.Description, Protein.ID, Gene)) |>
    dplyr::distinct(Protein.ID, .keep_all = TRUE) |>
    dplyr::arrange(Gene, Protein.Description)

  shiny::updateSelectizeInput(
    session,
    inputId = "fraction_protein",
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

  shiny::updateSelectizeInput(
    session,
    inputId = "time_protein",
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

  fraction_fasta <- reactive({
    req(input$fraction_protein)
    selected_protein <- stringr::str_extract(
      input$fraction_protein,
      "(?<=: )\\w+$"
    )
    get_protein_fasta(fasta, protein = selected_protein)
  })

  time_fasta <- reactive({
    req(input$time_protein)
    selected_protein <- stringr::str_extract(input$time_protein, "(?<=: )\\w+$")
    get_protein_fasta(fasta, protein = selected_protein)
  })

  output$fraction_text <- renderUI({
    req(input$fraction_protein)
    if (stringr::str_detect(input$fraction_protein, "020931560")) {
      return(NULL)
    } else {
      uni_id <- stringr::str_extract(input$fraction_protein, "(?<=: )\\w+$")
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

  output$time_text <- renderUI({
    req(input$time_protein)
    if (stringr::str_detect(input$time_protein, "020931560")) {
      return(NULL)
    } else {
      uni_id <- stringr::str_extract(input$time_protein, "(?<=: )\\w+$")
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
    req(input$fraction_protein)
    selected_protein <- stringr::str_extract(
      input$fraction_protein,
      "(?<=: )\\w+$"
    )

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
      ggtitle(label = input$fraction_protein) +
      labs(x = "Amino Acid Residue Position (N- to C-Term)")

    # Optionally add FASTA annotation if available and short enough
    if (
      nrow(fraction_fasta()) > 0 &&
        length(fraction_fasta()$Protein.Position) < 2500
    ) {
      p <- p +
        geom_text(
          data = fraction_fasta(),
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
          data = fraction_fasta(),
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
    p <- p +
      geom_text(
        data = data,
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
    req(input$fraction_protein)
    selected_protein <- stringr::str_extract(
      input$fraction_protein,
      "(?<=: )\\w+$"
    )

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
      ggtitle(label = input$fraction_protein) +
      labs(x = "Amino Acid Residue Position (N- to C-Term)")

    # Optionally add FASTA annotation if available and short enough
    if (
      nrow(fraction_fasta()) > 0 &&
        length(fraction_fasta()$Protein.Position) < 2500
    ) {
      p <- p +
        geom_text(
          data = fraction_fasta(),
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
          data = fraction_fasta(),
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

    p <- p +
      geom_text(
        data = data,
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

  output$time_01 <- plotly::renderPlotly({
    req(input$time_protein)
    selected_protein <- stringr::str_extract(input$time_protein, "(?<=: )\\w+$")

    # Build data for Sarcoplasmic fraction
    purge <- get_protein_data(
      all_total_peptides,
      "01",
      "Purge",
      selected_protein
    )
    sarco <- get_protein_data(
      all_total_peptides,
      "01",
      "Sarco",
      selected_protein
    )

    data <- dplyr::bind_rows(purge, sarco)

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
            Fraction,
            .keep_all = TRUE
          ),
        aes(x = Peptide.Position, y = 0, xend = Peptide.Position, yend = 0.1)
      ) +
      geom_point(
        data = data |>
          dplyr::distinct(
            Peptide.Sequence,
            Protein.Start,
            Fraction,
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
          Fraction,
          levels = c("Sarco", "Purge"),
          labels = c(
            "Sarcoplasmic Extract",
            "Purge Extract"
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
      ggtitle(label = input$time_protein) +
      labs(x = "Amino Acid Residue Position (N- to C-Term)")

    # Optionally add FASTA annotation if available and short enough
    if (
      nrow(time_fasta()) > 0 &&
        length(time_fasta()$Protein.Position) < 2500
    ) {
      p <- p +
        geom_text(
          data = time_fasta(),
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
          data = time_fasta(),
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

    p <- p +
      geom_text(
        data = data,
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

  output$time_07 <- plotly::renderPlotly({
    req(input$time_protein)
    selected_protein <- stringr::str_extract(input$time_protein, "(?<=: )\\w+$")

    # Build data for Sarcoplasmic fraction
    purge <- get_protein_data(
      all_total_peptides,
      "07",
      "Purge",
      selected_protein
    )
    sarco <- get_protein_data(
      all_total_peptides,
      "07",
      "Sarco",
      selected_protein
    )

    data <- dplyr::bind_rows(purge, sarco)

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
            Fraction,
            .keep_all = TRUE
          ),
        aes(x = Peptide.Position, y = 0, xend = Peptide.Position, yend = 0.1)
      ) +
      geom_point(
        data = data |>
          dplyr::distinct(
            Peptide.Sequence,
            Protein.Start,
            Fraction,
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
          Fraction,
          levels = c("Sarco", "Purge"),
          labels = c(
            "Sarcoplasmic Extract",
            "Purge Extract"
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
      ggtitle(label = input$time_protein) +
      labs(x = "Amino Acid Residue Position (N- to C-Term)")

    # Optionally add FASTA annotation if available and short enough
    if (
      nrow(time_fasta()) > 0 &&
        length(time_fasta()$Protein.Position) < 2500
    ) {
      p <- p +
        geom_text(
          data = time_fasta(),
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
          data = time_fasta(),
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

    p <- p +
      geom_text(
        data = data,
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

  output$time_14 <- plotly::renderPlotly({
    req(input$time_protein)
    selected_protein <- stringr::str_extract(input$time_protein, "(?<=: )\\w+$")

    # Build data for Sarcoplasmic fraction
    purge <- get_protein_data(
      all_total_peptides,
      "14",
      "Purge",
      selected_protein
    )
    sarco <- get_protein_data(
      all_total_peptides,
      "14",
      "Sarco",
      selected_protein
    )

    data <- dplyr::bind_rows(purge, sarco)

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
            Fraction,
            .keep_all = TRUE
          ),
        aes(x = Peptide.Position, y = 0, xend = Peptide.Position, yend = 0.1)
      ) +
      geom_point(
        data = data |>
          dplyr::distinct(
            Peptide.Sequence,
            Protein.Start,
            Fraction,
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
          Fraction,
          levels = c("Sarco", "Purge"),
          labels = c(
            "Sarcoplasmic Extract",
            "Purge Extract"
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
      ggtitle(label = input$time_protein) +
      labs(x = "Amino Acid Residue Position (N- to C-Term)")

    # Optionally add FASTA annotation if available and short enough
    if (
      nrow(time_fasta()) > 0 &&
        length(time_fasta()$Protein.Position) < 2500
    ) {
      p <- p +
        geom_text(
          data = time_fasta(),
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
          data = time_fasta(),
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

    p <- p +
      geom_text(
        data = data,
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
