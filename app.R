library(bslib)
library(ggplot2)
library(dplyr)
library(nanoparquet)
library(plotly)
library(tidyr)
library(shiny)

# Read in the data
purge_total_peptides <- nanoparquet::read_parquet(
  "data/2025_03_31_Purge_Combined_Total_Peptides.parquet"
)
sarco_total_peptides <- nanoparquet::read_parquet(
  "data/2025_03_31_Sarco_Combined_Total_Peptides.parquet"
)
all_total_peptides <- dplyr::bind_rows(
  sarco_total_peptides,
  purge_total_peptides
)
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


fraction_sidebar <- bslib::sidebar(
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


time_sidebar <- bslib::sidebar(
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
  titlePanel(
    "PepViewR: Tryptic Peptide Visualization from Proteomic Experiments"
  ),
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
        card_body(
          p(
            "Bioinformatic tools help generate meaningful outputs of large datasets from genomic, transcriptomic, proteomic, metabolomic, and other types of omic-based methods.
            These tools are almost exclusively trained and developed on data from living cellular systems.
            The meat science research community has employed these tools to draw conclusions from postmortem skeletal muscle.
            One could argue that the outputs from bioinformatic tools are helpful for interpreting these results in meat science; however, many of these tools provide results and outputs that are uninformative, given the experimental samples were in a postmortem context.
            This project and Shiny application aimed to illustrate an alternative approach that visualizes and contextualizes results from a meat science research experiment.
            This is ultimately an effort to demonstrate the need for creative alternatives to out-of-the-box bioinformatic tools to draw conclusions from samples in a postmortem context."
          )
        )
      ),
      layout_columns(
        card(
          card_header("Experimental Background"),
          card_body(
            p(
              "Pork loin samples (n = 50) were selected at 1 d postmortem from a commercial plant to be divertent in pH. A",
              tags$a(
                href = "https://www.iastatedigitalpress.com/mmb/article/id/18426/",
                "recent paper (Jess et al., 2025)"
              ),
              "describes the collection and differences in meat quality and soluble desmin degradation products through Western blots.
              A workflow diagram (right) broadly illustrates the protein extraction methods used for this experiment. 
              More details are available in the",
              tags$a(href = "", "abstract"),
              "submitted to the Reciprocal Meat Conference.
              Briefly, 2 separate proteomic experiments utilized the same population of pork loins at 1, 7, and 14 days postmortem. 
              The two experiments utilized either:"
            ),
            tags$ol(
              tags$li(
                "Proteins from the muscle tissue soluble in a low-ionic strength buffer, and"
              ),
              tags$li(
                "Proteins from the muscle exudate (purge)"
              )
            ),
            p(
              "Proteins were digest into tryptic peptides, labeled with Tandem Mass Tag 10plex reagents, and analyzed using a Q-Exactive Mass Spectrometer (Thermo Scientific), similar to prior work",
              tags$a(
                href = "https://academic.oup.com/jas/article/doi/10.1093/jas/skae355/7905115",
                "Johnson et al., 2024)."
              ),
              "The .raw files were converted to .mzML and analyzed with FragPipe (v. 22.0) with a Sus scrofa UniProt reference .fasta file.
              Data were summarized with MSstatsTMT (v. 2.14.2), and this Shiny app for viewing and interacting with the identified peptide data was developed in R Statistical Software (R Core Team, 2024)."
            )
          )
        ),
        card(
          card_header("Workflow Diagram"),
          card_body(
            imageOutput("illustration", height = "auto"),
            max_height = 400
          )
        ),
      )
    ),
    # Protein Fraction Tab
    nav_panel(
      h5("Protein Fraction"),
      bslib::layout_sidebar(
        sidebar = fraction_sidebar,

        navset_card_underline(
          id = "fraction_tab",
          nav_panel(
            "Muscle Exudate Extract",
            card(
              plotly::plotlyOutput("frac_purge"),
              max_height = 800
            )
          ),
          nav_panel(
            "Soluble Muscle Extract",
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

      layout_columns(
        card(
          card_header("Logan Johnson"),
          card_body(
            div(
              class = "row",
              div(
                class = "col-md-6",
                imageOutput("lgj", height = "auto")
              ),
              div(
                class = "col-md-6",
                p(
                  "Logan is currently a Postdoctoral Research Associate at Iowa State University in the Lonergan Lab.
                In July 2025, Logan will begin as an Assistant Professor of Food Science at Oklahoma State University in the",
                  tags$a(
                    href = "https://agriculture.okstate.edu/departments-programs/afs/",
                    "Department of Animal and Food Sciences."
                  ),
                )
              )
            )
          )
        ),
        card(
          card_header("Elisabeth Huff-Lonergan"),
          card_body(
            div(
              class = "row",
              div(
                class = "col-md-6",
                imageOutput("ehl", height = "auto")
              ),
              div(
                class = "col-md-6",
                p(
                  tags$a(
                    href = "https://www.ans.iastate.edu/people/elisabeth-huff-lonergan",
                    "Elisabeth"
                  ),
                  "is currently a University Professor at Iowa State University in the Department of Animal Science. \n
                  She serves as the Editor-In-Chief for the",
                  tags$a(
                    href = "https://academic.oup.com/jas",
                    "Journal of Animal Science"
                  ),
                  "and was previously the Editor-In-Chief for",
                  tags$a(
                    href = "https://www.iastatedigitalpress.com/mmb/",
                    "Meat and Muscle Biology."
                  )
                ),
              )
            )
          )
        )
      ),
      layout_columns(
        card(
          card_header("Steven Lonergan"),
          card_body(
            div(
              class = "row",
              div(
                class = "col-md-6",
                imageOutput("sml", height = "auto")
              ),
              div(
                class = "col-md-6",
                p(
                  tags$a(
                    href = "https://www.ans.iastate.edu/people/steven-lonergan",
                    "Steven"
                  ),
                  "is currently a Morrill Professor at Iowa State University in the Department of Animal Science.
                  He was recently elected to serve as the 2025-2026 President of the Board of Directors for the American Meat Science Association!
                  He serves as a Section Editor for the",
                  tags$a(
                    href = "https://academic.oup.com/jas",
                    "Journal of Animal Science."
                  )
                ),
              )
            )
          )
        ),
        card(
          card_header("Funding Acknowledgement"),
          card_body(
            p(
              "This work has been funded (in part) thanks to the following grants:"
            ),
            tags$ul(
              tags$li(
                "the United States Department of Agriculture–Agriculture and Food Research Initiative project 2019-67017- 29181"
              ),
              tags$li(
                "the National Science Foundation under Grant No. DGE-1828942, and"
              ),
              tags$li("the Iowa Pork Producers Association")
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  output$illustration <- shiny::renderImage(
    {
      list(
        src = "www/graphical_workflow.png",
        contentType = "image/png",
        height = 300,
        alt = "Protein sample preparation graphic"
      )
    },
    deleteFile = FALSE
  )

  output$lgj <- shiny::renderImage(
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
  output$ehl <- shiny::renderImage(
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
  output$sml <- shiny::renderImage(
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

  fraction_fasta <- shiny::reactive({
    req(input$fraction_protein)
    selected_protein <- stringr::str_extract(
      input$fraction_protein,
      "(?<=: )\\w+$"
    )
    get_protein_fasta(fasta, protein = selected_protein)
  })

  time_fasta <- shiny::reactive({
    req(input$time_protein)
    selected_protein <- stringr::str_extract(input$time_protein, "(?<=: )\\w+$")
    get_protein_fasta(fasta, protein = selected_protein)
  })

  output$fraction_text <- shiny::renderUI({
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

  output$time_text <- shiny::renderUI({
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

    # Build data for Low-Ionic Strength fraction
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

    # Build data for Low-Ionic Strength fraction
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
            "Low-Ionic Strength Soluble Protein Extract",
            "Muscle Exudate Protein Extract"
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

    # Build data for Low-Ionic Strength fraction
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
            "Low-Ionic Strength Soluble Protein Extract",
            "Muscle Exudate Protein Extract"
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

    # Build data for Low-Ionic Strength Soluble fraction
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
            "Low-Ionic Strength Soluble Protein Extract",
            "Muscle Exudate Protein Extract"
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
