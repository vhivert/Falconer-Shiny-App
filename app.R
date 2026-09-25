###################################################################################
#
# FALCONER SHINY APP
#
# Description : The aim of this App is to show how the combination of gene action
# and allele frequencies at causal loci translate to genetic variance and genetic
# variance components for a complex trait. Although the theory underlying the App
# is more than a century old, it is highly relevant in the current era of
# genome-wide association studies (GWAS). It highlights the specific definition of
# the effect size estimates by GWAS and the variation it generates in the population,
# i.e. how locus-specific effects lead to individual differences. In addition, it can
# also be used to demonstrate how within and between locus interactions (dominance
# and epistasis, respectively) usually do not lead to a large amount of non-additive
# variance relative to additive variance, and therefore that these interactions
# usually do not explain individual differences in a population.
#
# The three models described mainly illustrate the Chapters 7 and 8 of Falconer
# and Mackay (1996) and the Chapter 5 of Lynch and Walsh (1998).
#
# Authors: Valentin Hivert, Naomi Wray and Peter Visscher
# Date: 03 Feb 2021 (v1.0), revised Sep 2026 (v1.1)
#
# Version 1.1
#
# Files:
#   App.R           user interface and server logic (this file)
#   genetics.R    model computations (population mean, average effects, variance components)
#   readme.R      content of the README dialog
#   www/            model figures (AD_model.png, AAA_model.png)
#   
# Citation:
# Hivert V, Wray NR, Visscher PM (2021) Gene action, genetic variation, and GWAS: A user-friendly web tool. PLOS Genetics 17(5): e1009548. https://doi.org/10.1371/journal.pgen.1009548
###################################################################################

## Libraries
library(shiny)
library(plot3D)
library(htmlTable)
library(pBrackets)
library(rhandsontable)

## Model computations and README content
source(file.path("./genetics.R"), local = TRUE)
source(file.path("./readme.R"),   local = TRUE)


### Display helpers ###############################################################

fmt <- function(x) format(x, digits = 2, nsmall = 2)

## One line of the variance summary, e.g. "Additive variance (VA) : 1.23 (locus A = ..., locus B = ...)"
var_line <- function(label, total, loci = NULL) {
  txt <- paste0(label, " : ", fmt(total))
  if (!is.null(loci)) {
    txt <- paste0(txt, " (locus A = ", fmt(loci[1]), ", locus B = ", fmt(loci[2]), ")")
  }
  HTML(txt)
}

REF_FALCONER <- "Falconer, D.S., and Mackay T.F.C. (1996). Introduction to quantitative Genetics, Ed. 4th. Longmans Green, Harlow, Essex."
REF_MAKI     <- "M&auml;ki-Tanila A., Hill W.G. (2014). Influence of gene interaction on complex trait variation with multi-locus models. <i>Genetics</i>, 198(1):355-367."
REF_LYNCH    <- "Lynch, M. and Walsh, B. (1998). Genetics and Analysis of Quantitative Traits. Sinauer Associates."

LAB_P <- bquote("Frequency of " ~ A[1] ~ (italic(p)))
LAB_Q <- bquote("Frequency of " ~ B[1] ~ (italic(q)))
GENO_AXIS_LABELS <- c(expression(A[2] * A[2] ~ (0)), expression(A[1] * A[2] ~ (1)), expression(A[1] * A[1] ~ (2)))

###################################################################################

ui <- fluidPage(
  tags$head(tags$style(HTML("
    .info-box     { text-align: left; background-color: #f2f2f2; font-size: 110%; }
    .model-figure { text-align: center; background-color: #f2f2f2; }
  "))),

  titlePanel(title = "The Falconer ShinyApp"),
  HTML("<b>Choose a model to display:</b>"),
  navbarPage(tags$b("Model"), id = "model",
             tabPanel("Single-locus Additive and Dominance", value = "AD"),
             tabPanel("Two-locus Additive and Additive-by-Additive", value = "AA"),
             tabPanel("General two-locus model", value = "Perso")
  ),
  conditionalPanel(
    condition = "input.model == 'AD'",
    div(class = "info-box",
        HTML("This model uses the notation of Falconer and Mackay (1996) Chapters 7 & 8, and uses their equations to generate the population specific variance components given the allele frequency <i>p</i> and the genotypic values <i>a</i> and <i>d</i>.<br><br>"))
  ),
  conditionalPanel(
    condition = "input.model == 'AA'",
    div(class = "info-box",
        HTML("This model uses the notation of Falconer and Mackay (1996) Chapters 7 & 8 as well as M&auml;ki-Tanila and Hill (2014). It uses their equations to generate the population specific variance components given the allele frequencies <i>p</i> (locus A) and <i>q</i> (locus B), the genotypic values of each locus <i>a</i><sub>A</sub> and <i>a</i><sub>B</sub>, as well as the additive-by-additive effect <i>a</i><sub>AB</sub>.<br><br>"))
  ),
  conditionalPanel(
    condition = "input.model == 'Perso'",
    div(class = "info-box",
        HTML("This model uses the general least square model described in Lynch and Walsh (1998) Chapter 5 to derive the different variance components given the genotypic values input by the user as well as the chosen allele frequencies <i>p</i> (locus A) and <i>q</i> (locus B).<br><br>"))
  ),
  div(class = "info-box",
      HTML("To obtain help and details about the application, please click on the README button: "),
      actionButton("preview", "README", icon = icon("circle-question"))),
  HTML("<br>"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("sliderP", HTML("A<sub>1</sub> allele frequency (<i>p</i>)"),
                  min = 0.001, max = 0.999, value = 0.3, step = 0.01),

      conditionalPanel(
        condition = "input.model == 'AA' || input.model == 'Perso'",
        sliderInput("sliderQ", HTML("B<sub>1</sub> allele frequency (<i>q</i>)"),
                    min = 0.001, max = 0.999, value = 0.3, step = 0.01)
      ),

      conditionalPanel(
        condition = "input.model == 'AD' || input.model == 'AA'",
        uiOutput("SliderA1")
      ),

      conditionalPanel(
        condition = "input.model == 'AD'",
        sliderInput("sliderD", HTML("Genotypic value <i>d</i>"),
                    min = -10.0, max = 10.0, value = 4.0, step = 1)
      ),

      conditionalPanel(
        condition = "input.model == 'AA'",
        sliderInput("sliderA2", HTML("Locus B genotypic value <i>a</i><sub>B</sub>"),
                    min = -10.0, max = 10.0, value = 4.0),
        sliderInput("sliderAA", HTML("Additive-by-Additive interaction effect (<i>a</i><sub>AB</sub>)"),
                    min = -10.0, max = 10.0, value = 2.0)
      ),

      conditionalPanel(
        condition = "input.model == 'Perso'",
        HTML("<b>Input your genotypic values in the table:</b><br><br>"),
        rHandsontableOutput("hot"),
        HTML("<br>")
      )
    ),
    mainPanel(
      conditionalPanel(
        condition = "input.model == 'AD'",
        div(class = "model-figure",
            img(src = 'AD_model.png', width = "50%"), br(),
            tags$b("Arbitrarily assigned genotypic values (i.e., trait means per genotype class)"), br(), br())
      ),
      conditionalPanel(
        condition = "input.model == 'AA'",
        div(class = "model-figure",
            img(src = 'AAA_model.png', width = "40%"), br(),
            tags$b("Arbitrarily assigned genotypic values (i.e., trait means per genotype class)"), br(), br())
      ),

      splitLayout(
        align = "left",
        cellWidths = c("40%", "60%"),
        verticalLayout(htmlOutput("table"),
                       HTML("<b>Variance components</b>"),
                       htmlOutput("varianceA"),
                       conditionalPanel(
                         condition = "input.model == 'AD' || input.model == 'Perso'",
                         htmlOutput("varianceD")
                       ),
                       conditionalPanel(
                         condition = "input.model == 'AA' || input.model == 'Perso'",
                         htmlOutput("varianceAA")
                       ),
                       conditionalPanel(
                         condition = "input.model == 'Perso'",
                         htmlOutput("varianceAD"),
                         htmlOutput("varianceDD")
                       ),
                       htmlOutput("varianceG"),
                       htmlOutput("varianceRatio")
        ),
        plotOutput("plot"),
        cellArgs = list(style = "vertical-align: top;")
      ),

      fluidRow(
        column(4, HTML("")),
        column(8,
               conditionalPanel(
                 condition = "input.model == 'AA' || input.model == 'Perso'",
                 HTML("Graphical representation of genotypic values (closed circles) at two biallelic loci A and B. The horizontal scale shows the number of A<sub>1</sub> alleles in the genotype. The different genotypes at locus B are depicted in different colors. For each of the locus B genotypes, the linear regression line between the number of A<sub>1</sub> alleles and the genotypic value is drawn.")
               ),
               conditionalPanel(
                 condition = "input.model == 'AD'",
                 HTML("Reproduction of the Figure 7.2 of Falconer and Mackay (1996). Graphical representation of genotypic (closed blue circles) and breeding (open blue circles) values,
             of the genotypes for a locus with two alleles A<sub>1</sub> and A<sub>2</sub> at frequencies <i>p</i> and 1-<i>p</i>.
             Horizontal scale: number of A<sub>1</sub> alleles in the genotype. Vertical scales of values are: on the left, arbitrarily assigned values (see figure at the top of the page);
             on the right, deviation from the population mean (black cross). Each point size is weighted by its genotype frequency. A linear regression line between the number of A<sub>1</sub> alleles and the genotypic values is fitted by weighted least squares.")
               )
        )
      ),

      conditionalPanel(
        condition = "input.model == 'AA'",
        plotOutput(outputId = "contourplot")
      ),
      conditionalPanel(
        condition = "input.model == 'AD'",
        plotOutput(outputId = "plotDom")
      ),
      conditionalPanel(
        condition = "input.model == 'AA'",
        HTML("Additive (<i>V<sub>A</sub></i>), additive-by-additive (<i>V<sub>AA</sub></i>) and proportion of genotypic variance explained by additive variance (<i>V<sub>A</sub>/V<sub>G</sub></i>) as a function of the allele frequencies <i>p</i> and <i>q</i>. Current setting of <i>p</i> and <i>q</i> is depicted with a white cross.")
      ),
      conditionalPanel(
        condition = "input.model == 'AD'",
        HTML("Distributions of additive (<i>V<sub>A</sub></i>), dominance (<i>V<sub>D</sub></i>) and total (<i>V<sub>G</sub></i>) genetic variance on the left, and proportion of genotypic variance explained by additive variance (<i>V<sub>A</sub>/V<sub>G</sub></i>) as a function of the allele frequency <i>p</i> on the right. The current user input allele frequency <i>p</i> is depicted by a vertical red solid line in the left panel and by a red cross in the right panel.")
      )
    )
  ),
  hr(),
  uiOutput("References"),
  hr(),
  HTML("<br><b>Authors:</b> The Falconer ShinyApp was written by Valentin Hivert, based on the previous versions that had input from Luke Lloyd-Jones, Alex Holloway and Matt Robinson.<br><br><b>Contact:</b> v.hivert@imb.uq.edu.au")
)

server <- function(input, output, session) {

  ### README dialog (also shown at start-up) #############
  observeEvent(input$preview, {
    showModal(modalDialog(
      title = "README",
      readme_content(),
      easyClose = TRUE,
      size = "l",
      footer = tagList(modalButton("OK"))
    ))
  }, ignoreNULL = FALSE)

  ### Genotypic values of the general two-locus model (editable table) #############
  gv <- reactive({
    if (is.null(input$hot)) return(GV_default)
    m <- as.matrix(hot_to_r(input$hot))
    dimnames(m) <- dimnames(GV_default)
    m
  })

  # Rendered once; later edits stay in the widget and flow back through input$hot
  output$hot <- renderRHandsontable({
    rhandsontable(isolate(gv()))
  })

  ### Slider for a / aA: in renderUI so that its HTML label follows the model #############
  # The current value is isolated so the slider is only rebuilt when the model changes.
  output$SliderA1 <- renderUI({
    a <- isolate(input$sliderA1)
    if (is.null(a)) a <- 4
    label <- if (input$model == "AD") "Genotypic value <i>a</i> " else "Locus A genotypic value <i>a</i><sub>A</sub>"
    sliderInput("sliderA1", label = HTML(label), min = -10, max = 10, step = 1, value = a)
  })

  ### Single source of truth: population mean, average effects and variance components #############
  res <- reactive({
    switch(input$model,
      AD = {
        req(input$sliderA1, input$sliderD)
        model_AD(p = input$sliderP, a = input$sliderA1, d = input$sliderD)
      },
      AA = {
        req(input$sliderA1, input$sliderQ, input$sliderA2, input$sliderAA)
        model_AA(p = input$sliderP, q = input$sliderQ,
                 aA = input$sliderA1, aB = input$sliderA2, aAB = input$sliderAA)
      },
      Perso = {
        req(input$sliderQ)
        GV <- gv()
        validate(need(is.numeric(GV) && all(is.finite(GV)),
                      "Please enter a numeric genotypic value in every cell of the table."))
        model_general(GV = GV, p = input$sliderP, q = input$sliderQ)
      })
  })

  ### Variance components #############
  output$varianceA <- renderUI({
    r <- res()
    var_line("Additive variance (<b><i>V<sub>A</sub></i></b>)", r$Va, r$Va_loc)
  })
  output$varianceD <- renderUI({
    r <- res()
    var_line("Dominance variance (<b><i>V<sub>D</sub></i></b>)", r$Vd, r$Vd_loc)
  })
  output$varianceAA <- renderUI(var_line("Additive-by-Additive variance (<b><i>V<sub>AA</sub></i></b>)", res()$Vaa))
  output$varianceAD <- renderUI(var_line("Additive-by-Dominance variance (<b><i>V<sub>AD</sub></i></b>)", res()$Vad))
  output$varianceDD <- renderUI(var_line("Dominance-by-Dominance variance (<b><i>V<sub>DD</sub></i></b>)", res()$Vdd))

  output$varianceG <- renderUI({
    decomposition <- switch(res()$model,
                            AA    = "V<sub>A</sub> + V<sub>AA</sub>",
                            AD    = "V<sub>A</sub> + V<sub>D</sub>",
                            Perso = "V<sub>A</sub> + V<sub>D</sub> + V<sub>I</sub>")
    var_line(paste0("Genotypic variance (<b><i>V<sub>G</sub> = ", decomposition, "</i></b>)"), res()$Vg)
  })

  output$varianceRatio <- renderUI({
    r <- res()
    ratio <- if (is.finite(r$Va / r$Vg)) fmt(r$Va / r$Vg) else "&ndash;"   # V_G = 0 when all genotypic values are equal
    HTML(paste0("<b><i>V<sub>A</sub> / V<sub>G</sub></i></b> : ", ratio))
  })

  ### References #############
  output$References <- renderUI({
    refs <- switch(input$model,
                   AD    = REF_FALCONER,
                   AA    = c(REF_FALCONER, REF_MAKI),
                   Perso = c(REF_FALCONER, REF_MAKI, REF_LYNCH))
    HTML(paste0("<b>References</b><br>", paste(refs, collapse = "<br><br>")))
  })

  ### Main plot #############
  output$plot <- renderPlot({
    r <- res()

    if (r$model %in% c("AA", "Perso")) {
      # Genotypic values against A1 dosage, one line per locus B genotype
      par(mar = c(5.1, 4.1, 0.5, 2.1))
      x    <- c(-1, 0, 1)                     # locus A genotype coding
      GV   <- r$GV
      w    <- freq_2locus(r$p, r$q)           # rows = locus B genotype, cols = locus A genotype
      cols <- c("black", "red", "blue")       # B2B2, B1B2, B1B1
      fits <- lapply(1:3, function(k) lm(GV[k, ] ~ x, weights = w[k, ]))
      ylim <- range(GV, sapply(fits, fitted))

      plot(x, GV[1, ], cex = (1 + w[1, ])^2, ylim = ylim, xaxt = "n",
           xlab = "Genotype (effect allele counts)", ylab = "Genotypic value",
           pch = 16, cex.lab = 1.3, col = cols[1])
      axis(1, at = x, labels = GENO_AXIS_LABELS, cex.lab = 1.5)
      legend(x = "topleft", legend = c(expression(B[1] * B[1]), expression(B[1] * B[2]), expression(B[2] * B[2])),
             col = rev(cols), lty = 1, lwd = 2, title = expression(bold("Locus B genotype")), bty = "n")
      for (k in 1:3) {
        if (k > 1) points(x = x, y = GV[k, ], pch = 16, col = cols[k], cex = (1 + w[k, ])^2, lwd = 2)
        abline(fits[[k]], col = cols[k])
      }

    } else {
      # Single-locus model: reproduction of Figure 7.2 of Falconer and Mackay (1996)
      par(mar = c(5.1, 4.1, 0.5, 5.1))
      a <- r$a; d <- r$d; p <- r$p; mu <- r$mu
      alpha <- r$alpha[1]
      x    <- c(0, 1, 2)
      y    <- c(-a, d, a)                     # genotypic values
      w    <- freq_1locus(p)
      fit  <- lm(y ~ x, weights = w)
      pred <- as.vector(predict(fit))
      bv   <- c(-2 * p, 1 - 2 * p, 2 * (1 - p)) * alpha   # breeding values (deviations from M)

      plot(x, y, cex = (1 + w)^2, ylim = range(-10, 10, pred), xaxt = "n",
           xlab = "Genotype (effect allele counts)", ylab = "Value", pch = 16, col = "blue", cex.lab = 1.3)
      axis(1, at = x, labels = GENO_AXIS_LABELS)
      points(x = x, y = pred, col = "blue", cex = 2)
      abline(fit, col = "blue")
      axis(4, at = c(pred, mu), labels = round(c(bv, 0), digits = 2))
      mtext("Deviation from population mean", side = 4, line = 3, cex = 1.3)
      # The population mean lies on the weighted regression line at the mean dosage 2p
      points(x = 2 * p, y = mu, pch = 3, col = "black", cex = 2, lwd = 4)
      legend(x = if (a > 1) "topleft" else "bottomleft", legend = c("Genotypic value", "Additive (breeding) value"),
             col = c("blue", "blue"), pch = c(16, 1), bty = "n")

      if (d != 0) {
        segments(x0 = x, y0 = pred, y1 = y, lty = 3)   # dominance deviations
      }

      if (alpha != 0) {
        lines(x = c(0, 1), y = c(pred[1], pred[1]), lty = 2)
        lines(x = c(1, 1), y = c(pred[1], pred[2]), lty = 2)
        brackets(x1 = 1, x2 = 1, y1 = pred[1], y2 = pred[2], type = 1, h = if (alpha > 0) -0.1 else 0.1)
        text(x = 1.1, y = (pred[1] + pred[2]) / 2, bquote(alpha == .(round(alpha, 2))), srt = 0, cex = 1, pos = 4)
      }
    }
  })

  ### Variance as a function of p (single-locus model) #############
  output$plotDom <- renderPlot({
    r <- res()
    req(r$model == "AD")
    p     <- seq(0, 1, 0.01)
    alpha <- r$a + r$d * (1 - 2 * p)
    Va    <- 2 * p * (1 - p) * alpha^2
    Vd    <- (2 * p * (1 - p) * r$d)^2

    par(mar = c(5.1, 5.1, 2.1, 2.1), mfrow = c(1, 2))
    plot(x = p, y = Va + Vd, type = "l", col = "black", lwd = 2, lty = 4,
         xlab = LAB_P, ylab = "Genetic variance", cex.lab = 1.8, cex.axis = 1.8)
    lines(x = p, y = Va, col = "black", lwd = 2)
    lines(x = p, y = Vd, col = "black", lwd = 2, lty = 3)
    abline(v = r$p, col = "red", lwd = 2)
    legend(x = "topright", legend = c(expression(V[A]), expression(V[D]), expression(V[G])),
           lty = c(1, 3, 4), lwd = 2, bty = "n", cex = 1.5)

    plot(x = p, y = Va / (Va + Vd), type = "l", col = "black", lwd = 2,
         xlab = LAB_P, ylab = bquote(italic(V[A] / V[G])), cex.lab = 1.8, cex.axis = 1.8)
    points(x = r$p, y = r$Va / r$Vg, pch = 3, col = "red", cex = 4, lwd = 4)
  })

  ### Variance maps as a function of p and q (two-locus AA model) #############
  output$contourplot <- renderPlot({
    r <- res()
    req(r$model == "AA")
    aA  <- input$sliderA1
    aB  <- input$sliderA2
    aAB <- input$sliderAA
    pq  <- seq(0, 1, 0.01)
    VA  <- outer(pq, pq, function(p, q) 2 * p * (1 - p) * (aA + aAB * (2 * q - 1))^2 +
                                        2 * q * (1 - q) * (aB + aAB * (2 * p - 1))^2)
    VAA <- outer(pq, pq, function(p, q) 4 * p * (1 - p) * q * (1 - q) * aAB^2)

    draw_map <- function(z, clab) {
      image2D(z, pq, pq, contour = TRUE, rasterImage = TRUE, clab = clab, xlab = LAB_P, ylab = LAB_Q,
              cex.lab = 1.8, cex.axis = 1.5, colkey = list(cex.axis = 1.5, cex.clab = 1.8))
      points(x = r$p, y = r$q, pch = 3, col = "white", cex = 4, lwd = 4)
    }

    par(mar = c(5.1, 5.1, 4.1, 2.1), mfrow = c(1, 3))
    draw_map(VA,                expression(bold(V[A])))
    draw_map(VAA,               expression(bold(V[AA])))
    draw_map(VA / (VA + VAA),   bquote(bold(frac(V[A], V[G]))))
  })

  ### Main table #############
  output$table <- renderUI({
    r <- res()
    p <- r$p

    if (r$model == "AD") {
      q <- 1 - p
      a <- r$a; d <- r$d; alpha <- r$alpha[1]
      rowNames <- c("Frequencies", "Assigned values", "Genotypic value", "Additive (breeding) value", "Dominance deviation")
      colNames <- c("A<sub>2</sub>A<sub>2</sub>", "A<sub>1</sub>A<sub>2</sub>", "A<sub>1</sub>A<sub>1</sub>")
      data <- round(c(q^2, 2 * p * q, p^2,
                      -a, d, a,
                      -2 * p * (a + q * d), a * (q - p) + d * (1 - 2 * p * q), 2 * q * (a - p * d),
                      -2 * p * alpha, (q - p) * alpha, 2 * q * alpha,
                      -2 * p^2 * d, 2 * p * q * d, -2 * q^2 * d), digits = 2)

      HTML(htmlTable(matrix(data, ncol = 3, byrow = TRUE),
                     header = colNames,
                     rnames = rowNames,
                     rgroup = c("", "Deviations from population mean:"),
                     n.rgroup = c(2, 3),
                     cgroup = c("Genotypes"),
                     n.cgroup = c(3),
                     caption = paste0("<br>Population mean: M = ", round(r$mu, digits = 3),
                                      "<br>Average effect of gene-substitution: &#120572; = ", alpha)))
    } else {
      q <- r$q
      geno.freq <- round(freq_2locus(p, q), digits = 2)   # rows = locus B, cols = locus A (same layout as r$GV)
      data <- matrix(paste0(geno.freq, " (", r$GV, ")"), nrow = 3, ncol = 3)

      HTML(htmlTable(data,
                     header = c("A<sub>2</sub>A<sub>2</sub>", "A<sub>1</sub>A<sub>2</sub>", "A<sub>1</sub>A<sub>1</sub>"),
                     rnames = c("<b>B<sub>2</sub>B<sub>2</sub></b>", "<b>B<sub>1</sub>B<sub>2</sub></b>", "<b>B<sub>1</sub>B<sub>1</sub></b>"),
                     rowlabel = paste0("q=", q, "\\p=", p),
                     align.header = c('l', rep('c', 3)),
                     caption = paste0("<br>Population mean : M = ", round(r$mu, digits = 2),
                                      "<br>Locus A average effect of gene-substitution: &#120572;<sub>A</sub> = ", round(r$alpha[1], digits = 2),
                                      "<br>Locus B average effect of gene-substitution: &#120572;<sub>B</sub> = ", round(r$alpha[2], digits = 2),
                                      "<br><br>Genotype frequencies (genotypic values) are :")))
    }
  })
}

#############################################
# Run the app
#
shinyApp(ui = ui, server = server)
