# Falconer ShinyApp — v1.1 changes

This release is a code revision of v1.0 (Feb 2021). The user interface, the three models and what they are meant to show are unchanged. It fixes one display bug, one performance bug and several errors in the README text, and restructures the code so it is easier to maintain.


## New file layout

| File | Content |
|---|---|
| `app.R` | UI and server only (renamed from `App.R`) |
| `R/genetics.R` | *new*: model functions `model_AD()`, `model_AA()`, `model_general()`, `Compute_GeneticVariances()`, frequency helpers, default genotypic values |
| `README.md` | updated: same text as the in-app README dialog, with the corrections in section 4 |
| `R/readme.R` | *new*: `readme_content()`, the HTML of the README dialog (moved out of the server) |
| `www/` | unchanged (`AD_model.png`, `AAA_model.png`) |

`app.R` loads the two new files with `source()`. Recent Shiny versions also load `R/` automatically; loading them twice does no harm.

---

## 1. Bug fixes

### 1.1 Genotype frequencies transposed in the two-locus table (AA and general models)
- **Before:** `geno.freq <- freq.1 %*% t(freq.2)` put locus A in the rows. The table shows locus B genotypes in the rows and locus A in the columns, so each genotypic value was shown next to the frequency of the transposed genotype. At the default p = q = 0.3 this can't be seen. For example, with p = 0.1 and q = 0.5, the B₂B₂/A₁A₁ cell showed 0.20 instead of 0.00 (0.0025).
- **After:** `freq_2locus(p, q) = outer(freq_1locus(q), freq_1locus(p))` (rows = B, columns = A). The table, the plot weights and `Compute_GeneticVariances()` all use this one function. The variance components were never affected: they already used the correct orientation.

### 1.2 The a / a_A slider was rebuilt on every move
- **Before:** `output$SliderA1` read `input$sliderA1`, so each move re-rendered the whole slider widget, and an `init$idx` flag was needed for the first render.
- **After:** the current value is read with `isolate()` (defaulting to 4), so the slider is only rebuilt when the model changes (to update its label). `init` has been removed.

### 1.3 Unused `shinyalert` dependency
`useShinyalert()` was called, but no alert was ever used (the README uses `showModal()`). The function is deprecated in shinyalert ≥ 3.0. Removed along with the package. `htmltools` is also no longer loaded explicitly: every function used from it is re-exported by `shiny`.

### 1.4 Invalid input in the general model's table
If a cell of the handsontable was cleared or made non-numeric, `lm()` and `Compute_GeneticVariances()` failed, and the variance text kept showing the previous values. The app now shows the message *"Please enter a numeric genotypic value in every cell of the table."* in place of the outputs until the table is valid.

### 1.5 V_A / V_G when V_G = 0
With a = d = 0 (or all genotypic values equal), the ratio showed `NaN`. It now shows "–".

### 1.6 README icon
`icon("question-circle")` produced a deprecation message with Font Awesome 6, the version shipped with recent Shiny. It is now `icon("circle-question")`, which looks the same.

### 1.7 Main file renamed `App.R` → `app.R`
Shiny looks for `app.R`. The capitalized name only worked on case-insensitive file systems, so it could fail on Linux servers.

### 1.8 Population-mean cross in the single-locus plot when α = 0
- **Before:** the black cross was placed at `(mu - intercept) / slope`, which is undefined when α = 0, so the cross disappeared.
- **After:** it is placed at x = 2p, the mean allele dosage. The weighted regression line always passes through (2p, M), so the position is identical whenever α ≠ 0, and the cross now also shows when α = 0.

---

## 2. Reactivity restructure

- **Before:** all results (μ, α, variance components) were computed inside `output$plot` and written into `v <- reactiveValues(...)`. The table, the variance text and the red cross in `plotDom` read `v`. So:
  - they only updated after the main plot re-rendered;
  - they could render twice per slider move (once with old values);
  - α and the variances were computed separately in `plot`, `plotDom`, `contourplot` and `table`.
- **After:** a single `res <- reactive({...})` calls the matching model function and returns one list with the same fields for every model (`mu`, `alpha`, `Va`, `Va_loc`, `Vd`, `Vd_loc`, `Vaa`, `Vad`, `Vdd`, `Vg`, `GV`, …). Every output reads `res()`; the render functions only draw.
- `plotDom` and `contourplot` now start with `req(r$model == "AD")` / `req(r$model == "AA")`, so they never draw with another model's results during a tab switch.
- The handsontable logic is now a small reactive, `gv()`. The previous `observe()` with `values[["previous"]]` (never read) is gone, and the table is rendered once rather than after every edit.

## 3. Code cleanup (no change in output)

- `Compute_GeneticVariances()`:
  - Removed dead code: the unused `VAA/sqrt(...)` line, the unused `aa <- sum(...)`, the `rm()` calls and the leftover "Summary" comments.
  - Replaced the one-line `do.call(rbind, lapply(...))` for the conditional means of allele pairs with an explicit double loop.
  - Uses `<-` and `TRUE` consistently, and comments give the orientation of each matrix.
  - The maths is unchanged.
- The five near-identical variance renderers now use the `var_line()` / `fmt()` helpers. The displayed text is unchanged.
- References are built from three constants instead of three copied HTML strings.
- Main plot (two-locus models): the three copies of the weights / `cex` / `lm` / `abline` code are replaced by a loop over locus B genotypes with `freq_2locus()` weights.
- Main plot (single-locus model):
  - The fitted model is no longer named `lm` (which masked the function).
  - The three `lines()` calls for dominance deviations are replaced by one `segments()` call.
  - The α label is rounded to two decimals.
- Contour plot: the 101 × 101 double loop is replaced by `outer()`, and the three `image2D()` calls share a `draw_map()` helper.
- The unused `Geno_freq` computation in the server's general-model branch is removed.
- Uses `library()` instead of `require()`, so a missing package stops with a clear error.
- Uses `%in%` / `||` instead of `|` inside `if()`.
- Inline styles repeated three times per `div` are replaced by two CSS classes (`.info-box`, `.model-figure`) with the same values.
- The two `conditionalPanel`s for V_AD and V_DD (same condition) are merged into one.
- `paste(x, " (", gv, ")")` gave `"0.09  ( 18 )"`; it now gives `"0.09 (18)"`.

## 4. Text corrections

### README dialog: formulas
The app uses the −a / 0 / a coding at both loci. Three README formulas used the 0 / a / 2a (Mäki-Tanila & Hill) parametrization or were missing a square. The README now gives the equations for the −a / 0 / a coding. The code was correct in every case and has not changed.

| Quantity | v1.0 README | v1.1 README (matches code) |
|---|---|---|
| α_A | a_A + 2q·a_AB | a_A + (2q−1)·a_AB |
| α_B | a_B + 2p·a_AB | a_B + (2p−1)·a_AB |
| V_A | 2p(1−p)[a_A+2q·a_AB] + 2q(1−q)[a_B+2p·a_AB] | 2p(1−p)[a_A+(2q−1)a_AB]² + 2q(1−q)[a_B+(2p−1)a_AB]² |
| V_AA | 4p(1−p)q(1−q)·a_AB | 4p(1−p)q(1−q)·a_AB² |
| Freq. A₂A₂B₁B₁ (genotype table) | (1−p)q² | (1−p)²q² |

The new α and V_A expressions were checked numerically against the least-squares decomposition of `Compute_GeneticVariances()` (see tests, check 2).

### Typos and grammar
- **README:**
  - "Welcome in" → "Welcome to"
  - "allele frequencies of the 9 genotypes" → "frequencies of the 9 genotypes" (the table gives genotype frequencies)
  - duplicated comma ("at a given locus, , i.e.")
  - unclosed parenthesis after "the regression coefficient β"
  - "the variance simply become" → "becomes"
  - "and contain both a term" → "contains"
  - "a term due pairwise" → "due to the pairwise"
  - "partitioned in five components" → "into"
- **Model description (AA):** "aditive" → "additive". The AD description now shows *p* in italics, consistent with the other two.
- **Figure captions:**
  - "The horizontal scale show" → "shows"
  - "number of A₁ allele" → "alleles"
  - "freqencies" → "frequencies"
  - "weighted least square" → "least squares"
  - small punctuation fixes in the Figure 7.2 caption
- **Code comments:** "Conditionnal" → "Conditional"
- "Mäki" is written as `M&auml;ki` so it displays correctly whatever the file encoding.
