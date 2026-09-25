###################################################################################
#
# FALCONER SHINY APP - model computations
#
# Plain (non-reactive) functions used by App.R. Each model_*() function returns a
# list with the same fields so that every output in the server can read a single
# reactive result:
#   model, p, q, mu, alpha (locus A, locus B), Va, Va_loc, Vd, Vd_loc,
#   Vaa, Vad, Vdd, Vg, GV (3x3 genotypic values; rows = locus B, cols = locus A)
#
# Genotype order everywhere: (A2A2, A1A2, A1A1) and (B2B2, B1B2, B1B1).
#
###################################################################################

## Genotype labels (unicode subscripts) used for the 3x3 genotypic value matrices
GENO_A <- c("A₂A₂", "A₁A₂", "A₁A₁")
GENO_B <- c("B₂B₂", "B₁B₂", "B₁B₁")

## Default genotypic values for the general two-locus model
## (example values from Lynch and Walsh 1998). Rows = locus B, columns = locus A.
GV_default <- matrix(c(18.0, 40.9,  61.1,
                       54.6, 47.6,  66.5,
                       47.8, 83.6, 101.7),
                     ncol = 3, byrow = TRUE, dimnames = list(GENO_B, GENO_A))

## Hardy-Weinberg genotype frequencies (A2A2, A1A2, A1A1) for allele frequency p of A1
freq_1locus <- function(p) c((1 - p)^2, 2 * p * (1 - p), p^2)

## Two-locus genotype frequencies under linkage equilibrium (rows = locus B, cols = locus A)
freq_2locus <- function(p, q) outer(freq_1locus(q), freq_1locus(p))


### Single-locus additive and dominance model (Falconer & Mackay 1996, Ch. 7-8) ####
# Genotypic values -a, d, a for A2A2, A1A2, A1A1.
model_AD <- function(p, a, d) {
  H     <- 2 * p * (1 - p)
  alpha <- a + d * (1 - 2 * p)           # average effect of gene substitution
  mu    <- a * (2 * p - 1) + H * d       # population mean
  Va    <- H * alpha^2
  Vd    <- (H * d)^2
  list(model = "AD", p = p, q = NA, a = a, d = d,
       mu = mu, alpha = c(alpha, NA),
       Va = Va, Va_loc = NULL, Vd = Vd, Vd_loc = NULL,
       Vaa = 0, Vad = 0, Vdd = 0, Vg = Va + Vd, GV = NULL)
}


### Two-locus additive and additive-by-additive model ##############################
# Falconer coding (-1, 0, 1) at both loci: G = aA*xA + aB*xB + aAB*xA*xB
# (re-parametrization of Maki-Tanila and Hill 2014).
model_AA <- function(p, q, aA, aB, aAB) {
  x  <- c(-1, 0, 1)
  HA <- 2 * p * (1 - p)
  HB <- 2 * q * (1 - q)
  alphaA <- aA + aAB * (2 * q - 1)
  alphaB <- aB + aAB * (2 * p - 1)
  mu     <- aA * (2 * p - 1) + aB * (2 * q - 1) + aAB * (1 - 2 * p - 2 * q + 4 * p * q)
  Va_loc <- c(Locus_A = HA * alphaA^2, Locus_B = HB * alphaB^2)
  Vaa    <- HA * HB * aAB^2
  GV <- outer(x, x, function(xB, xA) aA * xA + aB * xB + aAB * xA * xB)  # rows = B, cols = A
  dimnames(GV) <- list(GENO_B, GENO_A)
  list(model = "AA", p = p, q = q,
       mu = mu, alpha = c(alphaA, alphaB),
       Va = sum(Va_loc), Va_loc = Va_loc, Vd = 0, Vd_loc = c(0, 0),
       Vaa = Vaa, Vad = 0, Vdd = 0, Vg = sum(Va_loc) + Vaa, GV = GV)
}


### General two-locus model (least squares, Lynch and Walsh 1998, Ch. 5) ###########
model_general <- function(GV, p, q) {
  res <- Compute_GeneticVariances(GV = GV, p = p, q = q)
  list(model = "Perso", p = p, q = q,
       mu = res$mu, alpha = unname(res$alpha),
       Va = res$var$VA[[2]], Va_loc = res$var$VA[[1]],
       Vd = res$var$VD[[2]], Vd_loc = res$var$VD[[1]],
       Vaa = res$var$VAA, Vad = res$var$VAD, Vdd = res$var$VDD,
       Vg = res$var$VG, GV = GV)
}


### Derive the components of genetic variance for any 3x3 genotypic value matrix ###
# GV: matrix of genotypic values, rows = locus B (B2B2, B1B2, B1B1),
#     columns = locus A (A2A2, A1A2, A1A1).
# p, q: frequencies of alleles A1 and B1.
Compute_GeneticVariances <- function(GV, p, q) {

  freqA     <- freq_1locus(p)
  freqB     <- freq_1locus(q)
  Geno_freq <- freqB %*% t(freqA)

  # Population mean
  M <- sum(GV * Geno_freq)

  ############################################################
  # Additive effects
  CondQ.freq <- c(1 - q, q) %*% t(freqA)   # frequencies conditional on one B allele
  CondP.freq <- freqB %*% t(c(1 - p, p))   # frequencies conditional on one A allele
  Cond.Mean  <- matrix(0, nrow = 2, ncol = 2)
  for (i in 1:2) {
    Cond.Mean[1, i] <- sum(GV[, c(i, i + 1)] * CondP.freq)  # conditional mean, locus A
    Cond.Mean[2, i] <- sum(GV[c(i, i + 1), ] * CondQ.freq)  # conditional mean, locus B
  }
  Alpha <- Cond.Mean - M   # average effects (row 1 = locus A, row 2 = locus B; col 1 = allele 2, col 2 = allele 1)
  Expected.Geno <- cbind(2 * Alpha[, 1], rowSums(Alpha), 2 * Alpha[, 2])

  ############################################################
  # Dominance effects (row 1 = locus A, row 2 = locus B; one column per genotype)
  D_effects <- matrix(NA, nrow = 2, ncol = 3)
  D_effects[1, ] <- freqB %*% GV    - M - Expected.Geno[1, ]
  D_effects[2, ] <- freqA %*% t(GV) - M - Expected.Geno[2, ]

  ############################################################
  # Additive-by-additive effects
  # Conditional mean of each pair of alleles (rows = B allele, cols = A allele)
  w_alleles <- c(1 - q, q) %*% t(c(1 - p, p))
  Cond.Mean.AB <- matrix(0, nrow = 2, ncol = 2)
  for (i in 1:2) {
    for (j in 1:2) {
      Cond.Mean.AB[i, j] <- sum(GV[c(i, i + 1), c(j, j + 1)] * w_alleles)
    }
  }
  AA_effects <- Cond.Mean.AB - M - rbind(Alpha[2, 1] + Alpha[1, ], Alpha[2, 2] + Alpha[1, ])

  # Additive-by-dominance effects
  AD_effects <- matrix(NA, nrow = 4, ncol = 3)
  for (i in 1:2) {
    # Conditional on one A allele (rows 1-2) and on one B allele (rows 3-4)
    AD_effects[i, ]     <- colSums(c(1 - p, p) * t(GV[, c(i, i + 1)])) - M - D_effects[2, ] - Expected.Geno[2, ] -
                           Alpha[1, i] - c(2 * AA_effects[1, i], sum(AA_effects[, i]), 2 * AA_effects[2, i])
    AD_effects[i + 2, ] <- colSums(c(1 - q, q) * GV[c(i, i + 1), ]) - M - D_effects[1, ] - Expected.Geno[1, ] -
                           Alpha[2, i] - c(2 * AA_effects[i, 1], sum(AA_effects[i, ]), 2 * AA_effects[i, 2])
  }

  ############################################################
  # Dominance-by-dominance effects (residual of all other terms)
  BV      <- t(outer(Expected.Geno[1, ], Expected.Geno[2, ], FUN = "+"))   # breeding values
  D.Geno  <- t(outer(D_effects[1, ], D_effects[2, ], FUN = "+"))
  AA.Geno <- t(apply(AA_effects, MARGIN = 1, FUN = function(x) c(2 * x[1], sum(x), 2 * x[2])))
  AA.Geno <- rbind(2 * AA.Geno[1, ], AA.Geno[1, ] + AA.Geno[2, ], 2 * AA.Geno[2, ])
  AD.Geno <- t(rbind(2 * AD_effects[1, ], colSums(AD_effects[1:2, ]), 2 * AD_effects[2, ])) +
               rbind(2 * AD_effects[3, ], colSums(AD_effects[3:4, ]), 2 * AD_effects[4, ])

  DD_effects <- GV - M - BV - D.Geno - AA.Geno - AD.Geno

  ############################################################
  # VARIANCE COMPONENTS

  # Total genotypic variance
  VG <- sum(Geno_freq * GV^2) - M^2

  # Additive variance
  VA_locus <- 2 * rowSums(matrix(c(1 - p, p, 1 - q, q), ncol = 2, byrow = TRUE) * Alpha^2)
  names(VA_locus) <- c("Locus_A", "Locus_B")
  VA <- sum(VA_locus)

  # Dominance variance
  VD_locus <- rowSums(rbind(freqA, freqB) * D_effects^2)
  names(VD_locus) <- c("Locus_A", "Locus_B")
  VD <- sum(VD_locus)

  # Additive-by-additive variance
  VAA <- 4 * weighted.mean(AA_effects^2,
                           w = matrix(c((1 - p) * (1 - q), (1 - q) * p, q * (1 - p), q * p), 2, 2, byrow = TRUE))

  # Additive-by-dominance variance
  VAD <- 4 * weighted.mean(AD_effects^2,
                           w = matrix(c((1 - p) * freqB, p * freqB, (1 - q) * freqA, q * freqA), ncol = 3, byrow = TRUE))

  # Dominance-by-dominance variance
  VDD <- weighted.mean(DD_effects^2, w = Geno_freq)

  res.var <- list(VA = list(VA_locus, VA), VD = list(VD_locus, VD),
                  VAA = VAA, VAD = VAD, VDD = VDD, VG = VG)
  list(mu = M, alpha = Alpha[, 2] - Alpha[, 1], var = res.var)
}
