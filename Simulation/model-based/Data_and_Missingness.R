


generateData <- function(n, scenario = "sce1", dgp = "DGP1") {
  # Coefficients for different scenarios and DGPs
  coefficients <- list(
    DGP1 = list(
      sce1 = list(a0 = -1.30, b0 = -1.90, b1 = 0.40, c0 = 0.40, c1 = 0.70, d0 = -0.60, d1 = -0.70, e0 = -0.50,
                  f0 = -1.60, f1 = 1.00, f2 = 0.60, f3 = -2.00, f4 = -1.20, f5 = 1.20, f6 = 0.70, f7 = 0.50, 
                  f8 = 0.10, f9 = 0.50, f10 = -0.20, f11 = -0.50, f12 = 0.80,
                  g0 = -0.80, g1 = 0.10, g2 = 0.40, g3 = 0.70, g4 = 0.20, g5 = 0.30, g6 = 0.20, g7 = -0.50, 
                  g8 = 1.00, g9 = 0.20, g10 = 0.20, g11 = 0.40, g12 = -0.20, g13 = -1.20, g14 = -1.00, g15 = -0.20, 
                  g16 = -0.40, g17 = 1.70),
      sce2 = list(a0 = -1.30, b0 = -1.90, b1 = 0.40, c0 = 0.40, c1 = 0.70, d0 = -0.60, d1 = -0.70, e0 = -0.50,
                  f0 = -1.80, f1 = 1.00, f2 = 0.60, f3 = -2.00, f4 = -1.20, f5 = 1.20, f6 = 0.70, f7 = 1.30, 
                  f8 = 0.30, f9 = 1.20, f10 = -2.00, f11 = -1.90, f12 = 2.00,
                  g0 = -0.80, g1 = 0.10, g2 = 0.40, g3 = 0.70, g4 = 0.20, g5 = 0.30, g6 = 0.20, g7 = -0.50, 
                  g8 = 1.00, g9 = 0.20, g10 = 0.20, g11 = 0.40, g12 = -0.20, g13 = -1.20, g14 = -1.00, g15 = -0.20, 
                  g16 = -0.40, g17 = 1.70),
      sce3 = list(a0 = -1.30, b0 = -1.90, b1 = 0.40, c0 = 0.40, c1 = 0.70, d0 = -0.60, d1 = -0.70, e0 = -0.50,
                  f0 = -2.85, f1 = 1.00, f2 = 0.60, f3 = -2.00, f4 = -1.20, f5 = 1.20, f6 = 0.70, f7 = 3.60, 
                  f8 = 0.60, f9 = 2.80, f10 = -4.40, f11 = -4.20, f12 = 4.40,
                  g0 = -0.80, g1 = 0.10, g2 = 0.40, g3 = 0.70, g4 = 0.20, g5 = 0.30, g6 = 0.20, g7 = -0.50, 
                  g8 = 1.00, g9 = 0.20, g10 = 0.20, g11 = 0.40, g12 = -0.20, g13 = -1.20, g14 = -1.00, g15 = -0.20, 
                  g16 = -0.40, g17 = 1.70)
    ),
    DGP2 = list(
      sce1 = list(a0 = -0.40, b0 = -0.90, b1 = -0.50, c0 = 0.40, c1 = 0.70, d0 = 3.00, d1 = -0.10, e0 = 1.00,
                  f0 = -2.95, f1 = 0.40, f2 = 0.80, f3 = 0.50, f4 = 0.20, f5 = 0.40, f6 = 0.70, f7 = -0.35, 
                  f8 = 0.25, f9 = 0.20, f10 = -0.25, f11 = 0.20, f12 = -0.30,
                  g0 = 0.20, g1 = -0.30, g2 = -0.40, g3 = 0.20, g4 = 0.20, g5 = 0.10, g6 = 0.20, g7 = -0.30, 
                  g8 = 0.20, g9 = -0.10, g10 = -0.40, g11 = 0.10, g12 = -0.10, g13 = -0.10, g14 = 0.10, g15 = -0.10, 
                  g16 = 0.10, g17 = 0.10),
      sce2 = list(a0 = -0.40, b0 = -0.90, b1 = -0.50, c0 = 0.40, c1 = 0.70, d0 = 3.00, d1 = -0.10, e0 = 1.00,
                  f0 = -2.95, f1 = 0.40, f2 = 0.80, f3 = 0.50, f4 = 0.20, f5 = 0.40, f6 = 0.70, f7 = -0.70, 
                  f8 = 0.50, f9 = 0.40, f10 = -0.50, f11 = 0.40, f12 = -0.55,
                  g0 = 0.20, g1 = -0.30, g2 = -0.40, g3 = 0.20, g4 = 0.20, g5 = 0.10, g6 = 0.20, g7 = -0.30, 
                  g8 = 0.20, g9 = -0.10, g10 = -0.40, g11 = 0.10, g12 = -0.10, g13 = -0.10, g14 = 0.10, g15 = -0.10, 
                  g16 = 0.10, g17 = 0.10),
      sce3 = list(a0 = -0.40, b0 = -0.90, b1 = -0.50, c0 = 0.40, c1 = 0.70, d0 = 3.00, d1 = -0.10, e0 = 1.00,
                  f0 = -3.25, f1 = 0.40, f2 = 0.80, f3 = 0.50, f4 = 0.20, f5 = 0.40, f6 = 0.70, f7 = -1.30, 
                  f8 = 1.00, f9 = 0.80, f10 = -1.00, f11 = 0.80, f12 = -0.95,
                  g0 = 0.20, g1 = -0.30, g2 = -0.40, g3 = 0.20, g4 = 0.20, g5 = 0.10, g6 = 0.20, g7 = -0.30, 
                  g8 = 0.20, g9 = -0.10, g10 = -0.40, g11 = 0.10, g12 = -0.10, g13 = -0.10, g14 = 0.10, g15 = -0.10, 
                  g16 = 0.10, g17 = 0.10)
    ),
    DGP3 = list(
      sce1 = list(a0 = 0.40, b0 = 0.85, b1 = 0.30, c0 = -0.35, c1 = 0.30, d0 = 3.00, d1 = -0.10, e0 = 1.00,
                  f0 = -2.80, f1 = 0.40, f2 = 0.80, f3 = 0.50, f4 = 0.20, f5 = 0.40, f6 = 0.70, f7 = -0.30, 
                  f8 = 0.25, f9 = 0.20, f10 = -0.25, f11 = 0.15, f12 = -0.25,
                  g0 = 0.30, g1 = -0.30, g2 = -0.40, g3 = 0.20, g4 = 0.20, g5 = 0.10, g6 = 0.20, g7 = -0.30, 
                  g8 = 0.20, g9 = -0.10, g10 = -0.40, g11 = 0.10, g12 = -0.10, g13 = -0.10, g14 = 0.10, g15 = -0.10, 
                  g16 = 0.10, g17 = 0.10),
      sce2 = list(a0 = 0.40, b0 = 0.85, b1 = 0.30, c0 = -0.35, c1 = 0.30, d0 = 3.00, d1 = -0.10, e0 = 1.00,
                  f0 = -2.50, f1 = 0.40, f2 = 0.80, f3 = 0.50, f4 = 0.20, f5 = 0.40, f6 = 0.70, f7 = -0.60, 
                  f8 = 0.50, f9 = 0.40, f10 = -0.50, f11 = 0.30, f12 = -0.50,
                  g0 = 0.30, g1 = -0.30, g2 = -0.40, g3 = 0.20, g4 = 0.20, g5 = 0.10, g6 = 0.20, g7 = -0.30, 
                  g8 = 0.20, g9 = -0.10, g10 = -0.40, g11 = 0.10, g12 = -0.10, g13 = -0.10, g14 = 0.10, g15 = -0.10, 
                  g16 = 0.10, g17 = 0.10),
      sce3 = list(a0 = 0.40, b0 = 0.85, b1 = 0.30, c0 = -0.35, c1 = 0.30, d0 = 3.00, d1 = -0.10, e0 = 1.00,
                  f0 = -2.80, f1 = 0.40, f2 = 0.80, f3 = 0.50, f4 = 0.20, f5 = 0.40, f6 = 0.70, f7 = -1.20, 
                  f8 = 1.00, f9 = 0.80, f10 = -1.00, f11 = 0.60, f12 = -0.85,
                  g0 = 0.30, g1 = -0.30, g2 = -0.40, g3 = 0.20, g4 = 0.20, g5 = 0.10, g6 = 0.20, g7 = -0.30, 
                  g8 = 0.20, g9 = -0.10, g10 = -0.40, g11 = 0.10, g12 = -0.10, g13 = -0.10, g14 = 0.10, g15 = -0.10, 
                  g16 = 0.10, g17 = 0.10)
    ),
    DGP4 = list(
      sce1 = list(a0 = 0.40, b0 = 0.90, b1 = -0.50, c0_1 = -0.20, c0_2 = 0.70, c0_3 = -0.30, c0_4 = -1.00,
                  c1_1 = 0.70, c1_2 = -0.30, c1_3 = -0.80, c1_4 = -0.60, d0 = 3.00, d1 = -0.10, e0 = 1.00,
                  f0 = -3.10, f1 = 0.40, f2 = 0.80, f3 = 0.50, f4 = 0.30, f5 = 0.60, f6 = 0.20, f7 = 0.40, f8 = 0.70,
                  f9 = -0.20, f10 = 0.15, f11 = 0.15, f12 = 0.15, f13 = -0.10, f14 = -0.15, f15 = -0.20, f16 = 0.05,   
                  f17 = -0.10, f18 = -0.20, f19 = -0.05, f20 = -0.10,
                  g0 = 0.15, g1 = -0.30, g2 = -0.40, g3 = 0.20, g4 = -0.30, g5 = 0.35, g6 = 0.20, g7 = 0.10, g8 = 0.20, 
                  g9 = -0.30, g10 = 0.40, g11 = -0.50, g12 = 0.20, g13 = -0.10, g14 = -0.40, g15 = 0.20, 
                  g16 = -0.40, g17 = 0.10, g18 = -0.30, g19 = -0.20, g20 = -0.10, g21 = -0.10, g22 = 0.10, 
                  g23 = -0.10, g24 = 0.10, g25 = 0.10),
      sce2 = list(a0 = 0.40, b0 = 0.90, b1 = -0.50, c0_1 = -0.20, c0_2 = 0.70, c0_3 = -0.30, c0_4 = -1.00,
                  c1_1 = 0.70, c1_2 = -0.30, c1_3 = -0.80, c1_4 = -0.60, d0 = 3.00, d1 = -0.10, e0 = 1.00,
                  f0 = -2.40, f1 = 0.40, f2 = 0.80, f3 = 0.50, f4 = 0.30, f5 = 0.60, f6 = 0.20, f7 = 0.40, f8 = 0.70,
                  f9 = -0.40, f10 = 0.30, f11 = 0.30, f12 = 0.30, f13 = -0.20, f14 = -0.30, f15 = -0.50, f16 = 0.10,   
                  f17 = -0.20, f18 = -0.60, f19 = -0.40, f20 = -0.30,
                  g0 = 0.15, g1 = -0.30, g2 = -0.40, g3 = 0.20, g4 = -0.30, g5 = 0.35, g6 = 0.20, g7 = 0.10, g8 = 0.20, 
                  g9 = -0.30, g10 = 0.40, g11 = -0.50, g12 = 0.20, g13 = -0.10, g14 = -0.40, g15 = 0.20, 
                  g16 = -0.40, g17 = 0.10, g18 = -0.30, g19 = -0.20, g20 = -0.10, g21 = -0.10, g22 = 0.10, 
                  g23 = -0.10, g24 = 0.10, g25 = 0.10),
      sce3 = list(a0 = 0.40, b0 = 0.90, b1 = -0.50, c0_1 = -0.20, c0_2 = 0.70, c0_3 = -0.30, c0_4 = -1.00,
                  c1_1 = 0.70, c1_2 = -0.30, c1_3 = -0.80, c1_4 = -0.60, d0 = 3.00, d1 = -0.10, e0 = 1.00,
                  f0 = -2.55, f1 = 0.40, f2 = 0.80, f3 = 0.50, f4 = 0.30, f5 = 0.60, f6 = 0.20, f7 = 0.40, f8 = 0.70,
                  f9 = -0.80, f10 = 0.60, f11 = 0.60, f12 = 0.50, f13 = -0.40, f14 = -0.50, f15 = -0.80, f16 = 0.20,
                  f17 = -0.40, f18 = -0.90, f19 = -0.70, f20 = -0.50,
                  g0 = 0.15, g1 = -0.30, g2 = -0.40, g3 = 0.20, g4 = -0.30, g5 = 0.35, g6 = 0.20, g7 = 0.10, g8 = 0.20, 
                  g9 = -0.30, g10 = 0.40, g11 = -0.50, g12 = 0.20, g13 = -0.10, g14 = -0.40, g15 = 0.20, 
                  g16 = -0.40, g17 = 0.10, g18 = -0.30, g19 = -0.20, g20 = -0.10, g21 = -0.10, g22 = 0.10, 
                  g23 = -0.10, g24 = 0.10, g25 = 0.10)
    ),
    DGP5 = list(
      sce1 = list(a0 = 0.40, b0 = 0.90, b1 = -0.50, c0_1 = -0.20, c0_2 = 0.70, c0_3 = -0.30, c0_4 = -1.00,
                  c1_1 = 0.70, c1_2 = -0.30, c1_3 = -0.80, c1_4 = -0.60, d0 = 3.00, d1 = -0.10, e0 = 1.00, s1 = 0.10, r1 = 0.10,
                  f0 = -3.45, f1 = 0.40, f2 = 0.80, f3 = 0.50, f4 = 0.30, f5 = 0.60, f6 = 0.20, f7 = 0.40, f8 = 0.40, f9 = 0.70,
                  f10 = -0.30, f11 = 0.20, f12 = 0.25, f13 = 0.30, f14 = -0.30, f15 = -0.25, f16 = 0.35, f17 = 0.10,   
                  f18 = -0.20, f19 = -0.30, f20 = -0.15, f21 = -0.15, f22 = -0.20, f23 = -0.15, f24 = 0.10,
                  g0 = -1.3, g1 = -0.30, g2 = -0.40, g3 = 0.20, g4 = -0.30, g5 = 0.35, g6 = 0.20, g7 = 0.10, g8 = -0.20, 
                  g9 = 0.20, g10 = -0.30, g11 = 0.40, g12 = -0.50, g13 = 0.20, g14 = -0.10, g15 = -0.40, g16 = 0.20, 
                  g17 = -0.40, g18 = 0.10, g19 = -0.30, g20 = -0.20, g21 = -0.10, g22 = -0.10, g23 = 0.10, 
                  g24 = -0.10, g25 = 0.10, g26 = -0.10, g27 = 0.10, g28 = 0.10, g29 = 0.10),
      sce2 = list(a0 = 0.40, b0 = 0.90, b1 = -0.50, c0_1 = -0.20, c0_2 = 0.70, c0_3 = -0.30, c0_4 = -1.00,
                  c1_1 = 0.70, c1_2 = -0.30, c1_3 = -0.80, c1_4 = -0.60, d0 = 3.00, d1 = -0.10, e0 = 1.00, s1 = 0.10, r1 = 0.10,
                  f0 = -3.00, f1 = 0.40, f2 = 0.80, f3 = 0.50, f4 = 0.30, f5 = 0.60, f6 = 0.20, f7 = 0.40, f8 = 0.40, f9 = 0.70,
                  f10 = -0.60, f11 = 0.40, f12 = 0.50, f13 = 0.60, f14 = -0.60, f15 = -0.50, f16 = 0.70, f17 = 0.20,   
                  f18 = -0.40, f19 = -0.60, f20 = -0.30, f21 = -0.30, f22 = -0.40, f23 = -0.30, f24 = 0.20,
                  g0 = -1.3, g1 = -0.30, g2 = -0.40, g3 = 0.20, g4 = -0.30, g5 = 0.35, g6 = 0.20, g7 = 0.10, g8 = -0.20, 
                  g9 = 0.20, g10 = -0.30, g11 = 0.40, g12 = -0.50, g13 = 0.20, g14 = -0.10, g15 = -0.40, g16 = 0.20, 
                  g17 = -0.40, g18 = 0.10, g19 = -0.30, g20 = -0.20, g21 = -0.10, g22 = -0.10, g23 = 0.10, 
                  g24 = -0.10, g25 = 0.10, g26 = -0.10, g27 = 0.10, g28 = 0.10, g29 = 0.10),
      sce3 = list(a0 = 0.40, b0 = 0.90, b1 = -0.50, c0_1 = -0.20, c0_2 = 0.70, c0_3 = -0.30, c0_4 = -1.00,
                  c1_1 = 0.70, c1_2 = -0.30, c1_3 = -0.80, c1_4 = -0.60, d0 = 3.00, d1 = -0.10, e0 = 1.00, s1 = 0.10, r1 = 0.10,
                  f0 = -2.80, f1 = 0.40, f2 = 0.80, f3 = 0.50, f4 = 0.30, f5 = 0.60, f6 = 0.20, f7 = 0.40, f8 = 0.40, f9 = 0.70,
                  f10 = -1.10, f11 = 0.80, f12 = 0.90, f13 = 1.20, f14 = -1.20, f15 = -0.90, f16 = 1.05, f17 = 0.40,   
                  f18 = -0.70, f19 = -1.10, f20 = -0.50, f21 = -0.50, f22 = -0.70, f23 = -0.50, f24 = 0.45,
                  g0 = -1.3, g1 = -0.30, g2 = -0.40, g3 = 0.20, g4 = -0.30, g5 = 0.35, g6 = 0.20, g7 = 0.10, g8 = -0.20, 
                  g9 = 0.20, g10 = -0.30, g11 = 0.40, g12 = -0.50, g13 = 0.20, g14 = -0.10, g15 = -0.40, g16 = 0.20, 
                  g17 = -0.40, g18 = 0.10, g19 = -0.30, g20 = -0.20, g21 = -0.10, g22 = -0.10, g23 = 0.10, 
                  g24 = -0.10, g25 = 0.10, g26 = -0.10, g27 = 0.10, g28 = 0.10, g29 = 0.10)
    )
  )
  
  # Select the coefficients based on the scenario and DGP
  coefs <- coefficients[[dgp]][[scenario]]
  
  # Generate covariates based on DGP
  if (dgp == "DGP1") {
    B <- rnorm(n, mean = 0, sd = 1)
    Z1 <- rbinom(n, size = 1, prob = plogis(coefs$a0))
    Z2 <- rbinom(n, size = 1, prob = plogis(coefs$b0 + coefs$b1 * B))
    Z3 <- rbinom(n, size = 1, prob = plogis(coefs$c0 + coefs$c1 * B))
    Z4 <- rbinom(n, size = 1, prob = plogis(coefs$d0 + coefs$d1 * B))
    Z5 <- rbinom(n, size = 1, prob = plogis(coefs$e0))
  } else if (dgp == "DGP2") {
    B <- rnorm(n, mean = 0, sd = 1)
    Z1 <- rbinom(n, size = 1, prob = plogis(coefs$a0))
    Z2 <- rbinom(n, size = 1, prob = plogis(coefs$b0 + coefs$b1 * B))
    Z3 <- rbinom(n, size = 1, prob = plogis(coefs$c0 + coefs$c1 * B))
    Z4 <- rnorm(n, mean = coefs$d0 + coefs$d1 * B, sd = 1)
    Z5 <- rnorm(n, mean = coefs$e0, sd = 2)
  } else if (dgp == "DGP3") {
    B <- rnorm(n, mean = 0, sd = 1)
    rho <- matrix(c(1, 0.3, -0.3, 0.3, 0.3, -0.3,
                    0.3, 1, 0.7, 0.3, 0.3, 0.3,
                    -0.3, 0.7, 1, 0.3, 0.3, 0.3,
                    0.3, 0.3, 0.3, 1, 0.7, 0.3,
                    0.3, 0.3, 0.3, 0.7, 1, -0.3,
                    -0.3, 0.3, 0.3, 0.3, -0.3, 1), nrow = 6, ncol = 6, byrow = TRUE)
    lower_triangular <- rho[lower.tri(rho)]
    gaussian_cop <- normalCopula(param = lower_triangular, dim = 6, dispstr = "un")
    sim_data_gaussian <- rCopula(n, gaussian_cop)
    Z1 <- ifelse(sim_data_gaussian[, 1] <= plogis(coefs$a0), 0, 1)
    Z2 <- ifelse(sim_data_gaussian[, 2] <= plogis(coefs$b0 + coefs$b1 * B), 0, 1)
    Z3 <- ifelse(sim_data_gaussian[, 3] <= plogis(coefs$c0 + coefs$c1 * B), 0, 1)
    Z4 <- qnorm(sim_data_gaussian[, 4], mean = coefs$d0 + coefs$d1 * B, sd = 1)
    Z5 <- qnorm(sim_data_gaussian[, 5], mean = coefs$e0, sd = 2)
  } else if (dgp == "DGP4") {
    B <- rnorm(n, mean = 0, sd = 1)
    rho <- matrix(c(1, 0.3, -0.3, 0.3, 0.3, -0.3,
                    0.3, 1, 0.7, 0.3, 0.3, 0.3,
                    -0.3, 0.7, 1, 0.3, 0.3, 0.3,
                    0.3, 0.3, 0.3, 1, 0.7, 0.3,
                    0.3, 0.3, 0.3, 0.7, 1, -0.3,
                    -0.3, 0.3, 0.3, 0.3, -0.3, 1), nrow = 6, ncol = 6, byrow = TRUE)
    lower_triangular <- rho[lower.tri(rho)]
    gaussian_cop <- normalCopula(param = lower_triangular, dim = 6, dispstr = "un")
    sim_data_gaussian <- rCopula(n, gaussian_cop)
    Z1 <- ifelse(sim_data_gaussian[, 1] <= plogis(coefs$a0), 0, 1)
    Z2 <- ifelse(sim_data_gaussian[, 2] <= plogis(coefs$b0 + coefs$b1 * B), 0, 1)
    
    # Categorical variable Z3 with 4 categories and dependency on B
    num_categories <- 4
    softmax <- function(x) exp(x) / sum(exp(x))
    cat_probs <- matrix(NA, nrow = n, ncol = num_categories)
    c0 <- c(coefs$c0_1, coefs$c0_2, coefs$c0_3, coefs$c0_4)
    c1 <- c(coefs$c1_1, coefs$c1_2, coefs$c1_3, coefs$c1_4)
    
    for (i in 1:n) {
      cat_probs[i, ] <- softmax(c0 + c1 * B[i])
    }
    cat_cum_probs <- t(apply(cat_probs, 1, cumsum))
    generate_categorical <- function(cum_probs, sample) {
      for (z in 1:length(cum_probs)) {
        if (sample <= cum_probs[z]) {
          return(z)
        }
      }
      return(length(cum_probs) + 1)
    }
    Z3 <- sapply(1:n, function(z) generate_categorical(cat_cum_probs[z, ], sim_data_gaussian[z, 3]))
    Z3 <- as.factor(Z3)
    
    Z4 <- qnorm(sim_data_gaussian[, 4], mean = coefs$d0 + coefs$d1 * B, sd = 1)
    Z5 <- qnorm(sim_data_gaussian[, 5], mean = coefs$e0, sd = 2)
  } else if (dgp == "DGP5") {
    B <- rnorm(n, mean = 0, sd = 1)
    rho <- matrix(c(1, 0.3, -0.3, 0.3, 0.3, -0.3,
                    0.3, 1, 0.7, 0.3, 0.3, 0.3,
                    -0.3, 0.7, 1, 0.3, 0.3, 0.3,
                    0.3, 0.3, 0.3, 1, 0.7, 0.3,
                    0.3, 0.3, 0.3, 0.7, 1, -0.3,
                    -0.3, 0.3, 0.3, 0.3, -0.3, 1), nrow = 6, ncol = 6, byrow = TRUE)
    lower_triangular <- rho[lower.tri(rho)]
    gaussian_cop <- normalCopula(param = lower_triangular, dim = 6, dispstr = "un")
    sim_data_gaussian <- rCopula(n, gaussian_cop)
    Z1 <- ifelse(sim_data_gaussian[, 1] <= plogis(coefs$a0), 0, 1)
    Z2 <- ifelse(sim_data_gaussian[, 2] <= plogis(coefs$b0 + coefs$b1 * B), 0, 1)
    
    # Categorical variable Z3 with 4 categories and dependency on B
    num_categories <- 4
    softmax <- function(x) exp(x) / sum(exp(x))
    cat_probs <- matrix(NA, nrow = n, ncol = num_categories)
    c0 <- c(coefs$c0_1, coefs$c0_2, coefs$c0_3, coefs$c0_4)
    c1 <- c(coefs$c1_1, coefs$c1_2, coefs$c1_3, coefs$c1_4)
    
    for (i in 1:n) {
      cat_probs[i, ] <- softmax(c0 + c1 * B[i])
    }
    cat_cum_probs <- t(apply(cat_probs, 1, cumsum))
    generate_categorical <- function(cum_probs, sample) {
      for (z in 1:length(cum_probs)) {
        if (sample <= cum_probs[z]) {
          return(z)
        }
      }
      return(length(cum_probs) + 1)
    }
    Z3 <- sapply(1:n, function(z) generate_categorical(cat_cum_probs[z, ], sim_data_gaussian[z, 3]))
    Z3 <- as.factor(Z3)
    
    Z4 <- qnorm(sim_data_gaussian[, 4], mean = coefs$d0 + coefs$d1 * B, sd = 1)
    Z5 <- qnorm(sim_data_gaussian[, 5], mean = coefs$e0, sd = 2)
    
    generate_truncated_normal <- function(input, lower, upper) {
      result <- input
      result[result < lower] <- lower
      result[result > upper] <- upper
      return(result)
    }
    t_B <- generate_truncated_normal(B, -0.99, 0.99)
    shape <- 2 + coefs$s1 * t_B
    rate <- 1 + coefs$r1 * t_B
    Z6 <- qgamma(sim_data_gaussian[, 6], shape = shape, rate = rate)
  }
  
  # Create exposure A
  if (dgp == "DGP4") {
    Z3_dummies <- model.matrix(~ Z3 - 1)
    A <- rbinom(n, size = 1, prob = plogis(coefs$f0 + coefs$f1 * Z1 + coefs$f2 * Z2 + coefs$f3 * Z3_dummies[,1] + coefs$f4 * Z3_dummies[,2] + coefs$f5 * Z3_dummies[,3] + 
                                             coefs$f6 * Z4 + coefs$f7 * Z5 + coefs$f8 * B + 
                                             coefs$f9 * Z1 * Z3_dummies[,1] + coefs$f10 * Z1 * Z3_dummies[,2] + coefs$f11 * Z1 * Z3_dummies[,3] + 
                                             coefs$f12 * Z1 * Z4 + coefs$f13 * Z1 * Z5 + 
                                             coefs$f14 * Z3_dummies[,1] * Z4 + coefs$f15 * Z3_dummies[,2] * Z4 + coefs$f16 * Z3_dummies[,3] * Z4 + 
                                             coefs$f17 * Z3_dummies[,1] * Z5 + coefs$f18 * Z3_dummies[,2] * Z5 + coefs$f19 * Z3_dummies[,3] * Z5 + 
                                             coefs$f20 * Z4 * Z5))
  } else if (dgp == "DGP5") {
    Z3_dummies <- model.matrix(~ Z3 - 1)
    A <- rbinom(n, size = 1, prob = plogis(coefs$f0 + coefs$f1 * Z1 + coefs$f2 * Z2 + coefs$f3 * Z3_dummies[,1] + coefs$f4 * Z3_dummies[,2] + coefs$f5 * Z3_dummies[,3] + 
                                             coefs$f6 * Z4 + coefs$f7 * Z5 + coefs$f8 * Z6 + coefs$f9 * B + 
                                             coefs$f10 * Z1 * Z3_dummies[,1] + coefs$f11 * Z1 * Z3_dummies[,2] + coefs$f12 * Z1 * Z3_dummies[,3] + 
                                             coefs$f13 * Z1 * Z4 + coefs$f14 * Z1 * Z5 + 
                                             coefs$f15 * Z3_dummies[,1] * Z4 + coefs$f16 * Z3_dummies[,2] * Z4 + coefs$f17 * Z3_dummies[,3] * Z4 + 
                                             coefs$f18 * Z3_dummies[,1] * Z5 + coefs$f19 * Z3_dummies[,2] * Z5 + coefs$f20 * Z3_dummies[,3] * Z5 + 
                                             coefs$f21 * Z4 * Z5 + 
                                             coefs$f22 * Z1 * Z6 + coefs$f23 * Z4 * Z6 + coefs$f24 * Z5 * Z6))
  } else {
    A <- rbinom(n, size = 1, prob = plogis(coefs$f0 + coefs$f1 * Z1 + coefs$f2 * Z2 + coefs$f3 * Z3 + coefs$f4 * Z4 + coefs$f5 * Z5 + coefs$f6 * B + 
                                             coefs$f7 * Z1 * Z3 + coefs$f8 * Z1 * Z4 + coefs$f9 * Z1 * Z5 + 
                                             coefs$f10 * Z3 * Z4 + coefs$f11 * Z3 * Z5 + coefs$f12 * Z4 * Z5))
  }
  
  # Create Outcome Y
  if (dgp == "DGP4") {
    Y <- rnorm(n, mean = (coefs$g0 + coefs$g1 * Z1 + coefs$g2 * Z2 + coefs$g3 * Z3_dummies[,1] + coefs$g4 * Z3_dummies[,2] + coefs$g5 * Z3_dummies[,3] + 
                            coefs$g6 * Z4 + coefs$g7 * Z5 + coefs$g8 * A +
                            coefs$g9 * Z1 * Z3_dummies[,1] + coefs$g10 * Z1 * Z3_dummies[,2] + coefs$g11 * Z1 * Z3_dummies[,3] + 
                            coefs$g12 * Z1 * Z4 + coefs$g13 * Z1 * Z5 +
                            coefs$g14 * Z3_dummies[,1] * Z4 + coefs$g15 * Z3_dummies[,2] * Z4 + coefs$g16 * Z3_dummies[,3] * Z4 + 
                            coefs$g17 * Z3_dummies[,1] * Z5 + coefs$g18 * Z3_dummies[,2] * Z5 + coefs$g19 * Z3_dummies[,3] * Z5 + 
                            coefs$g20 * Z4 * Z5 + 
                            coefs$g21 * Z1 * Z4 * Z2 + coefs$g22 * Z1 * Z5 * Z2 + coefs$g23 * Z1 * Z4 * Z5 + coefs$g24 * Z4 * Z5 * Z2 + coefs$g25 * Z1 * Z4 * Z5 * Z2), sd = 1)
  } else if (dgp == "DGP5") {
    Y <- rnorm(n, mean = (coefs$g0 + coefs$g1 * Z1 + coefs$g2 * Z2 + coefs$g3 * Z3_dummies[,1] + coefs$g4 * Z3_dummies[,2] + coefs$g5 * Z3_dummies[,3] + 
                            coefs$g6 * Z4 + coefs$g7 * Z5 + coefs$g8 * Z6 + coefs$g9 * A +
                            coefs$g10 * Z1 * Z3_dummies[,1] + coefs$g11 * Z1 * Z3_dummies[,2] + coefs$g12 * Z1 * Z3_dummies[,3] + 
                            coefs$g13 * Z1 * Z4 + coefs$g14 * Z1 * Z5 + 
                            coefs$g15 * Z3_dummies[,1] * Z4 + coefs$g16 * Z3_dummies[,2] * Z4 + coefs$g17 * Z3_dummies[,3] * Z4 + 
                            coefs$g18 * Z3_dummies[,1] * Z5 + coefs$g19 * Z3_dummies[,2] * Z5 + coefs$g20 * Z3_dummies[,3] * Z5 + 
                            coefs$g21 * Z4 * Z5 + 
                            coefs$g22 * Z1 * Z6 + coefs$g23 * Z4 * Z6 + coefs$g24 * Z5 * Z6 +
                            coefs$g25 * Z1 * Z4 * Z6 + coefs$g26 * Z1 * Z5 * Z6 + coefs$g27 * Z1 * Z4 * Z5 + coefs$g28 * Z4 * Z5 * Z6 + coefs$g29 * Z1 * Z4 * Z5 * Z6), sd = 1)
  } else {
    Y <- rnorm(n, mean = (coefs$g0 + coefs$g1 * Z1 + coefs$g2 * Z2 + coefs$g3 * Z3 + coefs$g4 * Z4 + coefs$g5 * Z5 + coefs$g6 * A +
                            coefs$g7 * Z1 * Z3 + coefs$g8 * Z1 * Z4 + coefs$g9 * Z1 * Z5 + 
                            coefs$g10 * Z3 * Z4 + coefs$g11 * Z3 * Z5 + coefs$g12 * Z4 * Z5 + 
                            coefs$g13 * Z1 * Z2 * Z4 + coefs$g14 * Z1 * Z2 * Z5 + coefs$g15 * Z1 * Z4 * Z5 + coefs$g16 * Z2 * Z4 * Z5 + coefs$g17 * Z1 * Z2 * Z4 * Z5), sd = 1)
  }
  
  if (dgp == "DGP5") {
    # Return data.frame
    return(data.frame(B, Z1, Z2, Z3, Z4, Z5, Z6, A, Y))
  } else {
    return(data.frame(B, Z1, Z2, Z3, Z4, Z5, A, Y)) 
  }
}


# Simulate datasets for different scenarios and DGPs
simulateDatasets <- function(scenario, dgp, n_datasets = 1000, n = 2000) {
  data_list <- vector("list", length = n_datasets)
  for (i in 1:n_datasets) {
    data_list[[i]] <- generateData(n = n, scenario = scenario, dgp = dgp)
  }
  data_list
}



