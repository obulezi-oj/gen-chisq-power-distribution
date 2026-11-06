# Generalized Chi-Squared Power Distribution ($\text{Gen-}\chi^2(n, \alpha, \sigma)$)

[![R Language](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![License: MIT](https.img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains R scripts for implementing, simulating, and evaluating the **Generalized Chi-Squared Power Distribution** ($\text{Gen-}\chi^2(n, \alpha, \sigma)$). The code supports theoretical density and hazard function plotting, parameter estimation via Maximum Likelihood Estimation (MLE) and MCMC Bayesian inference, Monte Carlo simulation studies, and empirical model benchmarks against Weibull and Chi distributions.

---

## 📁 Repository Structure

| File | Description |
| :--- | :--- |
| `plot of the chi-square script.R` | Computes PDF, SF, and HRF curves for various parameter settings ($n, \alpha, \sigma$) and renders synchronization plots. |
| `MLE and Bayes simulation script.R` | Monte Carlo & MCMC Metropolis-Hastings simulation comparing MLE and Bayesian estimation metrics (Bias, MSE, 95% Coverage Probability). |
| `real data script for gen-chi-squared power.R` | Fits the proposed model and benchmark distributions (Weibull, Chi) on real data, exporting diagnostic panels and Goodness-of-Fit tables. |

---

## 🧮 Model Overview

The probability density function (PDF) of the Generalized $\chi^2$ Power distribution for a random variable $Q > 0$ is defined as:

$$f_Q(q; n, \alpha, \sigma) = \frac{2^{1 - n/2}}{\Gamma(n/2)(1 + \alpha)\sigma^{\frac{n}{1+\alpha}}} \, q^{\frac{n}{1+\alpha} - 1} \exp\left\{ -\frac{1}{2} \left(\frac{q}{\sigma}\right)^{\frac{2}{1+\alpha}} \right\}$$

where:
- $n > 0$ is the degrees-of-freedom parameter.
- $\alpha > -1$ is the power shape parameter.
- $\sigma > 0$ is the scale parameter.

---

## 🚀 Getting Started

### Prerequisites

Ensure you have R (version 4.0 or higher) installed along with the required CRAN packages:

```R
install.packages(c("ggplot2", "gridExtra", "fitdistrplus", "stats", "goftest"))
