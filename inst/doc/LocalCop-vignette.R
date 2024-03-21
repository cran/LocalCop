params <-
list(do_calc = FALSE)

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.width = 7
  )
library(kableExtra)
library(readr)

## ----libs, message=FALSE------------------------------------------------------
# load required packages
library(LocalCop)   # for conditional copula modeling
library(VineCopula) # for simulating copula data 
library(dplyr)      # for data manipulations
library(tidyr)      # and
library(ggplot2)    # plotting

## ----calib, echo = FALSE------------------------------------------------------
tab <- read_csv("copula_table.csv", show_col_types = FALSE)
tab_cap <- "Copula families implemented in **LocalCop**."
kableExtra::kbl(tab, booktabs = TRUE, caption = tab_cap)

## ----dgm, message=FALSE, fig.height=3, fig.cap="Conditional Kendall $\\tau$ for Gumbel copula under DGM."----
# simulate copula data given a covariate
set.seed(2024) 
family <- 4  # Gumbel Copula 
n <- 1000    # number of observations
x <- sort(runif(n))  # covariate values
eta_fun <- function(x) sin(6*pi*x) # calibration function

# simulate data
eta_true <- eta_fun(x)
par_true <- BiCopEta2Par(family = family, eta = eta_true) # copula parameter
udata <- VineCopula::BiCopSim(n, family = family, par = par_true$par)

# plot tau(x)
tibble(
  x = seq(0, 1, len = 100),
) %>%
  mutate(
    tau = BiCopEta2Tau(family, eta = eta_fun(x))
  ) %>%
  ggplot(aes(x = x, y = tau)) +
  geom_line() +
  ylim(c(0, 1)) + 
  xlab(expression(x)) + ylab(expression(tau(x))) + 
  theme_bw()

## ----local1, message=FALSE----------------------------------------------------
x0 <- 0.1
band <- 0.1
degree <- 1
kernel <- KernEpa # Epanichov kernel (default value)
fit1 <- CondiCopLocFit(
  u1 = udata[,1], u2 = udata[,2], x = x,
  x0 = x0,
  family = family, 
  degree = degree,
  kernel = kernel,
  band = band
)
fit1

## ----localseq, message=FALSE, fig.height=3.5, warning=FALSE, fig.cap="Conditional Kendall $\\tau$ estimates under the Gumbel copula using bandwidth $h=0.1$."----
x0 <- seq(0, 1, by=0.02)
fitseq <- CondiCopLocFit(
  u1 = udata[,1], u2 = udata[,2], x = x, 
  x0 = x0,
  family = family, 
  degree = degree,
  kernel = kernel,
  band = band
)

# plot true vs fitted tau
legend_names <- c(expression(tau(x)),
                  expression(hat(tau)(x)))
tibble(
  x = x0,
  True = BiCopEta2Tau(family, eta = eta_fun(x0)),
  Fitted = BiCopEta2Tau(fitseq$eta, family= family)
) %>%
  pivot_longer(True:Fitted,
               names_to = "Type", values_to = "y") %>%
  mutate(
    Type = factor(Type, levels = c("True", "Fitted"))
  ) %>%
  ggplot(aes(x = x, y = y, group = Type)) +
  geom_line(aes(color = Type, linetype = Type)) + 
  geom_point(aes(shape = Type, color = Type)) +
  scale_color_manual(values = c("black", "red"), labels = legend_names) +
  scale_shape_manual(values = c(NA, 16), labels = legend_names) +
  scale_linetype_manual(values = c("solid", NA), labels = legend_names) +
  ylim(c(0, 1)) + 
  xlab(expression(x)) + ylab(expression("Kendall "*tau)) +
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

## ----select1-precalc----------------------------------------------------------
bandset <- c(0.1, 0.2, 0.4, 0.8, 1) # bandwidth set
famset <- c(1, 2, 3, 4, 5) # family set
n_loo <- 100  # number of leave-one-out observations in CV likelihood calculation

## ----select1-calc, eval = params$do_calc, message=FALSE-----------------------
#  system.time({
#    cvselect1 <- CondiCopSelect(
#      u1 = udata[,1], u2 = udata[,2], x = x,
#      family = famset,
#      degree = degree,
#      kernel = kernel,
#      band = bandset,
#      xind = n_loo
#    )
#  })

## ----select1-save, eval = params$do_calc, echo = FALSE------------------------
#  saveRDS(cvselect1, file = "cvselect1.rds")

## ----select1-load, eval = !params$do_calc, echo = FALSE-----------------------
cvselect1 <- readRDS("cvselect1.rds")

## ----select1------------------------------------------------------------------
cv_res1 <- cvselect1$cv
cv_res1

## ----select2-calc, eval = params$do_calc, message=FALSE-----------------------
#  system.time({
#    cvselect2 <- CondiCopSelect(
#      u1 = udata[,1], u2 = udata[,2], x = x,
#      family = famset,
#      degree = degree,
#      kernel = kernel,
#      band = bandset,
#      xind = nrow(udata)
#    )
#  })

## ----select2-save, eval = params$do_calc, echo = FALSE------------------------
#  saveRDS(cvselect2, file = "cvselect2.rds")

## ----select2-load, eval = !params$do_calc, echo = FALSE-----------------------
cvselect2 <- readRDS("cvselect2.rds")

## ----select2------------------------------------------------------------------
cv_res2 <- cvselect2$cv
cv_res2

## ----cvplot, message=FALSE, fig.height=3.5, fig.cap="Cross-validated likelihood for copula and bandwidth selection based on a subset and the full sample."----
fam_names <- c("Gaussian", "Student", "Clayton", "Gumbel", "Frank")
bind_rows(as_tibble(cvselect1$cv) %>%
          mutate(n = n_loo),
          as_tibble(cvselect2$cv) %>%
          mutate(n = nrow(udata))) %>%
  mutate(
    family = factor(family, levels = c(1,2,3,4,5),
                    labels = fam_names),
    Bandwidth = factor(band),
    n = factor(paste0("n = ", n))
  ) %>%
  ggplot(aes(x = family, y  = cv, fill = Bandwidth)) +
  geom_bar(stat="identity", position = "dodge") +
  facet_grid(. ~ n) +
  scale_fill_brewer(palette="Blues", direction=-1) +
  xlab("") + ylab("CV Likelihood") +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

## ----localseq2, message=FALSE, fig.align='center', fig.width=7, fig.height=5, fig.cap="Conditional Kendall $\\tau$ estimates under different copula families."----
x0 <- seq(0, 1, by=0.01)
tau_est <- sapply(1:5, function(fam_id) {
  fit <- CondiCopLocFit(
    u1 = udata[,1], u2 = udata[,2], x = x,
    x0 = x0,
    family = fam_id, 
    degree = degree,
    kernel = kernel,
    band = band)
  BiCopEta2Tau(fit$eta, family=fam_id)
})

colnames(tau_est) <- fam_names

as_tibble(tau_est) %>%
  mutate(
    x = x0,
    True = BiCopEta2Tau(family, eta = eta_fun(x0))
  ) %>%
  pivot_longer(!x, names_to = "Type", values_to = "tau") %>%
  mutate(
    Type = factor(Type, levels = c("True", fam_names))
  ) %>%
  ggplot(aes(x = x, y = tau, group = Type)) +
  geom_line(aes(col = Type, linewidth = Type)) +
  scale_color_manual(
    values = c("black", "red", "blue", "brown", "orange", "green4")
  ) + 
  scale_linewidth_manual(values = c(1, rep(.5, 5))) + 
  ylim(c(0, 1)) + 
  xlab(expression(x)) + ylab(expression(tau(x))) +
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
  )

