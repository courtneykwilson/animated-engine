---
Title: "Descriptive Epidemiology and Measures of Association: Risk Factors for Low Birth Weight"
Author: "Courtney Wilson"
Date: "29 July 2026"
output:
  html_document:
    toc: true
    toc_float: true
    theme: flatly
---


## 1. Background

Low birth weight (LBW), defined as birth weight below 2500g, is a key
public health indicator associated with increased infant morbidity and
mortality. This report uses the `birthwt` dataset (Hosmer & Lemeshow,
1989), collected at Baystate Medical Center, Springfield, MA, to explore
maternal risk factors associated with LBW using standard descriptive
epidemiology and measures of association.

**Research question:** Is maternal smoking during pregnancy associated
with low birth weight, and how does this association compare across
measures (risk ratio vs. odds ratio)?

## 2. Data preparation

```{r load-data}
data(birthwt, package = "MASS")

bw <- birthwt %>%
  mutate(
    low     = factor(low, levels = c(0, 1), labels = c("Normal weight", "Low birth weight")),
    smoke   = factor(smoke, levels = c(0, 1), labels = c("Non-smoker", "Smoker")),
    race    = factor(race, levels = c(1, 2, 3), labels = c("White", "Black", "Other")),
    ht      = factor(ht, levels = c(0, 1), labels = c("No hypertension", "Hypertension")),
    ui      = factor(ui, levels = c(0, 1), labels = c("No irritability", "Uterine irritability")),
    ptl_any = factor(ifelse(ptl > 0, 1, 0), levels = c(0, 1),
                      labels = c("No prior preterm labor", "Prior preterm labor"))
  )

glimpse(bw)
```

The dataset contains `r nrow(bw)` mothers, of whom `r sum(bw$low == "Low birth weight")`
(`r round(100 * mean(bw$low == "Low birth weight"), 1)`%) delivered a
low birth weight infant.

## 3. Table 1: descriptive summary by outcome status

```{r table1}
table1 <- bw %>%
  group_by(low) %>%
  summarise(
    n              = n(),
    mean_age       = round(mean(age), 1),
    sd_age         = round(sd(age), 1),
    mean_lwt       = round(mean(lwt), 1),
    pct_smoker     = round(100 * mean(smoke == "Smoker"), 1),
    pct_hypertension = round(100 * mean(ht == "Hypertension"), 1),
    pct_uterine_irr   = round(100 * mean(ui == "Uterine irritability"), 1),
    pct_prior_ptl     = round(100 * mean(ptl_any == "Prior preterm labor"), 1)
  )

kable(table1, caption = "Maternal characteristics by birth weight outcome")
```

## 4. Visualizing the exposure–outcome relationship

```{r plots, fig.width=7, fig.height=4}
ggplot(bw, aes(x = smoke, fill = low)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Proportion of low birth weight by smoking status",
       x = "Maternal smoking status", y = "Proportion", fill = "Outcome") +
  theme_minimal()

ggplot(bw, aes(x = race, y = bwt, fill = race)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Birth weight (g) by maternal race",
       x = "Race", y = "Birth weight (g)") +
  theme_minimal()
```

## 5. 2x2 table: smoking and low birth weight

```{r two-by-two}
tab <- table(Exposure = bw$smoke, Outcome = bw$low)
tab
```

## 6. Measures of association

Computing the risk ratio (RR) and odds ratio (OR) for low birth weight
comparing smokers to non-smokers, along with 95% confidence intervals
using standard log-scale formulas (Rothman, *Modern Epidemiology*).

```{r measures}
a <- tab["Smoker", "Low birth weight"]
b <- tab["Smoker", "Normal weight"]
c <- tab["Non-smoker", "Low birth weight"]
d <- tab["Non-smoker", "Normal weight"]

risk_exposed   <- a / (a + b)
risk_unexposed <- c / (c + d)
rr <- risk_exposed / risk_unexposed

se_log_rr <- sqrt((1 / a) - (1 / (a + b)) + (1 / c) - (1 / (c + d)))
rr_ci <- exp(log(rr) + c(-1, 1) * 1.96 * se_log_rr)

or <- (a * d) / (b * c)
se_log_or <- sqrt(1/a + 1/b + 1/c + 1/d)
or_ci <- exp(log(or) + c(-1, 1) * 1.96 * se_log_or)

results <- tibble::tibble(
  Measure = c("Risk in smokers", "Risk in non-smokers", "Risk Ratio", "Odds Ratio"),
  Estimate = round(c(risk_exposed, risk_unexposed, rr, or), 3),
  `Lower 95% CI` = c(NA, NA, round(rr_ci[1], 3), round(or_ci[1], 3)),
  `Upper 95% CI` = c(NA, NA, round(rr_ci[2], 3), round(or_ci[2], 3))
)

kable(results, caption = "Risk and odds ratios for low birth weight, smokers vs. non-smokers")
```

### Statistical test of association

```{r chisq}
chisq.test(tab)
fisher.test(tab)
```

## 7. Interpretation

Mothers who smoked during pregnancy had a higher proportion of low
birth weight deliveries than non-smokers. The risk ratio indicates the
relative increase in probability of LBW among smokers, while the odds
ratio commonly reported because it generalizes to case-control designs
and multivariable logistic regression is slightly larger in magnitude,
as expected when the outcome is not rare (the "rare disease assumption"
under which OR approximates RR does not hold well here, since LBW
prevalence is roughly 30% in this sample). The chi-square and Fisher's
exact tests both assess whether this association is unlikely to be due
to chance alone; the confidence intervals additionally quantify the
precision of the estimated effect. This descriptive analysis motivates
the multivariable logistic regression in the next report, which adjusts
for other maternal risk factors (age, weight, hypertension, uterine
irritability, prior preterm labor) that may confound the smoking–LBW
relationship.

## Session info

```{r session-info}
sessionInfo()
```
