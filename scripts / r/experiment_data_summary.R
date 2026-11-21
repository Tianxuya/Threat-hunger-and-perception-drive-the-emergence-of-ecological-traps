## ============================================================
## Full statistical analysis code
## Worm food-patch occupancy: count, proportion, and percentage
## ============================================================

## -----------------------------
## 0. Load packages
## -----------------------------

library(tidyverse)

# Rscript resolves this script independent of the current working directory.
script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_flag)) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_flag[[1]])))
} else {
  script_dir <- file.path(getwd(), "scripts", "r")
  if (!file.exists(file.path(script_dir, "experiment_patch_data.R"))) script_dir <- getwd()
}
source(file.path(script_dir, "experiment_patch_data.R"))

results_dir <- Sys.getenv("EXPERIMENT_RESULTS_DIR", unset = "")
if (!nzchar(results_dir)) results_dir <- file.path(script_dir, "..", "..", "results", "experiment")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(car)

library(emmeans)

## -----------------------------
## 1. Define raw data
## -----------------------------
## Each trial contains 10 worms.
## low, medium, high = number of worms in each food patch.

total_worms <- 10

dat <- get_patch_data()
dat$starvation <- factor(dat$starvation, levels = c("HS", "FS"),
                          labels = c("Low Hunger", "High Hunger"))

## -----------------------------
## 2. Prepare data
## -----------------------------

dat <- dat %>%
  mutate(
    starvation = factor(starvation, levels = c("Low Hunger", "High Hunger")),
    toxin = factor(toxin, levels = c("No toxin", "Toxin")),
    perception = factor(perception, levels = c("Local", "Nonlocal")),
    
    total_food = low + medium + high,
    elsewhere = total_worms - total_food,
    
    prop_food = total_food / total_worms,
    percent_food = 100 * total_food / total_worms,
    
    not_low = total_worms - low,
    not_medium = total_worms - medium,
    not_high = total_worms - high
  )

## Check data
print(dat)

dat %>%
  count(starvation, toxin, perception)

## -----------------------------
## 3. Summary table
## -----------------------------

summary_table <- dat %>%
  group_by(starvation, toxin, perception) %>%
  summarise(
    n_trials = n(),
    mean_total_food = mean(total_food),
    sd_total_food = sd(total_food),
    mean_prop_food = mean(prop_food),
    sd_prop_food = sd(prop_food),
    mean_percent_food = mean(percent_food),
    sd_percent_food = sd(percent_food),
    mean_low = mean(low),
    mean_medium = mean(medium),
    mean_high = mean(high),
    .groups = "drop"
  )

print(summary_table)

## ============================================================
## 4. Main correct analysis:
##    Quasibinomial GLM for total food-patch occupancy
## ============================================================
## Response = number of worms in food patches vs outside patches.
## This is the main model to report.

m_food <- glm(
  cbind(total_food, elsewhere) ~ starvation * perception * toxin,
  family = quasibinomial(link = "logit"),
  data = dat
)

summary(m_food)

food_anova <- Anova(m_food, type = 2, test.statistic = "F")
print(food_anova)

emm_food <- emmeans(
  m_food,
  ~ starvation * perception * toxin,
  type = "response"
)

print(emm_food)

## Hunger effect within each perception and toxin condition
hunger_pairs_food <- pairs(
  emmeans(m_food, ~ starvation | perception * toxin, type = "response"),
  adjust = "BH"
)

print(hunger_pairs_food)

## Perception effect within each hunger and toxin condition
perception_pairs_food <- pairs(
  emmeans(m_food, ~ perception | starvation * toxin, type = "response"),
  adjust = "BH"
)

print(perception_pairs_food)

## Toxin effect within each hunger and perception condition
toxin_pairs_food <- pairs(
  emmeans(m_food, ~ toxin | starvation * perception, type = "response"),
  adjust = "BH"
)

print(toxin_pairs_food)

## ============================================================
## 5. Simplified quasibinomial model without interactions
## ============================================================
## Use this if interaction terms are not significant.

m_food_simple <- glm(
  cbind(total_food, elsewhere) ~ starvation + perception + toxin,
  family = quasibinomial(link = "logit"),
  data = dat
)

summary(m_food_simple)

simple_anova <- Anova(m_food_simple, type = 2, test.statistic = "F")
print(simple_anova)

emm_food_simple <- emmeans(
  m_food_simple,
  ~ starvation + perception + toxin,
  type = "response"
)

print(emm_food_simple)

## Main-effect pairwise comparisons
simple_hunger_pairs <- pairs(
  emmeans(m_food_simple, ~ starvation, type = "response"),
  adjust = "BH"
)

simple_perception_pairs <- pairs(
  emmeans(m_food_simple, ~ perception, type = "response"),
  adjust = "BH"
)

simple_toxin_pairs <- pairs(
  emmeans(m_food_simple, ~ toxin, type = "response"),
  adjust = "BH"
)

print(simple_hunger_pairs)
print(simple_perception_pairs)
print(simple_toxin_pairs)

## ============================================================
## 6. Exploratory percentage-based analysis
## ============================================================
## This uses percentage values such as 30 instead of 0.30.
## This is only a robustness/exploratory check.
## Do not use it to pretend each trial has 100 worms.

m_percent <- lm(
  percent_food ~ starvation * perception * toxin,
  data = dat
)

summary(m_percent)
anova(m_percent)

percent_anova <- Anova(m_percent, type = 2)
print(percent_anova)

emm_percent <- emmeans(
  m_percent,
  ~ starvation * perception * toxin
)

print(emm_percent)

## Hunger comparisons using percentage model
percent_hunger_pairs <- pairs(
  emmeans(m_percent, ~ starvation | perception * toxin),
  adjust = "BH"
)

print(percent_hunger_pairs)

## Perception comparisons using percentage model
percent_perception_pairs <- pairs(
  emmeans(m_percent, ~ perception | starvation * toxin),
  adjust = "BH"
)

print(percent_perception_pairs)

## Toxin comparisons using percentage model
percent_toxin_pairs <- pairs(
  emmeans(m_percent, ~ toxin | starvation * perception),
  adjust = "BH"
)

print(percent_toxin_pairs)

## -----------------------------
## 6.1 Simplified percentage model
## -----------------------------

m_percent_simple <- lm(
  percent_food ~ starvation + perception + toxin,
  data = dat
)

summary(m_percent_simple)
anova(m_percent_simple)

percent_simple_anova <- Anova(m_percent_simple, type = 2)
print(percent_simple_anova)

emm_percent_simple <- emmeans(
  m_percent_simple,
  ~ starvation + perception + toxin
)

print(emm_percent_simple)

percent_simple_hunger_pairs <- pairs(
  emmeans(m_percent_simple, ~ starvation),
  adjust = "BH"
)

percent_simple_perception_pairs <- pairs(
  emmeans(m_percent_simple, ~ perception),
  adjust = "BH"
)

percent_simple_toxin_pairs <- pairs(
  emmeans(m_percent_simple, ~ toxin),
  adjust = "BH"
)

print(percent_simple_hunger_pairs)
print(percent_simple_perception_pairs)
print(percent_simple_toxin_pairs)

## ============================================================
## 7. Count-based robustness analysis
## ============================================================
## Since percent_food = total_food * 10,
## this should give the same p-values as the percentage linear model.

m_count_lm <- lm(
  total_food ~ starvation * perception * toxin,
  data = dat
)

summary(m_count_lm)
anova(m_count_lm)

count_lm_anova <- Anova(m_count_lm, type = 2)
print(count_lm_anova)

## ============================================================
## 8. Quasipoisson count model
## ============================================================
## This treats total_food as a count.
## It is a robustness check, but less appropriate than quasibinomial
## because total_food is bounded by 10.

m_count_qp <- glm(
  total_food ~ starvation * perception * toxin,
  family = quasipoisson(link = "log"),
  data = dat
)

summary(m_count_qp)

count_qp_anova <- Anova(m_count_qp, type = 2, test.statistic = "F")
print(count_qp_anova)

emm_count_qp <- emmeans(
  m_count_qp,
  ~ starvation * perception * toxin,
  type = "response"
)

print(emm_count_qp)

pairs(
  emmeans(m_count_qp, ~ starvation | perception * toxin, type = "response"),
  adjust = "BH"
)

pairs(
  emmeans(m_count_qp, ~ perception | starvation * toxin, type = "response"),
  adjust = "BH"
)

## ============================================================
## 9. Risk-taking analysis:
##    H--H patch selection
## ============================================================

m_high <- glm(
  cbind(high, not_high) ~ starvation * perception * toxin,
  family = quasibinomial(link = "logit"),
  data = dat
)

summary(m_high)

high_anova <- Anova(m_high, type = 2, test.statistic = "F")
print(high_anova)

emm_high <- emmeans(
  m_high,
  ~ starvation * perception * toxin,
  type = "response"
)

print(emm_high)

hunger_pairs_high <- pairs(
  emmeans(m_high, ~ starvation | perception * toxin, type = "response"),
  adjust = "BH"
)

perception_pairs_high <- pairs(
  emmeans(m_high, ~ perception | starvation * toxin, type = "response"),
  adjust = "BH"
)

print(hunger_pairs_high)
print(perception_pairs_high)

## ============================================================
## 10. Patch-specific analysis
## ============================================================

patch_long <- dat %>%
  select(rep, starvation, toxin, perception, low, medium, high) %>%
  pivot_longer(
    cols = c(low, medium, high),
    names_to = "patch",
    values_to = "count"
  ) %>%
  mutate(
    patch = factor(
      patch,
      levels = c("low", "medium", "high"),
      labels = c("L--L patch", "M--M patch", "H--H patch")
    ),
    not_patch = total_worms - count
  )

m_patch <- glm(
  cbind(count, not_patch) ~ patch * starvation * perception * toxin,
  family = quasibinomial(link = "logit"),
  data = patch_long
)

summary(m_patch)

patch_anova <- Anova(m_patch, type = 2, test.statistic = "F")
print(patch_anova)

emm_patch <- emmeans(
  m_patch,
  ~ patch | starvation * perception * toxin,
  type = "response"
)

print(emm_patch)

patch_pairs <- pairs(emm_patch, adjust = "BH")
print(patch_pairs)


## ============================================================
## 11. Toxin-only analysis:
##     Hunger × Perception interaction
## ============================================================

tox_dat <- dat %>%
  filter(toxin == "Toxin")

m_toxin_interaction <- glm(
  cbind(total_food, elsewhere) ~ starvation * perception,
  family = quasibinomial(link = "logit"),
  data = tox_dat
)

summary(m_toxin_interaction)

toxin_interaction_anova <- Anova(
  m_toxin_interaction,
  type = 2,
  test.statistic = "F"
)

print(toxin_interaction_anova)


## ------------------------------------------------------------
## Perception effect at each hunger level
## ------------------------------------------------------------

perception_under_toxin <- pairs(
  emmeans(
    m_toxin_interaction,
    ~ perception | starvation,
    type = "response"
  ),
  adjust = "BH"
)

print(perception_under_toxin)


## ------------------------------------------------------------
## Hunger effect at each perception level
## ------------------------------------------------------------

hunger_under_toxin <- pairs(
  emmeans(
    m_toxin_interaction,
    ~ starvation | perception,
    type = "response"
  ),
  adjust = "BH"
)

print(hunger_under_toxin)


## ============================================================
## 12. Save outputs
## ============================================================
local <- c(2, 1, 0, 1, 3)
nonlocal <- c(5, 3, 1, 5, 2)

n <- 10

elogit_local <-
  log((local + 0.5) / (n - local + 0.5))

elogit_nonlocal <-
  log((nonlocal + 0.5) / (n - nonlocal + 0.5))

elogit_local
elogit_nonlocal

t.test(
  elogit_nonlocal,
  elogit_local,
  alternative = "greater",
  var.equal = FALSE
)

## ============================================================
## 11. Save outputs
## ============================================================

write.csv(summary_table, file.path(results_dir, "summary_table_patch_occupancy.csv"), row.names = FALSE)

sink(file.path(results_dir, "statistical_results_all_models.txt"))

cat("============================================================\n")
cat("Summary table\n")
cat("============================================================\n")
print(summary_table)

cat("\n============================================================\n")
cat("Main quasibinomial model: total food-patch occupancy\n")
cat("============================================================\n")
print(summary(m_food))
print(food_anova)
print(emm_food)
print(hunger_pairs_food)
print(perception_pairs_food)
print(toxin_pairs_food)

cat("\n============================================================\n")
cat("Simplified quasibinomial model\n")
cat("============================================================\n")
print(summary(m_food_simple))
print(simple_anova)
print(emm_food_simple)
print(simple_hunger_pairs)
print(simple_perception_pairs)
print(simple_toxin_pairs)

cat("\n============================================================\n")
cat("Exploratory percentage model\n")
cat("============================================================\n")
print(summary(m_percent))
print(anova(m_percent))
print(percent_anova)
print(emm_percent)
print(percent_hunger_pairs)
print(percent_perception_pairs)
print(percent_toxin_pairs)

cat("\n============================================================\n")
cat("Simplified percentage model\n")
cat("============================================================\n")
print(summary(m_percent_simple))
print(anova(m_percent_simple))
print(percent_simple_anova)
print(emm_percent_simple)
print(percent_simple_hunger_pairs)
print(percent_simple_perception_pairs)
print(percent_simple_toxin_pairs)

cat("\n============================================================\n")
cat("Count linear model\n")
cat("============================================================\n")
print(summary(m_count_lm))
print(anova(m_count_lm))
print(count_lm_anova)

cat("\n============================================================\n")
cat("Quasipoisson count model\n")
cat("============================================================\n")
print(summary(m_count_qp))
print(count_qp_anova)
print(emm_count_qp)

cat("\n============================================================\n")
cat("High food patch occupancy model\n")
cat("============================================================\n")
print(summary(m_high))
print(high_anova)
print(emm_high)
print(hunger_pairs_high)
print(perception_pairs_high)

cat("\n============================================================\n")
cat("Patch-specific model\n")
cat("============================================================\n")
print(summary(m_patch))
print(patch_anova)
print(emm_patch)
print(patch_pairs)

cat("\n============================================================\n")
cat("Toxin-only analysis: Hunger x Perception\n")
cat("============================================================\n")

print(summary(m_toxin_interaction))
print(toxin_interaction_anova)

cat("\nPerception effect within each hunger level\n")
print(perception_under_toxin)

cat("\nHunger effect within each perception level\n")
print(hunger_under_toxin)


sink()
