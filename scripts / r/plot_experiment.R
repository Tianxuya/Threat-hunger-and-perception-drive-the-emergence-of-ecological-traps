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

figures_dir <- Sys.getenv("EXPERIMENT_FIGURES_DIR", unset = "")
if (!nzchar(figures_dir)) figures_dir <- file.path(script_dir, "..", "..", "figures", "experiment")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)


## ============================================================
## 1. Organize raw data
## ============================================================

total_worms <- 10

dat <- get_patch_data()

# Mutate factor classifications and compute metrics needed for validation checks
dat <- dat %>%
  mutate(
    starvation = factor(
      starvation, 
      levels = c("HS", "FS"), 
      labels = c("Low Hunger", "High Hunger")
    ),
    toxin = factor(toxin, levels = c("No toxin", "Toxin")),
    perception = factor(perception, levels = c("Local", "Nonlocal")),
    total_food = low + medium + high,
    elsewhere = total_worms - total_food
  )

## ============================================================
## 2. Convert to long format: Low, Medium, High food patches
## ============================================================

patch_long <- dat %>%
  select(starvation, toxin, perception, rep, low, medium, high) %>%
  pivot_longer(
    cols = c(low, medium, high),
    names_to = "patch",
    values_to = "abundance"
  ) %>%
  mutate(
    proportion = abundance / total_worms,
    patch = factor(
      patch,
      levels = c("low", "medium", "high"),
      labels = c("Low food patch", "Medium food patch", "High food patch")
    )
  )

# Calculate summary statistics for the bar plot columns
patch_sum <- patch_long %>%
  group_by(perception, toxin, starvation, patch) %>%
  summarise(
    mean_prop = mean(proportion),
    se = sd(proportion) / sqrt(n()),
    .groups = "drop"
  )

## ============================================================
## 3. Revised plot generation
## Rows = HS / FS
## X-axis = Local / Nonlocal
## Fill = Low / Medium / High food patch
## Columns = No toxin / Toxin
## ============================================================

p_patch_revised <- ggplot(
  patch_sum,
  aes(x = perception, y = mean_prop, fill = patch)
) +
  # Bar columns layer
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.25
  ) +
  # Standard Error Bars (bounded at 0)
  geom_errorbar(
    aes(
      ymin = pmax(mean_prop - se, 0),
      ymax = mean_prop + se
    ),
    position = position_dodge(width = 0.75),
    width = 0.18,
    linewidth = 0.4
  ) +
  # Raw replication data point overlay
  geom_point(
    data = patch_long,
    aes(x = perception, y = proportion, fill = patch),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.75
    ),
    shape = 21,
    color = "black",
    size = 1.8,
    alpha = 0.65,
    show.legend = FALSE
  ) +
  facet_grid(starvation ~ toxin) +
  scale_fill_manual(
    values = c(
      "Low food patch"    = "#c7e9c0",
      "Medium food patch" = "#74c476",
      "High food patch"   = "#238b45"
    )
  ) +
  # Set upper y-axis limit to 0.8 so data points/error bars don't get cut off
  scale_y_continuous(
    breaks = seq(0, 0.6, 0.2),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  coord_cartesian(ylim = c(0, 0.6)) +
  labs(
    x = NULL,
    y = "Proportion of abundance",
    fill = "Food patch"
  ) +
  theme_classic(base_size = 14, base_family = "Helvetica") +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 0, hjust = 0.5),
    legend.position  = "top",
    legend.title     = element_text(face = "bold"),
    panel.spacing    = unit(1, "lines")
  )

# Display final figure output
print(p_patch_revised)

## ============================================================
## 4. File Saving Outputs
## ============================================================

ggsave(
  file.path(figures_dir, "worm_patch_level_revised_grouping.png"),
  p_patch_revised,
  width = 10,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(figures_dir, "worm_patch_level_revised_grouping.pdf"),
  p_patch_revised,
  width = 10,
  height = 6.5,
  device = "pdf",
  bg = "white"
)

## ============================================================
## 5. Automated Validation Verification Checks
## ============================================================

stopifnot(nrow(dat) == 40)
stopifnot(all(dat$total_food <= total_worms))
stopifnot(all(dat$elsewhere >= 0))
