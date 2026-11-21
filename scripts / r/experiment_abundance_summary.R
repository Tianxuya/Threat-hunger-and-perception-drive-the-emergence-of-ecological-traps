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

# Mutate factor classifications and compute metrics
dat <- dat %>%
  mutate(
    starvation = factor(starvation, levels = c("HS", "FS"), labels = c("Low Hunger", "High Hunger")),
    toxin = factor(toxin, levels = c("No toxin", "Toxin")),
    perception = factor(perception, levels = c("Local", "Nonlocal"), labels = c("Local (MIA471)", "Nonlocal (MT1073)")),
    total_food = low + medium + high,
    elsewhere = total_worms - total_food,
    prop_inside = total_food / total_worms
  )

## ============================================================
## 2. Generate Summary Statistics (Mean & SE)
## ============================================================

inside_sum <- dat %>%
  group_by(starvation, toxin, perception) %>%
  summarise(
    mean_prop_inside = mean(prop_inside),
    se = sd(prop_inside) / sqrt(n()),
    .groups = "drop"
  )

## ============================================================
## 3. Generate Aligned Interaction Plot
## ============================================================

p_inside_plots <- ggplot(
  inside_sum,
  aes(x = starvation, y = mean_prop_inside, color = perception, group = perception)
) +
  # 1. Background raw individual data points to show distribution density
  geom_point(
    data = dat,
    aes(x = starvation, y = prop_inside, color = perception, group = perception),
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.25),
    size = 1.6,
    alpha = 0.45, # Slightly bumped up alpha so yellow individual points stay visible
    show.legend = FALSE
  ) +
  # 2. Standard Error tracking whiskers
  geom_errorbar(
    aes(
      ymin = pmax(mean_prop_inside - se, 0),
      ymax = pmin(mean_prop_inside + se, 1.0)
    ),
    width = 0.12,
    linewidth = 0.6,
    position = position_dodge(width = 0.25)
  ) +
  # 3. Connection path trendlines
  geom_line(
    linewidth = 1.0,
    position = position_dodge(width = 0.25)
  ) +
  # 4. Large crisp summary data points
  geom_point(
    aes(fill = perception),
    shape = 21,
    size = 4.2,
    color = "black",
    stroke = 0.8,
    position = position_dodge(width = 0.25)
  ) +
  # 5. Bold value data labels placed neatly over tracking metrics
  geom_text(
    aes(label = sprintf("%.2f", mean_prop_inside)),
    vjust = -1.8,
    fontface = "bold",
    size = 3.6,
    family = "Helvetica",
    position = position_dodge(width = 0.25),
    show.legend = FALSE
  ) +
  # Facet Side-by-Side: Columns = Toxin Status
  facet_wrap(~ toxin) +
  # Color Keys: Rich amber-yellow and deep royal blue theme
  scale_color_manual(values = c("Local (MIA471)" = "#ecb018", "Nonlocal (MT1073)" = "#1f77b4")) +
  scale_fill_manual(values = c("Local (MIA471)" = "#ecb018", "Nonlocal (MT1073)" = "#1f77b4")) +
  scale_y_continuous(
    breaks = seq(0, 0.7, 0.1),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  coord_cartesian(ylim = c(0, 0.7)) + 
  labs(
    x = "Hunger Level",
    y = "Proportion of abundance in food patches",
    color = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 14, base_family = "Helvetica") +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(face = "bold"),
    legend.position  = "top",
    legend.title     = element_text(face = "bold"),
    panel.spacing    = unit(1.5, "lines")
  )

# Output plot graphics layout
print(p_inside_plots)

## ============================================================
## 4. File Saving Outputs
## ============================================================

ggsave(
  file.path(figures_dir, "worm_abundance_inside_patches.png"),
  p_inside_plots,
  width = 10,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(figures_dir, "worm_abundance_inside_patches.pdf"),
  p_inside_plots,
  width = 10,
  height = 6.5,
  device = "pdf",
  bg = "white"
)

## ============================================================
## 5. Automated Integrity Checks
## ============================================================
stopifnot(nrow(dat) == 40)
stopifnot(all(dat$total_food <= total_worms))
stopifnot(all(dat$elsewhere >= 0))