# Recorded experimental counts: five replicates per treatment combination.
# Keep one copy so plots and statistical models use identical observations.
library(tidyverse)

get_patch_data <- function() {
  dat <- tribble(
    ~rep, ~starvation, ~toxin,     ~perception, ~low, ~medium, ~high,
    
    # Half starvation, no toxin
    1, "HS", "No toxin", "Local",    0, 2, 0,
    2, "HS", "No toxin", "Local",    1, 0, 0,
    3, "HS", "No toxin", "Local",    1, 1, 0,
    4, "HS", "No toxin", "Local",    2, 2, 0,
    5, "HS", "No toxin", "Local",    0, 2, 0,
    
    1, "HS", "No toxin", "Nonlocal", 1, 1, 1,
    2, "HS", "No toxin", "Nonlocal", 0, 1, 0,
    3, "HS", "No toxin", "Nonlocal", 2, 0, 1,
    4, "HS", "No toxin", "Nonlocal", 1, 0, 0,
    5, "HS", "No toxin", "Nonlocal", 1, 4, 0,
    
    # Half starvation, toxin
    1, "HS", "Toxin", "Local",       0, 1, 1,
    2, "HS", "Toxin", "Local",       0, 0, 1,
    3, "HS", "Toxin", "Local",       0, 0, 0,
    4, "HS", "Toxin", "Local",       0, 0, 1,
    5, "HS", "Toxin", "Local",       1, 1, 1,
    
    1, "HS", "Toxin", "Nonlocal",    4, 1, 0,
    2, "HS", "Toxin", "Nonlocal",    3, 0, 0,
    3, "HS", "Toxin", "Nonlocal",    1, 0, 0,
    4, "HS", "Toxin", "Nonlocal",    2, 2, 1,
    5, "HS", "Toxin", "Nonlocal",    2, 0, 0,
    
    # Full starvation, no toxin
    1, "FS", "No toxin", "Local",    1, 2, 3,
    2, "FS", "No toxin", "Local",    1, 1, 2,
    3, "FS", "No toxin", "Local",    0, 2, 2,
    4, "FS", "No toxin", "Local",    1, 0, 2,
    5, "FS", "No toxin", "Local",    1, 2, 1,
    
    1, "FS", "No toxin", "Nonlocal", 3, 0, 0,
    2, "FS", "No toxin", "Nonlocal", 1, 2, 1,
    3, "FS", "No toxin", "Nonlocal", 0, 0, 5,
    4, "FS", "No toxin", "Nonlocal", 1, 3, 1,
    5, "FS", "No toxin", "Nonlocal", 0, 3, 3,
    
    # Full starvation, toxin
    1, "FS", "Toxin", "Local",       0, 1, 0,
    2, "FS", "Toxin", "Local",       0, 1, 0,
    3, "FS", "Toxin", "Local",       2, 3, 0,
    4, "FS", "Toxin", "Local",       0, 2, 0,
    5, "FS", "Toxin", "Local",       3, 3, 1,
    
    1, "FS", "Toxin", "Nonlocal",    0, 4, 2,
    2, "FS", "Toxin", "Nonlocal",    1, 0, 2,
    3, "FS", "Toxin", "Nonlocal",    1, 2, 0,
    4, "FS", "Toxin", "Nonlocal",    3, 1, 0,
    5, "FS", "Toxin", "Nonlocal",    1, 4, 0
  )
  stopifnot(nrow(dat) == 40L, all(dat$low + dat$medium + dat$high <= 10L))
  dat
}
