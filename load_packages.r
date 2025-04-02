##################
# load R packages
##################

# List of required packages
packages <- c("fda", "tseries", "sandwich", "sde", "variables", "basefun", "polynom", "fracdiff", "LongMemoryTS", "arfima")

# Install missing packages
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(packages, library, character.only = TRUE)[[length(packages)]]
