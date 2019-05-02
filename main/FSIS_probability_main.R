options(scipen = 999)
options(stringsAsFactors = F)

# set working directory under FSIS folder
source("./shared/utility.R")

# read table
obs_dat <- read.csv("./data/shenzhen_variables.csv")
names(obs_dat) <- toupper(names(obs_dat))

# read variables
fs_obs_loc     = as.matrix(obs_dat[, var_nam])
fs_obs         = obs_dat[, y_nam]

# global parameters
fs_dim         = length(var_nam)
fs_map         = getFSMap(img_dir)


# check feature space distribution pattern





