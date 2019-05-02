options(scipen = 999)
options(stringsAsFactors = F)

# set working directory under FSIS folder
source("/shared/utility.R")

# hyperparameter
fs_obs_num     = 5      # -- number of nearest FS observations included
pixel_depth    = 256    # -- pixel depth used to form space
fs_dim         = 3      # -- dimensionality of the feature space
fs_block_size  = 5      # -- (1d) size of the high-dimensional block, should be an odd number
fs_kn_size     = 2      # -- number of nearest FS estimation for cdf
prob_precision = 3      # -- probility precision

# kernel function setting
density_kernel_function = "gaussian"   # -- kernel function for density estimation
vgm_model               = "sigm"       # -- kernel function for FS variogram model
fs_dist_method          = "euclidean"  # -- distance function

