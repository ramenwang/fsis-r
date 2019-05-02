# observation setting
obs_dir        = "./data/shenzhen_variables.csv"
var_nam        = c("B4", "SR536", "TVI", "SR32")
y_nam          = "VEG_CARBON"

# image setting
img_dir        = c("./data/B4.tif", "./data/SR36.tif", "./data/TVI.tif", "./dataSR32.tif")

# hyperparameter
fs_obs_num     = 5      # -- number of nearest FS observations included
pixel_depth    = 256    # -- pixel depth used to form space
fs_block_size  = 5      # -- (1d) size of the high-dimensional block, should be an odd number
fs_kn_size     = 2      # -- number of nearest FS estimation for cdf
prob_precision = 3      # -- probility precision

# kernel function setting
density_kernel_function = "gaussian"   # -- kernel function for density estimation
vgm_model               = "sigm"       # -- kernel function for FS variogram model
fs_dist_method          = "euclidean"  # -- distance function