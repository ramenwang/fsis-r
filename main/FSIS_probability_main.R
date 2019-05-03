options(scipen = 999)
options(stringsAsFactors = F)

# set working directory under FSIS folder
source("./main/arg.R")
source("./shared/utility.R")

# read table
obs_dat <- read.csv(obs_dir)
names(obs_dat) <- toupper(names(obs_dat))

# read variables
fs_obs_loc     = as.matrix(obs_dat[, var_nam])
fs_obs         = obs_dat[, y_nam]/100

# global parameters
fs_dim         = length(var_nam)
fs_map         = getFSMap(img_dir)
block_obj      = createBlock(fs_block_size, fs_dim, if_empty = F)

# dir setting
try(dir.create(out_dir))
try(dir.create(paste0(out_dir,"/sample_fs")))

# check feature space distribution pattern
for (i in 1:(fs_dim - 1)) {
    for (j in (i+1):fs_dim) {
        fs_obs_loc_2d <- fs_obs_loc[,c(i,j)]
        bitmap(file  = paste0(out_dir,"/sample_fs/", colnames(fs_obs_loc_2d)[1], 
                              "-",colnames(fs_obs_loc_2d)[2],".tiff"), height = 3,
               width = 3, units = "in", res = 900)
        plot2DFS(fs_obs, fs_obs_loc_2d, is_Classify)
        dev.off()
    }
}
cat("Please check <", paste0(out_dir,"/sample_fs"), "> for sample distribution in 2D feature space")

# fit feature space variogram
lsfit <- fitFSVGM(fs_obs, fs_obs_loc, vgm_model = "poly", ifplot = T, is_classify)
# split training and testing set
makeTrainingSet(fs_obs, fs_obs_loc, pct_train)
# test FSIS run on training and testing set
for (i in 1:sim_num) {
    FSIS_pred_single <- singleFSIS(fs_obs_train, fs_obs_loc_train, fs_obs_loc_test, 
                                   block_obj, fs_kn_size, density_kernel_function, 
                                   prob_precision, fs_obs_num, lsfit, vgm_model)
    if (i == 1) {
        FSIS_pred <- FSIS_pred_single
    } else {
        FSIS_pred <- cbind(FSIS_pred, FSIS_pred_single)
    }
    cat("Running on FSIS simulation: ", i,"\n")
}

# plot residule
plot(fs_obs_test, FSIS_pred_single, xlab = "Prediction", ylab = "Observation")

