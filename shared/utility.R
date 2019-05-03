################################################################################
####################### Utility Functions for FSIS #############################
####################### Author:   Qing Wang        #############################
####################### Date Updated: May 2nd 2019 #############################
################################################################################
library(parallelDist)
library(spatstat)
library(FNN)
library(raster)

# FUNCTION: drawFromCDF
# Draw random number form cdf estimated by feature-space observations based on 
# kernel function. The density to cdf function is calculated based on spatstat.
# see https://cran.r-project.org/web/packages/spatstat/spatstat.pdf
drawFromCDF <- function(fs_kn_pred, density_kernel_function, prob_precision) {
    ds      = density(fs_kn_pred, kernel = density_kernel_function, from = 0, 
                      to = 1)
    kcdf    = CDF(ds, warn = F)
    x_axis  = seq(0, 1, 0.1^prob_precision)
    prob    = round(kcdf(x_axis), prob_precision)
    r_prob  = round(runif(1), prob_precision)
    return(x_axis[which.min(abs(prob - r_prob))])
}

# FUNCTION: distKNN
# Calculate feature space indicator kriging parameters
# The distances in high-dimension feature space for k neareast neighbors are
# calculated based on a fast neighborhood searching algorithm -- Approximate 
# Near Neighbor, and see http://www.cs.umd.edu/~mount/ANN/
# The feature space distance matrix is calculated by parallel distance package
# see https://cran.r-project.org/web/packages/parallelDist/parallelDist.pdf
distKNN <- function(fs_loc, fs_obs_loc, fs_obs_num) {
    knn         = get.knnx(fs_obs_loc, fs_loc, k=fs_obs_num, algo="kd_tree")
    knn_loc     = fs_obs_loc[knn[[1]],]
    loc_matrix  = rbind(fs_loc, knn_loc)
    dist_matrix = as.matrix(parDist(loc_matrix, method = "euclidean"))[-1,]
    return(list(dist_matrix, knn[[1]]))
}

# FUNCTION: fitFSVGM
# Form a feature space variogram function for indicator kriging. This is an
# autofitting process based on least-square esitmation
fitFSVGM <- function(fs_obs, fs_obs_loc, vgm_model = "sigm", ifplot = T) {
    h = c(parDist(fs_obs_loc, method = "sigm"))
    obs_dif =  c(parDist(t(rbind(fs_obs, rep(0,length(fs_obs))))))
    # fitting model by using least-square estimation
    if (vgm_model == "poly") {
        lsfit = lm(obs_dif ~ h + I(h^2) + I(h^3))
    } else if (vgm_model == "sigm") {
        lsfit = glm(obs_dif ~ h, family=binomial)
    } else if (vgm_model == "relu") {
        lsfit = lm(obs_dif ~ h)
    } else {
        cat("Please select model from sph, gau or pwd, which stand for spheric
            Gaussian and power model")
        stop
    }
    # plot variogram model
    if (ifplot == T) {
        plot(h, obs_dif, ylab = "Gamma", xlab = "Feature-space distance", 
             col = "blue", pch = 4)
        peudo_h = seq(0, max(h), max(h)/100)
        if (vgm_model %in% c("poly","relu")) {
            lines(peudo_h, predict(lsfit, data.frame(h = peudo_h)), col = "red")
        } else {
            lines(peudo_h, predict(lsfit, data.frame(h = peudo_h), 
                                   type="response"), col = "red")
        }
    }
    return(lsfit)
}

# FUNCTION: predFSVGM
# Predict indicator variogram value based on fitted variogram model
predFSVGM <- function(lsfit, h, vgm_model) {
    if (vgm_model %in% c("poly","relu")) {
        output = predict(lsfit, data.frame(h = h))
    } else {
        output = predict(lsfit, data.frame(h = h), type="response")
    }
    output = ifelse(output > 1, 1, output)
    output = ifelse(output < 0, 0, output)
    return(output)
}

# FUNCTION: FSIK
# Calculate indicator Kriging for estimated location in feature space
# calculation see http://www.imm.dtu.dk/~alan/krexample.pdf
FSIK <- function(fs_loc, fs_obs_loc, fs_obs, fs_obs_num, lsfit, vgm_model) {
    dist_knn = distKNN(fs_loc, fs_obs_loc, fs_obs_num)
    D = dist_knn[[1]][,1]
    C = dist_knn[[1]][,-1]
    knn_obs = fs_obs[dist_knn[[2]]]
    C.dim = dim(C)
    # Compute spatial correlagram
    D = 1 - predFSVGM(lsfit, D, vgm_model)
    C = 1 - predFSVGM(lsfit, c(C), vgm_model)
    D = matrix(D, ncol = 1)
    C = matrix(C, ncol = C.dim[1])
    # Kriging calculation
    ones   = as.vector(rep(1, C.dim[1]))
    C.inv  = solve(C)
    lambda = (t(D) %*% C.inv %*% ones - 1) / (t(ones) %*% C.inv %*% ones)
    krig.w = C.inv %*% (D - lambda[1,1] * ones)
    return((t(knn_obs) %*% krig.w)[1])
}

# FUNCTION: blockFSIK
# Calculate indicator Kriging for multiple neighbors around the location to be
# estimated in feature space
blockFSIK <- function(fs_block_loc, fs_obs_loc, fs_obs, fs_obs_num, 
                       lsfit, vgm_model) {
    fs_kn_pred = array(dim = nrow(fs_block_loc))
    for (i in 1:fs_kn_size) {
        fs_kn_pred[i] <- FSIK(t(fs_block_loc[i,]), fs_obs_loc, fs_obs, 
                               fs_obs_num, lsfit, vgm_model)
    }
    return(fs_kn_pred)
}

# FUNCTION: createBlock
# Create a high-dimension block object later used in sampleBlock function.
# The return would be a FS location matrix with nrow = fs_block_size * fs_dim,
# and ncol = fs_dim, centered at FS location 0,0,...,0
createBlock <- function(fs_block_size, fs_dim, if_empty) {
    block_center = matrix(rep(0, fs_dim), ncol = fs_dim)
    block_ext    = (fs_block_size - 1) / 2
    block_obj    = matrix(nrow = fs_block_size ^ fs_dim, ncol = fs_dim)
    if (if_empty) return(block_obj)
    rep_num = rep(x = (block_center[1,i] - block_ext) : 
                      (block_center[1,i] + block_ext))
    for (i in 1:fs_dim) {
         ind_dup    = fs_block_size ^ (i-1)
         ary_dup    = fs_block_size ^ (fs_dim - i)
         ind_dupped = c(sapply(X = rep_num, function(x) {x = rep(x, ind_dup)}))
         ary_dupped = rep(ind_dupped, ary_dup)
         block_obj[,i] = ary_dupped
    }
    return(block_obj)
}

# FUNCTION: getFSMap
# Read imagery data and create a map of the corresponding feature space
getFSMap <- function(img_dir) {
    for (i in 1:length(img_dir)) {
        img_obj = readImageRS(img_dir[i])
        # make sure all the image inputs have the same extend, projects and resolution
        if (i > 2 ) {
            if (setdiff(img_obj[[2]], ext)) { 
                cat("ERROR: image extends are different")
                stop
            }
            if (setdiff(img_obj[[3]], img_res) != 0) {
                cat("ERROR: image resolutions are different")
                stop
            }
            if (img_obj[[4]] != prj) {
                cat("ERROR: image projections are different")
                stop
            }
            fs_map = cbind(fs_map, img_obj[[1]])
        }
        ext     <- img_obj[[2]]
        img_res <- img_obj[[3]]
        prj     <- img_obj[[4]]
        if (i == 1) {
            fs_map = img_obj[[1]]
        }
    }
    img_ext <<- ext
    img_res <<- img_res
    img_prj <<- prj
    return(unique(fs_map))
}

# FUNCTION: readImageRS
# Read single layer remotely sensed image as a matrix, also return projection,  
# extend, and resolution. Based on raster package: 
# https://cran.r-project.org/web/packages/raster/raster.pdf
readImageRS <- function(img_dir_sig) {
    img     <- raster(img_dir_sig)
    ext     <- extend(img)
    prj     <- proj4string(img)
    img_res <- res(img)
    img_mtx <- as.matrix(img)
    return(list(c(img_mtx), ext, img_res, prj))
}

# FUNCTION: plot2DFS
# Plot a 2D feature space for selected features
plot2DFS <- function(fs_obs, fs_obs_loc_2d, is_Classify = F){
    if (is_Classify == F) {
        colfunc <- colorRampPalette(c("#fee5d9", "#fcae91", "#fb6a4a", 
                                      "#de2d26", "#a50f15"))
        colfunc <- colfunc(5)[as.numeric(cut(fs_obs,breaks = 5))]
        lgd     <- sort(unique(cut(fs_obs,breaks = 5)))
        pch_lgd <- rep(15,5)
        col_lgd <- c("#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15")
    } else {
        colfunc <- colorRampPalette(c("#deebf7", "#3182bd"))
        colfunc <- colfunc(2)[as.numeric(cut(fs_obs,breaks = 2))] 
    }
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    plot(fs_obs_loc_2d[,1], fs_obs_loc_2d[,2],
         col = colfunc, pch = 15, xlab = colnames(fs_obs_loc_2d)[1],
         ylab = colnames(fs_obs_loc_2d)[2],
         main = paste0("2D Feature Space of ", colnames(fs_obs_loc_2d)[1], "-", colnames(fs_obs_loc_2d)[2]))
    legend(par('usr')[2], par('usr')[4], inset=c(-1,0), legend=lgd, 
           col = col_lgd, pch=pch_lgd, title="Value", bty='n')
}

# FUNCTION: makeTrainingSet
# Make training set and test set based on samples. This can be used before fully
# running the model for image classifiation or probability mapping
makeTrainingSet <- function(fs_obs, fs_obs_loc, pct_train) {
    random_index     = sample(1:length(fs_obs), round(pct_train*length(fs_obs)))
    fs_obs_train     <<- fs_obs[random_index]
    fs_obs_loc_train <<- fs_obs_loc[random_index,]
    fs_obs_test      <<- fs_obs[-random_index]
    fs_obs_loc_test  <<- fs_obs_loc[-random_index,]
}

# FUNCTION: singleFSIS
# Running FSIS for a single time
singleFSIS <- function(fs_obs, fs_obs_loc, fs_map, block_obj,
                       density_kernel_function, prob_precision, fs_obs_num) {
    fs_map = unique(fs_map)
    fs_path <- sample(1:nrow(fs_map), nrow(fs_map))
    for (i in fs_path) {
        fs_loc = fs_map[i,]
        block_loc <- t(fs_loc+t(block_obj))
        
    }
}