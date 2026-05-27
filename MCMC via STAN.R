# Lade Zeug-------------
{
  realDataAvailable <- TRUE
  
  library(SoundModeling)
  if (realDataAvailable)
    library(SoundModelingData)
  
  ntp_data_acc <- ntp_data[!ntp_data$seltsam_mbbm,]
  
  imso_matrix <- as.matrix(ntp_data_acc[,13:33])
  
  sides <- volume_to_side_length(ntp_data_acc$E_VolA, output = "both")
  
  stan_imso_list <- list(
    y = imso_matrix,
    N = nrow(imso_matrix),
    rho = 3500,
    l1 = sides$l1,
    l2 = sides$l2
  )
  
  
  stan_imso_list$y <- t(
    apply(
      stan_imso_list$y,
      1, linear_continuation, orientation = 3
    )
  )
  
  
  
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  rstan_options(threads_per_chain = 1)
}


# Prior helpmat-------------------




bp1 <- 86
bp2 <- 1000


x1 <- rep(1, 21)
x2 <- log(freq)
x3 <- ifelse(freq > bp1, 1,0) * log(freq/bp1) # segment 2
x4 <- ifelse(freq > bp2, 1,0) * log(freq/bp2) # segment 3

X3seg <- cbind(x1,x2,x3,x4)
X2seg <- cbind(x1,x2,x3)


beta_prior_help_3seg <- X3seg
beta_prior_help_2seg <- X2seg




# Fix 1 break---------


imso_list_mshift <- stan_imso_list
imso_list_mshift$max_shift <- 1
scofflen <- 21*imso_list_mshift$max_shift-imso_list_mshift$max_shift*(imso_list_mshift$max_shift+1)/2
imso_list_mshift$priorHelpMat <- beta_prior_help_2seg

#best2seg <- readRDS("Saves/ImSo tvAR/New Prior/ImSo_fix_1break_tvAR_3")
#best2segFit <- extract(best2seg)


## TVAR(1)-------



imso_fix_1b <- stan("STAN/ImSo_Fix_1break.stan", data = imso_list_mshift)

isf_ex1b <- extract(imso_fix_1b)

loo::loo(isf_ex1b$log_lik) # 27567.8 

saveRDS(imso_fix_1b, "Saves/ImSo_fix_1break_tvAR_1")


data_mean <- apply(stan_imso_list$y, 2, mean)
model_mean <- apply(isf_ex1b$mu, 2, mean)

plot(data_mean)
points(model_mean, col = "blue")



## While LOOP---------

last_loo <- Inf
next_shift <- TRUE

L <- 1L

loo_list <- list()
imso_list_mshift <- stan_imso_list

#imso_list_mshift$priorHelpMat <- beta_prior_help_2seg

initlist <- list(
  fe_betas = c(100,-12,5),
  b_point_diff = 50,
  shift_coeff = rep(0.7, scofflen),
  variances = rep(5,21)
)



while (next_shift) {
  imso_list_mshift$max_shift <- L
  scofflen <- 21*imso_list_mshift$max_shift-imso_list_mshift$max_shift*(imso_list_mshift$max_shift+1)/2
  
  
  initlist$shift_coeff <- rep(0, scofflen)
  initlist4 <- list(initlist,initlist,initlist,initlist)
  
  
  poster <- stan("STAN/ImSo_Fix_1break.stan",
                 data = imso_list_mshift,
                 init = initlist4,
                 save_warmup = FALSE)
  
  loo_list[[L]] <- loo::loo(extract(poster)$log_lik)
  print(check_hmc_diagnostics(poster))
  print(loo_list[[L]])
  
  loo_val <- loo_list[[L]]$estimates[3]
  loo_is_better <- loo_val <= last_loo
  
  if (loo_is_better) {
    best2b_poster <- poster
    last_loo <- loo_val
  }
  if (!loo_is_better) {
    saveRDS(
      best2b_poster,
      paste("Saves/ImSo_fix_1break_tvAR",
            L-1, sep = "_")
    )
    
    next_shift <- FALSE
  }
  
  saveRDS(loo_list, "Saves/LOOlist_Fix1b")
  
  L <- L + 1L
  
}

# Fix 2 breaks---------


imso_list_mshift <- stan_imso_list
imso_list_mshift$max_shift <- 3
scofflen <- 21*imso_list_mshift$max_shift-imso_list_mshift$max_shift*(imso_list_mshift$max_shift+1)/2


#best3seg <- readRDS("Saves/ImSo tvAR/New Prior/ImSo_fix_2breaks_tvAR_3")
#best3segFit <- extract(best3seg)



## TVAR(3)--------


imso_list_mshift <- stan_imso_list
imso_list_mshift$max_shift <- 3L
scofflen <- 21*imso_list_mshift$max_shift-imso_list_mshift$max_shift*(imso_list_mshift$max_shift+1)/2



initlist <- list(
  fe_betas = c(100,-12,5,-13),
  b_point_diff = 50,
  b_point_diff2 = 1330,
  shift_coeff = rep(0.7, scofflen),
  variances = rep(5,21)
)
initlist4 <- list(initlist,initlist,initlist,initlist)


imso_fix <- stan("STAN/ImSo_Fix_2breaks.stan",
                 data = imso_list_mshift,
                 init = initlist4,
                 save_warmup = FALSE)

isf_ex <- extract(imso_fix)

loo::loo(isf_ex$log_lik) # n_breaks = 2 => loo = 13893.0 

saveRDS(imso_fix, "Saves/ImSo tvAR/New Prior/ImSo_Fix_2breaks_tvAR_3")


# Ist Präzisionsmatrix * Kovarianzmatrix = Identitätsmatrix?
s <- 100
round(isf_ex$precisionMatrix[s,,] %*% isf_ex$cov_mat[s,,], 5)

solve(apply(isf_ex$betaPriorPrecision, 2:3, mean))


data_mean <- apply(stan_imso_list$y, 2, mean)
model_mean <- apply(isf_ex$mu, 2, mean)

plot(data_mean)
points(model_mean, col = "blue")


## While LOOP---------

last_loo <- Inf
next_shift <- TRUE

L <- 1L

loo_list <- list()
imso_list_mshift <- stan_imso_list



initlist <- list(
  fe_betas = c(100,-12,5,-13),
  b_point_diff = 50,
  b_point_diff2 = 1330,
  shift_coeff = rep(0.7, scofflen),
  variances = rep(5,21)
)



while (next_shift) {
  imso_list_mshift$max_shift <- L
  scofflen <- 21*imso_list_mshift$max_shift-imso_list_mshift$max_shift*(imso_list_mshift$max_shift+1)/2
  
  initlist$shift_coeff <- rep(0, scofflen)
  initlist4 <- list(initlist,initlist,initlist,initlist)
  
  
  poster <- stan("STAN/ImSo_Fix_2breaks.stan",
                 data = imso_list_mshift,
                 init = initlist4,
                 save_warmup = FALSE)
  
  loo_list[[L]] <- loo::loo(extract(poster)$log_lik)
  print(check_hmc_diagnostics(poster))
  print(loo_list[[L]])
  
  loo_val <- loo_list[[L]]$estimates[3]
  loo_is_better <- loo_val <= last_loo
  
  if (loo_is_better) {
    best2b_poster <- poster
    last_loo <- loo_val
  }
  if (!loo_is_better) {
    saveRDS(
      best2b_poster,
      paste("Saves/ImSo tvAR/New Prior/ImSo_Fix_2breaks_tvAR",
            L-1, sep = "_")
    )
    
    next_shift <- FALSE
  }
  
  saveRDS(loo_list, "Saves/ImSo tvAR/New Prior/LOOlist_Fix2b")
  
  L <- L + 1L
  
}


# Fix 3 breaks---------



## TVAR(1)------------


imso_list_mshift <- stan_imso_list
imso_list_mshift$max_shift <- 1L
scofflen <- 21*imso_list_mshift$max_shift-imso_list_mshift$max_shift*(imso_list_mshift$max_shift+1)/2


imso_fix <- stan("STAN/ImSo_Fix_3breaks.stan", data = imso_list_mshift)

isf_ex <- extract(imso_fix)


loo::loo(isf_ex$log_lik) # loo = 27524.2  

saveRDS(imso_fix, "Saves/ImSo tvAR/ImSo_Fix_3breaks_tvAR_1")

# Konvergenzprobleme!

data_mean <- apply(stan_imso_list$y, 2, mean)
model_mean <- apply(isf_ex$mu, 2, mean)

plot(data_mean, ylim = c(0,60))
points(model_mean, col = "blue")





# LOO Summary Fix----------


loolist1b <- readRDS("Saves/ImSo tvAR/New Prior/LOOlist_Fix1b")
loolist2b <- readRDS("Saves/ImSo tvAR/New Prior/LOOlist_Fix2b")

loolistfix <- c(loolist1b, loolist2b)
names(loolistfix) <- c(
  paste("K2L", 1:4, sep = ""),
  paste("K3L", 1:4, sep = "")
)


saveRDS(loolistfix,"Saves/ImSo tvAR/New Prior/LOOlist_Fix")


loo::loo_compare(loolistfix)

