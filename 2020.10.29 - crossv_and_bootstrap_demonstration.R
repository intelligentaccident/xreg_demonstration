# clear? 
rm(list = ls())

options(digits=6)
options(scipen=15)


# Set correct working directory
wd <- "C:/Data/EQVT/"

# Change working directory
oldwd <- getwd()

# set working directory. 
setwd(wd) 

dta <- read.csv('sample_data.csv', stringsAsFactors = F)

# # If xreg is not installed, remove comments and run:
# install.packages("devtools")
# library(devtools)
# install_github("intelligentaccident/xreg")
library(xreg)


# Import functions for easy handling of EQ-5D states
source('2020.07.08 - Import EQ-5D data functions.R')

# Make a vector of crossed vectors, separated by sep
ovec <- function(..., typecast = as.character, FUN = function(a, b) paste0(a, b), rev = F) {
  args <- list(...)
  larg <- length(args)
  for(narg in 1:larg) args[[narg]] <- typecast(args[[narg]])
  if(larg < 2) return(unlist(args))
  retv <- outer(args[[1L]], do.call(ovec, c(args = args[-1L], typecast = typecast, FUN = FUN )), FUN = FUN)
  if(rev) retv <- t(retv)
  return(as.vector(retv))
}


# Laziness
dim.snames <- c("mo", "sc", "ua", "pd", "ad")
DIM.SNAMES <- toupper(dim.snames)


# List wrapper to automatically name variables
nlist <- function(...) {
  dots <- list(...)
  mcn <- as.list(match.call())[-1]
  nms <- names(dots)
  names(dots) <- mcn
  if(!is.null(nms)) {
    names(dots)[tmp] <- nms[(tmp <- which(nms != ""))]
  }
  return(dots)
}

# Function for predicting values based on formula object and 
# either a named vector of coefficient values, 
# or a data.frame with coefficients, where rownames indicate name of coefficient.
predict.formula <- function(object, newdata, args, return_df = FALSE, pred_colname = "predict") {
  if(class(args) %in% c("matrix", "data.frame")) {
    return(apply(X = args, MARGIN = 2, FUN = function(x) predict.formula(object, newdata, x, return_df = F, pred_colname)))
  }
  if(!all(all.vars(object[[3]]) %in% c(names(newdata), names(args)))) stop("Variables in object not defined in newdata or args.")
  newdata[, pred_colname] <- with(as.list(args), with(newdata, eval(parse(text = object))))
  if(return_df) return(newdata)
  return(newdata[, pred_colname])
}





# 20 parameter formula, no intercept
f20 <- formula(c ~ d_mo2 * MO2 + d_sc2 * SC2 + d_ua2 * UA2 + d_pd2 * PD2 + d_ad2 * AD2 + 
                  d_mo3 * MO3 + d_sc3 * SC3 + d_ua3 * UA3 + d_pd3 * PD3 + d_ad3 * AD3 + 
                  d_mo4 * MO4 + d_sc4 * SC4 + d_ua4 * UA4 + d_pd4 * PD4 + d_ad4 * AD4 + 
                  d_mo5 * MO5 + d_sc5 * SC5 + d_ua5 * UA5 + d_pd5 * PD5 + d_ad5 * AD5)

# corresponding formula for sigma
f20h <- formula(sigma_est ~ 
                   d_mo2 * HMO2 + d_sc2 * HSC2 + d_ua2 * HUA2 + d_pd2 * HPD2 + d_ad2 * HAD2 + 
                   d_mo3 * HMO3 + d_sc3 * HSC3 + d_ua3 * HUA3 + d_pd3 * HPD3 + d_ad3 * HAD3 + 
                   d_mo4 * HMO4 + d_sc4 * HSC4 + d_ua4 * HUA4 + d_pd4 * HPD4 + d_ad4 * HAD4 + 
                   d_mo5 * HMO5 + d_sc5 * HSC5 + d_ua5 * HUA5 + d_pd5 * HPD5 + d_ad5 * HAD5)

# 8 parameter formula
f8 <- formula(c~ (MO * (d_mo2 * L2 + d_mo3 * (L3-L2) + d_mo4 * (L4-L3) + d_mo5 * (1-L4)) + 
                     SC * (d_sc2 * L2 + d_sc3 * (L3-L2) + d_sc4 * (L4-L3) + d_sc5 * (1-L4)) + 
                     UA * (d_ua2 * L2 + d_ua3 * (L3-L2) + d_ua4 * (L4-L3) + d_ua5 * (1-L4)) + 
                     PD * (d_pd2 * L2 + d_pd3 * (L3-L2) + d_pd4 * (L4-L3) + d_pd5 * (1-L4)) + 
                     AD * (d_ad2 * L2 + d_ad3 * (L3-L2) + d_ad4 * (L4-L3) + d_ad5 * (1-L4))))

# 8 parameter formula for sigma
f8h <- formula(sigma_est ~ (HMO * (d_mo2 * HL2 + d_mo3 * (HL3-HL2) + d_mo4 * (HL4-HL3) + d_mo5 * (1-HL4)) + 
                           HSC * (d_sc2 * HL2 + d_sc3 * (HL3-HL2) + d_sc4 * (HL4-HL3) + d_sc5 * (1-HL4)) + 
                           HUA * (d_ua2 * HL2 + d_ua3 * (HL3-HL2) + d_ua4 * (HL4-HL3) + d_ua5 * (1-HL4)) + 
                           HPD * (d_pd2 * HL2 + d_pd3 * (HL3-HL2) + d_pd4 * (HL4-HL3) + d_pd5 * (1-HL4)) + 
                           HAD * (d_ad2 * HL2 + d_ad3 * (HL3-HL2) + d_ad4 * (HL4-HL3) + d_ad5 * (1-HL4))))

# 8 parametrer formula for sigma with intercept
f8hi <- appendmodel(object = f8h, prepend = "INTERCEPT_SIGMA + ")

# 20 parametrer formula for sigma with intercept
f20hi <- appendmodel(object = f20h, prepend = "INTERCEPT_SIGMA + ")

# 8-parameter with intercept
f8i <- appendmodel(object = f8, prepend = 'INTERCEPT + ')

# 20-parameter with intercept
f20i <- appendmodel(object = f20, prepend = 'INTERCEPT + ')

f8_stv <- c(L2 = 0.25, L3 = 0.5, L4 = 0.75, MO = 0.2, SC = 0.2, UA = 0.2, PD = 0.2, AD = 0.2)
f8i_stv <- c(INTERCEPT = 0, f8_stv)

f20_stv <- rep(0.1, 20)
names(f20_stv) <- ovec(c("MO", "SC", "UA", "PD", "AD"), 2:5)
f20i_stv <- c(INTERCEPT = 0, f20_stv)


# homoscedastic sigma
ho_sigma <- formula(sigma_est ~ SIGMA)

# sigma as an affine transformation of Xb
het_sigma <- formula(sigma_est ~ HET_INTERCEPT + HET_SLOPE * Xb)


tto <- dta[dta$type == "TTO",]
dce <- dta[dta$type == "DCE",]






# create dataframe with unique states and their observed means, as well as likelihood-based censored means
meandta <- aggregate(x = tto[, c("c"), drop = F], 
                     by = as.list(tto[, c('block', "type", "state_A", "state_B", grep(pattern = '^d_', x = colnames(tto), value = T))]),
                     FUN = mean)

meandta$id <- 0


meandta$mval <- meandta$mvalc1 <- meandta$mvalc2 <- meandta$c

for(lnum in 1:NROW(meandta)) {
  st <- meandta$state_B[lnum]
  bl <- meandta$block[lnum]
  print(paste0("State: ", st, " ,block: ", bl))
  mval <- meandta$c[lnum]
  meandta$mvalc2[lnum] <- xreg(controlList = xregControl(formulas = c ~ VAL, start_values = c(VAL = mval), name = "TTO", censor_bounds = c(-Inf, 2)), dataList = list(TTO = tto[tto$state_B == st & tto$block == bl,]))$coef[["VAL"]]
  meandta$mvalc1[lnum] <- xreg(controlList = xregControl(formulas = c ~ VAL, start_values = c(VAL = mval), name = "TTO", censor_bounds = c(-Inf, 1)), dataList = list(TTO = tto[tto$state_B == st & tto$block == bl,]))$coef[["VAL"]]
}

##########################################################################################
# Various model specifications. These first two are used in bootstrap demonstration
##########################################################################################



# 8 parameter model control, homoscedastic, non-censored
par8_c <-xregControl(formulas = list(f8i, sigma_est ~ SIGMA), start_values = c(SIGMA = 0.5, f8i_stv), censor_bounds = c(-Inf, Inf), name = "TTO")

# 20 parameter model control, homoscedastic, non-censored
par20_c <-xregControl(formulas = list(f20i, sigma_est ~ SIGMA), start_values = c(SIGMA = 0.5, f20i_stv), censor_bounds = c(-Inf, Inf), name = "TTO")

# 8 parameter model
(par8_x <- xreg(controlList = par8_c, dataList = list(TTO = tto)))

# 20 parameter model
(par20_x <- xreg(controlList = par20_c, dataList = list(TTO = tto)))

##########################################################################################
# These are more advanced, and are not used in boostrap demonstration
# For boostrap demonstration, go to line saying  # BOOTSTRAP
##########################################################################################


# 8 parameter model control, homoscedastic, censored
par8c_c <-xregControl(formulas = list(f8i, sigma_est ~ SIGMA), start_values = par8_x$coef, censor_bounds = c(-Inf, 2), name = "TTO")

# 20 parameter model control, homoscedastic, censored
par20c_c <-xregControl(formulas = list(f20i, sigma_est ~ SIGMA), start_values = par20_x$coef, censor_bounds = c(-Inf, 2), name = "TTO")

# 8 parameter model censored
(par8c_x <- xreg(controlList = par8c_c, dataList = list(TTO = tto)))

# 20 parameter model
(par20c_x <- xreg(controlList = par20c_c, dataList = list(TTO = tto)))

# 8 parameter model control, linear heteroscedastic, censored
par8ch_c <-xregControl(formulas = list(f8i, sigma_est ~ SIGMA_INTERCEPT + SIGMA_SLOPE * Xb), start_values = par8_x$coef, censor_bounds = c(-Inf, 2), name = "TTO")

# 20 parameter model control, linear heteroscedastic, censored
par20ch_c <-xregControl(formulas = list(f20i, sigma_est ~ SIGMA_INTERCEPT + SIGMA_SLOPE * Xb), start_values = par20_x$coef, censor_bounds = c(-Inf, 2), name = "TTO")

# 8 parameter model censored, linear heteroscedastic
(par8ch_x <- xreg(controlList = par8ch_c, dataList = list(TTO = tto)))

# 20 parameter modelcensored, , linear heteroscedastic
(par20ch_x <- xreg(controlList = par20ch_c, dataList = list(TTO = tto)))

# 8 parameter with 8 parameter heteroscedasticity
par8ch8_c <- xregControl(formulas = list(f8i, f8hi), start_values = par8_x$coef, censor_bounds = c(-Inf, 2), name = "TTO")

# 8 parameter with 20 parameter heteroscedasticity
par8ch20_c <- xregControl(formulas = list(f8i, f20hi), start_values = par8_x$coef, censor_bounds = c(-Inf, 2), name = "TTO")

# 20 parameter with 8 parameter heteroscedasticity
par20ch8_c <- xregControl(formulas = list(f20i, f8hi), start_values = par20_x$coef, censor_bounds = c(-Inf, 2), name = "TTO")

# 20 parameter with 20 parameter heteroscedasticity
par20ch20_c <- xregControl(formulas = list(f20i, f20hi), start_values = par20_x$coef, censor_bounds = c(-Inf, 2), name = "TTO")

# 8 parameter model censored, 8 parameter heteroscedasticity
(par8ch8_x <- xreg(controlList = par8ch8_c, dataList = list(TTO = tto)))

# 20 parameter model censored, 8 parameter heteroscedasticity
(par20ch8_x <- xreg(controlList = par20ch8_c, dataList = list(TTO = tto)))

# 8 parameter model censored, 20 parameter heteroscedasticity
(par8ch20_x <- xreg(controlList = par8ch20_c, dataList = list(TTO = tto)))

# 20 parameter model censored, 20 parameter heteroscedasticity
(par20ch20_x <- xreg(controlList = par20ch20_c, dataList = list(TTO = tto)))


var_colnames <- ovec(c("", "c", "ch", "ch8", "ch20"), c("par8", "par20"))
xreg_names <- ovec(c("par8", "par20"), c("", "c", "ch", "ch8", "ch20"), "_x")

for(i in 1:10) {
  meandta[, var_colnames[i]] <- predict(object = get(x = xreg_names[i]), newdata = list(TTO = meandta), return_vector = T)
  
}

##########################################################################################
# Hybrid models run using the hyreg wrapper function
##########################################################################################

# Hybrid 20 parameter
(h20 <- hyreg(formula = f20, df = dta, datatype = "continuous", ul = 2))


# Hybrid 20 parameter with 20 parameter sigma
(h20h <- hyreg(formula = f20, df = dta, datatype = "continuous", ul = 2, hetcont = f20h, init = h20))

# Hybrid 20 parameter with linear sigma
(h20hl <- hyreg(formula = f20, df = dta, datatype = "continuous", ul = 2, hetcont = het_sigma, init = h20))

# Hybrid 8 parameter
(h8 <- hyreg(formula = f8, df = dta, datatype = "continuous", ul = 2))

# Hybrid 8 parameter with 8 parameter sigma
(h8h <- hyreg(formula = f8, df = dta, datatype = "continuous", ul = 2, hetcont = f8h, init = h8))

# Hybrid 8 parameter with linear sigma
(h8hl <- hyreg(formula = f8, df = dta, datatype = "continuous", ul = 2, hetcont = het_sigma, init = h8))


##########################################################################################
# Cross-validation demonstration
##########################################################################################


# list of models for crossvalidation
modlist <- list(h20 = h20,
                h20h = h20h,
                h20hl = h20hl,
                h8 = h8,
                h8h = h8h,
                h8hl = h8hl)

# Define crossvalidation groups
dta$crossblock <- as.character(interaction(dta[, c("type", "block")]))

# Find all unique ids
all_ids <- unique(dta$id)

# Sorted groups
crossgroups <- c(paste0("TTO.", 1:10), paste0("DCE.", 1:28))

# Crossvalidation test
crossv <- t(vapply(X = crossgroups, FUN.VALUE = rep(0, 6) , FUN = function(block) {
  # block <- unique(dta$crossblock)[i]
  
  # Identify ids and data to be excluded
  excl_ids <- unique(dta$id[dta$crossblock == block])
  excl_data <- dta[dta$id %in% excl_ids,]
  # Format data to fit with requirements for re-running xreg
  excl_dl <- list(Continuous = excl_data[excl_data$continuous == 1,], Dichotomous = excl_data[excl_data$continuous == 0,])
  
  
  # Identify ids and data to be included
  incl_ids <- all_ids[!all_ids %in% excl_ids]
  incl_data <- dta[dta$id %in% incl_ids,]
  incl_dl <- list(Continuous = incl_data[incl_data$continuous == 1,], Dichotomous = incl_data[incl_data$continuous == 0,])
  
  # Print progress
  cat(paste('\n', block, '\n'))
  return(unlist(lapply(X = modlist, FUN = function(mod) {
    # Fit model to included data using main model coefficients from full data as start values
    modf <- xreg(controlList = mod, dataList = incl_dl, start_values = mod$coef)
    # Fit model to excluded data with all parameters fixed
    fitf <- xreg(controlList = mod, dataList = excl_dl, fixed_values = modf$coef)
    # Extract total minus-log-likelihood
    return(-fitf$minima['total'])
  })))
  
}))

# Calculate sum
crossv <- rbind(crossv[1:10, ], SUMTTO = colSums(crossv[1:10, ]), crossv[11:38, ],SUMDCE = colSums(crossv[11:38, ]))

# Print resulting table
crossv


##########################################################################################
# BOOTSTRAP demonstration
# This example shows how to generate boostrap-based SE and 95% CI for coefficients and
# all 3125 health states using the vanilla 8- and 20-parameter models run on TTO data.
# For hybrid data, modifications are required.
##########################################################################################



# Datalist
dta <- list(TTO = tto)


# Set attributes for dta
attributes(dta) <- c(attributes(dta),
                     list(# List of unique ids by datatype
                        uids_list = (tmpl<-lapply(X = dta, FUN = function(x) unique(x$id))),
                        # List of number of unique ids by datatype
                        nids_list = lapply(X = tmpl, length),
                        # Unique ids over all datatypes
                        uids = (uids <- unique(unlist(tmpl))),
                        # Total unique ids
                        nids = (nids = length(uids)),
                        # data.frame containing all unique ids and whether they exist in each data type
                        uids_df = cbind(all = uids, as.data.frame(lapply(X = tmpl, FUN = function(x) uids %in% x))),
                        # List of row numbers corresponding to each unique id
                        id_rows_list = lapply(X = dta, FUN = function(dtat) structure(.Data = lapply(X = uids, function(thisid) which(dtat$id == thisid)), .Names = uids))))


# Set attributes for dta contents

for(dtype in names(dta)) {
  attributes(dta[[dtype]]) <- c(attributes(dta[[dtype]]), 
                                list(uids = unique(dta[[dtype]]$id),
                                   nids = attr(dta, 'nids_list')[[dtype]],
                                   id_rows = attr(dta, 'id_rows_list')[[dtype]]))
}

# Make data.frame with all 3125 EQ-5D-5L health states and required dummies
all_EQ5D5L_states <- data.frame(state_A = (as.matrix(expand.grid(mo_A = 1:5, sc_A = 1:5, ua_A = 1:5, pd_A = 1:5, ad_A = 1:5)) %*% 10^(0:4))[,1], state_B = 11111)

# Generate all kinds of dummies
all_EQ5D5L_states <- cbind(all_EQ5D5L_states, diff_EQ_dummies(data = all_EQ5D5L_states, vector_a = "state_A", vector_b = "state_B", vars_vector = "d_"))

# Keep only relevant columns. Not required.
all_EQ5D5L_states <- all_EQ5D5L_states[, colnames(all_EQ5D5L_states) %in% colnames(dta$TTO)]

# Add id column, to prevent possible issues while running predictions
all_EQ5D5L_states$id <- 0







# Function for extracting bootstrap subsample. 
# Should work for hybrid data, as long as id_rows attribute is set for each item in data list.
# Returns list with data.frame(s) containing subsample
extract_subsample <- function(sample_vector, dta_list = dta) {
  return(lapply(X = dta_list, FUN = function(this_dta) {
    cur_rows <- attr(this_dta, which = "id_rows")[sample_vector]
    cur_dta <- this_dta[unlist(cur_rows),]
    cur_dta$old_id <- cur_dta$id
    cur_dta$id <- unlist(Map(function(x, i) x[] <- rep(i, length(x)), x = cur_rows, i = 1:length(cur_rows)))
    return(cur_dta)
  }))
}





# Models over which bootstraps will be run
bootmodels <- list(par8_x = par8_x, 
                   par20_x = par20_x)

# Number of bootstrap subsamples. This should be a large number, e.g. 1000 or 10000
# Currently set to just 10, so that the code below runs without having to wait a very long
# time to test whether it works.

nboot <- 10

# Make bootstrap sample matrix. 
# For simplicity of coding, the samples are between 1 and nids, not from ids directly.
# set.seed prior to random sampling in order to make results exactly reproducible.
set.seed(1)
boot_matrix <- matrix(data = sample(x = 1:nids, 
                                    size = nboot*nids, 
                                    replace = T), 
                      nrow = nids, 
                      ncol = nboot)

# Helper function to get names within lapply
nfrom <- function(x) structure(.Data = as.list(names(x)), .Names = names(x))


# This function should be modified to fit needs. In the example, it extracts only fitted coefficients
# If return_fits is set to TRUE, the full xreg object will be added for each run, which takes a lot of memory
run_bootstrap <- function(bootmods = bootmodels, bootmat = boot_matrix, return_fits = F, fix_coefs = NULL, ...) {
  dot_args <- list(...)
  lapply(X = nfrom(bootmods), FUN = function(this_model_name){
    this_model <- bootmods[[this_model_name]]
    cat(paste0("\n", this_model_name, '\n'))
    list(model_name = this_model_name,
         orig_model = this_model,
         control_list = this_model$controlList,
         start_values = this_model$coef,
         fits = lapply(X = 1:NCOL(boot_matrix), FUN = function(i) {
           cat(paste0("  subsample ", i, "/", NCOL(boot_matrix),'\r'))
           this_dta <- extract_subsample(boot_matrix[,i])
           
           if(length(fix_coefs)) {
             if('xreg' %in% class(this_model)) fix_coefs <- this_model$coef[fix_coefs]
             if('xregControlList' %in% class(this_model)) fix_coefs <- this_model[[1]]$start_values[fix_coefs]
             fix_coefs <- fix_coefs[!is.na(fix_coefs)]
             
           }
           argl <- c(list(controlList = this_model$controlList, dataList = this_dta, start_values = this_model$coef, fixed_values = fix_coefs), dot_args)
           # print(names(argl))
           cur_fit <- do.call(what = "xreg", args = argl)
           return(list(fit = if(return_fits) cur_fit else NULL,
                       coef = cur_fit$coef,
                       subsample = boot_matrix[,i]))    
         }))
  })
}

testall <- run_bootstrap()

# Make a list with data.frames with coefficients for each model. 
# Formula added as an attribute to each model for convenience.

coefficient_list <- lapply(X = testall, FUN = function(this_model) {
  structure(.Data = as.data.frame(lapply(this_model$fits, function(this_fit) this_fit$coef)), 
            .Names = 1:length(this_model$fits), 
            formula = this_model$control_list$TTO$formulas[[1]])
})

##########################################################################################
# Fix to add 20-parameter expansion of 8-parameter model coefficients. 

# Get relevant matrix, transposed so that coefficient names refer to columns
tmp <- as.data.frame(t(coefficient_list$par8_x))
# Get level parameters
incr_lvl <- cbind(tmp[, c("L2", "L3", "L4")], L5 = 1)
# Modify level parameters to incremental
incr_lvl[, 2:4] <- incr_lvl[,2:4]-incr_lvl[, 1:3]

# Calculate 20-parameter variants
form_20 <- as.data.frame(unlist(apply(X = tmp[, DIM.SNAMES], MARGIN = 2, FUN = function(x) incr_lvl*x), recursive = F))
# Set correct names
colnames(form_20) <- ovec(DIM.SNAMES, 2:5, rev = T)
# Add to existing data, transpose back
tmp <- as.data.frame(t(cbind(tmp, form_20)))
# Keep formula attribute
attr(tmp, 'formula') <- attr(coefficient_list$par8_x, 'formula')
# Copy into list
coefficient_list$par8_x <- tmp
# Remove temporary objects
rm(tmp, form_20, incr_lvl)

##########################################################################################

# Make a list with data.frames with predicted values for all 3125 states for each model

predicted_values_list <- lapply(X = coefficient_list, FUN = function(this_model) {
  this_formula <- attr(this_model, "formula")
  outv <- predict.formula(object = this_formula, newdata = all_EQ5D5L_states, args = this_model)
  rownames(outv) <- all_EQ5D5L_states$state_A
  return(outv)
})

# Lazyness
quant_fun <- function(x) c(MEAN = mean(x), SE = sd(x), quantile(x, probs = c("min" = 0, "2.5%" = 0.025, "25%" = 0.25, median = 0.5, "75%" = 0.75, "97.5%" = 0.975, "max" = 1)))
quant_funs <- function(x, MARGIN = 1) as.data.frame(t(apply(X = x, MARGIN = MARGIN, quant_fun)))

# Get mean, SE, and quantiles for each coefficient
(coefficient_statistics <- lapply(X = coefficient_list, FUN = quant_funs))

# Same, but for states, not automatically sent to output because the tables are very long
predicted_values_statistics <- lapply(X = predicted_values_list, FUN = quant_funs)

# Output of statistics for some states of potential interest, for 8-par model
predicted_values_statistics$par8_x[rownames(predicted_values_statistics$par8_x) %in% c("11111", "11112", "11121", "11211", "12111", "21111", "22222", "33333", "44444", "55555"),]

##########################################################################################
# More complex models
##########################################################################################

# Random intercept currently only supported over variable "id", by using the p function cont_r_normal (cont_normal is default)

# 8 parameter model control, linear heteroscedastic, censored, random intercept over variable "id"
par8chr_c <-xregControl(formulas = list(f8i, sigma_est ~ SIGMA_INTERCEPT + SIGMA_SLOPE * Xb), start_values = par8ch_x$coef+0.05, censor_bounds = c(-Inf, 2), name = "TTO", p_fun = cont_r_normal)
# 8 parameter model censored, linear heteroscedastic
(par8chr_x <- xreg(controlList = par8chr_c, dataList = list(TTO = tto)))

#compare same model with and without random intercept:

rbind(cbind(random_estimate = par8chr_x$coef, random_Std.Error = par8chr_x$full_coef[, "Std. Error"], fix_estimate = par8ch_x$coef[names(par8chr_x$coef)], fix_Std.Error = par8ch_x$full_coef[names(par8chr_x$coef), "Std. Error"]),
      logLik = c(-par8chr_x$minima[['total']], NA, -par8ch_x$minima[['total']], NA))

# For custom hybrids, use combined control objects and lists of data. Here is the model above, fitted with regular logit for DCE
# First, make a DCE control object
c8DCE <- xregControl(formulas = appendmodel(object = f8, prepend = "DCE_INTERCEPT + "), p_fun = dich_logistic, name = "DCE")

# Combined control object
c8chr_hybrid <- c(par8chr_c, c8DCE)

# Run model
(x8chr_hybrid <- xreg(controlList = c8chr_hybrid, dataList = list(TTO = tto, DCE = dce)))


