

tariff <- function(data, dim_names = c("mo", "sc", "ua", "pd", "ad") , colname = "state", country = "UK", return_all = FALSE, variant = "5L"){
  data_columns <- ncol(data) +1
  dft <- data[, dim_names]
  dft <- state_from_dims(dft)
  dft <- dummies_state(dft, variant = variant)
  country <- tolower(country)
  if("uk" %in% country) {
    data$uk <- with(dft, 1 
         - (c > 0) * 0.081
         - 0.069 * mo2
         - 0.104 * sc2
         - 0.036 * ua2
         - 0.123 * pd2
         - 0.071 * ad2
         - 0.314 * mo3
         - 0.214 * sc3
         - 0.094 * ua3
         - 0.386 * pd3
         - 0.236 * ad3
         - 0.269 * n3)
    
  }
  if("england" %in% country) {
    data$england <- with(dft, 1 - (
      mo2 * 0.058 + 
        mo3 * 0.076 + 
        mo4 * 0.207 + 
        mo5 * 0.274 + 
        sc2 * 0.05 + 
        sc3 * 0.08 + 
        sc4 * 0.164 + 
        sc5 * 0.203 + 
        ua2 * 0.05 + 
        ua3 * 0.063 + 
        ua4 * 0.162 + 
        ua5 * 0.184 + 
        pd2 * 0.063 + 
        pd3 * 0.084 + 
        pd4 * 0.276 + 
        pd5 * 0.335 + 
        ad2 * 0.078 + 
        ad3 * 0.104 + 
        ad4 * 0.285 + 
        ad5 * 0.289
      
    ))
  }
  if("us" %in% country) {
    data$us <- with(dft,
         1 - (
           mo2 * 0.146016 + 
           sc2 * 0.1753425 + 
           ua2 * 0.1397295 + 
           pd2 * 0.1728907 + 
           ad2 * 0.156223 + 
           mo3 * 0.557685 + 
           sc3 * 0.4711896 + 
           ua3 * 0.3742594 + 
           pd3 * 0.5371011 + 
           ad3 * 0.4501876 + 
           i3 * 0.0106868 + 
           i2_sq * -0.1215579 + 
           i3_sq * -0.0147963 + 
           d1 * -0.1395949)
         )
  }
  if("dk" %in% country) {
    data$dk <- with(dft,
                    1 - 
                    n2 * 0.114 - 
                    mo2 * 0.053 -
                    mo3 * 0.411 - 
                    sc2 * 0.063 - 
                    sc3 * 0.192 - 
                    ua2 * 0.048 - 
                    ua3 * 0.144 - 
                    pd2 * 0.062 - 
                    pd3 * 0.396 - 
                    ad2 * 0.068 - 
                    ad3 * 0.367)
  }
  if("eu_vas" %in% country) {
    data$eu_vas <- with(dft, 1 - 
                      n2 * 0.1279 -
                      n3 * 0.2288 -
                      mo2 * 0.0659 -
                      mo3 * 0.1829 -
                      sc2 * 0.1173 -
                      sc3 * 0.1559 - 
                      ua2 * 0.0265 -
                      ua3 * 0.0860 -
                      pd2 * 0.0930 -
                      pd3 * 0.1637 -
                      ad2 * 0.0891 -
                      ad3 * 0.1290
                    )
  }
  
  if(return_all) {
    return(data)
  } else {
    if(ncol(data)>data_columns) {
      
    } else {
      return(data[, ncol(data)])
    }
    
  }
}

# Function for splitting state variables into dimensions
split_state <- function(data, dim_names = c("mo", "sc", "ua", "pd", "ad"), colname = "state", include_col = FALSE) {
  tmp <- outer(data[,colname], 10^c(4:0), function(a,b) a%/%b %%10)
  colnames(tmp) <- dim_names
  
  # add state column if include_col == TRUE
  if(include_col == TRUE) {
    ret <- cbind(data[, colname], ret)
    colnames(ret)[1] <- colname
  }
  
  return(data.frame(tmp))
}


make_diff_dummies <- function(  data, 
                                dim_names = c("mo","sc","ua","pd","ad"), 
                                colname = "state",
                                prepend = "",
                                append = "",
                                variant = "5L") {
  if(variant == "3L") {
    maxval = 3
  } else if (variant == "5L") {
    maxval = 5
  } else {
    maxval <- max(as.vector(split))  
  }
  
  tmp_df <- data
  
  for(i in maxval:2) {
    
    now_cols <- paste0(dim_names, i)
    prev_cols <- paste0(dim_names, i+1)
    if(i == maxval) new_df <- tmp_df[, now_cols]
    else new_df <- cbind(tmp_df[, prev_cols] + tmp_df[, now_cols], new_df)
    colnames(new_df)[1:length(dim_names)] <- paste0(prepend, now_cols, append)
  }
  tmp_df[, colnames(new_df)] <- new_df
  return(tmp_df)
}

# Function for creating dummy variable matrix
dummies_state <- function(
  data, 
  dim_names = c("mo","sc","ua","pd","ad"), 
  colname = "state", 
  include_dimensions = TRUE, 
  vars_vector = c("d_", "l", "i", "n", "d", "i_sq", "d_sq", "misery", "c", "c_m1", "c_m1_sq", "m1", "c_sq", "c_sqm1", "cr", "crm1", "sum_dim_diff_sq", "variance"),
  keep_columns = NULL,
  variant = "5L",
  prepend = "",
  append = "")  {
  
  all_colnames <- colnames(data)
  
  names_vector <- vars_vector
  names(names_vector) <- vars_vector
  
  rown <- rownames(data)
  
  # create or get dimension matrix
  if(all(dim_names %in% colnames(data))) {
    split <- data[, dim_names]
  } else {
    split <- split_state(data, dim_names, colname)  
  }
  if(variant == "3L") {
    maxval = 3
  } else if (variant == "5L") {
    maxval = 5
  } else {
    maxval <- max(as.vector(split), na.rm = TRUE)  
  }
  
  # determine maximum value (for handling both -3L and -5L)

  # creating dummy names vector
  dummy_names <- as.vector(outer(2:maxval, dim_names, (function(a, b) paste(b, a, sep=""))))
  
  if(all(dummy_names %in% colnames(data))) {
    ret <- data[, dummy_names]
  } else {
    # adding data in all 5 columns to make sure that any missing levels are included
    rownames(split) <- NULL

    split <- rbind(as.matrix(split), matrix(rep(1:maxval, 5), ncol=5))
    
    rownames(split) <- 1:NROW(split)
    # create dummy variable matrix
    ret <- model.matrix(~ as.factor(mo) + as.factor(sc) + as.factor(ua) + as.factor(pd) + as.factor(ad), data.frame(split))
    
    ret <- ret[match(rownames(split), rownames(ret)),]
    
    # removing extra column and rows
    ret <- ret[1:(nrow(ret)-maxval), 2:ncol(ret)]
    
    split <- split[1:(nrow(split)-maxval),]
    # setting column names
    colnames(ret) = dummy_names
  }

  # add dimension variables if include_dimensions is TRUE
  if(include_dimensions == TRUE) {
    ret <- cbind(split, ret)
  }
  
  # add level variables if include_levels is TRUE
  if(length(vars_vector) > 0) {
    ret <- cbind(ret, levels_state(cbind(data[, colname], ret), 
                                   dim_names, 
                                   vars_vector,
                                   colname,
                                   FALSE, 
                                   variant = variant))
  }
  
  # keep columns from data
  if(!is.null(keep_columns)) {
    ret <- cbind(data[, keep_columns], ret)
    colnames(ret)[1:length(keep_columns)] <- keep_columns
  }
  rownames(ret) <- rown
  
  tmpcolnamesnum <- !colnames(ret) %in% all_colnames 
  colnames(ret)[tmpcolnamesnum] <- paste0(prepend, colnames(ret)[tmpcolnamesnum], append)
  
  
  # returning
  return(data.frame(ret))
}

# Function for splitting state variables into l2, l3... level variables. 
# If presence_dummies == TRUE, dummy variables representing the presence or absence of problems at level n will be added. 
levels_state <- function(data, 
                         dim_names = c("mo","sc","ua","pd","ad"), 
                         vars_vector = c("d_", "l", "i", "n", "d", "i_sq", "d_sq", "misery", "c", "c_m1", "c_m1_sq", "m1", "c_sq", "c_sqm1", "cr", "crm1", "sum_dim_diff_sq", "variance"), 
                         colname = "state", 
                         include_col = FALSE,
                         variant = "5L",
                         ...) {
  
  names_vector <- vars_vector
  names(names_vector) <- vars_vector
  
  # get dummies
  dummy_matrix <- dummies_state(data, dim_names, colname, TRUE, NULL, ...)
  # determine max value for dimensions (-3L or -5L)
  if(variant == "3L") {
    maxval = 3
  } else if (variant == "5L") {
    maxval = 5
  }
  # create initial column to enable cbind
  ret <- rep(NA, nrow(dummy_matrix))
  
  # incremental dummies
  if("d_" %in% vars_vector) {
    colprefix = names_vector["d_"]
    for(dim in dim_names) {
      for(i in 2:maxval) {
        ret <- cbind(ret, rowSums(dummy_matrix[, paste0(dim, i:maxval), F]))
        colnames(ret)[ncol(ret)] <- paste0(colprefix, dim, i)
      }
    }
  }
  
  
  # generate columns with number of dimensions at level x
  if("l" %in% vars_vector) {
    colprefix = names_vector["l"]
    for(i in 2:maxval) {
      ret <- cbind(ret, rowSums(dummy_matrix[,paste0(dim_names, i), F]))
      colnames(ret)[ncol(ret)] <- paste(colprefix, i, sep = "_")
    }  
  }
  
  
  # generate columns representing the presence of problems at level x
  if("n" %in% vars_vector) {
    current_name = names_vector["n"]
   for(i in 2:maxval) {
     ret <- cbind(ret, (rowSums(dummy_matrix[,paste(dim_names, i, sep="")]) > 0)*1)
     colnames(ret)[ncol(ret)] <- paste(current_name, i, sep = "")
   } 
   if(maxval == 5) {
     ret <- cbind(ret, (rowSums(ret[, paste(current_name, 2:3, sep="")])>0)*1, (rowSums(ret[, paste(current_name, 4:5, sep = "")])>0)*1)
     colnames(ret)[((ncol(ret)-1):ncol(ret))] <- paste(current_name, c(23, 45), sep="")
   }
  }
  
  if("i" %in% vars_vector) {
    current_name = names_vector["i"]
    for(i in 2:maxval) {
      tmp <- rowSums(dummy_matrix[,paste(dim_names, i, sep="")])-1
      tmp[tmp < 0] <- 0
      ret <- cbind(ret, tmp)
      colnames(ret)[ncol(ret)] <- paste(current_name, i, sep = "")
    }
  }

  if("i_sq" %in% vars_vector) {
    current_name = names_vector["i_sq"]
    for(i in 2:maxval) {
      tmp <- rowSums(dummy_matrix[,paste(dim_names, i, sep="")])-1
      tmp[tmp < 0] <- 0
      tmp <- tmp^2
      ret <- cbind(ret, tmp)
      colnames(ret)[ncol(ret)] <- paste("i", i, "_sq", sep = "")
    }
  }
  
  if("d" %in% vars_vector) {
    current_name <- names_vector["d"]
    for(i in 2:maxval) {
      tmp <- rowSums(dummy_matrix[, as.vector(outer(dim_names, i:maxval, function(a, b) paste(a, b, sep = "")))])-1
      tmp[tmp < 0] <- 0
      ret <- cbind(ret, tmp)
      colnames(ret)[ncol(ret)] <- paste(current_name, i-1, sep = "")
    }
  }

  
  if("d_sq" %in% vars_vector) {
    current_name <- names_vector["d_sq"]
    for(i in 2:maxval) {
      tmp <- rowSums(dummy_matrix[, as.vector(outer(dim_names, i:maxval, function(a, b) paste(a, b, sep = "")))])-1
      tmp[tmp < 0] <- 0
      tmp <- tmp^2
      ret <- cbind(ret, tmp)
      colnames(ret)[ncol(ret)] <- paste(current_name, i-1, sep = "")
    }
  }
  
  if("misery" %in% vars_vector) {
    tmp <- rowSums(dummy_matrix[, dim_names])
    ret <- cbind(ret, tmp)
    colnames(ret)[ncol(ret)] <- names_vector["misery"]
  }
  
  if("c" %in% vars_vector) {
    tmp <- rowSums(dummy_matrix[, dim_names])
    tmp <- tmp -5
    ret <- cbind(ret, tmp)
    colnames(ret)[ncol(ret)] <- names_vector["c"]
  }  

  if("c_m1" %in% vars_vector) {
    tmp <- rowSums(dummy_matrix[, dim_names])
    tmp <- tmp -6
    tmp[tmp < 0] <- 0
    ret <- cbind(ret, tmp)
    colnames(ret)[ncol(ret)] <- names_vector["c_m1"]
  }  

  if("c_m1_sq" %in% vars_vector) {
    tmp <- rowSums(dummy_matrix[, dim_names])
    tmp <- tmp -6
    tmp[tmp < 0] <- 0
    tmp <- tmp^2
    ret <- cbind(ret, tmp)
    colnames(ret)[ncol(ret)] <- names_vector["c_m1_sq"]
  }
  
  if("m1" %in% vars_vector) {
    tmp <- (rowSums(dummy_matrix[, dim_names]) == 6)*1
    ret <- cbind(ret, tmp)
    colnames(ret)[ncol(ret)] <- names_vector["m1"]
  }

  if("c_sq" %in% vars_vector) {
    tmp <- rowSums(dummy_matrix[, dim_names])
    tmp <- (tmp -5)^2
    ret <- cbind(ret, tmp)
    colnames(ret)[ncol(ret)] <- names_vector["c_sq"]
  }   
  if("c_sqm1" %in% vars_vector) {
    tmp <- rowSums(dummy_matrix[, dim_names])
    tmp <- (tmp -5)^2 -1
    ret <- cbind(ret, tmp)
    colnames(ret)[ncol(ret)] <- names_vector["c_sqm1"]
  }  
  if("cr" %in% vars_vector) {
    tmp <- rowSums(dummy_matrix[, dim_names])
    tmp <- (tmp -5)^0.5
    ret <- cbind(ret, tmp)
    colnames(ret)[ncol(ret)] <- names_vector["cr"]
  }  
  if("crm1" %in% vars_vector) {
    tmp <- rowSums(dummy_matrix[, dim_names])
    tmp <- (tmp -5)^0.5 - 1
    ret <- cbind(ret, tmp)
    colnames(ret)[ncol(ret)] <- names_vector["crm1"]
  }  
  if("sum_dim_diff_sq" %in% vars_vector) {
    tmp <- rowSums((dummy_matrix[, dim_names]-1)^2)
    ret <- cbind(ret, tmp)
    colnames(ret)[ncol(ret)] <- names_vector["sum_dim_diff_sq"]
  }  

  if("variance" %in% vars_vector) {
    tmp <- apply(dummy_matrix[, dim_names], MARGIN = 1, FUN = var)
    ret <- cbind(ret, tmp)
    colnames(ret)[ncol(ret)] <- names_vector["variance"]
  }
  
  ret <- ret[, 2:ncol(ret)]
  
  # add state column if include_col == TRUE
  if(include_col == TRUE) {
    ret <- cbind(data[, colname], ret)
    colnames(ret)[1] <- colname
  }
  
  return(data.frame(ret))
}


# Create state variable from dimension variables
state_from_dims <- function(data, dim_names = c("mo","sc","ua","pd","ad"), colname = "state", stateonly = FALSE) {
  data <- data.frame(data)
  if(!colname %in% colnames(data) & all(dim_names %in% colnames(data))) {
    data[, colname] <- as.integer(do.call(paste, c(as.data.frame(data[,dim_names]), sep="")))
  }
  if(stateonly) {
    ret <- data.frame(data[, colname])
    colnames(ret) <- colname
    return(ret)  
  }else{
    return(data.frame(data))
  }
}

TTO_import <- function(filename, included_variables, dim.snames = c("mo", "sc", "ua", "pd", "ad"), removestates = c(0, 21121, 35554, 15411)) {
  
  # read data from csv
  TTO_data <- read.csv3(filename)

  # drop variables not listed in included_variables
  TTO_data <- TTO_data[, colnames(TTO_data) %in% included_variables]

  # create state if missing
  if(!"state" %in% colnames(TTO_data) & all(dim.snames %in% colnames(TTO_data))) {
    TTO_data$state <- as.integer(do.call(paste, c(as.data.frame(TTO_data[,dim.snames]), sep="")))
  }
  
  # create mo, sc, ua, pd, ad if missing and state is present
  TTO_data <- state_from_dims(TTO_data)

  # check if data appears to be -3L or -5L. True if max observed level is >3. 
  EQ5L <- max(as.vector(TTO_data[, c("mo", "sc", "ua", "pd", "ad")])) >3
  
  # remove rows corresponding to practice states
  #removestates = c(0, 21121, 35554, 15411)
  TTO_data <- TTO_data[!TTO_data$state %in% removestates,]

  # create pseudonyme id if id is missing
  if(!"id" %in% colnames(TTO_data)) {
    TTO_data$id <- 1:NROW(TTO_data)
    if(NROW(TTO_data)> 100)  TTO_data$id <- floor((as.integer(rownames(TTO_data))-1)/10)+1  
  }
  
  # recreate block information if missing
  if(!"block" %in% colnames(TTO_data)) {
    # Hack to generate blocks on the basis of one state from each block
    blocks <- data.frame(state = c(11221, 11235, 12514, 34515, 35245, 45144, 51451, 54231, 12121, 12543, 23514, 32443, 34155, 43542, 45133, 52215, 11421, 12244, 13313, 25122, 31525, 45233, 52455, 55233, 12344, 12513, 14554, 21112, 44125, 44345, 53221, 54342, 14113, 15151, 21315, 24443, 31524, 43315, 52431, 54153, 11212, 12112, 21345, 23152, 34244, 43514, 44553, 55424, 11425, 13122, 22434, 24553, 35332, 42115, 45413, 51152, 12334, 21334, 23242, 24342, 32314, 33253, 53412, 55225, 11414, 21444, 25222, 25331, 31514, 35143, 53243, 53244, 11122, 13224, 24445, 34232, 35311, 42321, 43555, 52335), blocks = rep(1:10, each = 8))
    
    #blocks <- c(11221,12121,11421,21112,14113,12112,13122,21334,11414,11122)
    
    TTO_data$block <- (blocks[match(TTO_data$state, blocks$state), "blocks"] -> tmpblocks)[!is.na(tmpblocks)->tmpincl][match(TTO_data$id, TTO_data[tmpincl, "id"])]
    rm(tmpblocks, tmpincl, blocks)
    
  }

  # If variable "included" is found in the CSV, only records witha 1 are kept.
  if("included" %in% colnames(TTO_data)){
    TTO_data <- TTO_data[TTO_data$included == 1 & !is.na(TTO_data$included),]
    TTO_data <- subset(TTO_data, select= -included)
  }

  # misery index
  TTO_data$misery <- rowSums(TTO_data[, dim.snames])

  TTO_data <- TTO_data[!is.na(rowSums(TTO_data[,-which(colnames(TTO_data) == "block")])) & TTO_data$tto <= 20 & TTO_data$tto >= -1,]

  # Determine if the tto-values are between 0 and 20 or between -1 and 0, and make disutilities
  temp <- TTO_data$tto
  if(sum(temp>=-1 & temp <= 1)>sum(temp>=0 & temp<= 20)) {
    TTO_data$tto<- -(TTO_data$tto-1)
  } else {
    if(cor(TTO_data$tto, TTO_data$misery)<0) {
      TTO_data$tto <- -(TTO_data$tto-20)/10
    }
  }
  
  #remove any values outside the appropriate scale.
  TTO_data <- TTO_data[TTO_data$tto >= 0 & TTO_data$tto <=2,]
  return(TTO_data)
}


genLatin <- function(n) array((rep(sample(n), each=n) + sample(n))%%n+1, dim = c(n, n))[sample(n), sample(n)]

add_pseudo_block <- function(data, pseudoblockname = "pseudo_block", blockname = "block")
{
  data[, pseudoblockname] <- NA
  data[data$c == 1, pseudoblockname] <- 1
  data[data$c == 20, pseudoblockname] <- 10
  data[is.na(data[, pseudoblockname]),] <- balanced_subsample(data[is.na(data[, pseudoblockname]),], "state", blockname, 1:10, pseudoblockname)
  data[, pseudoblockname] <- genLatin(10)[as.matrix(data[, c(blockname, pseudoblockname)])]
  data
}

diff_EQ_dummies <- function(data, 
                            vector_a = NULL, 
                            vector_b = NULL,
                            dim_names = c("mo","sc","ua","pd","ad"),
                            vars_vector = c("d_", "l", "i", "n", "d", "i_sq", "d_sq", "misery", "c", "c_m1", "c_m1_sq", "m1", "c_sq", "c_sqm1", "cr", "crm1", "sum_dim_diff_sq", "variance"),
                            keep_columns = NULL,
                            append = NULL,
                            prepend = NULL,
                            posdiff = T) {
  
  
  if(is.null(vector_a) & is.null(vector_b) & NCOL(data) ==2) {
    vector_a <- colnames(data)[1]
    vector_a <- colnames(data)[2]
  }
  
  all_colnames <- colnames(data)
  rown <- rownames(data)
  
  df_a <- dummies_state(data[, vector_a, F],colname = vector_a, vars_vector = vars_vector)
  df_b <- dummies_state(data[, vector_b, F],colname = vector_b, vars_vector = vars_vector)
  ret <- df_a-df_b
  
  colnames(ret) <- paste0(prepend, colnames(ret), append)
  
  if(posdiff) {
    
    ret_a <- ret
    ret_a[ret_a <0] <- 0
    ret_b <- -ret
    ret_b[ret_b <0] <- 0
    
    colnames(ret_a) <- paste0(colnames(ret), "_A")
    colnames(ret_b) <- paste0(colnames(ret), "_B")
    
    ret <- cbind(ret, ret_a, ret_b)
  }
  
  
  # keep columns from data
  if(!is.null(keep_columns)) {
    ret <- cbind(data[, keep_columns], ret)
    colnames(ret)[1:length(keep_columns)] <- keep_columns
  }
  rownames(ret) <- rown
  
  
  
  return(ret)
  
}


