#' Randomize blocks
#'
#' This function creates a set of randomized blocks for future use. It assumes
#' only two blocks with a randomized trial order and block order designated by 
#' the operator.
#'
#' @param num_rands Number of desired randomizations
#' @param Blocks Dataframes of the desired blocks to be randomized, input as c(1,2)
#' @return A list that contains a dataframe for each randomization in the order specified
#' @export
rand_blocks <- function(num_rands, Blocks){
  list_randomizations <- list()
  for (x in 1:num_rands) {
    rand <- rbind(Blocks[[1]][sample(nrow(Blocks[[1]]), nrow(Blocks[[1]])), ], 
                  Blocks[[2]][sample(nrow(Blocks[[2]]), nrow(Blocks[[2]])), ])
    rownames(rand) <- c()
    rand$Trial <- as.numeric(rownames(rand))
    list_randomizations[[x]] <- rand
  }
  return(list_randomizations)
}

#' Update beliefs
#'
#' This function runs the input randomizations through the belief_update function
#' from David Kleinschmidt's beliefupdatr package. It assumes two perceptual categories
#' that vary along a single acoustic parameter. These categories are assumed to form
#' normal distributions.
#'
#' @param randomizations Takes a list of randomized trials.
#' @param priorParamsG A series of parameters associated with the G distribution.
#' Entered as c(1,2,3,4). Where 1 = mean, 2 = variance, 3 = mean psuedocount, and
#' 4 = variance psuedocount.
#' @param priorParamsK A series of parameters associated with the K distribution.
#' Entered as c(1,2,3,4). Where 1 = mean, 2 = variance, 3 = mean psuedocount, and
#' 4 = variance psuedocount.
#' @return A list of the randomizations in the same order as presented with an
#' additional column corresponding to the beliefs updated at each trial.
#' @export
updateBeliefs <- function(randomizations, priorParamsG, priorParamsK){
  list_beliefs <- list()
  for (x in 1:length(randomizations)) {
    list_beliefs[[x]] <- as.data.frame(beliefupdatr::belief_update(randomizations[[x]], 'rawVOT', 'Label',
                                                                   list(G = beliefupdatr::nix2_params(priorParamsG[1], priorParamsG[2],
                                                                                                      priorParamsG[3], priorParamsG[4]),
                                                                        K = beliefupdatr::nix2_params(priorParamsK[1], priorParamsK[2],
                                                                                                      priorParamsK[3], priorParamsK[4]))))
  }
  return(list_beliefs)
}

#' Extract mu and sigma
#'
#' This function extracts the mu and sigma associated with each randomization.
#'
#' @param beliefs Takes the output of the updateBeliefs function.
#' @return A list of the mu and sigma of the respective distributions at each
#' trial for each randomization.
#' @export
extractMuSigma <- function(beliefs){
  list_params <- list()
  for (x in 1:length(beliefs)) {
    BUpdates_var <- as.data.frame(lapply(unlist(beliefs[[x]]$beliefs, recursive = FALSE), `[[`, "sigma2"))
    G_updates <- dplyr::select(BUpdates_var, tidyselect::contains("G"))
    G_var_long <- as.data.frame(matrix(t(G_updates), ncol=1))
    colnames(G_var_long) <- c("G_var")
    K_updates <- dplyr::select(BUpdates_var, tidyselect::contains("K"))
    K_var_long <- as.data.frame(matrix(t(K_updates), ncol=1))
    colnames(K_var_long) <- c("K_var")
    
    BUpdates_bound <- as.data.frame(lapply(unlist(beliefs[[x]]$beliefs, recursive = FALSE), `[[`, "mu"))
    G_updates_b <- dplyr::select(BUpdates_bound, tidyselect::contains("G"))
    G_mu_long <- as.data.frame(matrix(t(G_updates_b), ncol=1))
    colnames(G_mu_long) <- c("G_mu")
    K_updates_b <- dplyr::select(BUpdates_bound, tidyselect::contains("K"))
    K_mu_long <- as.data.frame(matrix(t(K_updates_b), ncol=1))
    colnames(K_mu_long) <- c("K_mu")
    
    Param_updates <- cbind(G_mu_long, K_mu_long, G_var_long, K_var_long, 
                           beliefs[[x]]$Distribution, beliefs[[x]]$rawVOT,
                           beliefs[[x]]$Trial, beliefs[[x]]$Label)
    colnames(Param_updates) <- c("G_mu", "K_mu", "G_var", "K_var", 
                                 "Distribution", "VOT", "trial_number", 
                                 "Response")
    Param_updates <- dplyr::mutate(Param_updates, G_sd = sqrt(G_var), K_sd = sqrt(K_var))
    
    list_params[[x]] <- Param_updates
  }
  return(list_params)
}

#' Select specific trials
#'
#' This function allows the user to select specific trials of interest.
#'
#' @param paramUpdates Takes the output of the extractMuSigma function.
#' @param trials A list of the trials of interest input as c(1,2,3,etc.).
#' @return A list of the mu and sigma at each designated trial for each
#' randomization.
#' @export
selectTrials <- function(paramUpdates, trials){
  list_trials <- list()
  for (x in 1:length(paramUpdates)) {
    trialParams <- dplyr::filter(paramUpdates[[x]],
                          trial_number %in% trials)
    list_trials[[x]] <- trialParams
  }
  return(list_trials)
}

#' Extract slope
#'
#' This function allows the user extract slopes of the internal category
#' distributions. It first calculates the distributions of the internal
#' categories along the designated VOT space. Then uses these distributions
#' to create an indentification function using a simplied Bayes rule. Then
#' takes the implicit differential along the ID function and pulls the slope
#' closest to the provided boundary percentile.
#'
#' @param trialParameters Takes the output of the selectTrials function.
#' @param VOTspace A list of the parameters for the VOT space. Input as
#' c(1,2,3), where 1 = minimum, 2 = maximum, and 3 = step count.
#' @param boundary The designated percentile at with to extract the slope.
#' This is typically done at the boundary, designated by 0.5. The function
#' pulls the slopes closest to the designated percentile, which may not 
#' always be exact.
#' @return A dataframe of the extracted slopes for each randomization at 
#' each trial. The dataframe inclues the associated VOT, percentile, trial
#' number and randomization. Also returns three lists of dataframes that 
#' contains to the probabilities generated by the simplified Bayes rule at
#' the trials indicated in the selectTrials function.
#' @export
slopeExtract <- function(trialParameters, VOTspace, boundary){
  Probs.1 <<- list()
  Probs.2 <<- list()
  Probs.3 <<- list()
  linspace <- seq(VOTspace[1], VOTspace[2], length.out = VOTspace[3])
  slopes <- as.data.frame(matrix(ncol = ncol(trialParameters[[1]]), nrow = 0))
  for (x in 1:length(trialParameters)) {
    for (r in 1:nrow(trialParameters[[x]])) {
      K = list()
      G = list()
      bayes_curvesK = list()
      for (i in linspace){
        G <- c(G, dnorm(i, 
                        mean = trialParameters[[x]][r,]$G_mu, 
                        sd = trialParameters[[x]][r,]$G_sd))
      }
      
      for (i in linspace){
        K <- c(K, dnorm(i, 
                        mean = trialParameters[[x]][r,]$K_mu, 
                        sd = trialParameters[[x]][r,]$K_sd))
      }
      
      for (i in 1:length(K)){
        bayes_curvesK <- c(bayes_curvesK, K[[i]] / (K[[i]] + G[[i]]))
      }
      df_bayes <- do.call(rbind, Map(data.frame, Linspace=linspace,
                                     as_K=bayes_curvesK))
      slopes_ <- diff(df_bayes$as_K)/diff(df_bayes$Linspace)
      # add NA to end of slopes in order to match lengths
      slopes_ <- append(slopes_, NA)
      df_bayes <- do.call(rbind, Map(data.frame, Linspace=linspace,
                                     as_K=bayes_curvesK, slope=slopes_))
      df_bayes$ID <- as.character(x)
      df_bayes$Trial <- as.numeric(trialParameters[[x]][r,]$trial_number)
      if (r==1) {
        Probs.1[[x]] <<- df_bayes 
      } else if (r==2) {
        Probs.2[[x]] <<- df_bayes
      } else {
        Probs.3[[x]] <<- df_bayes
      }
      # get slope closest to 0.5 boundary
      slopes <- rbind(slopes, dplyr::filter(df_bayes, abs(as_K - boundary) == min(abs(as_K - boundary))))
    }
  }
  return(slopes)
}

#' Standard error of the mean
#'
#' This function calculates the standard error of the mean 
#' given a series of numbers.
#'
#' @param x A series of number with which to calculate the standard error.
#' @return The standard error of the mean.
#' @export
se <- function(x) {
  sd(x)/sqrt(length(x))
}
