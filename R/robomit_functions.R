################################################################################
# R functions of the robomit R package                                         #
# -----------------------------------------------------------------------------#
# 2020, ETH Zurich                                                             #
# developed by Sergei Schaub                                                   #
# e-mail: seschaub#@ethz.ch                                                    #
# first version: August, 7, 2020                                               #
# lst update: August, 12, 2020                                                 #
################################################################################


################################################################################
#------------------------------- 0. settings
################################################################################
#' @importFrom  plm plm pdata.frame
#' @import ggplot2
#' @import dplyr
#' @import broom
#' @import tidyr
#' @import tibble
#' @importFrom  stats as.formula lm sd var

utils::globalVariables(c("Beta", "..density..", "Rmax","Delta","na.exclude","dnorm"))


################################################################################
#------------------------------- 1. functions related to Oster (2019)
################################################################################

##################--------------------------------------------------------------
#------------------------------- 1.1 o_delta - function 1
##################--------------------------------------------------------------
#' @title Delta*
#'
#' @description Estimates delta*, i.e. the degree of selection on unobservables relative to observables that would be necessary to explain away the result, following Oster (2019).
#' @usage o_delta(y, x, con, id = "none", time = "none", beta = 0, R2max, type, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent variable of interest (treatment variable; as string).
#' @param con Name of the other control variables. Provided as string in the format: "w + z +...".
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect models.
#' @param time Name of the time variable (e.g. year or month; as string). Only applicable for fixed effect models.
#' @param beta Beta for which delta* should be estimated (default is beta = 0).
#' @param R2max Max R-square for which beta* should be estimated.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param data Data.
#' @details Estimates delta*, i.e. the degree of selection on unobservables relative to observables that would be necessary to explain away the result, following Oster (2019). The function supports linear cross sectional (see \emph{lm} objects in R) and panel fixed effect (see \emph{plm} objects in R) models.
#' @return Returns tibble object. Including delta* and various other information.
#' @references Oster, E. (2019). Unobservable selection and coefficient stability: Theory and evidence. Journal of Business & Economic Statistics, 37, 187-204.
#' @examples
#' # load data, e.g. the in-build mtcars dataset
#' data("mtcars")
#' data_oster <- mtcars
#'
#' # preview of data
#' head(data_oster)
#'
#' # load robomit
#' require(robomit)
#'
#' # estimate delta*
#' o_delta(y = "mpg",               # define the dependent variable name
#'             x = "wt",            # define the main independent variable name
#'             con = "hp + qsec",   # other control variables
#'             beta = 0,            # define beta. This is usually set to 0
#'             R2max = 0.9,         # define the max R-square.
#'             type = "lm",         # define model type
#'             data = data_oster)   # define dataset
#' @export

o_delta <- function(y, x, con,  id = "none", time = "none", beta = 0, R2max, type, data) {


  # rename variables and create, if type is plm, a panel dataset
  ## lm
  if (type == "lm") {data <- data %>% dplyr::rename(y_var = y, x_var = x)}
  ## plm
  if (type == "plm") {data <- data %>% dplyr::rename(y_var = y, x_var = x, id_var = id, time_var = time)
  data_plm <- pdata.frame(data, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)}


  # define formulas of the different models
  model0_formula    <- as.formula("y_var ~ x_var")              # uncontrolled model
  model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) # controlled model
  aux_model_formula <- as.formula(paste("x_var ~",con))         # auxiliary model


  # define model to compute sigma_xx when model type is plm
  if (type == "plm") {sigma_xx_model_formula <- as.formula(paste("x_var ~ factor(id_var)"))}


  # run models
  ## lm
  if (type == "lm") {
    model0    <- lm(model0_formula,    data = data, na.action = na.exclude) # run uncontrolled model
    model1    <- lm(model1_formula,    data = data, na.action = na.exclude) # run controlled model
    aux_model <- lm(aux_model_formula, data = data, na.action = na.exclude) # run auxiliary model
  }
  ## plm
  if (type == "plm") {
    model0    <- plm(model0_formula,    data = data_plm, model = "within", na.action = na.exclude) # run uncontrolled model
    model1    <- plm(model1_formula,    data = data_plm, model = "within", na.action = na.exclude) # run controlled model
    aux_model <- plm(aux_model_formula, data = data_plm, model = "within", na.action = na.exclude) # run auxiliary model
    model_xx  <-  lm(sigma_xx_model_formula, data = data, na.action = na.exclude)                  # run model to obtain singa_xx
  }


  # variables based on model outputs
  if (type == "lm") {b0 = as.numeric(tidy(model0)[2,2])} else if (type == "plm") {b0 = as.numeric(tidy(model0)[1,2])}           # beta uncontrolled model
  if (type == "lm") {b1 = as.numeric(tidy(model1)[2,2])} else if (type == "plm") {b1 = as.numeric(tidy(model1)[1,2])}           # beta controlled model
  if (type == "lm") {R20 = summary(model0)$r.squared} else if (type == "plm") {R20 = as.numeric(summary(model0)$r.squared[1])}  # R-square uncontrolled model
  if (type == "lm") {R21 = summary(model1)$r.squared} else if (type == "plm") {R21 = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model
  sigma_yy = var(data$y_var, na.rm = T)                                                                                         # variance of dependent variable
  if (type == "lm") {sigma_xx = var(data$x_var, na.rm = T)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)}       # variance of independent variable
  t_x  = var(aux_model$residuals)


  # create some additional variables
  bt_m_b = b1 - beta
  rt_m_ro_t_syy = (R21-R20) * sigma_yy
  b0_m_b1 = b0 - b1
  rm_m_rt_t_syy = (R2max - R21) * sigma_yy


  # compute numerator, in steps
  num1 = bt_m_b * rt_m_ro_t_syy * t_x
  num2 = bt_m_b * sigma_xx * t_x * b0_m_b1^2
  num3 = 2 * bt_m_b^2 * (t_x * b0_m_b1 * sigma_xx)
  num4 = bt_m_b^3 * (t_x * sigma_xx - t_x^2)
  num  = num1 + num2 + num3 + num4


  # compute denominator, in steps
  den1 = rm_m_rt_t_syy * b0_m_b1 * sigma_xx
  den2 = bt_m_b * rm_m_rt_t_syy * (sigma_xx - t_x)
  den3 = bt_m_b^2 * (t_x * b0_m_b1 * sigma_xx)
  den4 = bt_m_b^3 * (t_x * sigma_xx - t_x^2)
  den  = den1 + den2 + den3 + den4


  # finally compute delta_star
  delta_star = num/den


  # create triblle object that contains results
  result_delta <- tribble(
    ~Name, ~Value,
    "Delta*",                   round(delta_star,6),
    "Uncontrolled Coefficient", b0,
    "Controlled Coefficient",   b1,
    "Uncontrolled R-square",    R20,
    "Controlled R-square",      R21,
    "Max R-square",             R2max,
    "Beta hat (defined)",       beta
  )

  # warning if max R-square is smaller than the R-square of the control model
  if (R21 > R2max)
    warning("The max R-square value is smaller than the R-square of the controlled model")


  # define return object
  return(result_delta)
}



##################--------------------------------------------------------------
#------------------------------- 1.2 o_delta_rsq - function 2
##################--------------------------------------------------------------
#' @title Deltas* over a range of max R-squares
#' @description Estimates deltas*, i.e. the degree of selection on unobservables relative to observables that would be necessary to explain away the result, following Oster (2019) over a range of max R-squares following Oster (2019).
#' @usage o_delta_rsq(y, x, con, id = "none", time = "none", beta = 0, type, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent variable of interest (treatment variable; as string).
#' @param con Name of the other control variables. Provided as string in the format: "w + z +...".
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect models.
#' @param time Name of the time variable (e.g. year or month; as string). Only applicable for fixed effect models.
#' @param beta Beta for which delta* should be estimated (default is beta = 0).
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param data Data.
#' @details Estimates deltas*, i.e. the degree of selection on unobservables relative to observables that would be necessary to explain away the result, following Oster (2019) over a range of max R-squares. The range of max R-squares starts from the R-square of the controlled model rounded up to the next 1/100 to 1. The function supports linear cross sectional (see \emph{lm} objects in R) and panel fixed effect (see \emph{plm} objects in R) models.
#' @return Returns tibble object. Including deltas* over a range of max R-squares.
#' @references Oster, E. (2019). Unobservable selection and coefficient stability: Theory and evidence. Journal of Business & Economic Statistics, 37, 187-204.
#' @examples
#' # load data, e.g. the in-build mtcars dataset
#' data("mtcars")
#' data_oster <- mtcars
#'
#' # preview of data
#' head(data_oster)
#'
#' # load robomit
#' require(robomit)
#'
#' # estimate deltas* over a range of max R-squares
#' o_delta_rsq(y = "mpg",               # define the dependent variable name
#'                 x = "wt",            # define the main independent variable name
#'                 con = "hp + qsec",   # other control variables
#'                 beta = 0,            # define beta. This is usually set to 0
#'                 type = "lm",         # define model type
#'                 data = data_oster)   # define dataset
#' @export

o_delta_rsq <- function(y, x, con, id = "none", time = "none", beta = 0, type, data) {


  # rename variables and create, if type is plm, a panel dataset
  ## lm
  if (type == "lm") {data <- data %>% dplyr::rename(y_var = y, x_var = x)}
  ## plm
  if (type == "plm") {data <- data %>% dplyr::rename(y_var = y, x_var = x, id_var = id, time_var = time)
  data_plm <- pdata.frame(data, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)}


  # define formulas of the different models
  model0_formula    <- as.formula("y_var ~ x_var")              # uncontrolled model
  model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) # controlled model
  aux_model_formula <- as.formula(paste("x_var ~",con))         # auxiliary model


  # define model to compute sigma_xx when model type is plm
  if (type == "plm") {sigma_xx_model_formula <- as.formula(paste("x_var ~ factor(id_var)"))}


  # run models
  ## lm
  if (type == "lm") {
    model0    <- lm(model0_formula,    data = data, na.action = na.exclude) # run uncontrolled model
    model1    <- lm(model1_formula,    data = data, na.action = na.exclude) # run controlled model
    aux_model <- lm(aux_model_formula, data = data, na.action = na.exclude) # run auxiliary model
  }
  ## plm
  if (type == "plm") {
    model0    <- plm(model0_formula,    data = data_plm, model = "within", na.action = na.exclude) # run uncontrolled model
    model1    <- plm(model1_formula,    data = data_plm, model = "within", na.action = na.exclude) # run controlled model
    aux_model <- plm(aux_model_formula, data = data_plm, model = "within", na.action = na.exclude) # run auxiliary model
    model_xx  <-  lm(sigma_xx_model_formula, data = data, na.action = na.exclude)                  # run model to obtain singa_xx
  }


  # variables based on model outputs
  if (type == "lm") {b0 = as.numeric(tidy(model0)[2,2])} else if (type == "plm") {b0 = as.numeric(tidy(model0)[1,2])}           # beta uncontrolled model
  if (type == "lm") {b1 = as.numeric(tidy(model1)[2,2])} else if (type == "plm") {b1 = as.numeric(tidy(model1)[1,2])}           # beta controlled model
  if (type == "lm") {R20 = summary(model0)$r.squared} else if (type == "plm") {R20 = as.numeric(summary(model0)$r.squared[1])}  # R-square uncontrolled model
  if (type == "lm") {R21 = summary(model1)$r.squared} else if (type == "plm") {R21 = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model
  sigma_yy = var(data$y_var, na.rm = T)                                                                                         # variance of dependent variable
  if (type == "lm") {sigma_xx = var(data$x_var, na.rm = T)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)}       # variance of independent variable
  t_x  = var(aux_model$residuals)


  # define max R-square range
  r_start <- ceiling(R21/0.01)*0.01 # start max R-square at the next highest 1/100
  n_runs <- length(seq(r_start:1,by = 0.01))
  delta__over_rsq_results <- tibble(
    "x" = 1:n_runs,
    "Rmax" = 0,
    "Delta*" = 0) # create tibble object to store results


  # run loop over different max R-squares
  for (i in 1:n_runs) {
    R2max = seq(r_start:1,by = 0.01)[i] # current max R-square


    # create some additional variables
    bt_m_b = b1 - beta
    rt_m_ro_t_syy = (R21-R20) * sigma_yy
    b0_m_b1 = b0 - b1
    rm_m_rt_t_syy = (R2max - R21) * sigma_yy


    # compute numerator, in steps
    num1 = bt_m_b * rt_m_ro_t_syy * t_x
    num2 = bt_m_b * sigma_xx * t_x * b0_m_b1^2
    num3 = 2 * bt_m_b^2 * (t_x * b0_m_b1 * sigma_xx)
    num4 = bt_m_b^3 * (t_x * sigma_xx - t_x^2)
    num  = num1 + num2 + num3 + num4


    # compute denominator, in steps
    den1 = rm_m_rt_t_syy * b0_m_b1 * sigma_xx
    den2 = bt_m_b * rm_m_rt_t_syy * (sigma_xx - t_x)
    den3 = bt_m_b^2 * (t_x * b0_m_b1 * sigma_xx)
    den4 = bt_m_b^3 * (t_x * sigma_xx - t_x^2)
    den  = den1 + den2 + den3 + den4


    # finally compute delta_star
    delta_star = num/den


    # store results
    delta__over_rsq_results[i,2] <- R2max
    delta__over_rsq_results[i,3] <- round(delta_star,6)
  }


  # define return object
  return(delta__over_rsq_results)
}




##################--------------------------------------------------------------
#------------------------------- 1.3 o_delta_rsq_viz - function 3
##################--------------------------------------------------------------
#' @title Visualization of deltas* over a range of max R-squares
#'
#' @description Estimates and visualizes deltas*, i.e. the degree of selection on unobservables relative to observables that would be necessary to explain away the result, following Oster (2019) over a range of max R-squares.
#' @usage o_delta_rsq_viz(y, x, con, id = "none", time = "none", beta = 0, type, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent variable of interest (treatment variable; as string).
#' @param con Name of the other control variables. Provided as string in the format: "w + z +...".
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect models.
#' @param time Name of the time variable (e.g. year or month; as string). Only applicable for fixed effect models.
#' @param beta Beta for which delta* should be estimated (default is beta = 0).
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param data Data.
#' @details Estimates and visualizes deltas*, i.e. the degree of selection on unobservables relative to observables that would be necessary to explain away the result, following Oster (2019) over a range of max R-squares. The range of max R-squares starts from the R-square of the controlled model rounded up to the next 1/100 to 1. The function supports linear cross sectional (see \emph{lm} objects in R) and panel fixed effect (see \emph{plm} objects in R) models.
#' @return Returns ggplot object. Including deltas* over a range of max R-squares.
#' @references Oster, E. (2019). Unobservable selection and coefficient stability: Theory and evidence. Journal of Business & Economic Statistics, 37, 187-204.
#' @examples
#' # load data, e.g. the in-build mtcars dataset
#' data("mtcars")
#' data_oster <- mtcars
#'
#' # preview of data
#' head(data_oster)
#'
#' # load robomit
#' require(robomit)
#'
#' # estimate and visualize deltas* over a range of max R-squares
#' o_delta_rsq_viz(y = "mpg",               # define the dependent variable name
#'                     x = "wt",            # define the main independent variable name
#'                     con = "hp + qsec",   # other control variables
#'                     beta = 0,            # define beta. This is usually set to 0
#'                     type = "lm",         # define model type
#'                     data = data_oster)   # define dataset
#' @export

o_delta_rsq_viz <- function(y, x, con, id = "none", time = "none", beta = 0, type, data) {


  # rename variables and create, if type is plm, a panel dataset
  ## lm
  if (type == "lm") {data <- data %>% dplyr::rename(y_var = y, x_var = x)}
  ## plm
  if (type == "plm") {data <- data %>% dplyr::rename(y_var = y, x_var = x, id_var = id, time_var = time)
  data_plm <- pdata.frame(data, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)}


  # define formulas of the different models
  model0_formula    <- as.formula("y_var ~ x_var")              # uncontrolled model
  model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) # controlled model
  aux_model_formula <- as.formula(paste("x_var ~",con))         # auxiliary model


  # define model to compute sigma_xx when model type is plm
  if (type == "plm") {sigma_xx_model_formula <- as.formula(paste("x_var ~ factor(id_var)"))}


  # run models
  ## lm
  if (type == "lm") {
    model0    <- lm(model0_formula,    data = data, na.action = na.exclude) # run uncontrolled model
    model1    <- lm(model1_formula,    data = data, na.action = na.exclude) # run controlled model
    aux_model <- lm(aux_model_formula, data = data, na.action = na.exclude) # run auxiliary model
  }
  ## plm
  if (type == "plm") {
    model0    <- plm(model0_formula,    data = data_plm, model = "within", na.action = na.exclude) # run uncontrolled model
    model1    <- plm(model1_formula,    data = data_plm, model = "within", na.action = na.exclude) # run controlled model
    aux_model <- plm(aux_model_formula, data = data_plm, model = "within", na.action = na.exclude) # run auxiliary model
    model_xx  <-  lm(sigma_xx_model_formula, data = data, na.action = na.exclude)                  # run model to obtain singa_xx
  }


  # variables based on model outputs
  if (type == "lm") {b0 = as.numeric(tidy(model0)[2,2])} else if (type == "plm") {b0 = as.numeric(tidy(model0)[1,2])}           # beta uncontrolled model
  if (type == "lm") {b1 = as.numeric(tidy(model1)[2,2])} else if (type == "plm") {b1 = as.numeric(tidy(model1)[1,2])}           # beta controlled model
  if (type == "lm") {R20 = summary(model0)$r.squared} else if (type == "plm") {R20 = as.numeric(summary(model0)$r.squared[1])}  # R-square uncontrolled model
  if (type == "lm") {R21 = summary(model1)$r.squared} else if (type == "plm") {R21 = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model
  sigma_yy = var(data$y_var, na.rm = T)                                                                                         # variance of dependent variable
  if (type == "lm") {sigma_xx = var(data$x_var, na.rm = T)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)}       # variance of independent variable
  t_x  = var(aux_model$residuals)


  # define max R-square range
  r_start <- ceiling(R21/0.01)*0.01 # sthart max R-square at the next highest 1/100
  n_runs <- length(seq(r_start:1,by = 0.01))
  delta_over_rsq_results <- tibble(
    x = 1:n_runs,
    Rmax = 0,
    Delta = 0) # create tibble object to store results


  # run loop over different max R-squares
  for (i in 1:n_runs) {
    R2max = seq(r_start:1,by = 0.01)[i] # current max R-square


    # create some additional variables
    bt_m_b = b1 - beta
    rt_m_ro_t_syy = (R21-R20) * sigma_yy
    b0_m_b1 = b0 - b1
    rm_m_rt_t_syy = (R2max - R21) * sigma_yy


    # compute numerator, in steps
    num1 = bt_m_b * rt_m_ro_t_syy * t_x
    num2 = bt_m_b * sigma_xx * t_x * b0_m_b1^2
    num3 = 2 * bt_m_b^2 * (t_x * b0_m_b1 * sigma_xx)
    num4 = bt_m_b^3 * (t_x * sigma_xx - t_x^2)
    num  = num1 + num2 + num3 + num4


    # compute denominator, in steps
    den1 = rm_m_rt_t_syy * b0_m_b1 * sigma_xx
    den2 = bt_m_b * rm_m_rt_t_syy * (sigma_xx - t_x)
    den3 = bt_m_b^2 * (t_x * b0_m_b1 * sigma_xx)
    den4 = bt_m_b^3 * (t_x * sigma_xx - t_x^2)
    den  = den1 + den2 + den3 + den4


    # finally compute delta_star
    delta_star = num/den


    # store results
    delta_over_rsq_results[i,2] <- R2max
    delta_over_rsq_results[i,3] <- round(delta_star,6)
  }


  # plot figure
  theme_set(theme_bw()) # set main design
  result_plot <- ggplot(data=delta_over_rsq_results, aes(x=Rmax, y=Delta)) +
    geom_line(size = 1.3)+
    scale_y_continuous(name = expression(delta^"*"))+
    scale_x_continuous(name = expression("max R"^2))+
    theme(axis.title = element_text( size=15),
          axis.text  = element_text( size=13))


  # warning if max R-square is smaller than the R-square of the control model
  if (R21 > R2max)
    warning("The max R-square value is smaller than the R-square of the controlled model")


  # define return object
  return(result_plot)
}



##################--------------------------------------------------------------
#------------------------------- 1.4 o_delta_boot - function 4
##################--------------------------------------------------------------
#' @title Bootstrapped deltas*
#'
#' @description Estimates bootstrapped deltas*, i.e. the degree of selection on unobservables relative to observables that would be necessary to explain away the result, following Oster (2019).
#' @usage o_delta_boot(y, x, con, id = "none", time = "none", beta = 0, R2max, sim, obs,
#' rep, type, useed = NA, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent variable of interest (treatment variable; as string).
#' @param con Name of the other control variables. Provided as string in the format: "w + z +...".
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect models.
#' @param time Name of the time variable (e.g. year or month; as string). Only applicable for fixed effect models.
#' @param beta Beta for which delta* should be estimated (default is beta = 0).
#' @param R2max Max R-square for which beta* should be estimated.
#' @param sim Number of simulations.
#' @param obs Number of draws per simulation.
#' @param rep Bootstrapping either with (= TRUE) or without (= FALSE) replacement.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param useed Seed number defined by user.
#' @param data Data.
#' @details Estimates bootstrapped deltas*, i.e. the degree of selection on unobservables relative to observables that would be necessary to explain away the result, following Oster (2019). Bootstrapping can either be done with or without replacement. The function supports linear cross sectional (see \emph{lm} objects in R) and panel fixed effect (see \emph{plm} objects in R) models.
#' @return Returns tibble object. Including bootstrapped deltas*.
#' @references Oster, E. (2019). Unobservable selection and coefficient stability: Theory and evidence. Journal of Business & Economic Statistics, 37, 187-204.
#' @examples
#' # load data, e.g. the in-build mtcars dataset
#' data("mtcars")
#' data_oster <- mtcars
#'
#' # preview of data
#' head(data_oster)
#'
#' # load robomit
#' require(robomit)
#'
#' # estimate bootstrapped deltas*
#' o_delta_boot(y = "mpg",              # define the dependent variable name
#'                  x = "wt",           # define the main independent variable name
#'                  con = "hp + qsec",  # other control variables
#'                  beta = 0,           # define beta. This is usually set to 0
#'                  R2max = 0.9,        # define the max R-square.
#'                  sim = 100,          # define number of simulations
#'                  obs = 30,           # define number of drawn observations per simulation
#'                  rep = FALSE,        # define if bootstrapping is with or without replacement
#'                  type = "lm",        # define model type
#'                  useed = 123,        # define seed
#'                  data = data_oster)  # define dataset
#' @export

o_delta_boot <- function(y, x, con, id = "none", time = "none", beta = 0, R2max, sim, obs, rep, type, useed = NA, data) {


  # rename variables and create
  ## lm
  if (type == "lm") {data <- data %>% dplyr::rename(y_var = y, x_var = x)}
  ## plm
  if (type == "plm") {data <- data %>% dplyr::rename(y_var = y, x_var = x, id_var = id, time_var = time)}


  # define formulas of the different models
  model0_formula    <- as.formula("y_var ~ x_var")              # uncontrolled model
  model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) # controlled model
  aux_model_formula <- as.formula(paste("x_var ~",con))         # auxiliary model


  # define model to compute sigma_xx when model type is plm
  if (type == "plm") {sigma_xx_model_formula <- as.formula(paste("x_var ~",con,"+ factor(id_var)"))}


  # create tibble object to store results
  simulation_results <- tibble(
    "x" = 1:sim,
    "Delta*" = 0)


  for (s in 1:sim) {


    # fun models
    ## lm
    if (type == "lm") {


      if (!is.na(useed)) {set.seed(useed+s) } # set seed (defined by user) that it varies with runs
      if (rep) {
        data_current <- sample_n(data, size = obs, replace = TRUE) # with replacement
      } else {
        data_current <- sample_n(data, size = obs, replace = FALSE) # without replacement
      }
      model0    <- lm(model0_formula,    data = data_current, na.action = na.exclude) # run uncontrolled model
      model1    <- lm(model1_formula,    data = data_current, na.action = na.exclude) # run controlled model
      aux_model <- lm(aux_model_formula, data = data_current, na.action = na.exclude) # run auxiliary model
    }
    ## plm
    if (type == "plm") {

      if (!is.na(useed)) {set.seed(useed+s) } # set seed (defined by user) that it varies with runs
      if (rep) {
        data_current     <- sample_n(data, size = obs, replace = TRUE) # with replacement
        data_current_plm <- pdata.frame(data_current, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)
      } else {
        data_current     <- sample_n(data, size = obs, replace = FALSE) # without replacement
        data_current_plm <- pdata.frame(data_current, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)
      }
      model0    <- plm(model0_formula,    data = data_current_plm, model = "within", na.action = na.exclude) # run uncontrolled model
      model1    <- plm(model1_formula,    data = data_current_plm, model = "within", na.action = na.exclude) # run controlled model
      aux_model <- plm(aux_model_formula, data = data_current_plm, model = "within", na.action = na.exclude) # run auxiliary model
      model_xx  <-  lm(sigma_xx_model_formula, data = data_current, na.action = na.exclude)                  # run model to obtain singa_xx
    }


    # variables based on model outputs
    if (type == "lm") {b0 = as.numeric(tidy(model0)[2,2])} else if (type == "plm") {b0 = as.numeric(tidy(model0)[1,2])}           # beta uncontrolled model
    if (type == "lm") {b1 = as.numeric(tidy(model1)[2,2])} else if (type == "plm") {b1 = as.numeric(tidy(model1)[1,2])}           # beta controlled model
    if (type == "lm") {R20 = summary(model0)$r.squared} else if (type == "plm") {R20 = as.numeric(summary(model0)$r.squared[1])}  # R-square uncontrolled model
    if (type == "lm") {R21 = summary(model1)$r.squared} else if (type == "plm") {R21 = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model
    sigma_yy = var(data_current$y_var, na.rm = T)                                                                                          # variance of dependent variable
    if (type == "lm") {sigma_xx = var(data_current$x_var, na.rm = T)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)}       # variance of independent variable
    t_x  = var(aux_model$residuals)


    # create some additional variables
    rt_m_ro_t_syy = (R21-R20) * sigma_yy
    b0_m_b1 = b0 - b1
    rm_m_rt_t_syy = (R2max - R21) * sigma_yy


    # create some additional variables
    bt_m_b = b1 - beta
    rt_m_ro_t_syy = (R21-R20) * sigma_yy
    b0_m_b1 = b0 - b1
    rm_m_rt_t_syy = (R2max - R21) * sigma_yy


    # compute numerator, in steps
    num1 = bt_m_b * rt_m_ro_t_syy * t_x
    num2 = bt_m_b * sigma_xx * t_x * b0_m_b1^2
    num3 = 2 * bt_m_b^2 * (t_x * b0_m_b1 * sigma_xx)
    num4 = bt_m_b^3 * (t_x * sigma_xx - t_x^2)
    num  = num1 + num2 + num3 + num4


    # compute denominator, in steps
    den1 = rm_m_rt_t_syy * b0_m_b1 * sigma_xx
    den2 = bt_m_b * rm_m_rt_t_syy * (sigma_xx - t_x)
    den3 = bt_m_b^2 * (t_x * b0_m_b1 * sigma_xx)
    den4 = bt_m_b^3 * (t_x * sigma_xx - t_x^2)
    den  = den1 + den2 + den3 + den4


    # finally compute delta_star
    delta_star = num/den


    # store results
    simulation_results[s,2] <- round(delta_star,6)
  }


  # warning if max R-square is smaller than the R-square of the control model
  if (R21 > R2max)
    warning("The max R-square value is smaller than the R-square of the controlled model")


  # warning if the bootstrapped observations are larger than the number of observations of the original dataset
  if (obs >= length(data$y_var))
    warning("Number of bootstrapped observation is larger than/equal to number of observation in the provided dataset")


  #define return object
  return(simulation_results)
}



##################--------------------------------------------------------------
#------------------------------- 1.5 o_delta_boot_viz - function 5
##################--------------------------------------------------------------
#' @title  Visualization of bootstrapped deltas*
#'
#' @description Estimates and visualizes bootstrapped deltas*, i.e. the degree of selection on unobservables relative to observables that would be necessary to explain away the result, following Oster (2019).
#' @usage o_delta_boot_viz(y, x, con, id = "none", time = "none", beta = 0, R2max, sim, obs, rep,
#' CI, type, norm = TRUE, bin, col = c("#08306b","#4292c6","#c6dbef"),
#' nL = FALSE, mL = TRUE, useed = NA, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent variable of interest (treatment variable; as string).
#' @param con Name of the other control variables. Provided as string in the format: "w + z +...".
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect models.
#' @param time Name of the time variable (e.g. year or month; as string). Only applicable for fixed effect models.
#' @param beta Beta for which delta* should be estimated (default is beta = 0).
#' @param R2max Max R-square for which beta* should be estimated.
#' @param sim Number of simulations.
#' @param obs Number of draws per simulation.
#' @param rep Bootstrapping either with (= TRUE) or without (= FALSE) replacement
#' @param CI  Confidence intervals, indicated as vector. Can be and/or 90,95,99.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param norm Option to include a normal distribution in the plot (default is norm = TURE).
#' @param bin Number of bins used for the histogram.
#' @param col Colors used to indicate different confidence interval levels (indicated as vector). Needs to be the same length as the variable CI. The default is a blue color range.
#' @param nL Option to include a red vertical line at 0 (default is nL = TRUE).
#' @param mL Option to include a vertical line at beta* mean (default is mL = TRUE).
#' @param useed Seed number defined by user.
#' @param data Data.
#' @details Estimates and visualizes bootstrapped deltas*, i.e. the degree of selection on unobservables relative to observables that would be necessary to explain away the result, following Oster (2019). Bootstrapping can either be done with or without replacement. The function supports linear cross sectional (see \emph{lm} objects in R) and panel fixed effect (see \emph{plm} objects in R) models.
#' @return Returns ggplot object. Including bootstrapped deltas*.
#' @references Oster, E. (2019). Unobservable selection and coefficient stability: Theory and evidence. Journal of Business & Economic Statistics, 37, 187-204.
#' @examples
#' # load data, e.g. the in-build mtcars dataset
#' data("mtcars")
#' data_oster <- mtcars
#'
#' # preview of data
#' head(data_oster)
#'
#' # load robomit
#' require(robomit)
#'
#' # estimate and visualize bootstrapped deltas*
#' o_delta_boot_viz(y = "mpg",               # define the dependent variable name
#'                      x = "wt",            # define the main independent variable name
#'                      con = "hp + qsec",   # other control variables
#'                      beta = 0,            # define beta. This is usually set to 0
#'                      R2max = 0.9,         # define the max R-square.
#'                      sim = 100,           # define number of simulations
#'                      obs = 30,            # define number of drawn observations per simulation
#'                      rep = FALSE,         # define if bootstrapping is with or without replacement
#'                      CI = c(90,95,99),    # define confidence intervals.
#'                      type = "lm",         # define model type
#'                      norm = TRUE,         # include normal distribution
#'                      bin = 200,           # set number of bins
#'                      useed = 123,         # define seed
#'                      data = data_oster)   # define dataset
#' @export


o_delta_boot_viz <- function(y, x, con, id = "none", time = "none", beta = 0, R2max, sim, obs, rep, CI, type, norm = TRUE, bin, col = c("#08306b","#4292c6","#c6dbef"), nL = FALSE, mL = TRUE, useed = NA, data) {


  # rename variables and create
  ## lm
  if (type == "lm") {data <- data %>% dplyr::rename(y_var = y, x_var = x)}
  ## plm
  if (type == "plm") {data <- data %>% dplyr::rename(y_var = y, x_var = x, id_var = id, time_var = time)}


  # define formulas of the different models
  model0_formula    <- as.formula("y_var ~ x_var")              # uncontrolled model
  model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) # controlled model
  aux_model_formula <- as.formula(paste("x_var ~",con))         # auxiliary model


  # define model to compute sigma_xx when model type is plm
  if (type == "plm") {sigma_xx_model_formula <- as.formula(paste("x_var ~",con,"+ factor(id_var)"))}


  # create tibble object to store results
  simulation_results <- tibble(
    x = 1:sim,
    Delta = 0)


  for (s in 1:sim) {
    # fun models
    ## lm
    if (type == "lm") {


      if (!is.na(useed)) {set.seed(useed+s) } # set seed (defined by user) that it varies with runs
      if (rep) {
        data_current <- sample_n(data, size = obs, replace = TRUE) # with replacement
      } else {
        data_current <- sample_n(data, size = obs, replace = FALSE) # without replacement
      }
      model0    <- lm(model0_formula,    data = data_current, na.action = na.exclude) # run uncontrolled model
      model1    <- lm(model1_formula,    data = data_current, na.action = na.exclude) # run controlled model
      aux_model <- lm(aux_model_formula, data = data_current, na.action = na.exclude) # run auxiliary model
    }
    ## plm
    if (type == "plm") {

      if (!is.na(useed)) {set.seed(useed+s) } # set seed (defined by user) that it varies with runs
      if (rep) {
        data_current     <- sample_n(data, size = obs, replace = TRUE) # with replacement
        data_current_plm <- pdata.frame(data_current, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)
      } else {
        data_current     <- sample_n(data, size = obs, replace = FALSE) # without replacement
        data_current_plm <- pdata.frame(data_current, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)
      }
      model0    <- plm(model0_formula,    data = data_current_plm, model = "within", na.action = na.exclude) # run uncontrolled model
      model1    <- plm(model1_formula,    data = data_current_plm, model = "within", na.action = na.exclude) # run controlled model
      aux_model <- plm(aux_model_formula, data = data_current_plm, model = "within", na.action = na.exclude) # run auxiliary model
      model_xx  <-  lm(sigma_xx_model_formula, data = data_current, na.action = na.exclude)                  # run model to obtain singa_xx
    }


    # variables based on model outputs
    if (type == "lm") {b0 = as.numeric(tidy(model0)[2,2])} else if (type == "plm") {b0 = as.numeric(tidy(model0)[1,2])}           # beta uncontrolled model
    if (type == "lm") {b1 = as.numeric(tidy(model1)[2,2])} else if (type == "plm") {b1 = as.numeric(tidy(model1)[1,2])}           # beta controlled model
    if (type == "lm") {R20 = summary(model0)$r.squared} else if (type == "plm") {R20 = as.numeric(summary(model0)$r.squared[1])}  # R-square uncontrolled model
    if (type == "lm") {R21 = summary(model1)$r.squared} else if (type == "plm") {R21 = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model
    sigma_yy = var(data_current$y_var, na.rm = T)                                                                                          # variance of dependent variable
    if (type == "lm") {sigma_xx = var(data_current$x_var, na.rm = T)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)}       # variance of independent variable
    t_x  = var(aux_model$residuals)


    # create some additional variables
    rt_m_ro_t_syy = (R21-R20) * sigma_yy
    b0_m_b1 = b0 - b1
    rm_m_rt_t_syy = (R2max - R21) * sigma_yy


    # create some additional variables
    bt_m_b = b1 - beta
    rt_m_ro_t_syy = (R21-R20) * sigma_yy
    b0_m_b1 = b0 - b1
    rm_m_rt_t_syy = (R2max - R21) * sigma_yy


    # compute numerator, in steps
    num1 = bt_m_b * rt_m_ro_t_syy * t_x
    num2 = bt_m_b * sigma_xx * t_x * b0_m_b1^2
    num3 = 2 * bt_m_b^2 * (t_x * b0_m_b1 * sigma_xx)
    num4 = bt_m_b^3 * (t_x * sigma_xx - t_x^2)
    num  = num1 + num2 + num3 + num4


    # compute demoninator, in steps
    den1 = rm_m_rt_t_syy * b0_m_b1 * sigma_xx
    den2 = bt_m_b * rm_m_rt_t_syy * (sigma_xx - t_x)
    den3 = bt_m_b^2 * (t_x * b0_m_b1 * sigma_xx)
    den4 = bt_m_b^3 * (t_x * sigma_xx - t_x^2)
    den  = den1 + den2 + den3 + den4


    # finally compute delta_star
    delta_star = num/den


    # store results
    simulation_results[s,2] <- round(delta_star,6)
  }


  # compute mean, max absolute distance to mean, and confidence intervals
  ## mean
  meanDelta <- mean(simulation_results$Delta, na.rm = F)
  # sd
  sdDelta <- sd(simulation_results$Delta, na.rm = F)
  ## confidence interval
  ### 90%
  CI90_a <- meanDelta - 1.645 * sdDelta
  CI90_b <- meanDelta + 1.645 * sdDelta
  ### 95%
  CI95_a <- meanDelta - 1.96 * sdDelta
  CI95_b <- meanDelta + 1.96 * sd(simulation_results$Delta, na.rm = F)
  ### 99%
  CI99_a <- meanDelta - 2.58 * sdDelta
  CI99_b <- meanDelta + 2.58 * sdDelta


  # define colors of confidence intervals
  color1 = col[1]
  color2 = col[2]
  color3 = col[3]
  ## 90%
  if (is.element(90, CI)) {
    color90 <- color1}
  ## 95%
  if (is.element(95, CI)) {
    color95 <- ifelse(length(CI) == 3,color2,
                      ifelse(length(CI) == 2 & is.element(90, CI),color2,
                             ifelse(length(CI) == 2 & is.element(99, CI),color1,color1)))}
  ## 99%
  if (is.element(99, CI)) {
    color99 <- ifelse(length(CI) == 3,color3,
                      ifelse(length(CI) == 2,color2,color1))}


  # plot figure
  theme_set(theme_bw()) # set main design
  delta_plot <- ggplot(simulation_results, aes(x = Delta)) +
  {if(is.element(99, CI))annotate("rect", ymin = 0, ymax = +Inf, xmin = CI99_a, xmax = CI99_b, fill = color99, alpha = 0.6)}+
  {if(is.element(95, CI))annotate("rect", ymin = 0, ymax = +Inf, xmin = CI95_a, xmax = CI95_b, fill = color95, alpha = 0.6)}+
  {if(is.element(90, CI))annotate("rect", ymin = 0, ymax = +Inf, xmin = CI90_a, xmax = CI90_b, fill = color90, alpha = 0.6)}+
    geom_histogram(aes(y = ..density..), alpha = 0.8, colour = "#737373", fill = "#737373", size = 0.1, bins = bin)+
    {if(norm)stat_function(fun = dnorm, args = list(mean = meanDelta, sd = sdDelta), size = 1.3)}+
    scale_y_continuous(name = "Density",expand = expansion(mult = c(0, .1)))+
    scale_x_continuous(name = expression(delta^"*"))+
    expand_limits(x = 0, y = 0)+
    theme(axis.title = element_text( size=15),
          axis.text  = element_text( size=13))+
          {if(mL)geom_vline(xintercept = meanDelta, size = 0.8, color = "#000000", alpha = 0.8)}+
          {if(nL)geom_vline(xintercept = 0, size = 0.8, color = "#e41a1c", alpha = 0.8)}+
    geom_hline(yintercept = 0)


  # warning if max R-square is smaller than the R-square of the control model
  if (R21 > R2max)
    warning("The max R-square value is smaller than the R-square of the controlled model")


  # warning if the bootstrapped observations are larger than the number of observations of the original dataset
  if (obs >= length(data$y_var))
    warning("Number of bootstrapped observation is larger than/equal to number of observation in the provided dataset")


  # define return object
  return(delta_plot)
}


##################--------------------------------------------------------------
#------------------------------- 1.6 o_delta_boot_inf - function 6
##################--------------------------------------------------------------
#' @title Bootstrapped mean delta* and confidence intervals
#'
#' @description Estimates mean and provides confidence intervals of bootstrapped deltas*, i.e. the degree of selection on unobservables relative to observables that would be necessary to explain away the result, following Oster (2019).
#' @usage o_delta_boot_inf(y, x, con, id = "none", time = "none", beta = 0, R2max, sim, obs, rep,
#' CI, type, useed = NA, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent variable of interest (treatment variable; as string).
#' @param con Name of the other control variables. Provided as string in the format: "w + z +...".
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect models.
#' @param time Name of the time variable (e.g. year or month; as string). Only applicable for fixed effect models.
#' @param beta Beta for which delta* should be estimated (default is beta = 0)..
#' @param R2max Max R-square for which beta* should be estimated.
#' @param sim Number of simulations.
#' @param obs Number of draws per simulation.
#' @param rep Bootstrapping either with (= TRUE) or without (= FALSE) replacement
#' @param CI  Confidence intervals, indicated as vector. Can be and/or 90,95,99.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param useed Seed number defined by user.
#' @param data Data.
#' @details Estimates mean and provides confidence intervals of bootstrapped deltas*, i.e. the degree of selection on unobservables relative to observables that would be necessary to explain away the result, following Oster (2019). Bootstrapping can either be done with or without replacement. The function supports linear cross sectional (see \emph{lm} objects in R) and panel fixed effect (see \emph{plm} objects in R) models.
#' @return Returns tibble object. Including bootstrapped deltas* and confidence intervals.
#' @references Oster, E. (2019). Unobservable selection and coefficient stability: Theory and evidence. Journal of Business & Economic Statistics, 37, 187-204.
#' @examples
#' # load data, e.g. the in-build mtcars dataset
#' data("mtcars")
#' data_oster <- mtcars
#'
#' # preview of data
#' head(data_oster)
#'
#' # load robomit
#' require(robomit)
#'
#' # estimate and visualize bootstrapped deltas*
#' o_delta_boot_inf(y = "mpg",                # define the dependent variable name
#'                      x = "wt",             # define the main independent variable name
#'                      con = "hp + qsec",    # other control variables
#'                      beta = 0,             # define beta. This is usually set to 0
#'                      R2max = 0.9,          # define the max R-square.
#'                      sim = 100,            # define number of simulations
#'                      obs = 30,             # define number of drawn observations per simulation
#'                      rep = FALSE,          # define if bootstrapping is with or without replacement
#'                      CI = c(90,95,99),     # define confidence intervals.
#'                      type = "lm",          # define model type
#'                      useed = 123,          # define seed
#'                      data = data_oster)    # define dataset
#' @export


o_delta_boot_inf <- function(y, x, con, id = "none", time = "none", beta = 0, R2max, sim, obs, rep, CI, type, useed = NA, data) {


  # rename variables and create
  ## lm
  if (type == "lm") {data <- data %>% dplyr::rename(y_var = y, x_var = x)}
  ## plm
  if (type == "plm") {data <- data %>% dplyr::rename(y_var = y, x_var = x, id_var = id, time_var = time)}


  # define formulas of the different models
  model0_formula    <- as.formula("y_var ~ x_var")              # uncontrolled model
  model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) # controlled model
  aux_model_formula <- as.formula(paste("x_var ~",con))         # auxiliary model


  # define model to compute sigma_xx when model type is plm
  if (type == "plm") {sigma_xx_model_formula <- as.formula(paste("x_var ~",con,"+ factor(id_var)"))}


  # create tibble object to store results
  simulation_results <- tibble(
    x = 1:sim,
    Delta = 0)


  for (s in 1:sim) {
    # fun models
    ## lm
    if (type == "lm") {


      if (!is.na(useed)) {set.seed(useed+s) } # set seed (defined by user) that it varies with runs
      if (rep) {
        data_current <- sample_n(data, size = obs, replace = TRUE) # with replacement
      } else {
        data_current <- sample_n(data, size = obs, replace = FALSE) # without replacement
      }
      model0    <- lm(model0_formula,    data = data_current, na.action = na.exclude) # run uncontrolled model
      model1    <- lm(model1_formula,    data = data_current, na.action = na.exclude) # run controlled model
      aux_model <- lm(aux_model_formula, data = data_current, na.action = na.exclude) # run auxiliary model
    }
    ## plm
    if (type == "plm") {

      if (!is.na(useed)) {set.seed(useed+s) } # set seed (defined by user) that it varies with runs
      if (rep) {
        data_current     <- sample_n(data, size = obs, replace = TRUE) # with replacement
        data_current_plm <- pdata.frame(data_current, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)
      } else {
        data_current     <- sample_n(data, size = obs, replace = FALSE) # without replacement
        data_current_plm <- pdata.frame(data_current, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)
      }
      model0    <- plm(model0_formula,    data = data_current_plm, model = "within", na.action = na.exclude) # run uncontrolled model
      model1    <- plm(model1_formula,    data = data_current_plm, model = "within", na.action = na.exclude) # run controlled model
      aux_model <- plm(aux_model_formula, data = data_current_plm, model = "within", na.action = na.exclude) # run auxiliary model
      model_xx  <-  lm(sigma_xx_model_formula, data = data_current, na.action = na.exclude)                  # run model to obtain singa_xx
    }


    # variables based on model outputs
    if (type == "lm") {b0 = as.numeric(tidy(model0)[2,2])} else if (type == "plm") {b0 = as.numeric(tidy(model0)[1,2])}           # beta uncontrolled model
    if (type == "lm") {b1 = as.numeric(tidy(model1)[2,2])} else if (type == "plm") {b1 = as.numeric(tidy(model1)[1,2])}           # beta controlled model
    if (type == "lm") {R20 = summary(model0)$r.squared} else if (type == "plm") {R20 = as.numeric(summary(model0)$r.squared[1])}  # R-square uncontrolled model
    if (type == "lm") {R21 = summary(model1)$r.squared} else if (type == "plm") {R21 = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model
    sigma_yy = var(data_current$y_var, na.rm = T)                                                                                          # variance of dependent variable
    if (type == "lm") {sigma_xx = var(data_current$x_var, na.rm = T)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)}       # variance of independent variable
    t_x  = var(aux_model$residuals)


    # create some additional variables
    rt_m_ro_t_syy = (R21-R20) * sigma_yy
    b0_m_b1 = b0 - b1
    rm_m_rt_t_syy = (R2max - R21) * sigma_yy


    # create some additional variables
    bt_m_b = b1 - beta
    rt_m_ro_t_syy = (R21-R20) * sigma_yy
    b0_m_b1 = b0 - b1
    rm_m_rt_t_syy = (R2max - R21) * sigma_yy


    # compute numerator, in steps
    num1 = bt_m_b * rt_m_ro_t_syy * t_x
    num2 = bt_m_b * sigma_xx * t_x * b0_m_b1^2
    num3 = 2 * bt_m_b^2 * (t_x * b0_m_b1 * sigma_xx)
    num4 = bt_m_b^3 * (t_x * sigma_xx - t_x^2)
    num  = num1 + num2 + num3 + num4


    # compute denominator, in steps
    den1 = rm_m_rt_t_syy * b0_m_b1 * sigma_xx
    den2 = bt_m_b * rm_m_rt_t_syy * (sigma_xx - t_x)
    den3 = bt_m_b^2 * (t_x * b0_m_b1 * sigma_xx)
    den4 = bt_m_b^3 * (t_x * sigma_xx - t_x^2)
    den  = den1 + den2 + den3 + den4


    # finally compute delta_star
    delta_star = num/den


    # store results
    simulation_results[s,2] <- round(delta_star,6)
  }


  # compute mean, max absolute distance to mean, and confidence intervals
  ## mean
  meanDelta <- mean(simulation_results$Delta, na.rm = F)
  # sd
  sdDelta <- sd(simulation_results$Delta, na.rm = F)
  ## max absolute distance
  abDist <- if (abs(min(simulation_results$Delta, na.rm = F) - meanDelta) > abs(max(simulation_results$Delta, na.rm = F) - meanDelta)) {
    abs(min(simulation_results$Delta, na.rm = F) - meanDelta) } else {
      abs(max(simulation_results$Delta, na.rm = F) - meanDelta) }
  ## confidence interval
  ### 90%
  CI90_a <- meanDelta - 1.645 * sdDelta
  CI90_b <- meanDelta + 1.645 * sdDelta
  ### 95%
  CI95_a <- meanDelta - 1.96 * sdDelta
  CI95_b <- meanDelta + 1.96 * sd(simulation_results$Delta, na.rm = F)
  ### 99%
  CI99_a <- meanDelta - 2.58 * sdDelta
  CI99_b <- meanDelta + 2.58 * sdDelta


  # create tibble object of results
  result_Delta_inf_aux1 <- tribble(
    ~Name,                                 ~Value,
    "Delta* (mean)",                         round(meanDelta,6))
  result_Delta_inf_aux2 <- tribble(
    ~Name,                                 ~Value,
    if(is.element(90, CI)) {"CI_90_low"},  round(CI90_a,6),
    if(is.element(90, CI)) {"CI_90_high"}, round(CI90_b,6))
  result_Delta_inf_aux3 <- tribble(
    ~Name,                                 ~Value,
    if(is.element(95, CI)) {"CI_95_low"},  round(CI95_a,6),
    if(is.element(95, CI)) {"CI_95_high"}, round(CI95_b,6))
  result_Delta_inf_aux4 <- tribble(
    ~Name,                                 ~Value,
    if(is.element(99, CI)) {"CI_99_low"},  round(CI99_a,6),
    if(is.element(99, CI)) {"CI_99_high"}, round(CI99_b,6))
  result_Delta_inf_aux5 <- tribble(
    ~Name,                                 ~Value,
    "Simulations",                         sim,
    "Observations",                        obs,
    "Max R-square",                        R2max,
    "Beta (defined)",                      beta,
    paste("Model:",type),                  NA,
    paste("Replacement:",rep),     NA)
  if (is.element(90, CI)) {result_Delta_inf_aux1 <- result_Delta_inf_aux1 %>% add_row(result_Delta_inf_aux2)}
  if (is.element(95, CI)) {result_Delta_inf_aux1 <- result_Delta_inf_aux1 %>% add_row(result_Delta_inf_aux3)}
  if (is.element(99, CI)) {result_Delta_inf_aux1 <- result_Delta_inf_aux1 %>% add_row(result_Delta_inf_aux4)}
  result_Delta_inf <- result_Delta_inf_aux1 %>% add_row(result_Delta_inf_aux5)


  # warning if max R-square is smaller than the R-square of the control model
  if (R21 > R2max)
    warning("The max R-square value is smaller than the R-square of the controlled model")


  # warning if the bootstrapped observations are larger than the number of observations of the original dataset
  if (obs >= length(data$y_var))
    warning("Number of bootstrapped observation is larger than/equal to number of observation in the provided dataset")


  # define return object
  return(result_Delta_inf)
}





##################--------------------------------------------------------------
#------------------------------- 1.7 o_beta - function 7
##################--------------------------------------------------------------
#' @title Beta*
#'
#' @description Estimates beta*, i.e. the estimated bias-adjusted treatment effects, following Oster (2019).
#' @usage o_beta(y, x, con, id = "none", time = "none", delta = 1, R2max, type, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent variable of interest (treatment variable; as string).
#' @param con Name of the other control variables. Provided as string in the format: "w + z +...".
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect models.
#' @param time Name of the time variable (e.g. year or month; as string). Only applicable for fixed effect models.
#' @param delta Delta for which beta* should be estimated (default is delta = 1).
#' @param R2max Max R-square for which beta* should be estimated.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param data Data.
#' @details Estimates beta*, i.e. the estimated bias-adjusted treatment effects, following Oster (2019). The function supports linear cross sectional (see \emph{lm} objects in R) and panel fixed effect (see \emph{plm} objects in R) models.
#' @return Returns tibble object. Including beta* and various other information.
#' @references Oster, E. (2019). Unobservable selection and coefficient stability: Theory and evidence. Journal of Business & Economic Statistics, 37, 187-204.
#' @examples
#' # load data, e.g. the in-build mtcars dataset
#' data("mtcars")
#' data_oster <- mtcars
#'
#' # preview of data
#' head(data_oster)
#'
#' # load robomit
#' require(robomit)
#'
#' # estimate beta*
#' o_beta(y = "mpg",               # define the dependent variable name
#'            x = "wt",            # define the main independent variable name
#'            con = "hp + qsec",   # other control variables
#'            delta = 1,           # define delta. This is usually set to 1
#'            R2max = 0.9,         # define the max R-square.
#'            type = "lm",         # define model type
#'            data = data_oster)   # define dataset
#' @export


o_beta <- function(y, x, con, id = "none", time = "none", delta = 1 ,R2max, type, data) {


  # rename variables and create, if type is plm, a panel dataset
  ## lm
  if (type == "lm") {data <- data %>% dplyr::rename(y_var = y, x_var = x)}
  ## plm
  if (type == "plm") {data <- data %>% dplyr::rename(y_var = y, x_var = x, id_var = id, time_var = time)
  data_plm <- pdata.frame(data, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)}


  # define formulas of the different models
  model0_formula    <- as.formula("y_var ~ x_var")              # uncontrolled model
  model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) # controlled model
  aux_model_formula <- as.formula(paste("x_var ~",con))         # auxiliary model


  # define model to compute sigma_xx when model type is plm
  if (type == "plm") {sigma_xx_model_formula <- as.formula(paste("x_var ~ factor(id_var)"))}


  # run models
  ## lm
  if (type == "lm") {
    model0    <- lm(model0_formula,    data = data, na.action = na.exclude) # run uncontrolled model
    model1    <- lm(model1_formula,    data = data, na.action = na.exclude) # run controlled model
    aux_model <- lm(aux_model_formula, data = data, na.action = na.exclude) # run auxiliary model
  }
  ## plm
  if (type == "plm") {
    model0    <- plm(model0_formula,    data = data_plm, model = "within", na.action = na.exclude) # run uncontrolled model
    model1    <- plm(model1_formula,    data = data_plm, model = "within", na.action = na.exclude) # run controlled model
    aux_model <- plm(aux_model_formula, data = data_plm, model = "within", na.action = na.exclude) # run auxiliary model
    model_xx  <-  lm(sigma_xx_model_formula, data = data, na.action = na.exclude)                  # run model to obtain singa_xx
  }


  # variables based on model outputs
  if (type == "lm") {b0 = as.numeric(tidy(model0)[2,2])} else if (type == "plm") {b0 = as.numeric(tidy(model0)[1,2])}           # beta uncontrolled model
  if (type == "lm") {b1 = as.numeric(tidy(model1)[2,2])} else if (type == "plm") {b1 = as.numeric(tidy(model1)[1,2])}           # beta controlled model
  if (type == "lm") {R20 = summary(model0)$r.squared} else if (type == "plm") {R20 = as.numeric(summary(model0)$r.squared[1])}  # R-square uncontrolled model
  if (type == "lm") {R21 = summary(model1)$r.squared} else if (type == "plm") {R21 = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model
  sigma_yy = var(data$y_var, na.rm = T)                                                                                         # variance of dependent variable
  if (type == "lm") {sigma_xx = var(data$x_var, na.rm = T)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)}       # variance of independent variable
  t_x  = var(aux_model$residuals)                                                                                               # variance of residuals of the auxiliary model


  # create some additional variables
  rt_m_ro_t_syy = (R21-R20) * sigma_yy
  b0_m_b1 = b0 - b1
  rm_m_rt_t_syy = (R2max - R21) * sigma_yy


  if (delta == 1) { # quadratic function: v = (- cap_theta +/- sqrt ( cap_theta + d1_1))/(d1_2) (see psacalc command in Stata)


    # compute different variables
    cap_theta = rm_m_rt_t_syy*(sigma_xx-t_x)-rt_m_ro_t_syy*t_x-sigma_xx*t_x*(b0_m_b1^2)
    d1_1 = 4*rm_m_rt_t_syy*(b0_m_b1^2)*(sigma_xx^2)*t_x
    d1_2 = -2*t_x*b0_m_b1*sigma_xx

    sol1 = (-1*cap_theta-sqrt((cap_theta)^2+d1_1))/(d1_2)
    sol2 = (-1*cap_theta+sqrt((cap_theta)^2+d1_1))/(d1_2)

    beta1 = b1 - sol1
    beta2 = b1 - sol2

    if ( (beta1-b1)^2 < (beta2-b1)^2) {
      betax = beta1
      altsol1 = beta2
    } else {
      betax = beta2
      altsol1 = beta1
    }

    # change to alternative solution if bias of b0 vs b1 is of different sign that bias of b1 - beta
    if ( sign(betax-b1)!=sign(b1-b0) ) {
      solc = betax
      betax = altsol1
      altsol1 = solc
    }


    # compute squared distance measure
    distx = (betax - b1)^2
    dist1 = (altsol1 - b1)^2


    # create tibble object that contains results
    result_beta <- tribble(
      ~Name,                           ~Value,
      "Beta*",                         round(betax,6),
      "(Beta*-Beta controlled)^2",     round(distx,6),
      "Alternative Solution 1",        round(altsol1,6),
      "(Beta[AS1]-Beta controlled)^2", round(dist1,6),
      "Uncontrolled Coefficient",      b0,
      "Controlled Coefficient",        b1,
      "Uncontrolled R-square",         R20,
      "Controlled R-square",           R21,
      "Max R-square",                  R2max,
      "Detla (defined)",               delta
    )


    # warning if max R-square is smaller than the R-square of the control model
    if (R21 > R2max)
      warning("The max R-square value is smaller than the R-square of the controlled model")


    #define return object
    return(result_beta)
  }


  if (delta != 1) { # solve for cubic roots if delta is not 1


    # compute different variables
    A = t_x*b0_m_b1*sigma_xx*(delta-2)/((delta-1)*(t_x*sigma_xx-t_x^2))
    B = (delta*rm_m_rt_t_syy*(sigma_xx-t_x)-rt_m_ro_t_syy*t_x-sigma_xx*t_x*(b0_m_b1^2))/((delta-1)*(t_x*sigma_xx-t_x^2))
    C = (rm_m_rt_t_syy*delta*b0_m_b1*sigma_xx)/((delta-1)*(t_x*sigma_xx-t_x^2))
    Q = (A^2-3*B)/9
    R = (2*A^3-9*A*B+27*C)/54
    D = R^2-Q^3
    discrim = R^2-Q^3


    if (discrim <0) {

      theta = acos(R/sqrt(Q^3))

      sol1 = -2*sqrt(Q)*cos(theta/3)-(A/3)
      sol2 = -2*sqrt(Q)*cos((theta+2*pi)/3)-(A/3)
      sol3 = -2*sqrt(Q)*cos((theta-2*pi)/3)-(A/3)

      sols = c(sol1,sol2,sol3)
      sols =  b1 - sols

      dists = sols - b1
      dists = dists^2


      # change to alternative solutions if first solution violates assumption 3
      for (i in 1:3) {
        if ( sign(sols[i]-b1)!=sign(b1-b0) ) {
          dists[i]=max(dists)+1
        }
      }


      dists2 <- sort(dists)          # sort matrix
      ind1 <- match(dists2[1],dists) # find index of the min value in dist
      ind2 <- match(dists2[2],dists) # find index of the min value in dist
      ind3 <- match(dists2[3],dists) # find index of the min value in dist


      betax = sols[ind1]
      altsol1 = sols[ind2]
      altsol2 = sols[ind3]


      distx= dists[ind1]
      dist1= dists[ind2]
      dist2= dists[ind3]


    } else {


      t1=-1*R+sqrt(D)
      t2=-1*R-sqrt(D)


      crt1=sign(t1) * abs(t1)^(1/3)
      crt2=sign(t2) * abs(t2)^(1/3)

      sol1 = crt1+crt2-(A/3)
      betax=b1-sol1
    }


    # create tibble object that contains results
    result_beta <- tribble(
      ~Name, ~Value,
      "Beta*",                         round(betax,6),
      "(Beta*-Beta controlled)^2",     round(distx,6),
      "Alternative Solution 1",        round(altsol1,6),
      "(Beta[AS1]-Beta controlled)^2", round(dist1,6),
      "Alternative Solution 2",        round(altsol2,6),
      "(Beta[AS2]-Beta controlled)^2", round(dist2,6),
      "Uncontrolled Coefficient",      b0,
      "Controlled Coefficient",        b1,
      "Uncontrolled R-square",         R20,
      "Controlled R-square",           R21,
      "Max R-square",                  R2max,
      "Detla (defined)",               delta
    )


    # warning if max R-square is smaller than the R-square of the control model
    if (R21 > R2max)
      warning("The max R-square value is smaller than the R-square of the controlled model")


    # define return object
    return(result_beta)
  }

}



##################--------------------------------------------------------------
#------------------------------- 1.8 o_beta_rsq - function 8
##################--------------------------------------------------------------
#' @title Betas* over a range of max R-squares
#' @description Estimates betas*, i.e. the estimated bias-adjusted treatment effects, following Oster (2019) over a range of max R-squares.
#' @usage o_beta_rsq(y, x, con, id = "none", time = "none", delta = 1, type, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent variable of interest (treatment variable; as string).
#' @param con Name of the other control variables. Provided as string in the format: "w + z +...".
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect models.
#' @param time Name of the time variable (e.g. year or month; as string). Only applicable for fixed effect models.
#' @param delta Delta for which beta* should be estimated (default is delta = 1).
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param data Data.
#' @details Estimates betas*, i.e. the estimated bias-adjusted treatment effects, following Oster (2019) over a range of max R-squares. The range of max R-squares starts from the R-square of the controlled model rounded up to the next 1/100 to 1. The function supports linear cross sectional (see \emph{lm} objects in R) and panel fixed effect (see \emph{plm} objects in R) models.
#' @return Returns tibble object. Including betas* over a range of max R-squares.
#' @references Oster, E. (2019). Unobservable selection and coefficient stability: Theory and evidence. Journal of Business & Economic Statistics, 37, 187-204.
#' @examples
#' # load data, e.g. the in-build mtcars dataset
#' data("mtcars")
#' data_oster <- mtcars
#'
#' # preview of data
#' head(data_oster)
#'
#' # load robomit
#' require(robomit)
#'
#' # estimate deltas* over a range of max R-squares
#' o_beta_rsq(y = "mpg",                 # define the dependent variable name
#'                 x = "wt",             # define the main independent variable name
#'                 con = "hp + qsec",    # other control variables
#'                 delta = 1,            # define beta. This is usually set to 1
#'                 type = "lm",          # define model type
#'                 data = data_oster)    # define dataset
#' @export

o_beta_rsq <- function(y, x, con, id = "none", time = "none", delta = 1, type, data) {



  # rename variables and create, if type is plm, a panel dataset
  ## lm
  if (type == "lm") {data <- data %>% dplyr::rename(y_var = y, x_var = x)}
  ## plm
  if (type == "plm") {data <- data %>% dplyr::rename(y_var = y, x_var = x, id_var = id, time_var = time)
  data_plm <- pdata.frame(data, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)}


  # define formulas of the different models
  model0_formula    <- as.formula("y_var ~ x_var")              # uncontrolled model
  model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) # controlled model
  aux_model_formula <- as.formula(paste("x_var ~",con))         # auxiliary model

  # define model to compute sigma_xx when model type is plm
  if (type == "plm") {sigma_xx_model_formula <- as.formula(paste("x_var ~ factor(id_var)"))}

  # run models
  ## lm
  if (type == "lm") {
    model0    <- lm(model0_formula,    data = data, na.action = na.exclude) # run uncontrolled model
    model1    <- lm(model1_formula,    data = data, na.action = na.exclude) # run controlled model
    aux_model <- lm(aux_model_formula, data = data, na.action = na.exclude) # run auxiliary model
  }
  ## plm
  if (type == "plm") {
    model0    <- plm(model0_formula,    data = data_plm, model = "within", na.action = na.exclude) # run uncontrolled model
    model1    <- plm(model1_formula,    data = data_plm, model = "within", na.action = na.exclude) # run controlled model
    aux_model <- plm(aux_model_formula, data = data_plm, model = "within", na.action = na.exclude) # run auxiliary model
    model_xx  <-  lm(sigma_xx_model_formula, data = data, na.action = na.exclude)                  # run model to obtain singa_xx
  }


  # variables based on model outputs
  if (type == "lm") {b0 = as.numeric(tidy(model0)[2,2])} else if (type == "plm") {b0 = as.numeric(tidy(model0)[1,2])}           # beta uncontrolled model
  if (type == "lm") {b1 = as.numeric(tidy(model1)[2,2])} else if (type == "plm") {b1 = as.numeric(tidy(model1)[1,2])}           # beta controlled model
  if (type == "lm") {R20 = summary(model0)$r.squared} else if (type == "plm") {R20 = as.numeric(summary(model0)$r.squared[1])}  # R-square uncontrolled model
  if (type == "lm") {R21 = summary(model1)$r.squared} else if (type == "plm") {R21 = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model
  sigma_yy = var(data$y_var, na.rm = T)                                                                                         # variance of dependent variable
  if (type == "lm") {sigma_xx = var(data$x_var, na.rm = T)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)}       # variance of independent variable
  t_x  = var(aux_model$residuals)


  # define max R-square range
  r_start <- ceiling(R21/0.01)*0.01 # start max R-square at the next highest 1/100
  n_runs <- length(seq(r_start:1,by = 0.01))
  beta__over_rsq_results <- tibble(
    "x" = 1:n_runs,
    "Rmax" = 0,
    "Beta*" = 0) # create tibble object to store results

  # run loop over different max R-squares
  for (i in 1:n_runs) {
    R2max = seq(r_start:1,by = 0.01)[i] # current max R-square

    # create some additional variables
    rt_m_ro_t_syy = (R21-R20) * sigma_yy
    b0_m_b1 = b0 - b1
    rm_m_rt_t_syy = (R2max - R21) * sigma_yy


    if (delta == 1) { # quadratic function: v = (- cap_theta +/- sqrt ( cap_theta + d1_1))/(d1_2) (see psacalc command in Stata)


      # compute different variables
      cap_theta = rm_m_rt_t_syy*(sigma_xx-t_x)-rt_m_ro_t_syy*t_x-sigma_xx*t_x*(b0_m_b1^2)
      d1_1 = 4*rm_m_rt_t_syy*(b0_m_b1^2)*(sigma_xx^2)*t_x
      d1_2 = -2*t_x*b0_m_b1*sigma_xx

      sol1 = (-1*cap_theta-sqrt((cap_theta)^2+d1_1))/(d1_2)
      sol2 = (-1*cap_theta+sqrt((cap_theta)^2+d1_1))/(d1_2)

      beta1 = b1 - sol1
      beta2 = b1 - sol2

      if ( (beta1-b1)^2 < (beta2-b1)^2) {
        betax = beta1
        altsol1 = beta2
      } else {
        betax = beta2
        altsol1 = beta1
      }

      # change to alternative solution if bias of b0 vs b1 is of different sign that bias of b1 - beta
      if ( sign(betax-b1)!=sign(b1-b0) ) {
        solc = betax
        betax = altsol1
        altsol1 = solc
      }


      # compute squared distance measure
      distx = (betax - b1)^2
      dist1 = (altsol1 - b1)^2
    }


    if (delta != 1) { # solve for cubic roots if delta is not 1


      # compute different variables
      A = t_x*b0_m_b1*sigma_xx*(delta-2)/((delta-1)*(t_x*sigma_xx-t_x^2))
      B = (delta*rm_m_rt_t_syy*(sigma_xx-t_x)-rt_m_ro_t_syy*t_x-sigma_xx*t_x*(b0_m_b1^2))/((delta-1)*(t_x*sigma_xx-t_x^2))
      C = (rm_m_rt_t_syy*delta*b0_m_b1*sigma_xx)/((delta-1)*(t_x*sigma_xx-t_x^2))
      Q = (A^2-3*B)/9
      R = (2*A^3-9*A*B+27*C)/54
      D = R^2-Q^3
      discrim = R^2-Q^3


      if (discrim <0) {

        theta = acos(R/sqrt(Q^3))

        sol1 = -2*sqrt(Q)*cos(theta/3)-(A/3)
        sol2 = -2*sqrt(Q)*cos((theta+2*pi)/3)-(A/3)
        sol3 = -2*sqrt(Q)*cos((theta-2*pi)/3)-(A/3)

        sols = c(sol1,sol2,sol3)
        sols =  b1 - sols

        dists = sols - b1
        dists = dists^2


        # change to alternative solutions if first solution violates assumption 3
        for (i in 1:3) {
          if ( sign(sols[i]-b1)!=sign(b1-b0) ) {
            dists[i]=max(dists)+1
          }
        }


        dists2 <- sort(dists)          # sort matrix
        ind1 <- match(dists2[1],dists) # find index of the min value in dist
        ind2 <- match(dists2[2],dists) # find index of the min value in dist
        ind3 <- match(dists2[3],dists) # find index of the min value in dist


        betax = sols[ind1]
        altsol1 = sols[ind2]
        altsol2 = sols[ind3]


        distx= dists[ind1]
        dist1= dists[ind2]
        dist2= dists[ind3]

      } else {


        t1=-1*R+sqrt(D)
        t2=-1*R-sqrt(D)


        crt1=sign(t1) * abs(t1)^(1/3)
        crt2=sign(t2) * abs(t2)^(1/3)


        sol1 = crt1+crt2-(A/3)
        betax=b1-sol1
      }
    }

    beta__over_rsq_results[i,2] <- R2max
    beta__over_rsq_results[i,3] <- round(betax,6)

  }


  # define return object
  return(beta__over_rsq_results)
}



##################--------------------------------------------------------------
#------------------------------- 1.9 o_beta_rsq_viz - function 9
##################--------------------------------------------------------------
#' @title Visualization of betas* over a range of max R-squares
#'
#' @description Estimates and visualizes betas*, i.e. the estimated bias-adjusted treatment effects, following Oster (2019) over a range of max R-squares.
#' @usage o_beta_rsq_viz(y, x, con, id = "none", time = "none", delta = 1, type, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent variable of interest (treatment variable; as string).
#' @param con Name of the other control variables. Provided as string in the format: "w + z +...".
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect models.
#' @param time Name of the time variable (e.g. year or month; as string). Only applicable for fixed effect models.
#' @param delta Delta for which beta* should be estimated (default is delta = 1).
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param data Data.
#' @details For details about the estimation see Oster, E. (2019). Unobservable selection and coefficient stability: Theory and evidence. Journal of Business & Economic Statistics, 37, 187-204.
#' @details Estimates and visualizes betas*, i.e. the estimated bias-adjusted treatment effects, following Oster (2019) over a range of max R-squares. The range of max R-squares starts from the R-square of the controlled model rounded up to the next 1/100 to 1. The function supports linear cross sectional (see \emph{lm} objects in R) and panel fixed effect (see \emph{plm} objects in R) models.
#' @return Returns ggplot object. Including betas* over a range of max R-squares.
#' @references Oster, E. (2019). Unobservable selection and coefficient stability: Theory and evidence. Journal of Business & Economic Statistics, 37, 187-204.
#' @examples
#' # load data, e.g. the in-build mtcars dataset
#' data("mtcars")
#' data_oster <- mtcars
#'
#' # preview of data
#' head(data_oster)
#'
#' # load robomit
#' require(robomit)
#'
#' # estimate and visualize betas* over a range of max R-squares
#' o_beta_rsq_viz(y = "mpg",                # define the dependent variable name
#'                    x = "wt",             # define the main independent variable name
#'                    con = "hp + qsec",    # other control variables
#'                    delta = 1,            # define delta This is usually set to 1
#'                    type = "lm",          # define model type
#'                    data = data_oster)    # define dataset
#' @export


o_beta_rsq_viz <- function(y, x, con, id = "none", time = "none", delta = 1, type, data) {

  # rename variables and create, if type is plm, a panel dataset
  ## lm
  if (type == "lm") {data <- data %>% dplyr::rename(y_var = y, x_var = x)}
  ## plm
  if (type == "plm") {data <- data %>% dplyr::rename(y_var = y, x_var = x, id_var = id, time_var = time)
  data_plm <- pdata.frame(data, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)}


  # define formulas of the different models
  model0_formula    <- as.formula("y_var ~ x_var")              # uncontrolled model
  model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) # controlled model
  aux_model_formula <- as.formula(paste("x_var ~",con))         # auxiliary model

  # define model to compute sigma_xx when model type is plm
  if (type == "plm") {sigma_xx_model_formula <- as.formula(paste("x_var ~",con,"+ factor(id_var)"))}

  # run models
  ## lm
  if (type == "lm") {
    model0    <- lm(model0_formula,    data = data, na.action = na.exclude) # run uncontrolled model
    model1    <- lm(model1_formula,    data = data, na.action = na.exclude) # run controlled model
    aux_model <- lm(aux_model_formula, data = data, na.action = na.exclude) # run auxiliary model
  }
  ## plm
  if (type == "plm") {
    model0    <- plm(model0_formula,    data = data_plm, model = "within", na.action = na.exclude) # run uncontrolled model
    model1    <- plm(model1_formula,    data = data_plm, model = "within", na.action = na.exclude) # run controlled model
    aux_model <- plm(aux_model_formula, data = data_plm, model = "within", na.action = na.exclude) # run auxiliary model
    model_xx  <-  lm(sigma_xx_model_formula, data = data, na.action = na.exclude)                  # run model to obtain singa_xx
  }


  # variables based on model outputs
  if (type == "lm") {b0 = as.numeric(tidy(model0)[2,2])} else if (type == "plm") {b0 = as.numeric(tidy(model0)[1,2])}           # beta uncontrolled model
  if (type == "lm") {b1 = as.numeric(tidy(model1)[2,2])} else if (type == "plm") {b1 = as.numeric(tidy(model1)[1,2])}           # beta controlled model
  if (type == "lm") {R20 = summary(model0)$r.squared} else if (type == "plm") {R20 = as.numeric(summary(model0)$r.squared[1])}  # R-square uncontrolled model
  if (type == "lm") {R21 = summary(model1)$r.squared} else if (type == "plm") {R21 = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model
  sigma_yy = var(data$y_var, na.rm = T)                                                                                         # variance of dependent variable
  if (type == "lm") {sigma_xx = var(data$x_var, na.rm = T)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)}       # variance of independent variable
  t_x  = var(aux_model$residuals)

  # define max R-square range
  r_start <- ceiling(R21/0.01)*0.01 # start max R-square at the next highest 1/100
  n_runs <- length(seq(r_start:1,by = 0.01))
  beta__over_rsq_results <- tibble(
    x = 1:n_runs,
    Rmax = 0,
    Beta = 0) # create tibble object to store results

  # run loop over different max R-squares
  for (i in 1:n_runs) {
    R2max = seq(r_start:1,by = 0.01)[i] # current max R-square

    # create some additional variables
    rt_m_ro_t_syy = (R21-R20) * sigma_yy
    b0_m_b1 = b0 - b1
    rm_m_rt_t_syy = (R2max - R21) * sigma_yy


    if (delta == 1) { # quadratic function: v = (- cap_theta +/- sqrt ( cap_theta + d1_1))/(d1_2) (see psacalc command in Stata)


      # compute different variables
      cap_theta = rm_m_rt_t_syy*(sigma_xx-t_x)-rt_m_ro_t_syy*t_x-sigma_xx*t_x*(b0_m_b1^2)
      d1_1 = 4*rm_m_rt_t_syy*(b0_m_b1^2)*(sigma_xx^2)*t_x
      d1_2 = -2*t_x*b0_m_b1*sigma_xx

      sol1 = (-1*cap_theta-sqrt((cap_theta)^2+d1_1))/(d1_2)
      sol2 = (-1*cap_theta+sqrt((cap_theta)^2+d1_1))/(d1_2)

      beta1 = b1 - sol1
      beta2 = b1 - sol2

      if ( (beta1-b1)^2 < (beta2-b1)^2) {
        betax = beta1
        altsol1 = beta2
      } else {
        betax = beta2
        altsol1 = beta1
      }

      # change to alternative solution if bias of b0 vs b1 is of different sign that bias of b1 - beta
      if ( sign(betax-b1)!=sign(b1-b0) ) {
        solc = betax
        betax = altsol1
        altsol1 = solc
      }


      distx = (betax - b1)^2
      dist1 = (altsol1 - b1)^2

    }


    if (delta != 1) { # solve for cubic roots if delta is not 1


      # compute different variables
      A = t_x*b0_m_b1*sigma_xx*(delta-2)/((delta-1)*(t_x*sigma_xx-t_x^2))
      B = (delta*rm_m_rt_t_syy*(sigma_xx-t_x)-rt_m_ro_t_syy*t_x-sigma_xx*t_x*(b0_m_b1^2))/((delta-1)*(t_x*sigma_xx-t_x^2))
      C = (rm_m_rt_t_syy*delta*b0_m_b1*sigma_xx)/((delta-1)*(t_x*sigma_xx-t_x^2))
      Q = (A^2-3*B)/9
      R = (2*A^3-9*A*B+27*C)/54
      D = R^2-Q^3
      discrim = R^2-Q^3


      if (discrim <0) {

        theta = acos(R/sqrt(Q^3))

        sol1 = -2*sqrt(Q)*cos(theta/3)-(A/3)
        sol2 = -2*sqrt(Q)*cos((theta+2*pi)/3)-(A/3)
        sol3 = -2*sqrt(Q)*cos((theta-2*pi)/3)-(A/3)

        sols = c(sol1,sol2,sol3)
        sols =  b1 - sols

        dists = sols - b1
        dists = dists^2


        # change to alternative solutions if first solution violates assumption 3
        for (i in 1:3) {
          if ( sign(sols[i]-b1)!=sign(b1-b0) ) {
            dists[i]=max(dists)+1
          }
        }


        dists2 <- sort(dists)          # sort matrix
        ind1 <- match(dists2[1],dists) # find index of the min value in dist
        ind2 <- match(dists2[2],dists) # find index of the min value in dist
        ind3 <- match(dists2[3],dists) # find index of the min value in dist


        betax = sols[ind1]
        altsol1 = sols[ind2]
        altsol2 = sols[ind3]


        distx= dists[ind1]
        dist1= dists[ind2]
        dist2= dists[ind3]
      } else {


        t1=-1*R+sqrt(D)
        t2=-1*R-sqrt(D)


        crt1=sign(t1) * abs(t1)^(1/3)
        crt2=sign(t2) * abs(t2)^(1/3)


        sol1 = crt1+crt2-(A/3)
        betax=b1-sol1
      }
    }

    beta__over_rsq_results[i,2] <- R2max
    beta__over_rsq_results[i,3] <- round(betax,6)
  }


  # plot figure
  theme_set(theme_bw()) # set main design
  result_plot <- ggplot(data=beta__over_rsq_results, aes(x=Rmax, y=Beta)) +
    geom_line(size = 1.3)+
    scale_y_continuous(name = expression(beta^"*"))+
    scale_x_continuous(name = expression("max R"^2))+
    theme(axis.title = element_text( size=15),
          axis.text  = element_text( size=13))


  # define return object
  return(result_plot)
}





##################--------------------------------------------------------------
#------------------------------- 1.10 o_beta_boot - function 10
##################--------------------------------------------------------------
#' @title Bootstrapped betas*
#'
#' @description Estimates bootstrapped betas*, i.e. the estimated bias-adjusted treatment effects, following Oster (2019).
#' @usage o_beta_boot(y, x, con, id = "none", time = "none", delta = 1, R2max, sim, obs, rep,
#' type, useed = NA, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent variable of interest (treatment variable; as string).
#' @param con Name of the other control variables. Provided as string in the format: "w + z +...".
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect models.
#' @param time Name of the time variable (e.g. year or month; as string). Only applicable for fixed effect models.
#' @param delta Delta for which beta* should be estimated (default is delta = 1).
#' @param R2max Max R-square for which beta* should be estimated.
#' @param sim Number of simulations.
#' @param obs Number of draws per simulation.
#' @param rep Bootstrapping either with (= TRUE) or without (= FALSE) replacement.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param useed Seed number defined by user.
#' @param data Data.
#' @details Estimates bootstrapped betas*, i.e. the estimated bias-adjusted treatment effects, following Oster (2019). Bootstrapping can either be done with or without replacement. The function supports linear cross sectional (see \emph{lm} objects in R) and panel fixed effect (see \emph{plm} objects in R) models.
#' @return Returns tibble object. Including bootstrapped betas*.
#' @references Oster, E. (2019). Unobservable selection and coefficient stability: Theory and evidence. Journal of Business & Economic Statistics, 37, 187-204.
#' @examples
#' # load data, e.g. the in-build mtcars dataset
#' data("mtcars")
#' data_oster <- mtcars
#'
#' # preview of data
#' head(data_oster)
#'
#' # load robomit
#' require(robomit)
#'
#' # estimate bootstrapped beta*
#' o_beta_boot(y = "mpg",                 # define the dependent variable name
#'                  x = "wt",             # define the main independent variable name
#'                  con = "hp + qsec",    # other control variables
#'                  delta = 1,            # define beta. This is usually set to 1
#'                  R2max = 0.9,          # define the max R-square.
#'                  sim = 100,            # define number of simulations
#'                  obs = 30,             # define number of drawn observations per simulation
#'                  rep = FALSE,          # define if bootstrapping is with or without replacement
#'                  type = "lm",          # define model type
#'                  useed = 123,          # define seed
#'                  data = data_oster)    # define dataset
#' @export


o_beta_boot <- function(y, x, con, id = "none", time = "none", delta = 1, R2max, sim, obs, rep, type, useed = NA, data) {


  # rename variables and create
  ## lm
  if (type == "lm") {data <- data %>% dplyr::rename(y_var = y, x_var = x)}
  ## plm
  if (type == "plm") {data <- data %>% dplyr::rename(y_var = y, x_var = x, id_var = id, time_var = time)}


  # define formulas of the different models
  model0_formula    <- as.formula("y_var ~ x_var")              # uncontrolled model
  model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) # controlled model
  aux_model_formula <- as.formula(paste("x_var ~",con))         # auxiliary model

  # define model to compute sigma_xx when model type is plm
  if (type == "plm") {sigma_xx_model_formula <- as.formula(paste("x_var ~",con,"+ factor(id_var)"))}


  # create tibble object to store results
  simulation_results <- tibble(
    "x" = 1:sim,
    "Beta*" = 0)


  for (s in 1:sim) {

    # fun models
    ## lm
    if (type == "lm") {

      if (!is.na(useed)) {set.seed(useed+s) } # set seed (defined by user) that it varies with runs
      if (rep) {
        data_current <- sample_n(data, size = obs, replace = TRUE) # with replacement
      } else {
        data_current <- sample_n(data, size = obs, replace = FALSE) # without replacement
      }
      model0    <- lm(model0_formula,    data = data_current, na.action = na.exclude) # run uncontrolled model
      model1    <- lm(model1_formula,    data = data_current, na.action = na.exclude) # run controlled model
      aux_model <- lm(aux_model_formula, data = data_current, na.action = na.exclude) # run auxiliary model
    }
    ## plm
    if (type == "plm") {

      if (!is.na(useed)) {set.seed(useed+s) } # set seed (defined by user) that it varies with runs
      if (rep) {
        data_current     <- sample_n(data, size = obs, replace = TRUE) # with replacement
        data_current_plm <- pdata.frame(data_current, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)
      } else {
        data_current     <- sample_n(data, size = obs, replace = FALSE) # without replacement
        data_current_plm <- pdata.frame(data_current, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)
      }
      model0    <- plm(model0_formula,    data = data_current_plm, model = "within", na.action = na.exclude) # run uncontrolled model
      model1    <- plm(model1_formula,    data = data_current_plm, model = "within", na.action = na.exclude) # run controlled model
      aux_model <- plm(aux_model_formula, data = data_current_plm, model = "within", na.action = na.exclude) # run auxiliary model
      model_xx  <-  lm(sigma_xx_model_formula, data = data_current, na.action = na.exclude)                  # run model to obtain singa_xx
    }

    # variables based on model outputs
    if (type == "lm") {b0 = as.numeric(tidy(model0)[2,2])} else if (type == "plm") {b0 = as.numeric(tidy(model0)[1,2])}           # beta uncontrolled model
    if (type == "lm") {b1 = as.numeric(tidy(model1)[2,2])} else if (type == "plm") {b1 = as.numeric(tidy(model1)[1,2])}           # beta controlled model
    if (type == "lm") {R20 = summary(model0)$r.squared} else if (type == "plm") {R20 = as.numeric(summary(model0)$r.squared[1])}  # R-square uncontrolled model
    if (type == "lm") {R21 = summary(model1)$r.squared} else if (type == "plm") {R21 = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model
    sigma_yy = var(data_current$y_var, na.rm = T)                                                                                          # variance of dependent variable
    if (type == "lm") {sigma_xx = var(data_current$x_var, na.rm = T)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)}       # variance of independent variable
    t_x  = var(aux_model$residuals)

    # create some additional variables
    rt_m_ro_t_syy = (R21-R20) * sigma_yy
    b0_m_b1 = b0 - b1
    rm_m_rt_t_syy = (R2max - R21) * sigma_yy



    if (delta == 1) { # quadratic function: v = (- cap_theta +/- sqrt ( cap_theta + d1_1))/(d1_2) (see psacalc command in Stata)


      # compute different variables
      cap_theta = rm_m_rt_t_syy*(sigma_xx-t_x)-rt_m_ro_t_syy*t_x-sigma_xx*t_x*(b0_m_b1^2)
      d1_1 = 4*rm_m_rt_t_syy*(b0_m_b1^2)*(sigma_xx^2)*t_x
      d1_2 = -2*t_x*b0_m_b1*sigma_xx

      sol1 = (-1*cap_theta-sqrt((cap_theta)^2+d1_1))/(d1_2)
      sol2 = (-1*cap_theta+sqrt((cap_theta)^2+d1_1))/(d1_2)

      beta1 = b1 - sol1
      beta2 = b1 - sol2

      if ( (beta1-b1)^2 < (beta2-b1)^2) {
        betax = beta1
        altsol1 = beta2
      } else {
        betax = beta2
        altsol1 = beta1
      }

      # change to alternative solution if bias of b0 vs b1 is of different sign that bias of b1 - beta
      if ( sign(betax-b1)!=sign(b1-b0) ) {
        solc = betax
        betax = altsol1
        altsol1 = solc
      }


      distx = (betax - b1)^2
      dist1 = (altsol1 - b1)^2



    }


    if (delta != 1) { # solve for cubic roots if delta is not 1


      # compute different variables
      A = t_x*b0_m_b1*sigma_xx*(delta-2)/((delta-1)*(t_x*sigma_xx-t_x^2))
      B = (delta*rm_m_rt_t_syy*(sigma_xx-t_x)-rt_m_ro_t_syy*t_x-sigma_xx*t_x*(b0_m_b1^2))/((delta-1)*(t_x*sigma_xx-t_x^2))
      C = (rm_m_rt_t_syy*delta*b0_m_b1*sigma_xx)/((delta-1)*(t_x*sigma_xx-t_x^2))
      Q = (A^2-3*B)/9
      R = (2*A^3-9*A*B+27*C)/54
      D = R^2-Q^3
      discrim = R^2-Q^3


      if (discrim <0) {

        theta = acos(R/sqrt(Q^3))

        sol1 = -2*sqrt(Q)*cos(theta/3)-(A/3)
        sol2 = -2*sqrt(Q)*cos((theta+2*pi)/3)-(A/3)
        sol3 = -2*sqrt(Q)*cos((theta-2*pi)/3)-(A/3)

        sols = c(sol1,sol2,sol3)
        sols =  b1 - sols

        dists = sols - b1
        dists = dists^2


        # change to alternative solutions if first solution violates assumption 3
        for (i in 1:3) {
          if ( sign(sols[i]-b1)!=sign(b1-b0) ) {
            dists[i]=max(dists)+1
          }
        }


        dists2 <- sort(dists)          # sort matrix
        ind1 <- match(dists2[1],dists) # find index of the min value in dist
        ind2 <- match(dists2[2],dists) # find index of the min value in dist
        ind3 <- match(dists2[3],dists) # find index of the min value in dist


        betax = sols[ind1]
        altsol1 = sols[ind2]
        altsol2 = sols[ind3]


        distx= dists[ind1]
        dist1= dists[ind2]
        dist2= dists[ind3]
      } else {


        t1=-1*R+sqrt(D)
        t2=-1*R-sqrt(D)


        crt1=sign(t1) * abs(t1)^(1/3)
        crt2=sign(t2) * abs(t2)^(1/3)


        sol1 = crt1+crt2-(A/3)
        betax=b1-sol1
      }
    }

    # store results
    simulation_results[s,2] <- round(betax,6)
  }

  # warning if max R-square is smaller than the R-square of the control model
  if (R21 > R2max)
    warning("The max R-square value is smaller than the R-square of the controlled model")


  # warning if the bootstrapped observations are larger than the number of observations of the original dataset
  if (obs >= length(data$y_var))
    warning("Number of bootstrapped observation is larger than/equal to number of observation in the provided dataset")


  # define return object
  return(simulation_results)
}



##################--------------------------------------------------------------
#------------------------------- 1.11 o_beta_boot_viz - function 11
##################--------------------------------------------------------------
#' @title Visualization of bootstrapped betas*
#'
#' @description Estimates and visualizes bootstrapped betas*, i.e. the estimated bias-adjusted treatment effects, following Oster (2019).
#' @usage o_beta_boot_viz(y, x, con, id = "none", time = "none", delta = 1, R2max, sim, obs, rep,
#' CI, type, norm = TRUE, bin, col = c("#08306b","#4292c6","#c6dbef"),
#' nL = TRUE, mL = TRUE, useed = NA, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent variable of interest (treatment variable; as string).
#' @param con Name of the other control variables. Provided as string in the format: "w + z +...".
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect models.
#' @param time Name of the time variable (e.g. year or month; as string). Only applicable for fixed effect models.
#' @param delta Delta for which beta* should be estimated (default is delta = 1).
#' @param R2max Max R-square for which beta* should be estimated.
#' @param sim Number of simulations.
#' @param obs Number of draws per simulation.
#' @param rep Bootstrapping either with (= TRUE) or without (= FALSE) replacement
#' @param CI  Confidence intervals, indicated as vector. Can be and/or 90,95,99.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param norm Option to include a normal distribution in the plot (default is norm = TURE).
#' @param bin Number of bins used for the histogram.
#' @param col Colors used to indicate different confidence interval levels (indicated as vector). Needs to be the same length as the variable CI. The default is a blue color range.
#' @param nL Option to include a red vertical line at 0 (default is nL = TRUE).
#' @param mL Option to include a vertical line at beta* mean (default is mL = TRUE).
#' @param useed Seed number defined by user.
#' @param data Data.
#' @details Estimates and visualizes bootstrapped betas*, i.e. the estimated bias-adjusted treatment effects, following Oster (2019). Bootstrapping can either be done with or without replacement. The function supports linear cross sectional (see \emph{lm} objects in R) and panel fixed effect (see \emph{plm} objects in R) models.
#' @return Returns ggplot object. Including bootstrapped betas*.
#' @references Oster, E. (2019). Unobservable selection and coefficient stability: Theory and evidence. Journal of Business & Economic Statistics, 37, 187-204.
#' @examples
#' # load data, e.g. the in-build mtcars dataset
#' data("mtcars")
#' data_oster <- mtcars
#'
#' # preview of data
#' head(data_oster)
#'
#' # load robomit
#' require(robomit)
#'
#' # estimate and visualize bootstrapped betas*
#' o_beta_boot_viz(y = "mpg",                # define the dependent variable name
#'                     x = "wt",             # define the main independent variable name
#'                     con = "hp + qsec",    # other control variables
#'                     delta = 1,            # define delta This is usually set to 1
#'                     R2max = 0.9,          # define the max R-square.
#'                     sim = 100,            # define number of simulations
#'                     obs = 30,             # define number of drawn observations per simulation
#'                     rep = FALSE,          # define if bootstrapping is with or without replacement
#'                     CI = c(90,95,99),     # define confidence intervals.
#'                     type = "lm",          # define model type
#'                     norm = TRUE,          # include normal distribution
#'                     bin = 200,            # set number of bins
#'                     useed = 123,          # define seed
#'                     data = data_oster)    # define dataset
#' @export


o_beta_boot_viz <- function(y, x, con, id = "none", time = "none", delta = 1, R2max, sim, obs, rep, CI, type, norm = TRUE, bin, col = c("#08306b","#4292c6","#c6dbef"), nL = TRUE, mL = TRUE, useed = NA, data) {



  # rename variables and create, if type is plm, a panel dataset
  ## lm
  if (type == "lm") {data <- data %>% dplyr::rename(y_var = y, x_var = x)}
  ## plm
  if (type == "plm") {data <- data %>% dplyr::rename(y_var = y, x_var = x, id_var = id, time_var = time)}


  # define formulas of the different models
  model0_formula    <- as.formula("y_var ~ x_var")              # uncontrolled model
  model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) # controlled model
  aux_model_formula <- as.formula(paste("x_var ~",con))         # auxiliary model

  # define model to compute sigma_xx when model type is plm
  if (type == "plm") {sigma_xx_model_formula <- as.formula(paste("x_var ~",con,"+ factor(id_var)"))}


  # create tibble object to store results
  simulation_results <- tibble(
    x = 1:sim,
    Beta = 0)


  for (s in 1:sim) {

    # fun models
    ## lm
    if (type == "lm") {

      if (!is.na(useed)) {set.seed(useed+s) } # set seed (defined by user) that it varies with runs
      if (rep) {
        data_current <- sample_n(data, size = obs, replace = TRUE) # with replacement
      } else {
        data_current <- sample_n(data, size = obs, replace = FALSE) # without replacement
      }
      model0    <- lm(model0_formula,    data = data_current, na.action = na.exclude) # run uncontrolled model
      model1    <- lm(model1_formula,    data = data_current, na.action = na.exclude) # run controlled model
      aux_model <- lm(aux_model_formula, data = data_current, na.action = na.exclude) # run auxiliary model
    }
    ## plm
    if (type == "plm") {

      if (!is.na(useed)) {set.seed(useed+s) } # set seed (defined by user) that it varies with runs
      if (rep) {
        data_current     <- sample_n(data, size = obs, replace = TRUE) # with replacement
        data_current_plm <- pdata.frame(data_current, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)
      } else {
        data_current     <- sample_n(data, size = obs, replace = FALSE) # without replacement
        data_current_plm <- pdata.frame(data_current, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)
      }
      model0    <- plm(model0_formula,    data = data_current_plm, model = "within", na.action = na.exclude) # run uncontrolled model
      model1    <- plm(model1_formula,    data = data_current_plm, model = "within", na.action = na.exclude) # run controlled model
      aux_model <- plm(aux_model_formula, data = data_current_plm, model = "within", na.action = na.exclude) # run auxiliary model
      model_xx  <-  lm(sigma_xx_model_formula, data = data_current, na.action = na.exclude)                  # run model to obtain singa_xx
    }


    # variables based on model outputs
    if (type == "lm") {b0 = as.numeric(tidy(model0)[2,2])} else if (type == "plm") {b0 = as.numeric(tidy(model0)[1,2])}           # beta uncontrolled model
    if (type == "lm") {b1 = as.numeric(tidy(model1)[2,2])} else if (type == "plm") {b1 = as.numeric(tidy(model1)[1,2])}           # beta controlled model
    if (type == "lm") {R20 = summary(model0)$r.squared} else if (type == "plm") {R20 = as.numeric(summary(model0)$r.squared[1])}  # R-square uncontrolled model
    if (type == "lm") {R21 = summary(model1)$r.squared} else if (type == "plm") {R21 = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model
    sigma_yy = var(data_current$y_var, na.rm = T)                                                                                          # variance of dependent variable
    if (type == "lm") {sigma_xx = var(data_current$x_var, na.rm = T)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)}       # variance of independent variable
    t_x  = var(aux_model$residuals)

    # create some additional variables
    rt_m_ro_t_syy = (R21-R20) * sigma_yy
    b0_m_b1 = b0 - b1
    rm_m_rt_t_syy = (R2max - R21) * sigma_yy



    if (delta == 1) { # quadratic function: v = (- cap_theta +/- sqrt ( cap_theta + d1_1))/(d1_2) (see psacalc command in Stata)


      # compute different variables
      cap_theta = rm_m_rt_t_syy*(sigma_xx-t_x)-rt_m_ro_t_syy*t_x-sigma_xx*t_x*(b0_m_b1^2)
      d1_1 = 4*rm_m_rt_t_syy*(b0_m_b1^2)*(sigma_xx^2)*t_x
      d1_2 = -2*t_x*b0_m_b1*sigma_xx

      sol1 = (-1*cap_theta-sqrt((cap_theta)^2+d1_1))/(d1_2)
      sol2 = (-1*cap_theta+sqrt((cap_theta)^2+d1_1))/(d1_2)

      beta1 = b1 - sol1
      beta2 = b1 - sol2

      if ( (beta1-b1)^2 < (beta2-b1)^2) {
        betax = beta1
        altsol1 = beta2
      } else {
        betax = beta2
        altsol1 = beta1
      }

      # change to alternative solution if bias of b0 vs b1 is of different sign that bias of b1 - beta
      if ( sign(betax-b1)!=sign(b1-b0) ) {
        solc = betax
        betax = altsol1
        altsol1 = solc
      }


      distx = (betax - b1)^2
      dist1 = (altsol1 - b1)^2
    }


    if (delta != 1) { # solve for cubic roots if delta is not 1


      # compute different variables
      A = t_x*b0_m_b1*sigma_xx*(delta-2)/((delta-1)*(t_x*sigma_xx-t_x^2))
      B = (delta*rm_m_rt_t_syy*(sigma_xx-t_x)-rt_m_ro_t_syy*t_x-sigma_xx*t_x*(b0_m_b1^2))/((delta-1)*(t_x*sigma_xx-t_x^2))
      C = (rm_m_rt_t_syy*delta*b0_m_b1*sigma_xx)/((delta-1)*(t_x*sigma_xx-t_x^2))
      Q = (A^2-3*B)/9
      R = (2*A^3-9*A*B+27*C)/54
      D = R^2-Q^3
      discrim = R^2-Q^3


      if (discrim <0) {

        theta = acos(R/sqrt(Q^3))

        sol1 = -2*sqrt(Q)*cos(theta/3)-(A/3)
        sol2 = -2*sqrt(Q)*cos((theta+2*pi)/3)-(A/3)
        sol3 = -2*sqrt(Q)*cos((theta-2*pi)/3)-(A/3)

        sols = c(sol1,sol2,sol3)
        sols =  b1 - sols

        dists = sols - b1
        dists = dists^2


        # change to alternative solutions if first solution violates assumption 3
        for (i in 1:3) {
          if ( sign(sols[i]-b1)!=sign(b1-b0) ) {
            dists[i]=max(dists)+1
          }
        }


        dists2 <- sort(dists)          # sort matrix
        ind1 <- match(dists2[1],dists) # find index of the min value in dist
        ind2 <- match(dists2[2],dists) # find index of the min value in dist
        ind3 <- match(dists2[3],dists) # find index of the min value in dist


        betax = sols[ind1]
        altsol1 = sols[ind2]
        altsol2 = sols[ind3]


        distx= dists[ind1]
        dist1= dists[ind2]
        dist2= dists[ind3]
      } else {


        t1=-1*R+sqrt(D)
        t2=-1*R-sqrt(D)


        crt1=sign(t1) * abs(t1)^(1/3)
        crt2=sign(t2) * abs(t2)^(1/3)


        sol1 = crt1+crt2-(A/3)
        betax=b1-sol1
      }
    }


    # store results
    simulation_results[s,2] <- round(betax,6)
  }


  # compute mean, max absolute distance to mean, and confidence intervals
  ## mean
  meanBeta <- mean(simulation_results$Beta, na.rm = F)
  # sd
  sdBeta <- sd(simulation_results$Beta, na.rm = F)
  ## confidence interval
  ### 90%
  CI90_a <- meanBeta - 1.645 * sdBeta
  CI90_b <- meanBeta + 1.645 * sdBeta
  ### 95%
  CI95_a <- meanBeta - 1.96 * sdBeta
  CI95_b <- meanBeta + 1.96 * sd(simulation_results$Beta, na.rm = F)
  ### 99%
  CI99_a <- meanBeta - 2.58 * sdBeta
  CI99_b <- meanBeta + 2.58 * sdBeta

  # define colors of confidence intervals
  color1 = col[1]
  color2 = col[2]
  color3 = col[3]

  ## 90%
  if (is.element(90, CI)) {
    color90 <- color1}

  ## 95%
  if (is.element(95, CI)) {
    color95 <- ifelse(length(CI) == 3,color2,
                      ifelse(length(CI) == 2 & is.element(90, CI),color2,
                             ifelse(length(CI) == 2 & is.element(99, CI),color1,color1)))}
  ## 99%
  if (is.element(99, CI)) {
    color99 <- ifelse(length(CI) == 3,color3,
                      ifelse(length(CI) == 2,color2,color1))}



  # plot figure
  theme_set(theme_bw()) # set main design
  beta_plot <- ggplot(simulation_results, aes(x = Beta)) +
  {if(is.element(99, CI))annotate("rect", ymin = 0, ymax = +Inf, xmin = CI99_a, xmax = CI99_b, fill = color99, alpha = 0.6)}+
  {if(is.element(95, CI))annotate("rect", ymin = 0, ymax = +Inf, xmin = CI95_a, xmax = CI95_b, fill = color95, alpha = 0.6)}+
  {if(is.element(90, CI))annotate("rect", ymin = 0, ymax = +Inf, xmin = CI90_a, xmax = CI90_b, fill = color90, alpha = 0.6)}+
    geom_histogram(aes(y = ..density..), alpha = 0.8, colour = "#737373", fill = "#737373", size = 0.1, bins = bin)+
    {if(norm)stat_function(fun = dnorm, args = list(mean = meanBeta, sd = sdBeta), size = 1.3)}+
    scale_y_continuous(name = "Density",expand = expansion(mult = c(0, .1)))+
    scale_x_continuous(name = expression(beta^"*"))+
    expand_limits(x = 0, y = 0)+
    theme(axis.title = element_text( size=15),
          axis.text  = element_text( size=13))+
          {if(mL)geom_vline(xintercept = meanBeta, size = 0.8, color = "#000000", alpha = 0.8)}+
          {if(nL)geom_vline(xintercept = 0, size = 0.8, color = "#e41a1c", alpha = 0.8)}+
    geom_hline(yintercept = 0)


  # warning if max R-square is smaller than the R-square of the control model
  if (R21 > R2max)
    warning("The max R-square value is smaller than the R-square of the controlled model")


  # warning if the bootstrapped observations are larger than the number of observations of the original dataset
  if (obs >= length(data$y_var))
    warning("Number of bootstrapped observation is larger than/equal to number of observation in the provided dataset")


  # define return object
  return(beta_plot)


}



##################--------------------------------------------------------------
#------------------------------- 1.12 o_beta_boot_inf - function 12
##################--------------------------------------------------------------
#' @title Bootstrapped mean beta* and confidence intervals
#'
#' @description Estimates and provides confidence intervals of bootstrapped betas*, i.e. the estimated bias-adjusted treatment effects, following Oster (2019).
#' @usage o_beta_boot_inf(y, x, con, id = "none", time = "none", delta = 1, R2max, sim, obs, rep,
#' CI, type, useed = NA, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent variable of interest (treatment variable; as string).
#' @param con Name of the other control variables. Provided as string in the format: "w + z +...".
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect models.
#' @param time Name of the time variable (e.g. year or month; as string). Only applicable for fixed effect models.
#' @param delta Delta for which beta* should be estimated (default is delta = 1).
#' @param R2max Max R-square for which beta* should be estimated.
#' @param sim Number of simulations.
#' @param obs Number of draws per simulation.
#' @param rep Bootstrapping either with (= TRUE) or without (= FALSE) replacement
#' @param CI  Confidence intervals, indicated as vector. Can be and/or 90,95,99.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param useed Seed number defined by user.
#' @param data Data.
#' @details Estimates mean and provides confidence intervals of bootstrapped betas*, i.e. the estimated bias-adjusted treatment effects, following Oster (2019). Bootstrapping can either be done with or without replacement. The function supports linear cross sectional (see \emph{lm} objects in R) and panel fixed effect (see \emph{plm} objects in R) models.
#' @return Returns tibble object. Including bootstrapped betas* and confidence intervals.
#' @references Oster, E. (2019). Unobservable selection and coefficient stability: Theory and evidence. Journal of Business & Economic Statistics, 37, 187-204.
#' @examples
#' # load data, e.g. the in-build mtcars dataset
#' data("mtcars")
#' data_oster <- mtcars
#'
#' # preview of data
#' head(data_oster)
#'
#' # load robomit
#' require(robomit)
#'
#' # estimate and visualize bootstrapped deltas*
#' o_beta_boot_inf(y = "mpg",                # define the dependent variable name
#'                     x = "wt",             # define the main independent variable name
#'                     con = "hp + qsec",    # other control variables
#'                     delta = 0,            # define delta. This is usually set to 1
#'                     R2max = 0.9,          # define the max R-square.
#'                     sim = 100,            # define number of simulations
#'                     obs = 30,             # define number of drawn observations per simulation
#'                     rep = FALSE,          # define if bootstrapping is with or without replacement
#'                     CI = c(90,95,99),     # define confidence intervals.
#'                     type = "lm",          # define model type
#'                     useed = 123,          # define seed
#'                     data = data_oster)    # define dataset
#' @export


o_beta_boot_inf <- function(y, x, con, id = "none", time = "none", delta = 1, R2max, sim, obs, rep, CI, type, useed = NA, data) {


  # rename variables and create, if type is plm, a panel dataset
  ## lm
  if (type == "lm") {data <- data %>% dplyr::rename(y_var = y, x_var = x)}
  ## plm
  if (type == "plm") {data <- data %>% dplyr::rename(y_var = y, x_var = x, id_var = id, time_var = time)
  data_plm <- pdata.frame(data, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)}


  # define formulas of the different models
  model0_formula    <- as.formula("y_var ~ x_var")              # uncontrolled model
  model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) # controlled model
  aux_model_formula <- as.formula(paste("x_var ~",con))         # auxiliary model

  # define model to compute sigma_xx when model type is plm
  if (type == "plm") {sigma_xx_model_formula <- as.formula(paste("x_var ~",con,"+ factor(id_var)"))}


  # create tibble object to store results
  simulation_results <- tibble(
    x = 1:sim,
    Beta = 0)


  for (s in 1:sim) {

    # fun models
    ## lm
    if (type == "lm") {

      if (!is.na(useed)) {set.seed(useed+s) } # set seed (defined by user) that it varies with runs
      if (rep) {
        data_current <- sample_n(data, size = obs, replace = TRUE) # with replacement
      } else {
        data_current <- sample_n(data, size = obs, replace = FALSE) # without replacement
      }
      model0    <- lm(model0_formula,    data = data_current, na.action = na.exclude) # run uncontrolled model
      model1    <- lm(model1_formula,    data = data_current, na.action = na.exclude) # run controlled model
      aux_model <- lm(aux_model_formula, data = data_current, na.action = na.exclude) # run auxiliary model
    }
    ## plm
    if (type == "plm") {

      if (!is.na(useed)) {set.seed(useed+s) } # set seed (defined by user) that it varies with runs
      if (rep) {
        data_current     <- sample_n(data, size = obs, replace = TRUE) # with replacement
        data_current_plm <- pdata.frame(data_current, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)
      } else {
        data_current     <- sample_n(data, size = obs, replace = FALSE) # without replacement
        data_current_plm <- pdata.frame(data_current, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)
      }
      model0    <- plm(model0_formula,    data = data_current_plm, model = "within", na.action = na.exclude) # run uncontrolled model
      model1    <- plm(model1_formula,    data = data_current_plm, model = "within", na.action = na.exclude) # run controlled model
      aux_model <- plm(aux_model_formula, data = data_current_plm, model = "within", na.action = na.exclude) # run auxiliary model
      model_xx  <-  lm(sigma_xx_model_formula, data = data_current, na.action = na.exclude)                  # run model to obtain singa_xx
    }


    # variables based on model outputs
    if (type == "lm") {b0 = as.numeric(tidy(model0)[2,2])} else if (type == "plm") {b0 = as.numeric(tidy(model0)[1,2])}           # beta uncontrolled model
    if (type == "lm") {b1 = as.numeric(tidy(model1)[2,2])} else if (type == "plm") {b1 = as.numeric(tidy(model1)[1,2])}           # beta controlled model
    if (type == "lm") {R20 = summary(model0)$r.squared} else if (type == "plm") {R20 = as.numeric(summary(model0)$r.squared[1])}  # R-square uncontrolled model
    if (type == "lm") {R21 = summary(model1)$r.squared} else if (type == "plm") {R21 = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model
    sigma_yy = var(data_current$y_var, na.rm = T)                                                                                          # variance of dependent variable
    if (type == "lm") {sigma_xx = var(data_current$x_var, na.rm = T)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)}       # variance of independent variable
    t_x  = var(aux_model$residuals)

    # create some additional variables
    rt_m_ro_t_syy = (R21-R20) * sigma_yy
    b0_m_b1 = b0 - b1
    rm_m_rt_t_syy = (R2max - R21) * sigma_yy



    if (delta == 1) { # quadratic function: v = (- cap_theta +/- sqrt ( cap_theta + d1_1))/(d1_2) (see psacalc command in Stata)


      # compute different variables
      cap_theta = rm_m_rt_t_syy*(sigma_xx-t_x)-rt_m_ro_t_syy*t_x-sigma_xx*t_x*(b0_m_b1^2)
      d1_1 = 4*rm_m_rt_t_syy*(b0_m_b1^2)*(sigma_xx^2)*t_x
      d1_2 = -2*t_x*b0_m_b1*sigma_xx

      sol1 = (-1*cap_theta-sqrt((cap_theta)^2+d1_1))/(d1_2)
      sol2 = (-1*cap_theta+sqrt((cap_theta)^2+d1_1))/(d1_2)

      beta1 = b1 - sol1
      beta2 = b1 - sol2

      if ( (beta1-b1)^2 < (beta2-b1)^2) {
        betax = beta1
        altsol1 = beta2
      } else {
        betax = beta2
        altsol1 = beta1
      }

      # change to alternative solution if bias of b0 vs b1 is of different sign that bias of b1 - beta
      if ( sign(betax-b1)!=sign(b1-b0) ) {
        solc = betax
        betax = altsol1
        altsol1 = solc
      }


      distx = (betax - b1)^2
      dist1 = (altsol1 - b1)^2
    }


    if (delta != 1) { # solve for cubic roots if delta is not 1

      # compute different variables
      A = t_x*b0_m_b1*sigma_xx*(delta-2)/((delta-1)*(t_x*sigma_xx-t_x^2))
      B = (delta*rm_m_rt_t_syy*(sigma_xx-t_x)-rt_m_ro_t_syy*t_x-sigma_xx*t_x*(b0_m_b1^2))/((delta-1)*(t_x*sigma_xx-t_x^2))
      C = (rm_m_rt_t_syy*delta*b0_m_b1*sigma_xx)/((delta-1)*(t_x*sigma_xx-t_x^2))
      Q = (A^2-3*B)/9
      R = (2*A^3-9*A*B+27*C)/54
      D = R^2-Q^3
      discrim = R^2-Q^3


      if (discrim <0) {

        theta = acos(R/sqrt(Q^3))

        sol1 = -2*sqrt(Q)*cos(theta/3)-(A/3)
        sol2 = -2*sqrt(Q)*cos((theta+2*pi)/3)-(A/3)
        sol3 = -2*sqrt(Q)*cos((theta-2*pi)/3)-(A/3)

        sols = c(sol1,sol2,sol3)
        sols =  b1 - sols

        dists = sols - b1
        dists = dists^2


        # change to alternative solutions if first solution violates assumption 3
        for (i in 1:3) {
          if ( sign(sols[i]-b1)!=sign(b1-b0) ) {
            dists[i]=max(dists)+1
          }
        }


        dists2 <- sort(dists)          # sort matrix
        ind1 <- match(dists2[1],dists) # find index of the min value in dist
        ind2 <- match(dists2[2],dists) # find index of the min value in dist
        ind3 <- match(dists2[3],dists) # find index of the min value in dist


        betax = sols[ind1]
        altsol1 = sols[ind2]
        altsol2 = sols[ind3]


        distx= dists[ind1]
        dist1= dists[ind2]
        dist2= dists[ind3]
      } else {


        t1=-1*R+sqrt(D)
        t2=-1*R-sqrt(D)


        crt1=sign(t1) * abs(t1)^(1/3)
        crt2=sign(t2) * abs(t2)^(1/3)


        sol1 = crt1+crt2-(A/3)
        betax=b1-sol1
      }
    }


    # store results
    simulation_results[s,2] <- round(betax,6)
  }


  # compute mean, max absolute distance to mean, and confidence intervals
  ## mean
  meanBeta <- mean(simulation_results$Beta, na.rm = F)
  # sd
  sdBeta <- sd(simulation_results$Beta, na.rm = F)
  ## max absolute distance
  abDist <- if (abs(min(simulation_results$Beta, na.rm = F) - meanBeta) > abs(max(simulation_results$Beta, na.rm = F) - meanBeta)) {
    abs(min(simulation_results$Beta, na.rm = F) - meanBeta) } else {
      abs(max(simulation_results$Beta, na.rm = F) - meanBeta) }
  ## confidence interval
  ### 90%
  CI90_a <- meanBeta - 1.645 * sdBeta
  CI90_b <- meanBeta + 1.645 * sdBeta
  ### 95%
  CI95_a <- meanBeta - 1.96 * sdBeta
  CI95_b <- meanBeta + 1.96 * sd(simulation_results$Beta, na.rm = F)
  ### 99%
  CI99_a <- meanBeta - 2.58 * sdBeta
  CI99_b <- meanBeta + 2.58 * sdBeta


  # create tibble object of results
  result_beta_inf_aux1 <- tribble(
    ~Name,                                 ~Value,
    "Beta* (mean)",                         round(meanBeta,6))
  result_beta_inf_aux2 <- tribble(
    ~Name,                                 ~Value,
    if(is.element(90, CI)) {"CI_90_low"},  round(CI90_a,6),
    if(is.element(90, CI)) {"CI_90_high"}, round(CI90_b,6))
  result_beta_inf_aux3 <- tribble(
    ~Name,                                 ~Value,
    if(is.element(95, CI)) {"CI_95_low"},  round(CI95_a,6),
    if(is.element(95, CI)) {"CI_95_high"}, round(CI95_b,6))
  result_beta_inf_aux4 <- tribble(
    ~Name,                                 ~Value,
    if(is.element(99, CI)) {"CI_99_low"},  round(CI99_a,6),
    if(is.element(99, CI)) {"CI_99_high"}, round(CI99_b,6))
  result_beta_inf_aux5 <- tribble(
    ~Name,                                 ~Value,
    "Simulations",                         sim,
    "Observations",                        obs,
    "Max R-square",                        R2max,
    "Delta (defined)",                     delta,
    paste("Model:",type),                  NA,
    paste("Replacement:",rep),     NA)
  if (is.element(90, CI)) {result_beta_inf_aux1 <- result_beta_inf_aux1 %>% add_row(result_beta_inf_aux2)}
  if (is.element(95, CI)) {result_beta_inf_aux1 <- result_beta_inf_aux1 %>% add_row(result_beta_inf_aux3)}
  if (is.element(99, CI)) {result_beta_inf_aux1 <- result_beta_inf_aux1 %>% add_row(result_beta_inf_aux4)}
  result_beta_inf <- result_beta_inf_aux1 %>% add_row(result_beta_inf_aux5)


  # warning if max R-square is smaller than the R-square of the control model
  if (R21 > R2max)
    warning("The max R-square value is smaller than the R-square of the controlled model")


  # warning if the bootstrapped observations are larger than the number of observations of the original dataset
  if (obs >= length(data$y_var))
    warning("Number of bootstrapped observation is larger than/equal to number of observation in the provided dataset")


  # define return object
  return(result_beta_inf)
}
