################################################################################
# R functions of the robomit R package                                         #
# -----------------------------------------------------------------------------#
# 2020, ETH Zurich                                                             #
# developed by Sergei Schaub                                                   #
# e-mail: seschaub#@ethz.ch                                                    #
# first version: August, 7, 2020                                               #
# last update: June, 20, 2021                                                  #
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

utils::globalVariables(c("Beta", "..density..", "Rmax","Delta","na.exclude","dnorm","delta","beta","Name","rmax","R21","Value","w_var","data_plm"))

# note that parts of the code is based (but not copied) from the psacalc command in Stata, which is licensed under GNU GENERAL PUBLIC LICENSE.

################################################################################
#------------------------------- 1. functions related to Oster (2019)
################################################################################


##################--------------------------------------------------------------
#------------------------------- 1.1 o_beta - function 1
##################--------------------------------------------------------------
#' @title beta*
#'
#' @description Estimates beta*, i.e., the bias-adjusted treatment effect (or correlation) (following Oster 2019).
#' @usage o_beta(y, x, con, m = "none", w = NULL, id = "none", time = "none", delta = 1,
#' R2max, type, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent treatment variable (i.e., variable of interest; as string).
#' @param con Name of related control variables. Provided as string in the format: "w + z +...".
#' @param m Name of unrelated control variables (m; see Oster 2019; as string; default is m = "none").
#' @param w weights (only for weighted estimations). Warning: For weighted panel models R can report different R-square than Stata, leading deviation between R and Stata results.
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect panel models.
#' @param time Name of the time id variable (e.g. year or month; as string). Only applicable for fixed effect panel models.
#' @param delta delta for which beta* should be estimated (default is delta = 1).
#' @param R2max Maximum R-square for which beta* should be estimated.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param data Dataset.
#' @details Estimates beta*, i.e., the bias-adjusted treatment effect (or correlation).
#' @return Returns tibble object, which includes beta* and various other information.
#' @references Oster, E. (2019) Unobservable Selection and Coefficient Stability: Theory and Evidence. Journal of Business & Economic Statistics, 37, 187-204.
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
#' o_beta(y = "mpg",           # dependent variable
#'        x = "wt",            # independent treatment variable
#'        con = "hp + qsec",   # related control variables
#'        delta = 1,           # delta
#'        R2max = 0.9,         # maximum R-square
#'        type = "lm",         # model type
#'        data = data_oster)   # dataset
#' @export


o_beta <- function(y, x, con, m = "none", w = NULL, id = "none", time = "none", delta = 1 ,R2max, type, data) {



  # rename variables and create, if type is plm, a panel dataset
  ## lm
  if (type == "lm") {
    if (is.null(w)) { data <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x))
    w_var <- w } else { # no weights are assigned
      data <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), w_var = all_of(w))}} # weights are assigned
  ## plm
  if (type == "plm") {
    if (is.null(w)) { data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), id_var = all_of(id), time_var = all_of(time))
    w_var <- w
    } else { # no weights are assigned
      data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), id_var = all_of(id), time_var = all_of(time), w_var = all_of(w))}
    data_plm <- pdata.frame(data_aux, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)}



  # define formulas of the different models
  ## models without unrelated variables (i.e. m; see Oster 2019)
  if(m == "none"){
    model0_formula    <- as.formula("y_var ~ x_var")                             # uncontrolled model
    model1_formula    <- as.formula(paste("y_var ~ x_var +",con))                # controlled model
    aux_model_formula <- as.formula(paste("x_var ~",con))                        # auxiliary model
    if (type == "plm") {
      sigma_xx_model_formula <- as.formula(paste("x_var ~ factor(id_var)"))}     # define model to compute sigma_xx when model type is plm
  } else {
    ## models with unrelated variables (i.e. m; see Oster 2019)
    model0_formula    <- as.formula(paste("y_var ~ x_var +", m))                 # uncontrolled model
    model1_formula    <- as.formula(paste("y_var ~ x_var +",con, "+", m))        # controlled model
    aux_model_formula <- as.formula(paste("x_var ~",con, "+", m))                # auxiliary model
    if (type == "lm") {
      sigma_xx_model_formula <- as.formula(paste("x_var ~ ",m))}                 # if lm model with m: model to obtain sigma_xx
    if (type == "plm") {
      sigma_xx_model_formula <- as.formula(paste("x_var ~ factor(id_var) +",m))} # if plm model with m: model to obtain sigma_xx
  }



  # run models
  ## lm
  if (type == "lm") {
    model0    <- lm(model0_formula, weights = w_var,    data = data, na.action = na.exclude)                          # uncontrolled model
    model1    <- lm(model1_formula, weights = w_var,    data = data, na.action = na.exclude)                          # controlled model
    aux_model <- lm(aux_model_formula, weights = w_var, data = data, na.action = na.exclude)                          # auxiliary model
    if(m != "none") {model_xx  <-  lm(sigma_xx_model_formula, weights = w_var, data = data, na.action = na.exclude)}  # model to obtain singa_xx
  }
  ## plm
  if (type == "plm") {
    model0    <- plm(model0_formula, weights = w_var,    data = data_plm, model = "within", na.action = na.exclude) # uncontrolled model
    model1    <- plm(model1_formula, weights = w_var,    data = data_plm, model = "within", na.action = na.exclude) # controlled model
    aux_model <- plm(aux_model_formula, weights = w_var, data = data_plm, model = "within", na.action = na.exclude) # auxiliary model
    model_xx  <-  lm(sigma_xx_model_formula, weights = w_var, data = data, na.action = na.exclude)                  # model to obtain singa_xx
  }



  # variables based on model outputs
  if (type == "lm") {b0 = as.numeric(tidy(model0)[2,2])} else if (type == "plm") {b0 = as.numeric(tidy(model0)[1,2])}                     # beta uncontrolled model
  if (type == "lm") {b1 = as.numeric(tidy(model1)[2,2])} else if (type == "plm") {b1 = as.numeric(tidy(model1)[1,2])}                     # beta controlled model
  if (type == "lm") {R20 = summary(model0)$r.squared} else if (type == "plm") {R20 = as.numeric(summary(model0)$r.squared[1])}            # R-square uncontrolled model
  if (type == "lm") {R21 = summary(model1)$r.squared} else if (type == "plm") {R21 = as.numeric(summary(model1)$r.squared[1])}            # R-square uncontrolled model
  sigma_yy = var(data$y_var, na.rm = T)                                                                                                   # variance of dependent variable
  if(m == "none"){
    if (type == "lm") {sigma_xx = var(data$x_var, na.rm = T)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)} } else {      # variance of independent variable without m
      if (type == "lm") {sigma_xx = var(model_xx$residuals)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)}  }             # variance of independent variable with m
  t_x  = var(aux_model$residuals)



  # create some additional variables
  rt_m_ro_t_syy = (R21-R20) * sigma_yy
  b0_m_b1 = b0 - b1
  rm_m_rt_t_syy = (R2max - R21) * sigma_yy



  if (delta == 1) { # if delta = 1


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
      "beta*",                         round(betax,6),
      "(beta*-beta controlled)^2",     round(distx,6),
      "Alternative Solution 1",        round(altsol1,6),
      "(beta[AS1]-beta controlled)^2", round(dist1,6),
      "Uncontrolled Coefficient",      b0,
      "Controlled Coefficient",        b1,
      "Uncontrolled R-square",         R20,
      "Controlled R-square",           R21,
      "Max R-square",                  R2max,
      "delta",                         delta
    )


    # warning if max R-square is smaller than the R-square of the control model
    if (R21 > R2max)
      warning("The max R-square value is smaller than the R-square of the controlled model")


    #define return object
    return(result_beta)
  }



  if (delta != 1) { # if delta != 0


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
        if ( sign(sols[i]-b1)!=sign(b1-b0) ) { dists[i]=max(dists)+1 }
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
      "beta*",                         round(betax,6),
      "(beta*-beta controlled)^2",     if (exists("distx")) {round(distx,6)} else {NA},
      "Alternative Solution 1",        if (exists("altsol1")) {round(altsol1,6)} else {NA},
      "(beta[AS1]-beta controlled)^2", if (exists("dist1")) {round(dist1,6)} else {NA},
      "Alternative Solution 2",        if (exists("altsol2")) {round(altsol2,6)} else {NA},
      "(beta[AS2]-beta controlled)^2", if (exists("dist2")) {round(dist2,6)} else {NA},
      "Uncontrolled Coefficient",      b0,
      "Controlled Coefficient",        b1,
      "Uncontrolled R-square",         R20,
      "Controlled R-square",           R21,
      "Max R-square",                  R2max,
      "delta",                         delta
    )


    # warning if max R-square is smaller than the R-square of the control model
    if (R21 > R2max)
      warning("The max R-square value is smaller than the R-square of the controlled model")


    # define return object
    return(result_beta)
  }

}



##################--------------------------------------------------------------
#------------------------------- 1.2 o_beta_rsq - function 2
##################--------------------------------------------------------------
#' @title beta*s over a range of maximum R-squares
#'
#' @description Estimates beta*s, i.e., the bias-adjusted treatment effects (or correlations) (following Oster 2019) over a range of maximum R-squares.
#' @usage o_beta_rsq(y, x, con, m = "none", w = NULL, id = "none", time = "none", delta = 1,
#' type, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent treatment variable (i.e., variable of interest; as string).
#' @param con Name of related control variables. Provided as string in the format: "w + z +...".
#' @param m Name of unrelated control variables (m; see Oster 2019; as string; default is m = "none").
#' @param w weights (only for weighted estimations). Warning: For weighted panel models R can report different R-square than Stata, leading deviation between R and Stata results.
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect panel models.
#' @param time Name of the time id variable (e.g. year or month; as string). Only applicable for fixed effect panel models.
#' @param delta delta for which beta*s should be estimated (default is delta = 1).
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param data Dataset.
#' @details Estimates beta*s, i.e., the bias-adjusted treatment effects (or correlations) (following Oster 2019) over a range of maximum R-squares. The range of maximum R-squares starts from the R-square of the controlled model rounded up to the next 1/100 to 1. The function supports linear cross-sectional (see \emph{lm} objects in R) and fixed effect panel (see \emph{plm} objects in R) models.
#' @return Returns tibble object, which includes beta*s over a range of maximum R-squares.
#' @references Oster, E. (2019). Unobservable Selection and Coefficient Stability: Theory and Evidence. Journal of Business & Economic Statistics, 37, 187-204.
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
#' # estimate delta*s over a range of maximum R-squares
#' o_beta_rsq(y = "mpg",            # dependent variable
#'            x = "wt",             # independent treatment variable
#'            con = "hp + qsec",    # related control variables
#'            delta = 1,            # delta
#'            type = "lm",          # model type
#'            data = data_oster)    # dataset
#' @export

o_beta_rsq <- function(y, x, con, m = "none", w = NULL, id = "none", time = "none", delta = 1, type, data) {



  # first define R-square of uncontrolled model
  ## rename variables and create, if type is plm, a panel dataset
  ### lm
  if (type == "lm") {
    if (is.null(w)) { data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x))
    w_var <- w } else { # no weights are assigned
      data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), w_var = all_of(w))}} # weights are assigned
  ### plm
  if (type == "plm") {
    if (is.null(w)) { data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), id_var = all_of(id), time_var = all_of(time))
    w_var <- w } else { # no weights are assigned
      data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), id_var = all_of(id), time_var = all_of(time), w_var = all_of(w))}
    data_plm_aux <- pdata.frame(data_aux, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)}

  ## define uncontrolled model
  ### models without variables m
  if(m == "none"){model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) } else { # controlled model
  ### models with
  model1_formula    <- as.formula(paste("y_var ~ x_var +",con, "+", m)) }                # controlled model

  ## run model
  ### lm
  if (type == "lm") { model1    <- lm(model1_formula, weights = w_var,    data = data_aux, na.action = na.exclude)}                                      # run controlled model
  ### plm
  if (type == "plm") {model1    <- plm(model1_formula, weights = w_var,    data = data_plm_aux, model = "within", na.action = na.exclude) }              # run controlled model

  ## get R-square
  if (type == "lm") {R21_aux = summary(model1)$r.squared} else if (type == "plm") {R21_aux = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model



  # define max R-square range
  r_start <- ceiling(R21_aux/0.01)*0.01      # start max R-square at the next highest 1/100
  n_runs <- length(seq(r_start:1,by = 0.01)) # number of different R-squares
  beta_over_rsq_results <- tibble(           # create tibble object to store results
    "x" = 1:n_runs,
    "Max R-square" = 0,
    "beta*" = 0)



  # run loop over different max R-squares
  for (i in 1:n_runs) {
    R2max = seq(r_start:1,by = 0.01)[i] # current max R-square

    # call o_delta function
    beta_results_aux <- o_beta(y = y, x = x, con = con, m = m, id = id, time = time, delta = delta, R2max = R2max, type = type, data = data)

    # store results
    beta_over_rsq_results[i,2] <- R2max
    beta_over_rsq_results[i,3] <- round(beta_results_aux[1,2],6)
  }


  # define return object
  return(beta_over_rsq_results)
}



##################--------------------------------------------------------------
#------------------------------- 1.3 o_beta_rsq_viz - function 3
##################--------------------------------------------------------------
#' @title Visualization of beta*s over a range of maximum R-squares
#'
#' @description Estimates and visualizes beta*s, i.e., the bias-adjusted treatment effects (or correlations) (following Oster 2019) over a range of maximum R-squares.
#' @usage o_beta_rsq_viz(y, x, con, m = "none", w = NULL, id = "none", time = "none", delta = 1,
#' type, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent treatment variable (i.e., variable of interest; as string).
#' @param con Name of related control variables. Provided as string in the format: "w + z +...".
#' @param m Name of unrelated control variables (m; see Oster 2019; as string; default is m = "none").
#' @param w weights (only for weighted estimations). Warning: For weighted panel models R can report different R-square than Stata, leading deviation between R and Stata results.
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect panel models.
#' @param time Name of the time id variable (e.g. year or month; as string). Only applicable for fixed effect panel models.
#' @param delta delta for which beta*s should be estimated (default is delta = 1).
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param data Dataset.
#' @details Estimates and visualizes beta*s, i.e., the bias-adjusted treatment effects (or correlations) (following Oster 2019) over a range of maximum R-squares. The range of maximum R-squares starts from the R-square of the controlled model rounded up to the next 1/100 to 1. The function supports linear cross-sectional (see \emph{lm} objects in R) and fixed effect panel (see \emph{plm} objects in R) models.
#' @return Returns ggplot2 object, which depicts beta*s over a range of maximum R-squares.
#' @references Oster, E. (2019). Unobservable Selection and Coefficient Stability: Theory and Evidence. Journal of Business & Economic Statistics, 37, 187-204.
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
#' # estimate and visualize beta*s over a range of maximum R-squares
#' o_beta_rsq_viz(y = "mpg",            # dependent variable
#'                x = "wt",             # independent treatment variable
#'                con = "hp + qsec",    # related control variables
#'                delta = 1,            # delta
#'                type = "lm",          # model type
#'                data = data_oster)    # dataset
#' @export


o_beta_rsq_viz <- function(y, x, con, m = "none", w = NULL, id = "none", time = "none", delta = 1, type, data) {



  # call o_beta_rsq
  beta_over_rsq_results <- o_beta_rsq(y = y, x = x, con = con, m = m, id = id, time = time, delta = delta, type = type, data = data)



  # rename variables
  beta_over_rsq_results <- beta_over_rsq_results %>% dplyr::rename(Rmax  = 2, beta = 3)



  # plot figure
  theme_set(theme_bw()) # set main design
  result_plot <- ggplot(data=beta_over_rsq_results, aes(x=Rmax, y=beta)) +
    geom_line(size = 1.3)+
    scale_y_continuous(name = expression(beta^"*"))+
    scale_x_continuous(name = expression("maximum R"^2))+
    theme(axis.title = element_text( size=15),
          axis.text  = element_text( size=13))



  # define return object
  return(result_plot)
}



##################--------------------------------------------------------------
#------------------------------- 1.4 o_beta_boot - function 4
##################--------------------------------------------------------------
#' @title Bootstrapped beta*s
#'
#' @description Estimates bootstrapped beta*s, i.e., the bias-adjusted treatment effects (or correlations) (following Oster 2019).
#' @usage o_beta_boot(y, x, con, m = "none", w = NULL, id = "none", time = "none", delta = 1,
#' R2max, sim, obs, rep, type, useed = NA, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent treatment variable (i.e., variable of interest; as string).
#' @param con Name of related control variables. Provided as string in the format: "w + z +...".
#' @param m Name of unrelated control variables (m; see Oster 2019; as string; default is m = "none").
#' @param w weights (only for weighted estimations). Warning: For weighted panel models R can report different R-square than Stata, leading deviation between R and Stata results.
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect panel models.
#' @param time Name of the time id variable (e.g. year or month; as string). Only applicable for fixed effect panel models.
#' @param delta delta for which beta*s should be estimated (default is delta = 1).
#' @param R2max Maximum R-square for which beta*s should be estimated.
#' @param sim Number of simulations.
#' @param obs Number of draws per simulation.
#' @param rep Bootstrapping either with (= TRUE) or without (= FALSE) replacement.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param useed User defined seed.
#' @param data Dataset.
#' @details Estimates bootstrapped beta*s, i.e., the bias-adjusted treatment effects (or correlations) (following Oster 2019). Bootstrapping can either be done with or without replacement. The function supports linear cross-sectional (see \emph{lm} objects in R) and fixed effect panel (see \emph{plm} objects in R) models.
#' @return Returns tibble object, which includes bootstrapped beta*s.
#' @references Oster, E. (2019). Unobservable Selection and Coefficient Stability: Theory and Evidence. Journal of Business & Economic Statistics, 37, 187-204.
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
#' # estimate bootstrapped beta*s
#' o_beta_boot(y = "mpg",            # dependent variable
#'             x = "wt",             # independent treatment variable
#'             con = "hp + qsec",    # related control variables
#'             delta = 1,            # delta
#'             R2max = 0.9,          # maximum R-square
#'             sim = 100,            # number of simulations
#'             obs = 30,             # draws per simulation
#'             rep = FALSE,          # bootstrapping with or without replacement
#'             type = "lm",          # model type
#'             useed = 123,          # seed
#'             data = data_oster)    # dataset
#' @export



o_beta_boot <- function(y, x, con, m = "none", w = NULL, id = "none", time = "none", delta = 1, R2max, sim, obs, rep, type, useed = NA, data) {




  # first define R-square of uncontrolled model
  ## rename variables and create, if type is plm, a panel dataset
  ### lm
  if (type == "lm") {
    if (is.null(w)) { data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x))
    w_var <- w } else { # no weights are assigned
      data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), w_var = all_of(w))}} # weights are assigned
  ### plm
  if (type == "plm") {
    if (is.null(w)) {
      data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), id_var = all_of(id), time_var = all_of(time))
      w_var <- w } else { # no weights are assigned
      data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), id_var = all_of(id), time_var = all_of(time), w_var = all_of(w))}
    data_plm_aux <- pdata.frame(data_aux, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)}

  ## define uncontrolled model
  ### models without variables m
  if(m == "none"){model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) } else { # controlled model
  ### models with
    model1_formula    <- as.formula(paste("y_var ~ x_var +",con, "+", m)) }              # controlled model

  ## run model
  ### lm
  if (type == "lm") { model1    <- lm(model1_formula, weights = w_var,    data = data_aux, na.action = na.exclude)}                                      # run controlled model
  ### plm
  if (type == "plm") {model1    <- plm(model1_formula, weights = w_var,    data = data_plm_aux, model = "within", na.action = na.exclude) }              # run controlled model
  ## get R-square
  if (type == "lm") {R21_aux = summary(model1)$r.squared} else if (type == "plm") {R21_aux = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model


  # keep original data
  data_original <- data


  # create tibble object to store results
  simulation_results <- tibble(
    "x" = 1:sim,
    "beta*" = 0)


  for (s in 1:sim) {

    if (!is.na(useed)) {set.seed(useed+s) } # set seed (defined by user) that it varies with runs

    # draw data
    if (rep) {
      data_current <- sample_n(data_original, size = obs, replace = TRUE)  # with replacement
    } else {
      data_current <- sample_n(data_original, size = obs, replace = FALSE) # without replacement
    }


    # call o_delta function
    beta_results_aux <- o_beta(y = y, x = x, con = con, m = m, id = id, time = time,
                               delta = delta, R2max = R2max, type = type, data = data_current)

    # store results
    simulation_results[s,2] <- round(beta_results_aux[1,2],6)
  }


  # warning if max R-square is smaller than the R-square of the control model
  if (R21_aux > R2max)
    warning("The max R-square value is smaller than the R-square of the controlled model")


  # warning if the bootstrapped observations are larger than the number of observations of the original dataset
  data_aux2 <- data %>% dplyr::rename(y_var = all_of(y))
  if (obs >= length(data_aux2$y_var))
    warning("Number of bootstrapped observation is larger than/equal to number of observation in the provided dataset")


  # define return object
  return(simulation_results)
}



##################--------------------------------------------------------------
#------------------------------- 1.5 o_beta_boot_viz - function 5
##################--------------------------------------------------------------
#' @title Visualization of bootstrapped beta*s
#'
#' @description Estimates and visualizes bootstrapped beta*s, i.e., the bias-adjusted treatment effects (or correlations) (following Oster 2019).
#' @usage o_beta_boot_viz(y, x, con, m = "none", w = NULL, id = "none", time = "none",
#' delta = 1, R2max, sim, obs, rep, CI, type, norm = TRUE, bin,
#' col = c("#08306b","#4292c6","#c6dbef"), nL = TRUE, mL = TRUE, useed = NA, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent treatment variable (i.e., variable of interest; as string).
#' @param con Name of related control variables. Provided as string in the format: "w + z +...".
#' @param m Name of unrelated control variables (m; see Oster 2019; as string; default is m = "none").
#' @param w weights (only for weighted estimations). Warning: For weighted panel models R can report different R-square than Stata, leading deviation between R and Stata results.
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect panel models.
#' @param time Name of the time id variable (e.g. year or month; as string). Only applicable for fixed effect panel models.
#' @param delta delta for which beta*s should be estimated (default is delta = 1).
#' @param R2max Maximum R-square for which beta*s should be estimated.
#' @param sim Number of simulations.
#' @param obs Number of draws per simulation.
#' @param rep Bootstrapping either with (= TRUE) or without (= FALSE) replacement
#' @param CI  Confidence intervals, indicated as vector. Can be and/or 90, 95, 99.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param norm Option to include a normal distribution in the plot (default is norm = TURE).
#' @param bin Number of bins used in the histogram.
#' @param col Colors used to indicate different confidence interval levels (indicated as vector). Needs to be the same length as the variable CI. The default is a blue color range.
#' @param nL Option to include a red vertical line at 0 (default is nL = TRUE).
#' @param mL Option to include a vertical line at mean of all beta*s (default is mL = TRUE).
#' @param useed User defined seed.
#' @param data Dataset.
#' @details Estimates and visualizes bootstrapped beta*s, i.e., the bias-adjusted treatment effects (or correlations) (following Oster 2019). Bootstrapping can either be done with or without replacement. The function supports linear cross-sectional (see \emph{lm} objects in R) and fixed effect panel (see \emph{plm} objects in R) models.
#' @return Returns ggplot2 object, which depicts the bootstrapped beta*s.
#' @references Oster, E. (2019). Unobservable Selection and Coefficient Stability: Theory and Evidence. Journal of Business & Economic Statistics, 37, 187-204.
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
#' # estimate and visualize bootstrapped beta*s
#' o_beta_boot_viz(y = "mpg",            # dependent variable
#'                 x = "wt",             # independent treatment variable
#'                 con = "hp + qsec",    # related control variables
#'                 delta = 1,            # delta
#'                 R2max = 0.9,          # maximum R-square
#'                 sim = 100,            # number of simulations
#'                 obs = 30,             # draws per simulation
#'                 rep = FALSE,          # bootstrapping with or without replacement
#'                 CI = c(90,95,99),     # confidence intervals
#'                 type = "lm",          # model type
#'                 norm = TRUE,          # normal distribution
#'                 bin = 200,            # number of bins
#'                 useed = 123,          # seed
#'                 data = data_oster)    # dataset
#' @export


o_beta_boot_viz <- function(y, x, con, m = "none", w = NULL, id = "none", time = "none", delta = 1, R2max, sim, obs, rep, CI, type, norm = TRUE, bin, col = c("#08306b","#4292c6","#c6dbef"), nL = TRUE, mL = TRUE, useed = NA, data) {


  # call o_beta_boot function and store results
  simulation_results <- o_beta_boot(y = y, x = x, con = con, m = m, id = id, time = time, delta = delta, R2max = R2max,
                                    sim = sim, obs = obs, rep = rep, type = type, useed = useed, data = data)
  simulation_results <- simulation_results %>% dplyr::rename("beta" = 2) # rename variable



  # compute mean, max absolute distance to mean, and confidence intervals
  ## mean
  meanBeta <- mean(simulation_results$beta, na.rm = F)
  # sd
  sdBeta <- sd(simulation_results$beta, na.rm = F)
  ## confidence interval
  ### 90%
  CI90_a <- meanBeta - 1.645 * sdBeta
  CI90_b <- meanBeta + 1.645 * sdBeta
  ### 95%
  CI95_a <- meanBeta - 1.96 * sdBeta
  CI95_b <- meanBeta + 1.96 * sdBeta
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
  beta_plot <- ggplot(simulation_results, aes(x = beta)) +
    {if(is.element(99, CI))annotate("rect", ymin = 0, ymax = +Inf, xmin = CI99_a, xmax = CI99_b, fill = color99, alpha = 0.6)}+
    {if(is.element(95, CI))annotate("rect", ymin = 0, ymax = +Inf, xmin = CI95_a, xmax = CI95_b, fill = color95, alpha = 0.6)}+
    {if(is.element(90, CI))annotate("rect", ymin = 0, ymax = +Inf, xmin = CI90_a, xmax = CI90_b, fill = color90, alpha = 0.6)}+
    geom_histogram(aes(y = ..density..), alpha = 0.8, colour = "#737373", fill = "#737373", size = 0.1, bins = bin)+
    {if(norm)stat_function(fun = dnorm, args = list(mean = meanBeta, sd = sdBeta), size = 1.3)}+
    scale_y_continuous(name = "Frequency",expand = expansion(mult = c(0, .1)))+
    scale_x_continuous(name = expression(beta^"*"))+
    expand_limits(x = 0, y = 0)+
    theme(axis.title = element_text( size=15),
          axis.text  = element_text( size=13))+
    {if(mL)geom_vline(xintercept = meanBeta, size = 0.8, color = "#000000", alpha = 0.8)}+
    {if(nL)geom_vline(xintercept = 0, size = 0.8, color = "#e41a1c", alpha = 0.8)}+
    geom_hline(yintercept = 0)


  # define return object
  return(beta_plot)


}



##################--------------------------------------------------------------
#------------------------------- 1.6 o_beta_boot_inf - function 6
##################--------------------------------------------------------------
#' @title Bootstrapped mean beta* and confidence intervals
#'
#' @description Provides the mean and confidence intervals of estimated bootstrapped beta*s, i.e., the bias-adjusted treatment effects (or correlations) (following Oster 2019).
#' @usage o_beta_boot_inf(y, x, con, m = "none", w = NULL, id = "none", time = "none",
#' delta = 1, R2max, sim, obs, rep, CI, type, useed = NA, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent treatment variable (i.e., variable of interest; as string).
#' @param con Name of related control variables. Provided as string in the format: "w + z +...".
#' @param m Name of unrelated control variables (m; see Oster 2019; as string; default is m = "none").
#' @param w weights (only for weighted estimations). Warning: For weighted panel models R can report different R-square than Stata, leading deviation between R and Stata results.
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect panel models.
#' @param time Name of the time id variable (e.g. year or month; as string). Only applicable for fixed effect panel models.
#' @param delta delta for which beta*s should be estimated (default is delta = 1).
#' @param R2max Maximum R-square for which beta*s should be estimated.
#' @param sim Number of simulations.
#' @param obs Number of draws per simulation.
#' @param rep Bootstrapping either with (= TRUE) or without (= FALSE) replacement
#' @param CI  Confidence intervals, indicated as vector. Can be and/or 90, 95, 99.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param useed User defined seed.
#' @param data Dataset.
#' @details Provides the mean and confidence intervals of estimated bootstrapped beta*s, i.e., the bias-adjusted treatment effects (or correlations) (following Oster 2019). Bootstrapping can either be done with or without replacement. The function supports linear cross-sectional (see \emph{lm} objects in R) and fixed effect panel (see \emph{plm} objects in R) models.
#' @return Returns tibble object, which includes the mean and confidence intervals of estimated bootstrapped beta*s.
#' @references Oster, E. (2019). Unobservable Selection and Coefficient Stability: Theory and Evidence. Journal of Business & Economic Statistics, 37, 187-204.
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
#' # compute the mean and confidence intervals of estimated bootstrapped beta*s
#' o_beta_boot_inf(y = "mpg",            # dependent variable
#'                 x = "wt",             # independent treatment variable
#'                 con = "hp + qsec",    # related control variables
#'                 delta = 1,            # delta
#'                 R2max = 0.9,          # maximum R-square
#'                 sim = 100,            # number of simulations
#'                 obs = 30,             # draws per simulation
#'                 rep = FALSE,          # bootstrapping with or without replacement
#'                 CI = c(90,95,99),     # confidence intervals
#'                 type = "lm",          # model type
#'                 useed = 123,          # seed
#'                 data = data_oster)    # dataset
#' @export


o_beta_boot_inf <- function(y, x, con, m = "none", w = NULL, id = "none", time = "none", delta = 1, R2max, sim, obs, rep, CI, type, useed = NA, data) {


  # call o_beta_boot function and store results
  simulation_results <- o_beta_boot(y = y, x = x, con = con, m = m, id = id, time = time, delta = delta, R2max = R2max,
                                    sim = sim, obs = obs, rep = rep, type = type, useed = useed, data = data)
  simulation_results <- simulation_results %>% dplyr::rename("beta" = 2) # rename variable



  # compute mean, max absolute distance to mean, and confidence intervals
  ## mean
  meanBeta <- mean(simulation_results$beta, na.rm = F)
  # sd
  sdBeta <- sd(simulation_results$beta, na.rm = F)
  ## confidence interval
  ### 90%
  CI90_a <- meanBeta - 1.645 * sdBeta
  CI90_b <- meanBeta + 1.645 * sdBeta
  ### 95%
  CI95_a <- meanBeta - 1.96 * sdBeta
  CI95_b <- meanBeta + 1.96 * sdBeta
  ### 99%
  CI99_a <- meanBeta - 2.58 * sdBeta
  CI99_b <- meanBeta + 2.58 * sdBeta

  # create tibble object of results
  result_beta_inf_aux1 <- tribble(
    ~Name,                                 ~Value,
    "beta* (mean)",                         round(meanBeta,6))
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
    "delta",                               delta,
    "Model:",                              NA,
    "Replacement:",                        NA)
  if (is.element(90, CI)) {result_beta_inf_aux1 <- result_beta_inf_aux1 %>% add_row(result_beta_inf_aux2)}
  if (is.element(95, CI)) {result_beta_inf_aux1 <- result_beta_inf_aux1 %>% add_row(result_beta_inf_aux3)}
  if (is.element(99, CI)) {result_beta_inf_aux1 <- result_beta_inf_aux1 %>% add_row(result_beta_inf_aux4)}
  result_beta_inf <- result_beta_inf_aux1 %>% add_row(result_beta_inf_aux5)%>%
    mutate(Value = as.character(Value),
           Value = ifelse(Name == "Model:", type,
                          ifelse(Name == "Replacement:", rep, Value)))


  # define return object
  return(result_beta_inf)
}


##################--------------------------------------------------------------
#------------------------------- 1.7 o_delta - function 7
##################--------------------------------------------------------------
#' @title delta*
#'
#' @description Estimates delta*, i.e., the degree of selection on unobservables relative to observables (with respect to the treatment variable) that would be necessary to eliminate the result (following Oster 2019).
#' @usage o_delta(y, x, con, m = "none", w = NULL, id = "none", time = "none", beta = 0, R2max,
#' type, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent treatment variable (i.e., variable of interest; as string).
#' @param con Name of related control variables. Provided as string in the format: "w + z +...".
#' @param m Name of unrelated control variables (m; see Oster 2019; as string; default is m = "none").
#' @param w weights (only for weighted estimations). Warning: For weighted panel models R can report different R-square than Stata, leading deviation between R and Stata results.
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect panel models.
#' @param time Name of the time id variable (e.g. year or month; as string). Only applicable for fixed effect panel models.
#' @param beta beta for which delta* should be estimated (default is beta = 0).
#' @param R2max Maximum R-square for which delta* should be estimated.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param data Dataset.
#' @details Estimates delta*, i.e., the degree of selection on unobservables relative to observables (with respect to the treatment variable) that would be necessary to eliminate the result (following Oster 2019). The function supports linear cross-sectional (see \emph{lm} objects in R) and fixed effect panel (see \emph{plm} objects in R) models.
#' @return Returns tibble object, which includes delta* and various other information.
#' @references Oster, E. (2019). Unobservable Selection and Coefficient Stability: Theory and Evidence. Journal of Business & Economic Statistics, 37, 187-204.
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
#' o_delta(y = "mpg",           # dependent variable
#'         x = "wt",            # independent treatment variable
#'         con = "hp + qsec",   # related control variables
#'         beta = 0,            # beta
#'         R2max = 0.9,         # maximum R-square
#'         type = "lm",         # model type
#'         data = data_oster)   # dataset
#' @export

o_delta <- function(y, x, con, m = "none", w = NULL, id = "none", time = "none", beta = 0, R2max, type, data) {


  # rename variables and create, if type is plm, a panel dataset
  ## lm
  if (type == "lm") {
    if (is.null(w)) { data <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x))
    w_var <- w } else { # no weights are assigned
      data <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), w_var = all_of(w))}} # weights are assigned
  ## plm
  if (type == "plm") {
    if (is.null(w)) { data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), id_var = all_of(id), time_var = all_of(time))
    w_var <- w } else { # no weights are assigned
      data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), id_var = all_of(id), time_var = all_of(time), w_var = all_of(w))}
    data_plm <- pdata.frame(data_aux, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)}


  # define formulas of the different models
  ## models without variables m
  if(m == "none"){
    model0_formula    <- as.formula("y_var ~ x_var")                                            # uncontrolled model
    model1_formula    <- as.formula(paste("y_var ~ x_var +",con))                               # controlled model
    aux_model_formula <- as.formula(paste("x_var ~",con))                                       # auxiliary model
    if (type == "plm") { sigma_xx_model_formula <- as.formula(paste("x_var ~ factor(id_var)"))} # define model to compute sigma_xx when model type is plm
  } else {
    ## models with variables m
    model0_formula    <- as.formula(paste("y_var ~ x_var +", m))                 # uncontrolled model
    model1_formula    <- as.formula(paste("y_var ~ x_var +",con, "+", m))        # controlled model
    aux_model_formula <- as.formula(paste("x_var ~",con, "+", m))                # auxiliary model
    if (type == "lm") {
      sigma_xx_model_formula <- as.formula(paste("x_var ~ ",m))}                 # if lm model: model to obtain sigma_xx with m
    if (type == "plm") {
      sigma_xx_model_formula <- as.formula(paste("x_var ~ factor(id_var) +",m))} # if plm model: model to obtain sigma_xx with m
  }



  # run models
  ## lm
  if (type == "lm") {
    model0    <- lm(model0_formula, weights = w_var,    data = data, na.action = na.exclude)                          # run uncontrolled model
    model1    <- lm(model1_formula, weights = w_var,    data = data, na.action = na.exclude)                          # run controlled model
    aux_model <- lm(aux_model_formula, weights = w_var, data = data, na.action = na.exclude)                          # run auxiliary model
    if(m != "none") {model_xx  <-  lm(sigma_xx_model_formula, weights = w_var, data = data, na.action = na.exclude)}  # run model to obtain singa_xx
  }
  ## plm
  if (type == "plm") {
    model0    <- plm(model0_formula, weights = w_var,    data = data_plm, model = "within", na.action = na.exclude) # run uncontrolled model
    model1    <- plm(model1_formula, weights = w_var,    data = data_plm, model = "within", na.action = na.exclude) # run controlled model
    aux_model <- plm(aux_model_formula, weights = w_var, data = data_plm, model = "within", na.action = na.exclude) # run auxiliary model
    model_xx  <-  lm(sigma_xx_model_formula, weights = w_var, data = data, na.action = na.exclude)                  # run model to obtain singa_xx
  }


  # variables based on model outputs
  if (type == "lm") {b0 = as.numeric(tidy(model0)[2,2])} else if (type == "plm") {b0 = as.numeric(tidy(model0)[1,2])}           # beta uncontrolled model
  if (type == "lm") {b1 = as.numeric(tidy(model1)[2,2])} else if (type == "plm") {b1 = as.numeric(tidy(model1)[1,2])}           # beta controlled model
  if (type == "lm") {R20 = summary(model0)$r.squared} else if (type == "plm") {R20 = as.numeric(summary(model0)$r.squared[1])}  # R-square uncontrolled model
  if (type == "lm") {R21 = summary(model1)$r.squared} else if (type == "plm") {R21 = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model
  sigma_yy = var(data$y_var, na.rm = T)                                                                                         # variance of dependent variable
  if(m == "none"){
    if (type == "lm") {sigma_xx = var(data$x_var, na.rm = T)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)} } else {      # variance of independent variable
      if (type == "lm") {sigma_xx = var(model_xx$residuals)} else if (type == "plm") {sigma_xx =  var(model_xx$residuals)}  }             # variance of independent variable
  t_x  = var(aux_model$residuals)


  # create some additional variables
  bt_m_b = b1 - beta
  rt_m_ro_t_syy = (R21-R20) * sigma_yy
  b0_m_b1 = b0 - b1
  rm_m_rt_t_syy = (R2max - R21) * sigma_yy


  # compute numerator
  num1 = bt_m_b * rt_m_ro_t_syy * t_x
  num2 = bt_m_b * sigma_xx * t_x * b0_m_b1^2
  num3 = 2 * bt_m_b^2 * (t_x * b0_m_b1 * sigma_xx)
  num4 = bt_m_b^3 * (t_x * sigma_xx - t_x^2)
  num  = num1 + num2 + num3 + num4


  # compute denominator
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
    "delta*",                   round(delta_star,6),
    "Uncontrolled Coefficient", b0,
    "Controlled Coefficient",   b1,
    "Uncontrolled R-square",    R20,
    "Controlled R-square",      R21,
    "Max R-square",             R2max,
    "beta hat",                 beta
  )

  # warning if max R-square is smaller than the R-square of the control model
  if (R21 > R2max)
    warning("The max R-square value is smaller than the R-square of the controlled model")


  # define return object
  return(result_delta)
}



##################--------------------------------------------------------------
#------------------------------- 1.8 o_delta_rsq - function 8
##################--------------------------------------------------------------
#' @title delta*s over a range of maximum R-squares
#' @description Estimates delta*s, i.e., the degree of selection on unobservables relative to observables (with respect to the treatment variable) that would be necessary to eliminate the result (following Oster 2019) over a range of maximum R-squares following Oster (2019).
#' @usage o_delta_rsq(y, x, con, m = "none", w = NULL, id = "none", time = "none", beta = 0,
#' type, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent treatment variable (i.e., variable of interest; as string).
#' @param con Name of related control variables. Provided as string in the format: "w + z +...".
#' @param m Name of unrelated control variables (m; see Oster 2019; as string; default is m = "none").
#' @param w weights (only for weighted estimations). Warning: For weighted panel models R can report different R-square than Stata, leading deviation between R and Stata results.
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect panel models.
#' @param time Name of the time id variable (e.g. year or month; as string). Only applicable for fixed effect panel models.
#' @param beta beta for which delta*s should be estimated (default is beta = 0).
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param data Dataset.
#' @details Estimates delta*s, i.e., the degree of selection on unobservables relative to observables (with respect to the treatment variable) that would be necessary to eliminate the result (following Oster 2019) over a range of maximum R-squares. The range of maximum R-squares starts from the R-square of the controlled model rounded up to the next 1/100 to 1. The function supports linear cross-sectional (see \emph{lm} objects in R) and fixed effect panel (see \emph{plm} objects in R) models.
#' @return Returns tibble object, which includes delta*s over a range of maximum R-squares.
#' @references Oster, E. (2019). Unobservable Selection and Coefficient Stability: Theory and Evidence. Journal of Business & Economic Statistics, 37, 187-204.
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
#' # estimate delta*s over a range of maximum R-squares
#' o_delta_rsq(y = "mpg",           # dependent variable
#'             x = "wt",            # independent treatment variable
#'             con = "hp + qsec",   # related control variables
#'             beta = 0,            # beta
#'             type = "lm",         # model type
#'             data = data_oster)   # dataset
#' @export

o_delta_rsq <- function(y, x, con, m = "none", w = NULL, id = "none", time = "none", beta = 0, type, data) {


  # first define R-square of uncontrolled model
  ## rename variables and create, if type is plm, a panel dataset
  ### lm
  if (type == "lm") {
    if (is.null(w)) { data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x))
    w_var <- w } else { # no weights are assigned
      data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), w_var = all_of(w))}} # weights are assigned
  ### plm
  if (type == "plm") {
    if (is.null(w)) { data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), id_var = all_of(id), time_var = all_of(time))
    w_var <- w } else { # no weights are assigned
      data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), id_var = all_of(id), time_var = all_of(time), w_var = all_of(w))}
    data_plm_aux <- pdata.frame(data_aux, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)}

  ## define uncontrolled model
  ### models without variables m
  if(m == "none"){model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) } else { # controlled model
  ### models with variables m
    model1_formula    <- as.formula(paste("y_var ~ x_var +",con, "+", m)) }              # controlled model

  ## run model
  ### lm
  if (type == "lm") { model1    <- lm(model1_formula, weights = w_var,    data = data_aux, na.action = na.exclude)}                                      # run controlled model
  ### plm
  if (type == "plm") {model1    <- plm(model1_formula, weights = w_var,    data = data_plm_aux, model = "within", na.action = na.exclude) }              # run controlled model
  ## get R-square
  if (type == "lm") {R21_aux = summary(model1)$r.squared} else if (type == "plm") {R21_aux = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model

  # define max R-square range
  r_start <- ceiling(R21_aux/0.01)*0.01       # start max R-square at the next highest 1/100
  n_runs <- length(seq(r_start:1,by = 0.01))  # number of different R-squares
  delta_over_rsq_results <- tibble(           # create tibble object to store results
    "x" = 1:n_runs,
    "Max R-square" = 0,
    "delta*" = 0)


  # run loop over different max R-squares
  for (i in 1:n_runs) {
    R2max = seq(r_start:1,by = 0.01)[i] # current max R-square

    # call o_delta function and store results
    delta_results_aux <- o_delta(y = y, x = x, con = con, m = m, id = id, time = time, beta = beta, R2max = R2max, type = type, data = data)

    # store results
    delta_over_rsq_results[i,2] <- R2max
    delta_over_rsq_results[i,3] <- round(delta_results_aux[1,2],6)
  }


  # define return object
  return(delta_over_rsq_results)
}




##################--------------------------------------------------------------
#------------------------------- 1.9 o_delta_rsq_viz - function 9
##################--------------------------------------------------------------
#' @title Visualization of delta*s over a range of maximum R-squares
#'
#' @description Estimates and visualizes delta*s, i.e., the degree of selection on unobservables relative to observables (with respect to the treatment variable) that would be necessary to eliminate the result (following Oster 2019) over a range of maximum R-squares.
#' @usage o_delta_rsq_viz(y, x, con, m = "none", w = NULL, id = "none", time = "none", beta = 0,
#' type, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent treatment variable (i.e., variable of interest; as string).
#' @param con Name of related control variables. Provided as string in the format: "w + z +...".
#' @param m Name of unrelated control variables (m; see Oster 2019; as string; default is m = "none").
#' @param w weights (only for weighted estimations). Warning: For weighted panel models R can report different R-square than Stata, leading deviation between R and Stata results.
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect panel models.
#' @param time Name of the time id variable (e.g. year or month; as string). Only applicable for fixed effect panel models.
#' @param beta beta for which delta*s should be estimated (default is beta = 0).
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param data Dataset.
#' @details Estimates and visualizes delta*s, i.e., the degree of selection on unobservables relative to observables (with respect to the treatment variable) that would be necessary to eliminate the result (following Oster 2019) over a range of maximum R-squares. The range of maximum R-squares starts from the R-square of the controlled model rounded up to the next 1/100 to 1. The function supports linear cross-sectional (see \emph{lm} objects in R) and fixed effect panel (see \emph{plm} objects in R) models.
#' @return Returns ggplot2 object, which depicts delta*s over a range of maximum R-squares.
#' @references Oster, E. (2019). Unobservable Selection and Coefficient Stability: Theory and Evidence. Journal of Business & Economic Statistics, 37, 187-204.
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
#' # estimate and visualize delta*s over a range of maximum R-squares
#' o_delta_rsq_viz(y = "mpg",           # dependent variable
#'                 x = "wt",            # independent treatment variable
#'                 con = "hp + qsec",   # related control variables
#'                 beta = 0,            # beta
#'                 type = "lm",         # model type
#'                 data = data_oster)   # dataset
#' @export

o_delta_rsq_viz <- function(y, x, con, m = "none", w = NULL, id = "none", time = "none", beta = 0, type, data) {


  # first define R-square of uncontrolled model
  ## rename variables and create, if type is plm, a panel dataset
  ### lm
  if (type == "lm") {
    if (is.null(w)) { data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x))
    w_var <- w } else { # no weights are assigned
      data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), w_var = all_of(w))}} # weights are assigned
  ### plm
  if (type == "plm") {
    if (is.null(w)) { data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), id_var = all_of(id), time_var = all_of(time))
    w_var <- w} else { # no weights are assigned
      data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), id_var = all_of(id), time_var = all_of(time), w_var = all_of(w))}
    data_plm_aux <- pdata.frame(data_aux, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)}

  ## define uncontrolled model
  ### models without variables m
  if(m == "none"){model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) } else { # controlled model
  ### models with variables m
    model1_formula    <- as.formula(paste("y_var ~ x_var +",con, "+", m)) }              # controlled model

  ## run model
  ### lm
  if (type == "lm") { model1    <- lm(model1_formula, weights = w_var,    data = data_aux, na.action = na.exclude)}                                      # run controlled model
  ### plm
  if (type == "plm") {model1    <- plm(model1_formula, weights = w_var,    data = data_plm_aux, model = "within", na.action = na.exclude) }              # run controlled model

  ## get R-square
  if (type == "lm") {R21_aux = summary(model1)$r.squared} else if (type == "plm") {R21_aux = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model



  # define max R-square range
  r_start <- ceiling(R21_aux/0.01)*0.01      # start max R-square at the next highest 1/100
  n_runs <- length(seq(r_start:1,by = 0.01)) # number of different R-squares
  delta_over_rsq_results <- tibble(          # create tibble object to store results
    x = 1:n_runs,
    rmax = 0,
    delta = 0)


  # run loop over different max R-squares
  for (i in 1:n_runs) {
    R2max = seq(r_start:1,by = 0.01)[i] # current max R-square

    # call o_delta function
    delta_results_aux <- o_delta(y = y, x = x, con = con, m = m, id = id, time = time, beta = beta, R2max = R2max, type = type, data = data)

    # store results
    delta_over_rsq_results[i,2] <- R2max
    delta_over_rsq_results[i,3] <- round(delta_results_aux[1,2],6)
  }


  # plot figure
  theme_set(theme_bw()) # set main design
  result_plot <- ggplot(data=delta_over_rsq_results, aes(x=rmax, y=delta)) +
    geom_line(size = 1.3)+
    scale_y_continuous(name = expression(delta^"*"))+
    scale_x_continuous(name = expression("maximum R"^2))+
    theme(axis.title = element_text( size=15),
          axis.text  = element_text( size=13))


  # define return object
  return(result_plot)
}



##################--------------------------------------------------------------
#------------------------------- 1.10 o_delta_boot - function 10
##################--------------------------------------------------------------
#' @title Bootstrapped delta*s
#'
#' @description Estimates bootstrapped delta*s, i.e., the degree of selection on unobservables relative to observables (with respect to the treatment variable) that would be necessary to eliminate the result (following Oster 2019).
#' @usage o_delta_boot(y, x, con, m = "none", w = NULL, id = "none", time = "none", beta = 0, R2max,
#' sim, obs, rep, type, useed = NA, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent treatment variable (i.e., variable of interest; as string).
#' @param con Name of related control variables. Provided as string in the format: "w + z +...".
#' @param m Name of unrelated control variables (m; see Oster 2019; as string; default is m = "none").
#' @param w weights (only for weighted estimations). Warning: For weighted panel models R can report different R-square than Stata, leading deviation between R and Stata results.
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect panel models.
#' @param time Name of the time id variable (e.g. year or month; as string). Only applicable for fixed effect panel models.
#' @param beta beta for which delta*s should be estimated (default is beta = 0).
#' @param R2max Maximum R-square for which delta*s should be estimated.
#' @param sim Number of simulations.
#' @param obs Number of draws per simulation.
#' @param rep Bootstrapping either with (= TRUE) or without (= FALSE) replacement.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param useed User defined seed.
#' @param data Dataset.
#' @details Estimates bootstrapped delta*s, i.e., the degree of selection on unobservables relative to observables (with respect to the treatment variable) that would be necessary to eliminate the result (following Oster 2019). Bootstrapping can either be done with or without replacement. The function supports linear cross-sectional (see \emph{lm} objects in R) and fixed effect panel (see \emph{plm} objects in R) models.
#' @return Returns tibble object, which includes bootstrapped delta*s.
#' @references Oster, E. (2019). Unobservable Selection and Coefficient Stability: Theory and Evidence. Journal of Business & Economic Statistics, 37, 187-204.
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
#' # estimate bootstrapped delta*s
#' o_delta_boot(y = "mpg",          # dependent variable
#'              x = "wt",           # independent treatment variable
#'              con = "hp + qsec",  # related control variables
#'              beta = 0,           # beta
#'              R2max = 0.9,        # maximum R-square
#'              sim = 100,          # number of simulations
#'              obs = 30,           # draws per simulation
#'              rep = FALSE,        # bootstrapping with or without replacement
#'              type = "lm",        # model type
#'              useed = 123,        # seed
#'              data = data_oster)  # dataset
#' @export

o_delta_boot <- function(y, x, con, m = "none", w = NULL, id = "none", time = "none", beta = 0, R2max, sim, obs, rep, type, useed = NA, data) {


  # first define R-square of uncontrolled model
  ## rename variables and create, if type is plm, a panel dataset
  ### lm
  if (type == "lm") {
    if (is.null(w)) { data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x))
    w_var <- w } else { # no weights are assigned
      data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), w_var = all_of(w))}} # weights are assigned
  ### plm
  if (type == "plm") {
    if (is.null(w)) { data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), id_var = all_of(id), time_var = all_of(time))
    w_var <- w } else { # no weights are assigned
      data_aux <- data %>% dplyr::rename(y_var = all_of(y), x_var = all_of(x), id_var = all_of(id), time_var = all_of(time), w_var = all_of(w))}
    data_plm_aux <- pdata.frame(data_aux, index=c("id_var","time_var"), drop.index=TRUE, row.names=TRUE)}

  ## define uncontrolled model
  ### models without variables m
  if(m == "none"){model1_formula    <- as.formula(paste("y_var ~ x_var +",con)) } else { # controlled model
  ### models with
    model1_formula    <- as.formula(paste("y_var ~ x_var +",con, "+", m)) }              # controlled model

  ## run model
  ### lm
  if (type == "lm") { model1    <- lm(model1_formula,  weights = w_var,   data = data_aux, na.action = na.exclude)}                                      # run controlled model
  ### plm
  if (type == "plm") {model1    <- plm(model1_formula, weights = w_var,    data = data_plm_aux, model = "within", na.action = na.exclude) }              # run controlled model

  ## get R-square
  if (type == "lm") {R21_aux = summary(model1)$r.squared} else if (type == "plm") {R21_aux = as.numeric(summary(model1)$r.squared[1])}  # R-square uncontrolled model


  # keep original data
  data_original <- data


  # create tibble object to store results
  simulation_results <- tibble(
    "x" = 1:sim,
    "delta*" = 0)


  for (s in 1:sim) {

    if (!is.na(useed)) {set.seed(useed+s) } # set seed (defined by user) that it varies with runs

    # draw data
    if (rep) {
      data_current <- sample_n(data_original, size = obs, replace = TRUE) # with replacement
    } else {
      data_current <- sample_n(data_original, size = obs, replace = FALSE) # without replacement
    }


    # call o_delta function
    delta_results_aux <- o_delta(y = y, x = x, con = con, m = m, id = id, time = time, beta = beta, R2max = R2max, type = type, data = data_current)


    # store results
    simulation_results[s,2] <- round(delta_results_aux[1,2],6)
  }


  # warning if max R-square is smaller than the R-square of the control model
  if (R21_aux > R2max)
    warning("The max R-square value is smaller than the R-square of the controlled model")


  # warning if the bootstrapped observations are larger than the number of observations of the original dataset
  data_aux2 <- data %>% dplyr::rename(y_var = all_of(y))
  if (obs >= length(data_aux2$y_var))
    warning("Number of bootstrapped observation is larger than/equal to number of observation in the provided dataset")


  #define return object
  return(simulation_results)
}



##################--------------------------------------------------------------
#------------------------------- 1.11 o_delta_boot_viz - function 11
##################--------------------------------------------------------------
#' @title  Visualization of bootstrapped delta*s
#'
#' @description Estimates and visualizes bootstrapped delta*s, i.e., the degree of selection on unobservables relative to observables (with respect to the treatment variable) that would be necessary to eliminate the result (following Oster 2019).
#' @usage o_delta_boot_viz(y, x, con, m = "none", w = NULL, id = "none", time = "none",
#' beta = 0, R2max, sim, obs, rep, CI, type, norm = TRUE, bin,
#' col = c("#08306b","#4292c6","#c6dbef"), nL = TRUE, mL = TRUE, useed = NA, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent treatment variable (i.e., variable of interest; as string).
#' @param con Name of related control variables. Provided as string in the format: "w + z +...".
#' @param m Name of unrelated control variables (m; see Oster 2019; as string; default is m = "none").
#' @param w weights (only for weighted estimations). Warning: For weighted panel models R can report different R-square than Stata, leading deviation between R and Stata results.
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect panel models.
#' @param time Name of the time id variable (e.g. year or month; as string). Only applicable for fixed effect panel models.
#' @param beta beta for which delta*s should be estimated (default is beta = 0).
#' @param R2max Maximum R-square for which delta*s should be estimated.
#' @param sim Number of simulations.
#' @param obs Number of draws per simulation.
#' @param rep Bootstrapping either with (= TRUE) or without (= FALSE) replacement
#' @param CI  Confidence intervals, indicated as vector. Can be and/or 90, 95, 99.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param norm Option to include a normal distribution in the plot (default is norm = TURE).
#' @param bin Number of bins used in the histogram.
#' @param col Colors used to indicate different confidence interval levels (indicated as vector). Needs to be the same length as the variable CI. The default is a blue color range.
#' @param nL Option to include a red vertical line at 0 (default is nL = TRUE).
#' @param mL Option to include a vertical line at beta* mean (default is mL = TRUE).
#' @param useed User defined seed.
#' @param data Dataset.
#' @details Estimates and visualizes bootstrapped delta*s, i.e., the degree of selection on unobservables relative to observables (with respect to the treatment variable) that would be necessary to eliminate the result (following Oster 2019). Bootstrapping can either be done with or without replacement. The function supports linear cross-sectional (see \emph{lm} objects in R) and fixed effect panel (see \emph{plm} objects in R) models.
#' @return Returns ggplot2 object, which depicts the bootstrapped delta*s.
#' @references Oster, E. (2019). Unobservable Selection and Coefficient Stability: Theory and Evidence. Journal of Business & Economic Statistics, 37, 187-204.
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
#' # estimate and visualize bootstrapped delta*s
#' o_delta_boot_viz(y = "mpg",           # dependent variable
#'                  x = "wt",            # independent treatment variable
#'                  con = "hp + qsec",   # related control variables
#'                  beta = 0,            # beta
#'                  R2max = 0.9,         # maximum R-square
#'                  sim = 100,           # number of simulations
#'                  obs = 30,            # draws per simulation
#'                  rep = FALSE,         # bootstrapping with or without replacement
#'                  CI = c(90,95,99),    # confidence intervals
#'                  type = "lm",         # model type
#'                  norm = TRUE,         # normal distribution
#'                  bin = 200,           # number of bins
#'                  useed = 123,         # seed
#'                  data = data_oster)   # dataset
#' @export


o_delta_boot_viz <- function(y, x, con, m = "none", w = NULL, id = "none", time = "none", beta = 0, R2max, sim, obs, rep, CI, type, norm = TRUE, bin, col = c("#08306b","#4292c6","#c6dbef"), nL = TRUE, mL = TRUE, useed = NA, data) {

  # call function o_delta_boot and store results
  simulation_results <- o_delta_boot(y = y, x = x, con = con, m = m, id = id, time = time, beta = beta, R2max = R2max,
                                     sim = sim, obs = obs, rep = rep, type = type, useed = useed, data = data)
  simulation_results <- simulation_results %>% dplyr::rename("delta" = 2) # rename variable


  # compute mean, max absolute distance to mean, and confidence intervals
  ## mean
  meanDelta <- mean(simulation_results$delta, na.rm = F)
  # sd
  sdDelta <- sd(simulation_results$delta, na.rm = F)
  ## confidence interval
  ### 90%
  CI90_a <- meanDelta - 1.645 * sdDelta
  CI90_b <- meanDelta + 1.645 * sdDelta
  ### 95%
  CI95_a <- meanDelta - 1.96 * sdDelta
  CI95_b <- meanDelta + 1.96 * sdDelta
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
  delta_plot <- ggplot(simulation_results, aes(x = delta)) +
  {if(is.element(99, CI))annotate("rect", ymin = 0, ymax = +Inf, xmin = CI99_a, xmax = CI99_b, fill = color99, alpha = 0.6)}+
  {if(is.element(95, CI))annotate("rect", ymin = 0, ymax = +Inf, xmin = CI95_a, xmax = CI95_b, fill = color95, alpha = 0.6)}+
  {if(is.element(90, CI))annotate("rect", ymin = 0, ymax = +Inf, xmin = CI90_a, xmax = CI90_b, fill = color90, alpha = 0.6)}+
    geom_histogram(aes(y = ..density..), alpha = 0.8, colour = "#737373", fill = "#737373", size = 0.1, bins = bin)+
    {if(norm)stat_function(fun = dnorm, args = list(mean = meanDelta, sd = sdDelta), size = 1.3)}+
    scale_y_continuous(name = "Frequency",expand = expansion(mult = c(0, .1)))+
    scale_x_continuous(name = expression(delta^"*"))+
    expand_limits(x = 0, y = 0)+
    theme(axis.title = element_text( size=15),
          axis.text  = element_text( size=13))+
          {if(mL)geom_vline(xintercept = meanDelta, size = 0.8, color = "#000000", alpha = 0.8)}+
          {if(nL)geom_vline(xintercept = 0, size = 0.8, color = "#e41a1c", alpha = 0.8)}+
    geom_hline(yintercept = 0)



  # define return object
  return(delta_plot)
}


##################--------------------------------------------------------------
#------------------------------- 1.12 o_delta_boot_inf - function 12
##################--------------------------------------------------------------
#' @title Bootstrapped mean delta* and confidence intervals
#'
#' @description Provides the mean and confidence intervals of bootstrapped delta*s, i.e., the degree of selection on unobservables relative to observables (with respect to the treatment variable) that would be necessary to eliminate the result (following Oster 2019).
#' @usage o_delta_boot_inf(y, x, con, m = "none", w = NULL, id = "none", time = "none",
#' beta = 0, R2max, sim, obs, rep, CI, type, useed = NA, data)
#' @param y Name of the dependent variable (as string).
#' @param x Name of the independent treatment variable (i.e., variable of interest; as string).
#' @param con Name of related control variables. Provided as string in the format: "w + z +...".
#' @param m Name of unrelated control variables (m; see Oster 2019; as string; default is m = "none").
#' @param w weights (only for weighted estimations). Warning: For weighted panel models R can report different R-square than Stata, leading deviation between R and Stata results.
#' @param id Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect panel models.
#' @param time Name of the time id variable (e.g. year or month; as string). Only applicable for fixed effect panel models.
#' @param beta beta for which delta*s should be estimated (default is beta = 0)..
#' @param R2max Maximum R-square for which delta*s should be estimated.
#' @param sim Number of simulations.
#' @param obs Number of draws per simulation.
#' @param rep Bootstrapping either with (= TRUE) or without (= FALSE) replacement
#' @param CI  Confidence intervals, indicated as vector. Can be and/or 90, 95, 99.
#' @param type Model type (either \emph{lm} or \emph{plm}; as string).
#' @param useed User defined seed.
#' @param data Dataset.
#' @details Provides the mean and confidence intervals of bootstrapped delta*s, i.e., the degree of selection on unobservables relative to observables (with respect to the treatment variable) that would be necessary to eliminate the result (following Oster 2019). Bootstrapping can either be done with or without replacement. The function supports linear cross-sectional (see \emph{lm} objects in R) and fixed effect panel (see \emph{plm} objects in R) models.
#' @return Returns tibble object, which includes the mean and confidence intervals of bootstrapped delta*s.
#' @references Oster, E. (2019). Unobservable Selection and Coefficient Stability: Theory and Evidence. Journal of Business & Economic Statistics, 37, 187-204.
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
#' # compute the mean and confidence intervals of estimated bootstrapped delta*s
#' o_delta_boot_inf(y = "mpg",            # dependent variable
#'                  x = "wt",             # independent treatment variable
#'                  con = "hp + qsec",    # related control variables
#'                  beta = 0,             # beta
#'                  R2max = 0.9,          # maximum R-square
#'                  sim = 100,            # number of simulations
#'                  obs = 30,             # draws per simulation
#'                  rep = FALSE,          # bootstrapping with or without replacement
#'                  CI = c(90,95,99),     # confidence intervals
#'                  type = "lm",          # model type
#'                  useed = 123,          # seed
#'                  data = data_oster)    # dataset
#' @export


o_delta_boot_inf <- function(y, x, con, m = "none", w = NULL, id = "none", time = "none", beta = 0, R2max, sim, obs, rep, CI, type, useed = NA, data) {

  # call function o_delta_boot and store results
  simulation_results <- o_delta_boot(y = y, x = x, con = con, m = m, id = id, time = time, beta = beta, R2max = R2max,
                                     sim = sim, obs = obs, rep = rep, type = type, useed = useed, data = data)
  simulation_results <- simulation_results %>% dplyr::rename("delta" = 2) # rename variable


  # compute mean, max absolute distance to mean, and confidence intervals
  ## mean
  meanDelta <- mean(simulation_results$delta, na.rm = F)
  # sd
  sdDelta <- sd(simulation_results$delta, na.rm = F)
  ## confidence interval
  ### 90%
  CI90_a <- meanDelta - 1.645 * sdDelta
  CI90_b <- meanDelta + 1.645 * sdDelta
  ### 95%
  CI95_a <- meanDelta - 1.96 * sdDelta
  CI95_b <- meanDelta + 1.96 * sdDelta
  ### 99%
  CI99_a <- meanDelta - 2.58 * sdDelta
  CI99_b <- meanDelta + 2.58 * sdDelta

  # create tibble object of results
  result_Delta_inf_aux1 <- tribble(
    ~Name,                                 ~Value,
    "delta* (mean)",                         round(meanDelta,6))
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

  empty_aux <- ""
  result_Delta_inf_aux5 <- tribble(
    ~Name,                                 ~Value,
    "Simulations",                         sim,
    "Observations",                        obs,
    "Max R-square",                        R2max,
    "Beta",                                beta,
    "Model:",                              NA,
    "Replacement:",                        NA)
  if (is.element(90, CI)) {result_Delta_inf_aux1 <- result_Delta_inf_aux1 %>% bind_rows(result_Delta_inf_aux2)}
  if (is.element(95, CI)) {result_Delta_inf_aux1 <- result_Delta_inf_aux1 %>% bind_rows(result_Delta_inf_aux3)}
  if (is.element(99, CI)) {result_Delta_inf_aux1 <- result_Delta_inf_aux1 %>% bind_rows(result_Delta_inf_aux4)}
  result_Delta_inf <- result_Delta_inf_aux1 %>% bind_rows(result_Delta_inf_aux5) %>%
    mutate(Value = as.character(Value),
           Value = ifelse(Name == "Model:", type,
                          ifelse(Name == "Replacement:", rep,Value)))



  # define return object
  return(result_Delta_inf)
}





