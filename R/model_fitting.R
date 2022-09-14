

#' Function: get_rawpt
#'
#' Fits models to each of a series of datasets, and extracts p values
#' and t statistics
#'
#' @param dat.all list of datasets (of length M)
#' @param d_m design/model
#' @param model.params.list list of model parameters
#'
#' @keywords internal
get_rawpt <- function(dat.all, d_m, model.params.list) {
    mods.out <- lapply(dat.all, function(m) make_model(m, d_m))
    mods <- lapply(mods.out, function(m){ return(m[['mod']]) })
    singular <- sapply(mods.out, function(m){ return(m[['singular']]) })
    failed.converge <- sapply(mods.out, function(m){ return(m[['failed.converge']]) })
    rawpt <- lapply(mods, function(x) get_pval_tstat(x, d_m, model.params.list))
    return(list(rawpt = rawpt, num.singular = sum(singular), num.failed.converge = sum(failed.converge)))
}


#' Function: get_pval_tstat	                                       
#'  
#' extracts p-value and t statistics from a given fitted model.
#'
#' @param mod model object
#' @param d_m design/model
#' @param model.params.list list of model parameters
#' 
#' @keywords internal

#'  Function: get_pval_tstat	                                       
#'  
#' extracts p-value and t statistics from a given model
#'
#' @param mod model object
#' @param model.params.list list of model parameters
#' 
#' @keywords internal
get_pval_tstat <- function(mod, d_m, model.params.list) {
    if(methods::is(mod, "lm")) {
        tstat <- summary(mod)$coefficients["T.x","t value"]
        pval <- summary(mod)$coefficients["T.x","Pr(>|t|)"]
    } else if(methods::is(mod, "lmerMod")) {
        df <- calc_df(d_m, model.params.list[['J']], model.params.list[['K']],
                      model.params.list[['nbar']],
                      numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1)
        tstat <- summary(mod)$coefficients["T.x","t value"]
        pval <- (1 - stats::pt(abs(tstat), df = df))*2
    } else if(methods::is(mod, "data.frame")) {
        # fixed effects models
        df <- calc_df(d_m, model.params.list[['J']], model.params.list[['K']],
                      model.params.list[['nbar']], numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1)
        tstat <- mod$ATE_hat[1]/mod$SE[1]
        pval <- 2*(1 - stats::pt(abs(tstat), df = df))
    } else
    {
        stop('Unknown model type')
    }
    return(list(tstat = tstat, pval = pval))
}



#' Function: make.model
#'
#' Function to generate and fit a model to a simulated dataset.
#'
#' @param dat.m a data frame for a single outcome
#' @param d_m design/model
#'
#' @keywords internal
make_model <- function(dat.m, d_m) {
    # dat.m = dat.all[[1]];
    
    singular <- FALSE
    failed.converge <- FALSE
    
    dat.m$S.id <- as.factor(dat.m$S.id)
    if(!is.null(dat.m$D.id)){ dat.m$D.id <- as.factor(dat.m$D.id) }
    
    if (d_m == "d1.1_m1c") {
        form <- stats::as.formula("Yobs ~ 1 + T.x + C.ijk")
        mod <- stats::lm(form, data = dat.m)
    } else if (d_m == "d2.1_m2fc") {
        form <- stats::as.formula("Yobs ~ 1 + T.x + C.ijk + S.id")
        mod <- stats::lm(form, data = dat.m)
    } else if (d_m == "d2.1_m2ff") {
        mod.out <- interacted_linear_estimators(
            Yobs = dat.m$Yobs, Z = dat.m$T.x,
            B = dat.m$S.id, data = dat.m,
            control_formula = "C.ijk", use.lmer = FALSE
        )
        singular <- mod.out[['singular']]
        mod <- mod.out[['mod']]
    } else if (d_m == "d2.1_m2fr") {
        form <- stats::as.formula(paste0(
            "Yobs ~ 1 + T.x + X.jk + C.ijk + (1 + T.x | S.id)")
        )
        mod <- suppressMessages(lme4::lmer(form, data = dat.m))
        singular <- lme4::isSingular(mod)
        failed.converge <- ifelse(!is.null(mod@optinfo$conv$lme4$code), TRUE, FALSE)
    } else if (d_m == "d3.1_m3rr2rr") {
        form <- stats::as.formula(paste0(
            "Yobs ~ 1 + T.x + V.k + X.jk + C.ijk + (1 + T.x | S.id) + (1 + T.x | D.id)")
        )
        mod <- suppressMessages(lme4::lmer(form, data = dat.m))
        singular <- lme4::isSingular(mod)
        failed.converge <- ifelse(!is.null(mod@optinfo$conv$lme4$code), TRUE, FALSE)
    } else if (d_m == "d2.2_m2rc") {
        form <- stats::as.formula(paste0(
            "Yobs ~ 1 + T.x + X.jk + C.ijk + (1 | S.id)")
        )
        mod <- suppressMessages(lme4::lmer(form, data = dat.m))
        singular <- lme4::isSingular(mod)
        failed.converge <- ifelse(!is.null(mod@optinfo$conv$lme4$code), TRUE, FALSE)
    } else if (d_m == "d3.3_m3rc2rc") {
        form <- stats::as.formula(paste0(
            "Yobs ~ 1 + T.x + V.k + X.jk + C.ijk + (1 | S.id) + (1 | D.id)")
        )
        mod <- suppressMessages(lme4::lmer(form, data = dat.m))
        singular <- lme4::isSingular(mod)
        failed.converge <- ifelse(!is.null(mod@optinfo$conv$lme4$code), TRUE, FALSE)
    } else if (d_m == "d3.2_m3ff2rc") {
        mod.out <- interacted_linear_estimators(
            Yobs = dat.m$Yobs, Z = dat.m$T.x,
            B = dat.m$D.id, data = dat.m,
            control_formula = "X.jk + C.ijk + (1 | S.id)",
            use.lmer = TRUE
        )
        singular <- mod.out[['singular']]
        mod <- mod.out[['mod']]
    }
    else if (d_m == "d3.2_m3fc2rc") {
        form <- stats::as.formula(paste0(
            "Yobs ~ 1 + T.x + X.jk + C.ijk + D.id + (1 | S.id)")
        )
        mod <- suppressMessages(lme4::lmer(form, data = dat.m))
        singular <- lme4::isSingular(mod)
        failed.converge <- ifelse(!is.null(mod@optinfo$conv$lme4$code), TRUE, FALSE)
    } else if (d_m == "d3.2_m3rr2rc") {
        form <- stats::as.formula(paste0(
            "Yobs ~ 1 + T.x + V.k + X.jk + C.ijk + (1 | S.id) + (1 + T.x | D.id)")
        )
        mod <- suppressMessages(lme4::lmer(form, data = dat.m))
        singular <- lme4::isSingular(mod)
        failed.converge <- ifelse(!is.null(mod@optinfo$conv$lme4$code), TRUE, FALSE)
    } else {
        stop(paste('Unknown d_m:', d_m)) 
    }
    
    return(list(mod = mod, singular = singular, failed.converge = failed.converge))
}


#' Interacted linear regression models
#' 
#' Code taken from:
#' https://github.com/lmiratrix/blkvar/blob_master/R/linear_model_method.R
#'
#' These linear models have block by treatment interaction terms.  The final ATE
#' estimates are then weighted average of the block (site) specific ATE
#' estimates.
#'
#' If siteID passed, it will weight the RA blocks within site and then average
#' these site estimates.
#'
#' SEs come from the overall variance-covariance matrix.
#'
#' @return Dataframe of the different versions of this estimator (person and
#'   site weighted)
#' @family linear model estimators
#' 
#' @keywords internal
interacted_linear_estimators <- function(Yobs, Z, B, siteID = NULL, data = NULL,
                                         control_formula = NULL, use.lmer = FALSE) {
    # siteID = NULL;
    # Yobs = dat.m$Yobs; Z = dat.m$T.x; B = dat.m$D.id; data = dat.m; control_formula = "X.jk + C.ijk + (1 | S.id)"; use.lmer = TRUE
    # Yobs = dat.m$Yobs; Z = dat.m$T.x; B = dat.m$S.id; data = dat.m; control_formula = "C.ijk"; use.lmer = FALSE
    
    # keep track of singularity
    singular <- FALSE
    
    # This code block takes the parameters of
    # Yobs, Z, B, siteID = NULL, data=NULL, ...
    # and makes a dataframe with canonical Yobs, Z, B, and siteID columns.
    d2 <- data
    d2$Yobs <- Yobs
    d2$Z <- Z
    d2$B <- B
    data <- d2
    rm(d2)
    
    
    data$B <- droplevels(as.factor(data$B))
    J <- length(unique(data$B))
    nj <- table(data$B)
    n <- nrow(data)
    
    formula <- stats::as.formula(sprintf(
        "%s ~ 0 + %s * %s - %s + %s", "Yobs", "Z", "B", "Z", control_formula))
    
    if(use.lmer)
    {
        M0.int <- lme4::lmer(formula, data = data)
        ids <- grep( "Z:", rownames(summary(M0.int)$coefficients))
    } else
    {
        M0.int <- stats::lm(formula, data = data)
        ids <- grep( "Z:", names(stats::coef(M0.int)))
    }
    
    if(length(ids) != J)
    {
        message('Proceeding with rank deficient model')
        singular <- TRUE
        J <- length(ids)
    }
    
    VC <- as.matrix(stats::vcov(M0.int))
    ATE_hats <- summary(M0.int)$coefficients[ids,1]
    wts <- rep(1 / J, J)
    
    # the block SEs from our linear model
    SE_hat <- diag(VC)[ids]
    
    ATE_hat_site <- stats::weighted.mean(ATE_hats, wts)
    
    # Calculate SE for ATE_hat_site
    SE_site <- sqrt(sum(wts ^ 2 * SE_hat))
    
    interactModels <- data.frame(method = c("FE-Int-Sites"),
                                 ATE_hat = c(ATE_hat_site),
                                 SE = c(SE_site), stringsAsFactors = FALSE)
    if (!is.null(control_formula)) {
        interactModels$method <- paste0(interactModels$method, "-adj")
    }
    return(list(mod = interactModels, singular = singular))
}
