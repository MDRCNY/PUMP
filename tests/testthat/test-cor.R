model.params.list <- list(
    M = 3                             # number of outcomes
    , J = 30                          # number of schools
    , K = 10                          # number of districts
    , nbar = 50                       # number of individuals per school
    , rho.default = 0.5               # default rho value (optional)
    ################################################## impact
    , MDES = 0.125                    # minimum detectable effect size      
    ################################################## level 3: districts
    , numCovar.3 = 0                  # number of district covariates
    , R2.3 = 0                        # percent of district variation
    # explained by district covariates
    , ICC.3 = 0.2                     # district intraclass correlation
    , omega.3 = 0                     # ratio of district effect size variability
    # to random effects variability
    ################################################## level 2: schools
    , numCovar.2 = 0                  # number of school covariates
    , R2.2 = 0                        # percent of school variation
    # explained by school covariates
    , ICC.2 = 0.2                     # school intraclass correlation	
    , omega.2 = 0                     # ratio of school effect size variability
    # to random effects variability
    ################################################## level 1: individuals
    , numCovar.1 = 0                  # number of individual covariates
    , R2.1 = 0                        # percent of indiv variation explained
    # by indiv covariates
)


test_that('Correlations checker runs', {
    
    skip_on_cran()
    
    sink("sink.txt")
    
    # note: this is not enough simulations to check if the tool works,
    # just if it runs!
    
    rho.default <- 0
    model.params.list$rho.default <- rho.default
    cor.tstat <- expect_warning(check_cor(
        d_m = 'd2.1_m2fr', model.params.list, Tbar = 0.5, n.sims = 10
    ))
    
    expect_true(length(cor.tstat) == 3)
    expect_true(all(!is.na(cor.tstat)))
    
    sink()
    file.remove("sink.txt")
    
})
