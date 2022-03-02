# simulation parameters
n.sims <- 20
Tbar <- 0.5

skip_on_cran()

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


test_that('Correlations are as expected for design d2.1_m2fr', {
    
    rho.default <- 0
    model.params.list$rho.default <- rho.default
    cor.tstat <- expect_warning(check_cor(
        d_m = 'd2.1_m2fr', model.params.list, Tbar = Tbar, n.sims = n.sims
    ))

    expect_equal(cor.tstat[1], rho.default, 0.1)
    expect_equal(cor.tstat[2], rho.default, 0.1)
    expect_equal(cor.tstat[3], rho.default, 0.1)

    rho.default <- 0.2
    model.params.list$rho.default <- rho.default
    cor.tstat <- expect_warning(check_cor(
        d_m = 'd2.1_m2fr', model.params.list, Tbar = Tbar, n.sims = n.sims
    ))

    expect_equal(cor.tstat[1], rho.default, 0.1)
    expect_equal(cor.tstat[2], rho.default, 0.1)
    expect_equal(cor.tstat[3], rho.default, 0.1)
    
    rho.default <- 0.5
    model.params.list$rho.default <- rho.default
    cor.tstat <- expect_warning(check_cor(
        d_m = 'd2.1_m2fr', model.params.list, Tbar = Tbar, n.sims = n.sims
    ))

    expect_equal(cor.tstat[1], rho.default, 0.1)
    expect_equal(cor.tstat[2], rho.default, 0.1)
    expect_equal(cor.tstat[3], rho.default, 0.1)

    rho.default <- 0.8
    model.params.list$rho.default <- rho.default
    cor.tstat <- expect_warning(check_cor(
        d_m = 'd2.1_m2fr', model.params.list, Tbar = Tbar, n.sims = n.sims
    ))

    expect_equal(cor.tstat[1], rho.default, 0.1)
    expect_equal(cor.tstat[2], rho.default, 0.1)
    expect_equal(cor.tstat[3], rho.default, 0.1)

    rho.default <- 0.99
    model.params.list$rho.default <- rho.default
    cor.tstat <- expect_warning(check_cor(
        d_m = 'd2.1_m2fr', model.params.list, Tbar = Tbar, n.sims = n.sims
    ))

    expect_equal(cor.tstat[1], rho.default, 0.1)
    expect_equal(cor.tstat[2], rho.default, 0.1)
    expect_equal(cor.tstat[3], rho.default, 0.1)
})