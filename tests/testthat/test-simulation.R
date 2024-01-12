



test_that("simulation function works", {
    
    M <- 3
    rho.default <- 0.5
    default.rho.matrix <- gen_corr_matrix(M = M, rho.scalar = rho.default)
    default.kappa.matrix <- matrix(0, M, M) 
    model.params.list <- list(
        M = 3                             # number of outcomes
        , J = 30                          # number of schools
        , K = 10                          # number of districts
        # (for two-level model, set K = 1)
        , nbar = 50                       # number of individuals per school
        , S.id = NULL                     # N-length vector of school assignments
        , D.id = NULL                     # N-length vector of district assignments
        ################################################## grand mean outcome and impact
        , Xi0 = 0                         # scalar grand mean outcome under no treatment
        , MDES = rep(0.125, M)            # minimum detectable effect size      
        ################################################## level 3: districts
        , numCovar.3 = 1                  # number of district covariates
        , R2.3 = rep(0.1, M)              # percent of district variation
        # explained by district covariates
        , rho.V = default.rho.matrix      # MxM correlation matrix of district covariates
        , ICC.3 = rep(0.2, M)             # district intraclass correlation
        , omega.3 = rep(0.1, M)           # ratio of district effect size variability
        # to random effects variability
        , rho.w0 = default.rho.matrix     # MxM matrix of correlations for district random effects
        , rho.w1 = default.rho.matrix     # MxM matrix of correlations for district impacts
        , kappa.w =  default.kappa.matrix # MxM matrix of correlations between district
        # random effects and impacts
        ################################################## level 2: schools
        , numCovar.2 = 1                  # number of school covariates
        , R2.2 = rep(0.1, M)              # percent of school variation
        # explained by school covariates
        , rho.X = default.rho.matrix      # MxM correlation matrix of school covariates
        , ICC.2 = rep(0.2, M)             # school intraclass correlation	
        , omega.2 = rep(0.1, M)           # ratio of school effect size variability
        # to random effects variability
        , rho.u0 = default.rho.matrix     # MxM matrix of correlations for school random effects
        , rho.u1 = default.rho.matrix     # MxM matrix of correlations for school impacts
        , kappa.u = default.kappa.matrix  # MxM matrix of correlations between school
        # random effects and impacts
        ################################################## level 1: individuals
        , numCovar.1 = 1                  # number of individual covariates
        , R2.1 = rep(0.1, M)              # percent of indiv variation explained
        # by indiv covariates
        , rho.C = default.rho.matrix      # MxM correlation matrix of individual covariates
        , rho.r = default.rho.matrix      # MxM matrix of correlations for individual residuals 
    )
    
    sim.data <- gen_sim_data(d_m = 'd3.3',
                             param.list = model.params.list,
                             Tbar = 0.5, return.as.dataframe = FALSE)
    expect_true( !is.null(sim.data))
    expect_equal(8, length(sim.data))
    expect_equal(3, ncol(sim.data$Yobs))
    expect_equal(15000, nrow(sim.data$Yobs))
                 
    
    pp <- pump_power( d_m = "d3.2_m3ff2rc",
                      MTP = "BF",
                      MDES = rep( 0.10, 3 ),
                      M = 3,
                      J = 3, # number of schools/block
                      K = 21, # number RA blocks
                      nbar = 258,
                      Tbar = 0.50, # prop Tx
                      alpha = 0.05, # significance level
                      numCovar.1 = 5, numCovar.2 = 3,
                      R2.1 = 0.1, R2.2 = 0.7,
                      ICC.2 = 0.05, ICC.3 = 0.4,
                      rho = 0.4, # how correlated outcomes are
                      tnum = 200
    )
    sim.data <- gen_sim_data(pump.object = pp, return.as.dataframe = FALSE)
    expect_true( !is.null(sim.data))
    expect_equal(8, length(sim.data))
    expect_equal(3, ncol(sim.data$Yobs))
    expect_equal(16254, nrow(sim.data$Yobs))
    
    
    sim.data <- gen_sim_data(pump.object = pp, return.as.dataframe = TRUE, include_POs = TRUE)
    expect_true( length( sim.data ) == 3 )
    head( sim.data[[2]] )    
    expect_true( all( c( "Y0", "Y1", "Yobs", "T.x" ) %in% names( sim.data[[3]] ) ) )
})



test_that( "Simulation function works with MDES", {
    
    
    pp <- pump_sample( "d2.1_m2fr", MDES = 0.2, 
                       typesample = "nbar", alpha = 0.05,
                       power.definition = "D1indiv", target.power = 0.8,
                       M = 5, rho = 0.8,
                       MTP = "BH",
                       J = 15, Tbar = 0.5 )
    #summary( pp )
    ## 
    sim.data <- gen_sim_data( pp, return.as.dataframe = TRUE )
    
    expect_true( is.data.frame( head( sim.data[[1]] ) ) )
    expect_true( nrow( sim.data[[1]] ) == pp$Sample.size * 15 )
                 
})




test_that( "simulation function works (single outcome)", {
    
    model.params.list <- list(
        M = 1                            # number of outcomes
        , J = 30                          # number of schools
        , K = 10                          # number of districts
        # (for two-level model, set K = 1)
        , nbar = 10                       # number of individuals per school
        , S.id = NULL                     # N-length vector of school assignments
        , D.id = NULL                     # N-length vector of district assignments
        ################################################## grand mean outcome and impact
        , Xi0 = 0                         # scalar grand mean outcome under no treatment
        , MDES = 0.125            # minimum detectable effect size      
        ################################################## level 3: districts
        , R2.3 = 0.1              # percent of district variation
        # explained by district covariates
        , ICC.3 = 0.2 # district intraclass correlation
        , omega.3 = 0.1           # ratio of district effect size variability
        # random effects and impacts
        ################################################## level 2: schools
        , R2.2 = 0.1              # percent of school variation
        # explained by school covariates
        , ICC.2 = 0.2             # school intraclass correlation 
        , omega.2 = 0.1           # ratio of school effect size variability
        # to random effects variability
        # random effects and impacts
        ################################################## level 1: individuals
        , R2.1 = 0.1    # percent of indiv variation explained
        # by indiv covariates
    )
    
    # d3.2_m3fc2rc
    sim.data <- gen_sim_data( d_m = "d3.1_m3rr2rr", model.params.list, Tbar = 0.5,
                              return.as.dataframe = TRUE, no.list=FALSE )
    expect_true( length( sim.data ) == 1 )
    expect_true( is.data.frame( sim.data[[1]] ) )
    #head( sim.data[[1]] )
    
    sim.data <- gen_sim_data( d_m = "d3.1_m3rr2rr", model.params.list, Tbar = 0.5,
                              return.as.dataframe = TRUE, no.list=TRUE )
    expect_true( is.data.frame( sim.data ) )
    
    sim.data <- gen_sim_data( d_m = "d3.1_m3rr2rr", model.params.list, Tbar = 0.5,
                              return.as.dataframe = TRUE, no.list=TRUE, include_POs = TRUE )
    expect_true( is.data.frame( sim.data ) )
    expect_true( all( c( "Y0", "Y1", "Yobs" ) %in% names( sim.data ) ) )
    expect_true( all( sim.data$Y1 * sim.data$T.x == sim.data$Yobs * sim.data$T.x ) )
    
    sim.data <- gen_base_sim_data( model.params.list )
    #head( sim.data )
    expect_true( is.data.frame( sim.data ) )
    expect_true( all( c( "Y0", "Y1" ) %in% names( sim.data ) ) )
        
  
})

