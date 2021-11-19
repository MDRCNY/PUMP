


test_that( "optimize bounded logistic elements work", {
    
    # simulated data
    x <- c( seq(-10,20,len=4), 100 )
    size = rev( round( seq( 3, 50, length.out = length(x) )^2 ) )
    
    pmint = 0.2
    pmaxt = 0.8
    b0 = 1
    b1 = 0.5
    
    par_true = c( beta0=b0, beta1=b1, pmin=pmint, pmax=pmaxt)
    size
    
    probs = pum:::bounded_logistic_curve( x, par_true )
    p2 = pmint + (pmaxt-pmint) / ( 1 + exp( -1 * (b0 + b1*x) ) )
    p2
    
    probs == p2
    expect_true( probs[1] < probs[5] )
    
    y <- rbinom(n=length(x), size=size, prob=probs ) / size
    
    
    #qplot( x, y )  + coord_cartesian(ylim=c(0,1))
    par_est = pum:::fit_bounded_logistic( x, y, size )
    par_est
    
    if ( FALSE ) {
        dt = tibble( x = x, y=y, size=size )
        ggplot( dt ) +
            geom_point( aes( x, y, size=sqrt(size) ) )  +
            coord_cartesian(ylim=c(0,1)) +
            stat_function( aes( col="true" ), fun = function(x) { pum:::bounded_logistic_curve( x, par=par_true ) } ) +
            stat_function( aes( col="est" ), fun = function(x) { pum:::bounded_logistic_curve( x, par=par_est ) } )
    }
    
    
    xover = pum:::find_crossover( 0.4, par_true )
    
    expect_true( pum:::bounded_logistic_curve( xover, par_true ) == 0.4 )
    
    xover = pum:::find_crossover( 0.54, par_true )
    
    expect_true( pum:::bounded_logistic_curve( xover, par_true ) == 0.54 )
    

})



test_that( "optimize_power solves", {
    
    set.seed( 3042424 )
    op_pow <- pum:::optimize_power(
        MTP = "Holm", nbar = 200,
        power.definition = "D1indiv",
        design = "d2.1_m2fc", search.type = "J",
        start.low = 56, start.high = 75,
        start.tnum = 200,
        M = 3,
        MDES = 0.05, target.power = 0.80, tol = 0.01,
        Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
        rho = 0.4, max.tnum = 400, final.tnum = 2000
    )
    op_pow
    
    #pum:::plot_power_search(op_pow)
    expect_true( ncol( op_pow ) == 7 )
    expect_true( all( op_pow$w <= 2000 ) )
    expect_true( max( op_pow$w ) == 2000 )
    
})
