


quiet_pump_sample <- purrr::quietly(pump_sample)

# Test all three sample sizes with given set of parameters.
test_sample_triad <- function( target_power, nbar, J, K, seed, ... ) {
    set.seed(seed)
    noK <- is.null( K )

    if ( ! noK ) {
        K1 <- quiet_pump_sample(
            typesample = 'K', target.power = target_power,
            J = J, nbar = nbar, ... )
    }

    J1 <- quiet_pump_sample(
        typesample = 'J', target.power = target_power,
        K = K, nbar = nbar, ... )

    nbar1 <- quiet_pump_sample(
        typesample = 'nbar', target.power = target_power,
        J = J, K = K, ... )

    if ( noK ) {
        list( J = J1$result$`Sample.size`, nbar = nbar1$result$`Sample.size`,
              Jrun = J1$result, nbarrun= nbar1$result,
              Jwarn = J1$warnings, nbarwarn = nbar1$warnings)

    } else {
        list( K = K1$result$`Sample.size`, J = J1$result$`Sample.size`, nbar = nbar1$result$`Sample.size`,
              Krun = K1$result, Jrun = J1$result, nbarrun= nbar1$result,
              Kwarn = K1$warnings, Jwarn = J1$warnings, nbarwarn = nbar1$warnings )
    }
}


warning_pattern <- function( vals ) {
    if ( is.null( vals$Kwarn ) ) {
        c( length( vals$nbarwarn ) > 0,
           length( vals$Jwarn ) > 0 )
    } else {
        c( length( vals$nbarwarn ) > 0,
           length( vals$Jwarn ) > 0,
           length( vals$Kwarn ) > 0 )
    }
}
