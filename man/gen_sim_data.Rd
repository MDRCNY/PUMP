% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{gen_sim_data}
\alias{gen_sim_data}
\title{Generate simulated multi-level data (simulation function)}
\usage{
gen_sim_data(
  d_m = NULL,
  param.list = NULL,
  Tbar = 0.5,
  pump.object = NULL,
  return.as.dataframe = TRUE,
  no.list = TRUE,
  include_POs = FALSE
)
}
\arguments{
\item{d_m}{string; a single context, which is a design and model
code. See pump_info() for list of choices.}

\item{param.list}{list; model parameters such as ICC, R2, etc. See
simulation vignette for details.}

\item{Tbar}{scalar; the proportion of samples that are assigned to
the treatment.}

\item{pump.object}{A pumpresult object.}

\item{return.as.dataframe}{TRUE means return list of dataframes,
one for each outcome.  FALSE means return components of the
covariates, etc., in a list.}

\item{no.list}{Only relevant if return.as.dataframe=TRUE.
no.list=TRUE means if M=1 return the dataframe, not a list of
length 1.  FALSE means return a list of length 1, even if there
is only 1 outcome.}

\item{include_POs}{Include columns for the potential outcomes in
addition to the observed outcome.}
}
\value{
list; potential outcomes, covariates, observed outcomes,
  and treatment assignment.
}
\description{
Generates simulated data for multi-level RCTs for
  pump-supported designs and models for both unobserved and observed
  potential outcomes.

  Takes in two options:

  - a pumpresult object OR

  - a list of necessary data-generating parameters 
  - the context (d_m) 
  - Tbar (proportion assigned to treatment)

  This function is beyond the main scope of calculating power, and
  is instead used for simulating data. For more info on use, see
  the simulation vignette.
}
\examples{

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
                  rho = 0.4,
                  tnum = 200
)
sim.data <- gen_sim_data(pump.object = pp)

}
