% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FEISTY_plotting.R
\name{plotSimulationShiny}
\alias{plotSimulationShiny}
\title{Plot 'General' Panel for Shiny App}
\usage{
plotSimulationShiny(sim)
}
\arguments{
\item{sim}{The data frame of FEISTY simulation results.}
}
\description{
This function is designed to give the figures in the `General` tab in the Shiny App. 
It can be called independently.
}
\details{
It makes a plot combo displayed in the `General` tab in the Shiny App, 
including \code{\link{plotBiomasstime}}, \code{\link{plotSpectra}}, and \code{\link{plotRates}}.
}
\seealso{
\code{\link{webFEISTY}} A shiny interface for visualizing FEISTY model results \cr
\code{\link{plotBiomasstime}} Time series of biomass \cr
\code{\link{plotSpectra}} Biomass spectra plot \cr
\code{\link{plotRates}} Plots for growth rate, mortality, and feeding level \cr
}
\author{
Yixin Zhao
}
