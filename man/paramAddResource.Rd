% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FEISTY_parms.R
\name{paramAddResource}
\alias{paramAddResource}
\title{Add Resource Parameters}
\usage{
paramAddResource (p, 
                         K, 
                         r=1, 
                         dynamics="chemostat", 
                         mc, 
                         mLower = NA, 
                         mUpper = NA, 
                         names=NA, 
                         u0=NA, 
                         ixpelR=NA, 
                         ixbenR=NA)
}
\arguments{
\item{p}{Parameter list to be updated.}

\item{K}{A vector of carrying capacities of all resources [g/m2].}

\item{r}{A vector containing each resource nudging rate (growth rate) [1/year], default 1.}

\item{dynamics}{The type of resource dynamics, either "chemostat" or "logistic".}

\item{mc}{A vector containing each resource geometric mean weight [g].}

\item{mLower}{A vector containing the lower limit of each resource weight [g]. Optional, depending on the size-based preference calculation function.}

\item{mUpper}{A vector containing the upper limit of each resource weight [g]. Optional, depending on the size-based preference calculation function.}

\item{names}{A character vector of each resource name (acronym). Optional, if not provided, default names are assigned, e.g., Resource_1 and Resource_2.}

\item{u0}{A vector of the initial concentration of each resource. If not provided, defaults to the value of \code{K}.}

\item{ixpelR}{An integer vector of indices of pelagic resources.}

\item{ixbenR}{An integer vector of indices of benthic resources.}
}
\value{
The updated parameter list \code{p}:
\itemize{
\item nResources, the total number of resources, which is the length of the input \code{K}.
\item groupnames, a character vector of names of resources and each fish functional group; value is updated later by other function(s)
\item stagenames, a character vector of names of resources and each fish stage (functional group); value is updated later by other function(s)
\item dynamics, from parameter input.
\item Rtype, resource growth strategy number. Default 1: chemostat, 2: logistic growth. If 'dynamics' is not specified, Rtype=1.
\item K, from parameter input.
\item r, from parameter input.
\item ixR, indices for all resources, start from 1.
\item mc, resource geometric mean weight, fish values are added later by other function(s).
\item mLower, the lower limit of each resource weight, fish values are added later by other function(s).
\item mUpper, the upper limit of each resource weight, fish values are added later by other function(s).
\item u0, from parameter input or same as \code{K}, fish values are added later by other function(s).
\item ixpelR, indices of pelagic resources.
\item ixbenR, indices of benthic resources.
}
}
\description{
This function updates the parameter list by adding resource-related parameters.
}
\details{
This function is designed to add parameters of resources to the parameter list.
This function only can be called once in each simulation. All resources are added to the parameter list at one time.
Generally, this function needs to be called after the function \code{\link{paramInit}}.
The returned parameter list can be used for further updates.
}
\examples{
# Initialize a paramter list
p = paramInit()
# Add five resources. The values are provided solely for illustrative purposes.
p = paramAddResource (p=p, K=c(50,100,150,200,250), r=c(1,1,1,1,1), dynamics="chemostat",
   mc= c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500), 1e-07*sqrt(500)),
   mLower = c(2e-06,0.001, 0.5e-03, 0.25, 1e-07),
   mUpper = c(0.001, 0.5, 125, 125, 5e-05),
   names=c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", "Phytoplankton"), 
   u0=c(100,100,200,200,200))

}
\seealso{
\code{\link{paramAddGroup}} Add parameters for one functional type
}
\author{
Ken H. Andersen, Karline Soetaert, Yixin Zhao
}
