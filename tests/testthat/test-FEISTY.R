#
# Compare the total biomass of each functional type in the last year of the simulation for the four setups.
#

# Baseline values are from R methods (USEdll=F)
totB_setupBasic <- c(totBiomass.smallPel = 0.5946745, 
                     totBiomass.largePel = 4.0246085, 
                     totBiomass.demersals = 2.9236847)
totB_setupBasic2 <- c(totBiomass.smallPel = 13.155201, 
                      totBiomass.largePel = 12.596843, 
                      totBiomass.demersals = 4.715339)

totB_setupVertical <- c(totBiomass.smallPel = 2.847370, 
                        totBiomass.mesoPel = 6.582644, 
                        totBiomass.largePel = 10.344241, 
                        totBiomass.midwPred = 1.862001,
                        totBiomass.demersals = 2.870958)

totB_setupVertical2 <- c(totBiomass.smallPel = 5.2442065, 
                         totBiomass.mesoPel = 8.8442285, 
                         totBiomass.largePel = 6.3042523, 
                         totBiomass.midwPred = 0.9813705,
                         totBiomass.demersals = 2.8422796)

test_that('setupBasic total biomass', {
  p <- setupBasic(szprod = 100,
                  lzprod = 100,
                  depth  = 100,
                  Tp     = 10,
                  Tb     = 8)
  
  # Test with USEdll = T and bCust = F
  sim <- simulateFEISTY(p = p, 
                        tEnd   = 500,
                        tStep  = 1,
                        yini   = p$u0,
                        USEdll = T, 
                        bCust  = F)
  expect_equal(sim$totBiomass[sim$nTime,], 
               totB_setupBasic,
               tolerance = 1E-2,
               info = "Failure with USEdll = T and bCust = F")
  
  # Test with USEdll = T and bCust = T
  sim <- simulateFEISTY(p = p, 
                        tEnd   = 500,
                        tStep  = 1,
                        yini   = p$u0,
                        USEdll = T, 
                        bCust  = T)
  expect_equal(sim$totBiomass[sim$nTime,], 
               totB_setupBasic,
               tolerance = 1E-2,
               info = "Failure with USEdll = T and bCust = T")
  
  # Test with USEdll = F
  sim <- simulateFEISTY(p = p, 
                        tEnd   = 500,
                        tStep  = 1,
                        yini   = p$u0,
                        USEdll = F)
  expect_equal(sim$totBiomass[sim$nTime,], 
               totB_setupBasic,
               tolerance = 1E-2,
               info = "Failure with USEdll = F")
})


test_that('setupBasic2 total biomass', {
  p <- setupBasic2(szprod = 100,
                   lzprod = 100,
                   bprodin = NA,
                   dfbot  = NA,
                   depth  = 100,
                   Tp     = 10,
                   Tb     = 8,
                   nStages = 9,
                   etaMature = 0.25,
                   Fmax = 0,
                   etaF = 0.05,
                   bET  = TRUE)
  
  # Test with USEdll = T and bCust = F
  sim <- simulateFEISTY(p = p,
                        tEnd   = 500,
                        tStep  = 1,
                        yini   = p$u0,
                        USEdll = T, 
                        bCust  = F)
  expect_equal(sim$totBiomass[sim$nTime,], 
               totB_setupBasic2,
               tolerance = 1E-2,
               info = "Failure with USEdll = T and bCust = F")
  
  # Test with USEdll = T and bCust = T
  sim <- simulateFEISTY(p = p, 
                        tEnd   = 500,
                        tStep  = 1,
                        yini   = p$u0,
                        USEdll = T, 
                        bCust  = T)
  expect_equal(sim$totBiomass[sim$nTime,], 
               totB_setupBasic2,
               tolerance = 1E-2,
               info = "Failure with USEdll = T and bCust = T")
  
  # Test with USEdll = F
  sim <- simulateFEISTY(p = p, 
                        tEnd   = 500,
                        tStep  = 1,
                        yini   = p$u0,
                        USEdll = F)
  expect_equal(sim$totBiomass[sim$nTime,], 
               totB_setupBasic2,
               tolerance = 1E-2,
               info = "Failure with USEdll = F")
})


test_that('setupVertical total biomass', {
  p <- setupVertical(szprod = 80,
                     lzprod = 80,
                     dfpho  = 350,
                     region = 4,
                     depth  = 800,
                     photic = 150)
  
  # Test with USEdll = T and bCust = F
  sim <- simulateFEISTY(p = p, 
                        tEnd   = 500,
                        tStep  = 1,
                        yini   = p$u0,
                        USEdll = T, 
                        bCust  = F)
  expect_equal(sim$totBiomass[sim$nTime,], 
               totB_setupVertical,
               tolerance = 1E-2,
               info = "Failure with USEdll = T and bCust = F")
  
  # Test with USEdll = T and bCust = T
  sim <- simulateFEISTY(p = p, 
                        tEnd   = 500,
                        tStep  = 1,
                        yini   = p$u0,
                        USEdll = T, 
                        bCust  = T)
  expect_equal(sim$totBiomass[sim$nTime,], 
               totB_setupVertical,
               tolerance = 1E-2,
               info = "Failure with USEdll = T and bCust = T")
  
  # Test with USEdll = F
  sim <- simulateFEISTY(p = p, 
                        tEnd   = 500,
                        tStep  = 1,
                        yini   = p$u0,
                        USEdll = F)
  expect_equal(sim$totBiomass[sim$nTime,], 
               totB_setupVertical,
               tolerance = 1E-2,
               info = "Failure with USEdll = F")
})


test_that('setupVertical2 total biomass', {
  p <- setupVertical2(szprod = 80,
                      lzprod = 80,
                      dfpho  = 350,
                      nStages = 9,
                      Tp = 10,
                      Tm = 10,
                      Tb = 10,
                      depth = 800,
                      photic = 150,
                      shelfdepth = 250,
                      visual = 1.5,
                      etaMature = 0.25,
                      Fmax = 0,
                      etaF = 0.05)
  
  # Test with USEdll = T and bCust = F
  sim <- simulateFEISTY(p = p, 
                        tEnd   = 500,
                        tStep  = 1,
                        yini   = p$u0,
                        USEdll = T, 
                        bCust  = F)
  expect_equal(sim$totBiomass[sim$nTime,], 
               totB_setupVertical2,
               tolerance = 1E-2,
               info = "Failure with USEdll = T and bCust = F")
  
  # Test with USEdll = T and bCust = T
  sim <- simulateFEISTY(p = p, 
                        tEnd   = 500,
                        tStep  = 1,
                        yini   = p$u0,
                        USEdll = T, 
                        bCust  = T)
  expect_equal(sim$totBiomass[sim$nTime,], 
               totB_setupVertical2,
               tolerance = 1E-2,
               info = "Failure with USEdll = T and bCust = T")
  
  # Test with USEdll = F
  sim <- simulateFEISTY(p = p, 
                        tEnd   = 500,
                        tStep  = 1,
                        yini   = p$u0,
                        USEdll = F)
  expect_equal(sim$totBiomass[sim$nTime,], 
               totB_setupVertical2,
               tolerance = 1E-2,
               info = "Failure with USEdll = F")
})