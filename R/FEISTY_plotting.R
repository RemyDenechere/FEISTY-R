#===============================================================================
# Plotting routines for FEISTY
#===============================================================================

# tools and generic sizes
lwd_def     <- 1 # line thickness for time series, rates and spectra plots
min_bio     <- 1E-2 # minimum biomass or SSB in time series plots
min_bspec   <- 1E-3 # minimum biomass in spectra plot
min_g       <- 1E-1 # minimum growth rate in plotRates growth
min_mort    <- 1E-2 # minimum mortality in plotRates mortality

#' Spawning stock biomass plot
#' 
#' Makes a time series of spawning stock biomass (SSB) for each functional group
#' 
#' @details X-axis is time. Y-axis is the SSB data in log10 scale.
#' 
#' @author P. Daniël van Denderen, Yixin Zhao
#'
#' @usage plotSSBtime(sim)
#' 
#' @param sim The data frame of FEISTY simulation results.
#' 
#' @examples 
#' sim=simulateFEISTY()
#' plotSSBtime(sim)
#' 
#' @aliases plotSSBtime
#' 
#' @seealso 
#' \code{\link{simulateFEISTY}} Run FEISTY model simulations \cr
#' \code{\link{calcSSB}} 	Spawning stock biomass calculation
#' 
#' @export
#' 

plotSSBtime = function(sim) {
  p <- sim$p
  fish <- c((p$nResources+1):(p$nResources+p$nGroups))
# extract existing fish types according to initial values
  inival=p$u0[-p$ixR][grep("_1", names(p$u0[-p$ixR]))]
  names(inival)=p$groupnames[fish]
  fishexistname=names(inival[inival!=0])
 
  series <- getTimeseries(sim)
  series <- subset(series, series$group %in% fishexistname)
  
# convert zeros of SSB to a very small number (only for plotting purposes)
  series$SSB <- ifelse(series$SSB == 0, 1E-16,series$SSB) 
  
  plots <- defaultplot() + 
    geom_line(data = series, 
              aes(x = t, y = SSB, group = group, color = group),
              linewidth=lwd_def) + 
    scale_color_manual(name = "Groups", 
                       values = p$my_palette[fishexistname],
                       breaks = fishexistname,
                       labels = p$my_names[fishexistname]) +
    annotation_logticks(sides = "l",size = 0.2,colour = "darkgrey") +
    coord_cartesian(ylim = c(min_bio,max(min_bio*100,max(series$SSB)*5))) + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("Time (yr)") + ylab(expression("SSB (gWW m"^"-2"*")"))+
    theme(legend.key = element_blank())
    
  return(plots)
}
  
#' Fishing Yield plot
#' 
#' Makes a time series of the fishing yield for each functional group
#' 
#' @details X-axis is time. Y-axis is the fishing yield data.
#' 
#' @author P. Daniël van Denderen, Yixin Zhao
#'
#' @usage plotYieldtime(sim)
#' 
#' @param sim The data frame of FEISTY simulation results.
#' 
#' @examples 
#' sim=simulateFEISTY()
#' plotYieldtime(sim)
#' 
#' @aliases plotYieldtime
#' 
#' @seealso 
#' \code{\link{simulateFEISTY}} Run FEISTY model simulations \cr
#' \code{\link{calcYield}} Fishing yield calculation
#' 
#' @export
#' 

plotYieldtime = function(sim) {
  p <- sim$p
  fish <- c((p$nResources+1):(p$nResources+p$nGroups))
# extract existing fish types according to initial values
  inival=p$u0[-p$ixR][grep("_1", names(p$u0[-p$ixR]))]
  names(inival)=p$groupnames[fish]
  fishexistname=names(inival[inival!=0])
  
  series <- getTimeseries(sim)
  series <- subset(series, series$group %in% fishexistname)
  
  plots <- defaultplot() + 
    geom_line(data = series, 
              aes(x = t, y = yield, group = group, color = group),
              linewidth=lwd_def) + 
    scale_color_manual(name = "Groups", 
                       values = p$my_palette[fishexistname],
                       breaks = fishexistname,
                       labels = p$my_names[fishexistname]) +
    annotation_logticks(sides = "l",size = 0.2,colour = "darkgrey") +
    coord_cartesian(ylim = c(0,max(1,max(series$yield)*1.2))) + 
    xlab("Time (yr)") + ylab(expression("Yield (gWW m"^"-2"*" yr"^"-1"*")"))+
    theme(legend.key = element_blank())

  return(plots)
}

#' Total biomass plot
#' 
#' Makes a time series of total biomass for each fish functional group and zooplankton
#' and benthic resource
#' 
#' @details X-axis is time. Y-axis is the biomass data in log10 scale.
#' 
#' @author P. Daniël van Denderen, Yixin Zhao
#'
#' @usage plotBiomasstime(sim)
#' 
#' @param sim The data frame of FEISTY simulation results.
#' 
#' @examples 
#' sim=simulateFEISTY()
#' plotBiomasstime(sim)
#' 
#' @aliases plotBiomasstime
#' 
#' @seealso 
#' \code{\link{simulateFEISTY}} Run FEISTY model simulations
#' 
#' @export
#' 

plotBiomasstime = function(sim) {
  p <- sim$p
  fish <- c((p$nResources+1):(p$nResources+p$nGroups))
  all  <- c(p$ixR,fish)
  
# extract existing resources and fish types according to initial values
  inival=c( p$u0[p$ixR], p$u0[-p$ixR][grep("_1", names(p$u0[-p$ixR]))] )
  names(inival)=p$groupnames
  allexistname=names(inival[inival!=0])
  
  series <- getTimeseries(sim)
  series <- subset(series, series$group %in% allexistname)
  
# convert zeros of biomass to a very small number (only for plotting purposes)
  series$bio <- ifelse(series$bio == 0, 1E-16,series$bio) 
  
  plots <- defaultplot() + 
    geom_line(data = series, 
              aes(x = t, y = bio, group = group, color = group),
              linewidth=lwd_def) + 
    scale_color_manual(name = "Groups", 
                       values = p$my_palette[allexistname],
                       breaks = allexistname,
                       labels = p$my_names[allexistname]) +
    annotation_logticks(sides = "l",size = 0.2,colour = "darkgrey") +
    coord_cartesian(ylim = c(min_bio,max(min_bio*100,max(series$bio)*5))) + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("Time (yr)") + ylab(expression("Biomass (gWW m"^"-2"*")"))+
    theme(legend.key = element_blank())
  
  return(plots)
}

#' Biomass spectra plot
#' 
#' Makes a plot of the biomass of all functional groups over the size spectrum.
#' Data is averaged over the last 40\% simulation time.
#' 
#' @details The X-axis is individual weight on a log10 scale. The Y-axis is 
#' biomass on a log10 scale.
#' 
#' @author P. Daniël van Denderen, Yixin Zhao
#'
#' @usage plotSpectra(sim)
#' 
#' @param sim The data frame of FEISTY simulation results.
#' 
#' @examples 
#' sim=simulateFEISTY()
#' plotSpectra(sim)
#' 
#' @aliases plotSpectra
#' 
#' @seealso 
#' \code{\link{simulateFEISTY}} Run FEISTY model simulations
#' 
#' @export
#' 

plotSpectra = function(sim) {
  p <- sim$p
  Bpositive <- sim$B
  Bpositive[Bpositive<0] <- 0
  iTime <- sim$nTime
  fish <- c((p$nResources+1):(p$nResources+p$nGroups))
  i     <- which.max(sapply(p$ix,"length"))
  ix    <- p$ix[[i]]
  
  
# get all biomass spectra
  spec = matrix(nrow=p$nGroups, ncol=max(sapply(p$ix, length)),data=0)
  for (i in 1:p$nGroups) {
    spec[i,1:length(p$ix[[i]])] = colMeans(Bpositive[round(0.6*iTime):iTime,p$ix[[i]]-p$ixFish[1]+1])
  }
  spec <- rbind(colSums(spec,na.rm=T),spec)
  spec <- data.frame(group = rep(c("total",p$groupnames[fish]),each = max(sapply(p$ix, length))),
             bio   = as.vector(t(spec)),
             mc    = rep(p$mc[ix],(p$nGroups+1)))
# correct small fish size numbers (biomass are 0)
  spec <- subset(spec, !(spec$group %in% c("smallPel","mesoPel") & spec$bio ==0))

# extract existing fish types according to initial values
  inival=p$u0[-p$ixR][grep("_1", names(p$u0[-p$ixR]))]
  names(inival)=p$groupnames[fish]
  fishexistname=names(inival[inival!=0])
  
  spec <- subset(spec, spec$group %in% c("total",fishexistname))
  
# convert zeros of yield to a very small number (only for plotting purposes)
  spec$bio <- ifelse(spec$bio == 0, 1E-16,spec$bio)     
  
  plots <- defaultplot() +
    geom_line(data = spec, 
             aes(x = mc, y = bio, group = group, color = group),
             linewidth=lwd_def) + 
    scale_color_manual(name = "Groups", 
                       values = c("total" = "black", p$my_palette[fishexistname]),
                       breaks = c("total",fishexistname),
                       labels = c("total"= "Total", p$my_names[fishexistname])) +
    annotation_logticks(sides = "bl",size = 0.2,colour = "darkgrey") +
    coord_cartesian(ylim = c(min_bspec,max(min_bspec*100,max(spec$bio)*10))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("Weight (gWW)") + ylab(expression("Biomass (gWW m"^"-2"*")"))+
    theme(legend.key = element_blank())

  
  return(plots)
}

#' Food web plot
#' 
#' Makes a plot with all feeding interactions of fish functional groups and 
#' resources. The size of the biomass circles and feeding fluxes are divided by 
#' the maximum biomass or flux and re-scaled with the cube root.
#' The function works on the four prepared setups or revised 
#' versions based on these four setups. If customized setups by users have 
#' more resources and/or fish functional types, this plot function may need to be updated.
#' 
#' @details
#' To be added.
#' 
#' @author Daniel Ottmann Riera, P. Daniël van Denderen, Yixin Zhao
#'
#' @usage plotNetwork(sim)
#' 
#' @param sim The dataframe of FEISTY simulation results.
#' 
#' @examples 
#' sim=simulateFEISTY()
#' plotNetwork(sim)
#' 
#' @aliases plotNetwork
#' 
#' @seealso 
#' \code{\link{simulateFEISTY}} Run FEISTY model simulations
#' 
#' @export
#'

plotNetwork <- function(sim) {
  p <- sim$p
  u <- sim$u
  
  # Number of groups and biomass:
  ngroup <- p$ix[[length(p$ix)]][length(p$ix[[length(p$ix)]])] #resources (4) + fish 
  biomass <- u
  
  #Average of the biomass : 
  Bi <- colMeans(biomass[round(0.6*nrow(biomass), digits = 0):nrow(biomass),]) # mean value of the last 40% time 
  
  if (p$setup == "setupBasic"){
    
    # Set artificial depths to offset bubbles:
    Av_depth <- c(-1,-1,-4,-4,0,0,-2,-2,-2,-3,-3,-3)
    
    p$SpId <- c('smallPel','largePel', 'demersals')
    SpId <- c("smallZoo", "largeZoo", "benthos", "largeBenthos", 
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])))
    
    # Specify depth axis text:
    yaxis <- c("", "") # Leave blank because there is no depth in setup vertical
    
    # Set artificial depth:
    p$bottom <- -(min(Av_depth)) + 1
  }  
  
  if (p$setup == "setupBasic2"){
    
    # Set artificial depths to offset bubbles:
    Av_depth <- c(-1,-1,-4,-4,rep(0, length(p$ix[[1]])),rep(-2, length(p$ix[[2]])),rep(-3, length(p$ix[[3]])))
    
    p$SpId <- c('smallPel','largePel', 'demersals')
    SpId <- c("smallZoo", "largeZoo", "benthos", "largeBenthos", 
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])))
    
    # Specify depth axis text:
    yaxis <- c("", "") # Leave blank because there is no depth in setup vertical
    
    # Set artificial depth:
    p$bottom <- -(min(Av_depth)) + 1
  }  
  
  if (p$setup == "setupVertical" | p$setup == "setupVertical2"){
    
    #Calculate average depth day/night
    
    Av_depth_day <- 1 : p$nStages
    Av_depth_night <- 1 : p$nStages
    for (i in 1:p$nStages) {
      Av_depth_day[i] <- which.max(p$depthDay[ ,i])
      Av_depth_night[i] <- which.max(p$depthNight[ ,i]) 
    }
    
    Av_depth <- -(Av_depth_day + Av_depth_night) / 2
    
    # Change a bit for visualization:
    Av_depth[p$ix[[1]][1]:p$ix[[1]][length(p$ix[[1]])]] <- Av_depth[p$ix[[1]][1]:p$ix[[1]][length(p$ix[[1]])]] + 0.1 * p$bottom
    Av_depth[p$ix[[3]][1]:p$ix[[3]][length(p$ix[[3]])]] <- Av_depth[p$ix[[3]][1]:p$ix[[3]][length(p$ix[[3]])]] - 0.1 * p$bottom
    Av_depth[p$ix[[4]][1]:p$ix[[4]][length(p$ix[[4]])]] <- Av_depth[p$ix[[4]][1]:p$ix[[4]][length(p$ix[[4]])]] - 0.1 * p$bottom
    Av_depth[p$ix[[2]][1]:p$ix[[2]][length(p$ix[[2]])]] <- Av_depth[p$ix[[2]][1]:p$ix[[2]][length(p$ix[[2]])]] - 0.2 * p$bottom
    
    # Set color palette: 
    p$SpId <- c('smallPel','mesoPel','largePel', 'bathyPel', 'demersals')
    SpId <- c("smallZoo", "largeZoo", "benthos", "largeBenthos", 
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])),
              rep(p$SpId[4], length(p$ix[[4]])),
              rep(p$SpId[5], length(p$ix[[5]])))
    
    # Specify depth axis text:
    yaxis <- c("Surface   ", "     Bottom")
  }
  
  # Marker size depends on biomass following a cubic square transformation
  Msize <- Bi / max(Bi)
  Msize[Msize == 0] <- NA
  Msize <- Msize^(1/3)
  
  # Create line width: 
  Flux <- getFeeding(sim)
  Flux <- c(Flux) 
  threshold <- 0.05 
  indx <- which(Flux >= threshold) # takes the x highest values of theta
  
  # Set values of each coordinate (i.e. size and water column position) and put together:
  coord_1 <- data.frame(index = 1:p$nStages^2,
                        mc = rep(p$mc[1:p$nStages], p$nStages), 
                        depth = rep(Av_depth[1:p$nStages], p$nStages), 
                        SpId = rep(SpId, p$nStages),
                        Msize = rep(Msize, p$nStages), 
                        LineWdth = (Flux/max(Flux))^(1/3),
                        Alpha = (Flux/max(Flux))^(1/3))
  
  coord_2 <- data.frame(index = 1:p$nStages^2, # Notice repetition of ys grouped by "each" to change order
                        mc = rep(p$mc[1:p$nStages], each = p$nStages), 
                        depth = rep(Av_depth[1:p$nStages], each = p$nStages), 
                        SpId = rep(SpId, each = p$nStages),
                        Msize = rep(Msize, each = p$nStages),
                        LineWdth = (Flux/max(Flux))^(1/3),
                        Alpha = (Flux/max(Flux))^(1/3))
  
  # Combine in a data frame:
  df <- rbind(coord_1, coord_2)
  df <- df[order(-df$Msize),]   
  df <- subset(df,!(is.na(df$Msize))) 
  df2 <- subset(df,df$index %in% indx)
  df2 <- df2[order(-df2$Msize),]
  
  # Generate plot:
  plot <- defaultplot() +
    geom_line(data = df2, aes(x = mc, y = depth, group = index, color = SpId, alpha = Alpha),
              show.legend = F, linewidth = df2$LineWdth*3) +
    geom_point(data = df, aes(x = mc, y = depth, color = SpId, size = Msize), stroke = 0, shape = 16) +
    scale_color_manual(values = p$my_palette[attr(p$my_palette, "names") %in% df$SpId], 
                       labels = p$my_names[attr(p$my_palette, "names") %in% df$SpId]) +
    scale_radius(limits = c(0, NA), range = c(0, 10)) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_continuous(breaks = seq(0, round(-p$bottom - 1), by = -p$bottom), labels = yaxis) +
    annotation_logticks(sides = "b",size = 0.2,colour = "darkgrey") +
    labs(x ="Weight (gWW)", y = "", color = "Group") +
    guides(size = "none") +
    theme(legend.position = "bottom",legend.key = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(angle = 90, hjust = .5, margin = margin(r = 0), size=12),
          axis.ticks.y = element_blank())
  
  return(plot)
}

#' Diet plot
#' 
#' Makes a plot of prey groups in the stomach of each fish functional type and size class. 
#' The function only works on the four prepared setups or revised versions 
#' based on these four setups. If customized setups by users have more resources 
#' and/or fish functional types, this plot function needs to be updated.
#' 
#' @details
#' To be added.
#' 
#' @author Daniel Ottmann Riera, P. Daniël van Denderen, Yixin Zhao
#'
#' @usage plotDiet(sim)
#' 
#' @param sim The data frame of FEISTY simulation results.
#' 
#' @examples 
#' sim=simulateFEISTY()
#' plotDiet(sim)
#' 
#' @aliases plotDiet
#' 
#' @seealso 
#' \code{\link{simulateFEISTY}} Run FEISTY model simulations
#' 
#' @export
#'

plotDiet <- function(sim) {
  p <- sim$p
  u <- sim$u
  p$nstage <-lengths <- max(sapply(p$ix, length)) #maximum number of stages for one group
  fish <- c((p$nResources+1):(p$nResources+p$nGroups)) # get the fish
  biomass <- u 
  Bin <- round(0.6 * nrow(biomass), digits = 0) 
  biomassend   <- colMeans(biomass[Bin:nrow(biomass),]) # mean value of the last 40% time 
  biomassstage <- p$ixFish[length(p$ixFish)]
  biomasssmall <- p$nstage - round(2/3*p$nstage, digits = 0)
  Enc = p$V * (p$theta %*% biomassend)
  f   = Enc / (p$Cmax + Enc)
  f[is.na(f)] = 0  
  
  bom <- t(t(p$theta[5:biomassstage, ]) * colMeans(biomass[Bin:nrow(biomass),])) 
  fbom <- f[5:biomassstage] / rowSums(bom)
  output <- bom * fbom

  if(length(p$ix)==5){
    
    p$SpId <- c('smallPel','mesoPel','largePel', 'bathyPel', 'demersals')
    SpId <- c("smallZoo", "largeZoo", "benthos", "largeBenthos", 
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])),
              rep(p$SpId[4], length(p$ix[[4]])),
              rep(p$SpId[5], length(p$ix[[5]])))
    
    p$RSpName <- c("Small mesozooplankton", "Large mesozooplankton", "Benthos", "Small pelagics",
                   "Mesopelagics", "Large pelagics", "Bathypelagics", "Demersals")
    
  } else {
    
    p$SpId <- c('smallPel','largePel', 'demersals')
    SpId <- c("smallZoo", "largeZoo", "benthos", "largeBenthos", 
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])))
    
    p$RSpName <- c("Small mesozooplankton", "Large mesozooplankton", "Benthos", "Small pelagics",
                   "Large pelagics", "Demersals")   
    
  }
  
  # small pelagics: ----
  small_pel <- output[(p$ix[[1]][1] - length(p$ixR)):(p$ix[[1]][length(p$ix[[1]])] - length(p$ixR)), ] 
  small_pel <- t(rbind(small_pel, matrix(0, biomasssmall, biomassstage)))
  small_pel <- data.frame(val = c(small_pel), 
                          stage = rep(1:p$nstage, each = nrow(small_pel)), 
                          SpId = as.factor(rep(SpId, p$nstage)))
  
  p1 <- defaultplot_nl() + 
    geom_bar(data=small_pel, aes(x = stage, y = val, fill = SpId), stat = "identity") +
    scale_fill_manual(values = p$my_palette) + 
    scale_x_continuous(breaks = seq(0,p$nstage,length.out=4)) +
    labs(x = "", y = "Fraction in stomach", fill = "Prey group", title = p$RSpName[4])
  
  if (length(p$ix)==5){
    
    # Large pelagics: ----
    large_pel <- t(output[(p$ix[[3]][1] - length(p$ixR)):(p$ix[[3]][length(p$ix[[3]])] - length(p$ixR)), ])
    large_pel <- data.frame(val = c(large_pel), 
                            stage = rep(1:p$nstage, each = nrow(large_pel)), 
                            SpId = as.factor(rep(SpId, p$nstage)))
    
    p2 <- defaultplot_nl() +
      geom_bar(data = large_pel, aes(x = stage, y = val, fill = SpId), stat = "identity") +
      scale_fill_manual(values = p$my_palette) +
      scale_x_continuous(breaks = seq(0,p$nstage,length.out=4)) +
      labs(x ="size-class", y = "", fill = "Prey group", title = p$RSpName[6])
    
    
    # Demersal: ----
    demers <- t(output[(p$ix[[5]][1] - length(p$ixR)):(p$ix[[5]][length(p$ix[[5]])] - length(p$ixR)), ])
    demers <- data.frame(val = c(demers), 
                         stage = rep(1:p$nstage, each = nrow(demers)), 
                         SpId = as.factor(rep(SpId, p$nstage)))
    
    p3 <- defaultplot_nl() +
      geom_bar(data = demers, aes(x = stage, y = val, fill = SpId), stat = "identity") +
      scale_fill_manual(values = p$my_palette) +
      scale_x_continuous(breaks = seq(0,p$nstage,length.out=4)) +
      labs(x ="size-class", y = "", fill = "Prey group", title = p$RSpName[8])
    
    
    # Mesopelagics: ----
    meso_pel <- output[(p$ix[[2]][1] - length(p$ixR)):(p$ix[[2]][length(p$ix[[2]])] - length(p$ixR)), ]
    meso_pel <- t(rbind(meso_pel, matrix(0, biomasssmall, biomassstage)))
    meso_pel <- data.frame(val = c(meso_pel), 
                           stage = rep(1:p$nstage, each = nrow(meso_pel)), 
                           SpId = as.factor(rep(SpId, p$nstage)))
    
    p5 <- defaultplot_nl() +
      geom_bar(data = meso_pel, aes(x = stage, y = val, fill = SpId), stat = "identity") +
      scale_fill_manual(values = p$my_palette) +
      scale_x_continuous(breaks = seq(0,p$nstage,length.out=4)) +
      labs(x ="size-class", y = "Fraction in stomach", fill = "Prey group", title = p$RSpName[5])
    
    # Bathypelagics: ----
    bathy_pel <- t(output[(p$ix[[4]][1] - length(p$ixR)):(p$ix[[4]][length(p$ix[[4]])] - length(p$ixR)), ])
    bathy_pel <- data.frame(val = c(bathy_pel), 
                            stage = rep(1:p$nstage, each = nrow(bathy_pel)), 
                            SpId = as.factor(rep(SpId, p$nstage)))
    
    p6 <- defaultplot_nl() +
      geom_bar(data = bathy_pel, aes(x = stage, y = val, fill = SpId), stat = "identity") +
      scale_fill_manual(values = p$my_palette) +
      scale_x_continuous(breaks = seq(0,p$nstage,length.out=4)) +
      labs(x ="size-class", y = "", fill = "Prey group", title = p$RSpName[7])
    
  } else {
    
  # Large pelagics: ----
    large_pel <- t(output[(p$ix[[2]][1] - length(p$ixR)):(p$ix[[2]][length(p$ix[[2]])] - length(p$ixR)), ])
    large_pel <- data.frame(val = c(large_pel), 
                            stage = rep(1:p$nstage, each = nrow(large_pel)), 
                            SpId = as.factor(rep(SpId, p$nstage)))
    
    p2 <- defaultplot_nl() +
      geom_bar(data = large_pel,aes(x = stage, y = val, fill = SpId), stat = "identity") +
      scale_fill_manual(values = p$my_palette) +
      scale_x_continuous(breaks = seq(0,p$nstage,length.out=4)) +
      labs(x ="size-class", y = "", fill = "Prey group", title = p$RSpName[5])

    # Demersal: ----
    demers <- t(output[(p$ix[[3]][1] - length(p$ixR)):(p$ix[[3]][length(p$ix[[3]])] - length(p$ixR)), ])
    demers <- data.frame(val = c(demers), 
                         stage = rep(1:p$nstage, each = nrow(demers)), 
                         SpId = as.factor(rep(SpId, p$nstage)))
    
    p3 <- defaultplot_nl() +
      geom_bar(data = demers, aes(x = stage, y = val, fill = SpId), stat = "identity") +
      scale_fill_manual(values = p$my_palette) +
      scale_x_continuous(breaks = seq(0,p$nstage,length.out=4)) +
      labs(x ="size-class", y = "", fill = "Prey group", title = p$RSpName[6])
    } 
  
  # Legend: ----
  legend_colors <- data.frame(val = 0, stage = 0, SpId = unique(large_pel$SpId))
  scale_color_manual(name = "Groups", 
                     values = c("total" = "black", p$my_palette[p$groupnames[fish]]),
                     breaks = c("total",p$groupnames[fish]),
                     labels = c("total"= "Total", p$my_names[p$groupnames[fish]]))
  
  
  p7 <- ggplot(data = legend_colors) +
    geom_bar(aes(x = stage, y = val, fill = SpId), stat = "identity") +
    scale_fill_manual(name = "Prey group",
                      values = p$my_palette[attr(p$my_palette, "names") %in% legend_colors$SpId],
                      breaks = attr(p$my_palette, "names"),
                      labels = attr(p$my_palette, "names"))+
    theme_void() +
    theme(legend.position = c(.45,.4),
          legend.direction = "horizontal",
          legend.text = element_text(size = 12))+
    guides(fill = guide_legend(
      nrow = 3,  
      byrow = TRUE,  
      title.position = "top",  
    ))
  
  # Put all panels together:
  if (length(p$ix)==5){
    if(p$bottom > p$mesop) {
      plots <- plot_grid(p1,p2,p3,p5,p6,p7,nrow=2)
    }else{
      plots <- plot_grid(p1,p2,p3,p7,nrow=2)
    }
  } else {
    plots <- plot_grid(p1,p2,p3,p7,nrow=2)
  }
  return(plots)
}

#' Plots for growth rate, mortality, and feeding level 
#' 
#' Plots growth rate (1/year), mortality (1/year), and feeding level (dimensionless) over the size spectrum.
#' 
#' @details
#' Rate calculation:
#' \itemize{
#' \item Averaged values of rate data are estimated over a selected time period 
#' (default last 40% simulation time).
#' }
#' 
#' @author P. Daniël van Denderen, Yixin Zhao
#'
#' @usage plotRates(sim)
#' 
#' @param sim The data frame of FEISTY simulation results.
#' 
#' @examples 
#' sim=simulateFEISTY()
#' # averaged rate data of last 40%  simulation time
#' plotRates(sim=sim)
#' 
#' @aliases plotRates
#' 
#' @seealso 
#' \code{\link{simulateFEISTY}} Run FEISTY model simulations
#' 
#' @export
#' 

plotRates = function(sim) {
  p_rates <- getRates(sim)
  plots <- align_plots(p_rates[[1]],p_rates[[2]],p_rates[[3]], align = 'v')
  plots <- plot_grid(plots[[1]], plots[[2]], plots[[3]],ncol=1,rel_heights = c(1,1,1.6))
  return(plots)
}

#' Plot simulation results
#' 
#' Make a plot combo for the simulation results, including \code{\link{plotSSBtime}}, 
#' \code{\link{plotSpectra}}, \code{\link{plotNetwork}} and \code{\link{plotRates}}.
#' 
#' @details This function is designed to give a quick visualization of a simulation output. 
#' 
#' @author P. Daniël van Denderen, Rémy Denéchère, Yixin Zhao
#'
#' @usage plotSimulation(sim)
#' 
#' @param sim The dataframe of FEISTY simulation results.
#' 
#' @examples 
#' sim=simulateFEISTY()
#' plotSimulation(sim)
#' 
#' @aliases plotSimulation
#' 
#' @seealso 
#' \code{\link{webFEISTY}} A shiny interface for visualizing FEISTY model results \cr
#' \code{\link{simulateFEISTY}} Run FEISTY model simulations \cr
#' \code{\link{plotSSBtime}} Timeseries of spawning stock biomass \cr
#' \code{\link{plotSpectra}} Biomass spectra plot \cr
#' \code{\link{plotRates}} Plots for growth rate, mortality, and feeding level \cr
#' \code{\link{plotNetwork}} Plot with all feeding interactions
#' 
#' @export
#' 

plotSimulation = function(sim) {
  p_biomasstime     <- plotBiomasstime(sim)
  p_rates   <- getRates(sim)
  p_spectra <- plotSpectra(sim)
  p_network <- plotNetwork(sim)
  
  ## get legend from plotNetwork
  legend = get_legend_new(p_network)

  ## Turn legend off for later external plotting 
  # subpanel a
  
  # create a fake line for legend 
  text_loc <- data.frame(x=c(1,2,3,4,5,6),y=rep(1E-16,6),
                         group=c("total"))
  p_spectra     = p_spectra + labs(title = "a")+
    #fake line
    geom_line(data=text_loc,aes(x=x,y=y,linetype=group),linewidth=lwd_def) + 
    scale_linetype_manual(values=c(1),name = "") + guides(color = "none") +
    theme(legend.position=c(1,1),legend.justification=c("right"),
          legend.background = element_rect(colour = NA, fill = NA),
          axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          legend.spacing.x = unit(.1, 'cm'),
          legend.key = element_blank())
  
  # subpanel b-d
  p_rates[[1]]  = p_rates[[1]] + labs(title = "b")
  p_rates[[2]]  = p_rates[[2]] + labs(title = "c")
  p_rates[[3]]  = p_rates[[3]] + labs(title = "d") + 
                   theme(legend.position="none",
                        plot.background = element_rect(color = NA))
  
  # subpanel e
  p_network = p_network + labs(title = "e") + 
               theme(legend.position="none",
                    plot.background = element_rect(color = NA))
  
  # subpanel f - set new scale to make y-scale "cleaner"
  p_biomasstime     = p_biomasstime + labs(title = "f") +
               theme(legend.position="none",
                    plot.background = element_rect(color = NA)) #+
               # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=3), 
               #    labels = trans_format("log10", math_format(10^.x)))
              
  # create left and right and bottom 
  left  <- align_plots(p_spectra,p_rates[[1]],p_rates[[2]],p_rates[[3]],align="v")
  left  <- plot_grid(left[[1]],left[[2]],left[[3]],left[[4]],rel_heights = c(1,1.1,1.1,1.3),ncol=1)
  right <- align_plots(p_network,p_biomasstime,align = 'v')
  right <- plot_grid(right[[1]],right[[2]],ncol=1,rel_heights = c(3.2,1.3))
  plots <- plot_grid(left,right,nrow=1,rel_widths = c(1,1.4))  
  plots <- plot_grid(plots,legend,ncol=1,rel_heights = c(10,1))

  return(plots)
}

# Only for Shiny
#' @export
plotSimulationShiny = function(sim) {
  p_biomasstime   <- plotBiomasstime(sim)+theme(legend.position = "bottom")
  p_spectra       <- plotSpectra(sim)
  p_rates         <- getRates(sim)  
  
  ## get legend from plotNetwork
  legend = get_legend_new(p_biomasstime)
  
  ## Turn legend off for later external plotting 
  # subpanel a
  
  # create a fake line for legend 
  text_loc <- data.frame(x=c(1,2,3,4,5,6),y=rep(1E-16,6),
                         group=c("total"))
  p_spectra     = p_spectra + 
    #fake line
    geom_line(data=text_loc,aes(x=x,y=y,linetype=group),linewidth=lwd_def) + 
    scale_linetype_manual(values=c(1),name = "") + guides(color = "none") +
    theme(legend.position=c(1,1),legend.justification=c("right"),
          legend.background = element_rect(colour = NA, fill = NA),
          axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          legend.spacing.x = unit(.1, 'cm'),
          legend.key = element_blank())
  
  # subpanel b-d
  #p_rates[[1]]  = p_rates[[1]] + labs(title = "b")
  #p_rates[[2]]  = p_rates[[2]] + labs(title = "c")
  p_rates[[3]]  = p_rates[[3]]  + 
    theme(legend.position="none",
          plot.background = element_rect(color = NA))
  
  # subpanel f - set new scale to make y-scale "cleaner"
  p_biomasstime     = p_biomasstime  +
    theme(legend.position="none",
          plot.background = element_rect(color = NA)) #+
    # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=3), 
    #               labels = trans_format("log10", math_format(10^.x)))
  
  # create left and right and bottom 
  plots  <- align_plots(p_biomasstime,p_spectra,p_rates[[1]],p_rates[[2]],p_rates[[3]],align="v")
  plots  <- plot_grid(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],rel_heights = c(1.1,1,1,1,1),ncol=1)
  plots  <- plot_grid(plots,legend,ncol=1,rel_heights = c(10,1))
  
  return(plots)
}

# Support functions: -----
 
# creates empty ggplot as default
defaultplot = function() {
  ggplot() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill = NA))
}

# creates empty ggplot with no legend
defaultplot_nl = function() {
  defaultplot() + theme(legend.position = "none")
}

# creates dataframe for time series plots
getTimeseries = function(sim) {
  dat <- data.frame(t     = rep(sim$t,(sim$p$nResources + sim$p$nGroups)),
                    SSB   = c(as.vector(sim$R*NA), as.vector(sim$SSB)),
                    bio   = c(as.vector(sim$R), as.vector(sim$totBiomass)),
                    yield = c(as.vector(sim$R*0), as.vector(sim$yield)),
                    group = as.factor(rep(sim$p$groupnames, each = length(sim$t))))
  return(dat)
}

# creates feeding flux from prey to predator
getFeeding = function(sim) {
  p <- sim$p
  u <- sim$u
  
  # get last 40% of timeseries
  etaTime <- 0.4 
  ixTime  <- which(sim$t>=((1-etaTime)*sim$t[sim$nTime]))
  
  # Number of groups and biomass:
  ngroup <- p$ix[[length(p$ix)]][length(p$ix[[length(p$ix)]])] #resources (4) + fish 
  biomass <- u
  
  #Average of the biomass : 
  prey     <- colMeans(biomass[ixTime,]) # mean value of the last 40% time 
  predator <- colMeans(biomass[ixTime,]) # mean value of the last 40% time 
  
  # estimate encounters
  Enc <- p$V[1] * p$theta[1,] * prey 
  for (i in 2:length(predator)){
    Enc <- rbind(Enc, p$V[i] * p$theta[i,] * prey)
  }
  
  # estimate encounters per predator
  Encspecies = rowSums(Enc)
  
  # estimate the mortality generated
  mortpr <-  p$Cmax[1] * p$V[1] * p$theta[1,] / (p$Cmax[1]+ Encspecies[1]) * predator[1]
  for (i in 2:length(predator)){
    mortpr <- rbind(mortpr,p$Cmax[i] * p$V[i] * p$theta[i,] / (p$Cmax[i]+ Encspecies[i]) * predator[i])
  }
  
  # estimate the flux from prey to predator
  mortpr <- ifelse(is.na(mortpr),0,mortpr)
  flux_prey_to_pred <- t(t(mortpr)*prey)
  rownames(flux_prey_to_pred) <- paste("pred",colnames(p$theta),sep="_")
  colnames(flux_prey_to_pred) <- paste("prey",colnames(p$theta),sep="_")
  return(flux_prey_to_pred)
}

# creates individual outputs for combined rates plot
getRates = function(sim) {
  p <- sim$p
  
  # get last 40% of timeseries
  etaTime <- 0.4 
  ixTime  <- which(sim$t>=((1-etaTime)*sim$t[sim$nTime]))

    # estimate rates
  fish <- c((p$nResources+1):(p$nResources+p$nGroups))
  rates <- data.frame(group    = rep(p$groupnames[fish],sapply(p$ix,length)),
                      mc       = p$mc[p$ixFish],
                      g        = as.numeric(colMeans(sim$g[ixTime,])),
                      mortpred = as.numeric(colMeans(sim$mortpred[ixTime,p$ixFish])),
                      mortF    = as.numeric(p$mortF[p$ixFish]),
                      f        = as.numeric(colMeans(sim$f[ixTime,])))
  x_lim  <- range(p$mc[p$ixFish])
  
# extract existing fish types according to initial values
  inival=p$u0[-p$ixR][grep("_1", names(p$u0[-p$ixR]))]
  names(inival)=p$groupnames[fish]
  fishexistname=names(inival[inival!=0])
  
  rates <- subset(rates, rates$group %in% fishexistname)
  
# convert zeros of data to a very small number (only for plotting purposes)
  rates$g <- ifelse(rates$g == 0, 1E-16,rates$g)
  rates$mortpred <- ifelse(rates$mortpred == 0, 1E-16,rates$mortpred) 
  rates$mortF <- ifelse(rates$mortF == 0, 1E-16,rates$mortF) 
  rates$f <- ifelse(rates$f == 0, 1E-16,rates$f) 
  
  # growth rates: ----
  growth <- defaultplot_nl() +
    geom_line(data = rates, 
              aes(x = mc, y = g, group = group, color = group),
              linewidth=lwd_def) + 
    scale_color_manual(name = "Groups", 
                       values = p$my_palette[p$groupnames[fish]],
                       breaks = p$groupnames[fish],
                       labels = p$my_names[p$groupnames[fish]]) +
    coord_cartesian(ylim = c(min_g,max(min_g*100,max(rates$g)*5)),xlim=x_lim) + 
    annotation_logticks(sides = "bl",size = 0.2,colour = "darkgrey") +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))+
    xlab("") + ylab(expression("Growth rate (yr"^"-1"*")")) + 
    theme(axis.text.x=element_blank())
  
  # mortality rates: ----
  
  # add background mortality
  mort0_line  <- data.frame(mort0 = rep(p$mort0[p$ixFish][1],2), mc = range(p$mc[p$ixFish]))
  
  # convert zeros of fishing to a very small number (only for plotting purposes)
  rates$mortF <- ifelse(rates$mortF == 0, 1E-16,rates$mortF) 
  
  # add labels for legend 
  text_loc <- data.frame(x=c(1,2,3,4,5,6),y=rep(1E-5,6),
                         group=rep(c("backgr.","predation","fishing"),2))
  
  mort <- defaultplot() +
    geom_line(data = rates, 
              aes(x = mc, y = mortpred, group = group, color = group),
              linewidth=lwd_def) + 
    geom_line(data = rates, 
              aes(x = mc, y = mortF, group = group, color = group),
              linewidth=lwd_def,linetype="dashed") +
    scale_color_manual(name = "Groups", 
                       values = p$my_palette[p$groupnames[fish]],
                       breaks = p$groupnames[fish],
                       labels = p$my_names[p$groupnames[fish]]) +
    coord_cartesian(ylim = c(min_mort,max(min_mort*100,max(rates$mortpred)*5)),xlim=x_lim) + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))+
    annotation_logticks(sides = "bl",size = 0.2,colour = "darkgrey") + 
    xlab("") + ylab(expression("Mortality (yr"^"-1"*")")) +
    geom_line(data = mort0_line,
              aes(x = mc, y = mort0, color = "darkgrey"),
              linewidth=1,linetype = "dotted") +
    geom_line(data=text_loc,aes(x=x,y=y,linetype=group)) + 
    scale_linetype_manual(values=c(3,2,1),name = "") + guides(color = "none") +
    theme(legend.position=c(1,1),legend.justification=c("right"),
          legend.background = element_rect(colour = NA, fill = NA),
          legend.key = element_blank(),
          axis.text.x=element_blank(),legend.spacing.x = unit(.1, 'cm')) +
    guides(linetype = guide_legend(nrow = 1))
  
  # feeding levels: ----
  
  # estimate critical feeding level
  i     <- which.max(sapply(p$ix,"length"))
  ix    <- p$ix[[i]]
  fcrit <- data.frame(fc = p$metabolism[ix]/(p$epsAssim*p$Cmax[ix]), # Critical feeding level
                      mc = p$mc[ix])
  
  flevel <- defaultplot() +
    geom_line(data = rates, 
              aes(x = mc, y = f, group = group, color = group),
              linewidth=lwd_def)  +
    scale_color_manual(name = "Groups", 
                       values = p$my_palette[p$groupnames[fish]],
                       breaks = p$groupnames[fish],
                       labels = p$my_names[p$groupnames[fish]]) +
    coord_cartesian(ylim = c(0,1),xlim=x_lim) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "b",size = 0.2,colour = "darkgrey") + 
    xlab("Weight (gWW)") + ylab("Feeding level, f") +
    geom_line(data = fcrit,
              aes(x = mc, y = fc, color = "darkgrey"),
              linewidth=1,linetype = "dotted") + 
    theme(legend.position = "bottom",legend.key = element_blank()) +  
    guides(color = guide_legend(nrow = 2)) +
    annotate("text", x=fcrit$mc[2],y=fcrit$fc[1]-0.02,label="crit. f",
             color = "darkgrey")
  
  return(list(growth,mort,flevel)) 
}  

# Replace get_legend() from 'cowplot' package for temporary, 
# because get_legend() does not work with 'ggplot2' 3.5.0. Check the new version later.
get_legend_new <- function(plot, legend = NULL) {
  
  gt <- ggplotGrob(plot)
  
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  
  indices <- grep(pattern, gt$layout$name)
  
  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}