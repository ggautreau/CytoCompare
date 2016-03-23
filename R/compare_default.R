# title Internal - Comparison of two CELL, CLUSTER or gate profiles
#
# @description This function is a wrapper for the default comparison methods
#
# @param profile1.type a character specifying the type of the first profile (CELL, CLUSTER or GATE)
# @param profile2.type a character specifying the type of the second profile (CELL, CLUSTER or GATE)
# @param profile1.intensities a numeric vector containing the marker expressions for the first cell profile
# @param profile1.mean a numeric vector containing the marker expression means of the first cluster profile
# @param profile1.sd a numeric vector containing the marker expression standard deviations of the first cluster profile
# @param profile1.density a named list containing the marker expression densities (DENSITY objects) of the first cluster profile
# @param profile1.nbcells a numeric value containing the number of cells associated with the first cluster profile
# @param profile1.range a numeric matrix containing the marker expression ranges of the second gate profile
# @param profile2.intensities a numeric vector containing the marker expressions for the second cell profile
# @param profile2.mean a numeric vector containing the marker expression means of the second cluster profile
# @param profile2.sd a numeric vector containing the marker expression standard deviations of the second cluster profile
# @param profile2.density a named list containing the marker expression densities (DENSITY objects) of the second cluster profile
# @param profile2.nbcells a numeric value containing the number of cells associated with the second cluster profile
# @param profile2.range a numeric matrix containing the marker expression ranges of the second gate profile
# @param ... other parameters for the 'compare_default' methods
#
# @return a named list containing the marker measures, marker successes, the aggregated measure and the aggregated p-value
compare_default <- function(profile1.type,profile2.type,
                            profile1.intensities,profile1.mean,profile1.sd,profile1.density,profile1.nbcells,profile1.range,
                            profile2.intensities,profile2.mean,profile2.sd,profile2.density,profile2.nbcells,profile2.range,
                            ...){  
    if(profile1.type=="CELL" && profile2.type=="CELL"){
        res <- compare_default_CELL_CELL(profile1.intensities,profile2.intensities,...)
    }else if(profile1.type=="CLUSTER" && profile2.type=="CLUSTER"){
        res <- compare_default_CLUSTER_CLUSTER(profile1.mean,profile1.sd,profile1.density,profile1.nbcells,profile2.mean,profile2.sd,profile2.density,profile2.nbcells,...)
    }else if(profile1.type=="GATE" && profile2.type=="GATE"){
        res <- compare_default_GATE_GATE(profile1.range,profile2.range,...)
    }else if(profile1.type=="CELL" && profile2.type=="GATE"){
        res <- compare_default_CELL_GATE(profile1.intensities,profile2.range,...)
    }else if(profile1.type=="GATE" && profile2.type=="CELL"){
        res <- compare_default_CELL_GATE(profile2.intensities,profile1.range,...)
    }else if(profile1.type=="CELL" && profile2.type=="CLUSTER"){
        res <- compare_default_CELL_CLUSTER(profile1.intensities,profile2.mean,profile2.sd,profile2.density,profile2.nbcells,...)
    }else if(profile1.type=="CLUSTER" && profile2.type=="CELL"){
        res <- compare_default_CELL_CLUSTER(profile2.intensities,profile1.mean,profile1.sd,profile1.density,profile1.nbcells,...)
    }else if(profile1.type=="CLUSTER" && profile2.type=="GATE"){
        res <- compare_default_CLUSTER_GATE(profile1.mean,profile1.sd,profile1.density,profile1.nbcells,profile2.range,...)
    }else if(profile1.type=="GATE" && profile2.type=="CLUSTER"){
        res <- compare_default_CLUSTER_GATE(profile2.mean,profile2.sd,profile2.density,profile2.nbcells,profile1.range,...)
    }
    
    return(res)
}

# title Internal - Comparison of two cell profiles
#
# @description Performs a similarity comparison between two cell profiles.
#
# @details This function computes the marker similarity measures as the absolute difference betweens the expression markers of the two cell profiles.
#
# @param profile1.intensities a numeric vector containing the marker expressions for the first cell profile
# @param profile2.intensities a numeric vector containing the marker expressions for the second cell profile
# @param weights a numeric vector containing the marker weights
# @param success.th a numeric value specifying a threshold below which a marker similarity measure will be considered as a similarity success or fail
# @param success.p a numeric value specifying the expected proportion of marker similarity successes
#
# @return a named list containing the marker similarity measures, marker similarity successes, the aggregated similarity measure and the aggregated similarity p-value
compare_default_CELL_CELL <- function(profile1.intensities,
                                        profile2.intensities,
                                        weights,
                                        success.th=0.20,success.p=0.75){
                              
    measures   <- abs(profile1.intensities - profile2.intensities)
    
    successes  <- (measures <= success.th)
    
    measure    <- weighted.mean(measures,weights,na.rm = TRUE)
    pvalue     <- aggreg.binom(successes,weights,success.p)
    
    res <- list(measure  = measure,
        pvalue           = pvalue,
        marker.measures  = measures,
        marker.successes = successes)
    return(res)
}


# title Internal - Comparison of two cluster profiles
#
# @description Performs a similarity comparison between two cluster profiles.
#
# @details This function computes the marker similarity measures as the Kolmogorov-Smirnov statistics between the marker expression densities of the two cell cluster profiles. 
#
# @param profile1.mean a numeric vector containing the marker expression means of the first cluster profile
# @param profile1.sd a numeric vector containing the marker expression standard deviations of the first cluster profile
# @param profile1.density a named list containing the marker expression densities (DENSITY objects) of the first cluster profile
# @param profile1.nbcells a numeric value containing the number of cells in the first cluster profile
# @param profile2.mean a numeric vector containing the marker expression means of the second cluster profile
# @param profile2.sd a numeric vector containing the marker expression standard deviations of the second cluster profile
# @param profile2.density a named list containing the marker expression densities (DENSITY objects) of the second cluster profile
# @param profile2.nbcells a numeric value containing the number of cells in the second cluster profile
# @param weights a numeric vector containing the marker weights
# @param success.th a numeric value specifying a threshold below which a marker similarity measure will be considered as a similarity success or fail
# @param success.p a numeric value specifying the expected proportion of marker similarity successes
# @param nbcells.th a numeric value specifying a threshold below which the marker expression density will be approximated by a normal distribution
#
# @return a named list containing the marker similarity measures, marker similarity successes, the aggregated similarity measure and the aggregated similarity p-value
compare_default_CLUSTER_CLUSTER <- function(profile1.mean,profile1.sd,profile1.density,profile1.nbcells,
                                            profile2.mean,profile2.sd,profile2.density,profile2.nbcells,
                                            weights,
                                            success.th=0.30,success.p=0.75,
                                            nbcells.th=50){
    
    if(missing(profile1.density) || profile1.nbcells <= nbcells.th){
        profile1.fun.cdf <- function(x,i){ pnorm(x,profile1.mean[i],profile1.sd[i]) }
        profile1.fun.q   <- function(x,i){ qnorm(x,profile1.mean[i],profile1.sd[i]) }
    }else{
        profile1.fun.cdf <- function(x,i){ cdf.density(profile1.density[[i]],x) }
        profile1.fun.q   <- function(x,i){ quantiles.density(profile1.density[[i]],x) }
    }
    if(missing(profile2.density) || profile2.nbcells <= nbcells.th){
        profile2.fun.cdf <- function(x,i){ pnorm(x,profile2.mean[i],profile2.sd[i]) }
        profile2.fun.q   <- function(x,i){ qnorm(x,profile2.mean[i],profile2.sd[i]) }
    }else{
        profile2.fun.cdf <- function(x,i){ cdf.density(profile2.density[[i]],x) }
        profile2.fun.q   <- function(x,i){ quantiles.density(profile2.density[[i]],x) }
    }
    
    if(profile1.nbcells>=1 && profile2.nbcells>=1){
        measures   <- KS.stat(profile1.fun.cdf = profile1.fun.cdf,
                            profile1.fun.q     = profile1.fun.q,
                            profile2.fun.cdf   = profile2.fun.cdf,
                            profile2.fun.q     = profile2.fun.q,
                            markers.nb         = length(weights))
    }else{
        measures   <- rep(Inf,length(weights))
    }
    
    successes <- (measures <= success.th)
    
    measure   <- weighted.mean(measures,weights,na.rm = TRUE)
    pvalue    <- aggreg.binom(successes,weights,success.p)
    
    res <- list(measure  = measure,
        pvalue           = pvalue,
        marker.measures  = measures,
        marker.successes = successes)
    return(res)
}


# title Internal - Comparison of two gate profiles
#
# @description Performs a similarity comparison between two gate profiles.
#
# @details This function computes the marker similarity measures as the Kolmogorov-Smirnov statistics between the marker expression densities of the two gate profiles.
#
# @param profile1.range a numeric matrix containing the marker expression ranges of the first gate profile
# @param profile2.range a numeric matrix containing the marker expression ranges of the second gate profile
# @param weights a numeric vector containing the marker weights
# @param success.th a numeric value specifying a threshold below which a marker similarity measure will be considered as a similarity success or fail
# @param success.p a numeric value specifying the expected proportion of marker similarity successes
#
# @return a named list containing the marker similarity measures, marker similarity successes, the aggregated similarity measure and the aggregated similarity p-value
compare_default_GATE_GATE <- function(profile1.range,
                                profile2.range,
                                weights,
                                success.th=0.30,success.p=0.75){
    
    profile1.fun.cdf <- function(x,i){ cdf.uniform(profile1.range[i,],x) }
    profile1.fun.q   <- function(x,i){ quantiles.uniform(profile1.range[i,],x) }
    
    profile2.fun.cdf <- function(x,i){ cdf.uniform(profile2.range[i,],x) }
    profile2.fun.q   <- function(x,i){ quantiles.uniform(profile2.range[i,],x) }
    
    measures <- KS.stat(profile1.fun.cdf = profile1.fun.cdf,
                        profile1.fun.q   = profile1.fun.q,
                        profile2.fun.cdf = profile2.fun.cdf,
                        profile2.fun.q   = profile2.fun.q,
                        markers.nb       = length(weights))
    
    successes <- (measures <= success.th)
    
    measure   <- weighted.mean(measures,weights,na.rm = TRUE)
    pvalue    <- aggreg.binom(successes,weights,success.p)
    
    res <- list(measure  = measure,
        pvalue           = pvalue,
        marker.measures  = measures,
        marker.successes = successes)
    return(res)
}


# title Internal - Comparison of a cell profile with a gate profile
#
# @description Performs a comparison between a cell profile and a gate profile.
#
# @details This function computes the marker inclusion measures as the minimal distance between the cell expression and the gate ranges. The measure is equals to zero if the cell marker expression is included in the gate marker range.
#
# @param profile1.intensities a numeric vector containing the marker expressions for the cell profile
# @param profile2.range a numeric matrix containing the marker expression ranges for the gate profile
# @param weights a numeric vector containing the marker weights
# @param success.th a numeric value specifying a threshold below which a marker inclusion measure will be considered as an inclusion success
# @param success.p a numeric value specifying the expected proportion of marker inclusion successes
#
# @return a named list containing the marker inclusion measures, marker inclusion successes, the aggregated inclusion measure and the aggregated inclusion p-value
compare_default_CELL_GATE <- function(profile1.intensities,
                                profile2.range,
                                weights,
                                success.th=0.10,success.p=0.75){
                              
    measures         <- pmin(abs(profile1.intensities-profile2.range[,1]),abs(profile1.intensities-profile2.range[,2]))
    ingate           <- profile2.range[,1] <= profile1.intensities & profile1.intensities <= profile2.range[,2]
    measures[ingate] <- 0
    
    successes       <- (measures <= success.th)
    
    measure         <- weighted.mean(measures,weights,na.rm = TRUE)
    pvalue          <- aggreg.binom(successes,weights,success.p)
    

    res <- list(measure  = measure,
        pvalue           = pvalue,
        marker.measures  = measures,
        marker.successes = successes)
    return(res)
}


# title Internal - Comparison of a cell profile with a cluster profile
#
# @description Performs an inclusion comparison between a cell profile and a cluster profile.
#
# @details This function computes the marker inclusion measures as the minimal distance between the cell expression and the cell cluster ranges. The measure is equals to zero if the cell marker expression is included in the cell cluster marker range.
#
# The quantiles that define the cell cluster range can be defined using the 'cluster.quantiles' parameter.
#
# @param profile1.intensities a numeric vector containing the marker expression for the cell profile
# @param profile2.mean a numeric vector containing the marker expression means for the cluster profile
# @param profile2.sd a numeric vector containing the marker expression standard deviations for the cluster profile
# @param profile2.density a named list containing the marker expression densities (DENSITY objects) for the cluster profile
# @param profile2.nbcells a numeric value providing the number of cells associated with the cluster profile
# @param weights a numeric vector containing the marker weights
# @param success.th a numeric value specifying a threshold below which a marker inclusion measure will be considered as an inclusion success
# @param success.p a numeric value specifying the expected proportion of marker inclusion successes
# @param cluster.quantiles a numeric value indicating the percentile that define the expression ranges of the cell cluster profile
#
# @return a named list containing the marker inclusion measures, marker inclusion successes, the aggregated inclusion measure and the aggregated inclusion p-value
compare_default_CELL_CLUSTER <- function(profile1.intensities,
                                    profile2.mean,profile2.sd,profile2.density,profile2.nbcells,
                                    weights,
                                    success.th=0.10,success.p=0.75,
                                    cluster.quantiles=c(0.10,0.90)){
                                 
    if(missing(profile2.density))
        profile2.fun.q <- function(x,i) qnorm(x,profile2.mean[i],profile2.sd[i])
    else
        profile2.fun.q <- function(x,i) quantiles.density(profile2.density[[i]],x)
    
    if(profile2.nbcells>=1){
        measures <- rep(NA,length(weights))
        for(i in 1:length(weights)){
            cluster_low    <- profile2.fun.q(cluster.quantiles[1],i)[[1]]
            cluster_high   <- profile2.fun.q(cluster.quantiles[2],i)[[2]]
            measures[i]    <- max(abs(profile1.intensities[i]-cluster_low),abs(profile1.intensities[i]-cluster_high))
            if(cluster_low <= profile1.intensities[i] & profile1.intensities[i] <= cluster_high){
               measures[i] <- 0
            }
        }
    }else{
        measures  <- rep(Inf,length(weights))
    }
    
    successes    <- (measures <= success.th)
    
    measure      <- weighted.mean(measures,weights,na.rm = TRUE)
    pvalue       <- aggreg.binom(successes,weights,success.p)
    
    res <- list(measure  = measure,
        pvalue           = pvalue,
        marker.measures  = measures,
        marker.successes = successes)
    return(res)
}


# title Internal - Comparison of a cluster profile with a gate profile
#
# @description Performs an inclusion comparison between a cluster profile and a gate profile.
#
# @details This function computes the marker inclusion measures as the maximal absolute distance between the cluster and gate lower bounds and the cluster and gate upper bounds. The measure is equals to zero if the cell cluster marker expression is included in the gate marker range.
#
# @param profile1.mean a numeric vector containing the marker expression means for the cluster profile
# @param profile1.sd a numeric vector containing the marker expression standard deviation for the cluster profile
# @param profile1.density a named list containing the marker expression densities (DENSITY objects) for the cluster profile
# @param profile1.nbcells a numeric value providing the number of cells associated with the cluster profile
# @param profile2.range a numeric matrix containing the marker expression ranges of the gate profile
# @param weights a numeric vector containing the marker weights
# @param success.th a numeric value specifying a threshold below which a marker inclusion measure will be considered as an inclusion success
# @param success.p a numeric value specifying the expected proportion of marker inclusion successes
# @param cluster.quantiles a numeric value indicating the quantiles that define the marker expression ranges of the cell cluster
#
# @return a named list containing the marker inclusion measures, marker inclusion successes, the aggregated inclusion measure and the aggregated inclusion p-value
compare_default_CLUSTER_GATE <- function(profile1.mean,profile1.sd,profile1.density,profile1.nbcells,
                                    profile2.range,
                                    weights,
                                    success.th=0.10,success.p=0.75,
                                    cluster.quantiles=c(0.10,0.90)){
    if(missing(profile1.density)){
        profile1.fun.cdf <- function(x,i){ pnorm(x,profile1.mean[i],profile1.sd[i]) }
        profile1.fun.q   <- function(x,i){ qnorm(x,profile1.mean[i],profile1.sd[i]) }
    }else{
        profile1.fun.cdf <- function(x,i){ cdf.density(profile1.density[[i]],x) }
        profile1.fun.q   <- function(x,i){ quantiles.density(profile1.density[[i]],x) }
    }
    
    if(profile1.nbcells>=1){
        measures <- rep(1,length(weights))
        for(i in 1:length(weights)){
            cluster_low    <- profile1.fun.q(cluster.quantiles[1],i)[[1]]
            cluster_high   <- profile1.fun.q(cluster.quantiles[2],i)[[2]]
            measures[i]    <- min(abs(profile2.range[i,"inf"]-cluster_low),abs(profile2.range[i,"sup"]-cluster_high))
            if(profile2.range[i,"inf"] <= cluster_low && cluster_high <= profile2.range[i,"sup"])
                measures[i] <- 0    
        }
    }else{
        measures   <- rep(Inf,length(weights))
    }
    
    successes <- (measures <= success.th)
      
    measure   <- weighted.mean(measures,weights,na.rm = TRUE)
    pvalue    <- aggreg.binom(successes,weights,success.p)
    
    res <- list(measure  = measure,
        pvalue           = pvalue,
        marker.measures  = measures,
        marker.successes = successes)
    return(res)
}


# title Internal - Weighted binomial test for computations of profile similarity or inclusion p-values
#
# @description Performs a weighted binomial test to compute the significance of profile similarities or profile inclusions. 
#
# @details In the default statistical approach, each marker similarity or inclusion measure below a specific threshold value models a success in a Bernoulli experiment. The computed p-value represents then the significance of the proportion of similar or included cell markers when comparing two cytometry profiles. Each success or fail is repeated the number of times indicated by the 'weights' parameter. 
# 
# @param successes a boolean vector containing the marker successes (TRUE) or fail (FALSE)
# @param weights an integer vector containing the marker weights 
# @param success.p a numeric value specifying the expected success probability (set at 0.75 by default)
#
# @return a numeric providing the p-value of the weighted binomial test
aggreg.binom <- function(successes,weights,success.p=0.75){
    successes <- successes[!is.na(successes)]
    weights   <- weights[!is.na(successes)]
    pvalue    <- pbinom(sum(successes*weights), sum(weights), success.p, lower.tail = FALSE)
    return(pvalue)
}


# title Internal - Computation the Kolmogorov-Smirnov statistics of cell markers
#
# @description Computes the Kolmogorov-Smirnov statistics for each marker of two cytometry profiles.
#
# @details The Kolmogorov-Smirnov statistics corresponds to the maximal absolute difference between the cumulative distributions of two empirical distribution functions. A large statistics indicates that the two probability distributions are different, while a small statistics indicates that they are similar.
#
# @param profile1.fun.cdf the cumulative density function of the first cell cluster object  
# @param profile1.fun.q the quantile function of the first cell cluster object (used to obtain the marker expression ranges)
# @param profile2.fun.cdf the cumulative density function of the second cell cluster object 
# @param profile2.fun.q the quantile function of the second cell cluster object (used to obtain the marker expression ranges)
# @param markers.nb a numeric indicating the number of markers in both objects
#
# @return a number vector containing the Kolmogorov-Smirnov statistics for each marker
KS.stat <- function(profile1.fun.cdf,profile1.fun.q,profile2.fun.cdf,profile2.fun.q,markers.nb){
    measures <- numeric(markers.nb)
    bounds  <- c(0.0001,0.9999)
    for(i in 1:markers.nb){
        lower_bound <- min(profile1.fun.q(0,i),profile2.fun.q(0,i))
        if(is.infinite(lower_bound))
            lower_bound <- min(profile1.fun.q(bounds[1],i),profile2.fun.q(bounds[1],i))
        upper_bound <- max(profile1.fun.q(1,i),profile2.fun.q(1,i))
        if(is.infinite(upper_bound))
            upper_bound <- max(profile1.fun.q(bounds[2],i),profile2.fun.q(bounds[2],i))
        fun <- function(x){
            return(abs(profile1.fun.cdf(x,i)-profile2.fun.cdf(x,i)))
        }
        opt        <- optimize(f=fun,maximum=TRUE,interval=c(lower_bound,upper_bound))
        measures[i] <- opt$objective
    }
    return(measures)
}
