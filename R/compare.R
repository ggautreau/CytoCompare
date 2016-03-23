#' @title Compare two cytometry profiles.
#'
#' @description Cytometry profiles contained in CELL, CLUSTER, or GATE objects can be compared using the `compare()` function. Comparison results are stored in a RES object. Comparisons can be performed between profiles of same types or between profiles of different types:\cr
#' * in the default statistical approach, if the comparisons are performed on profiles of same type then profiles will be compared to identify similar profiles\cr
#' * in the default statistical approach, if the comparisons are performed on profiles of different types then profiles will be compared to identify included profiles
#'
#' For each comparison of two cytometry profiles, a similarity or an inclusion measure is provided as well as a p-value asserting the statistical significance of the similarity or inclusion. A similarity or inclusion measure is calculated for each marker of the profiles to compare. Different measures are used depending of the types of profiles to compare (please refer to the details section). An aggregated similarity or inclusion measure is computed based on the weighted sum of marker measures. Then, each marker having a measure below a specific threshold models a similarity or inclusion success, and a weighted binomial test provides the significance of the proportion of similar or included markers.
#' 
#' Comparisons can be performed based on the whole set of common markers between the two profiles, or based on a subset of markers specified by the user. Moreover, markers can be weighted in the comparison procedure, via a MWEIGHTS object.
#'
#' If only one object is provided to the `compare()` function then the comparisons will be performed between all profiles of this object. If two objects are provided to the `compare()` function then the comparisons will be performed between all possible pairs of profiles between these two objects.
#'
#' Importantly, users can define their own function to perform the statistical comparisons of the profiles, using the `method` parameter. 
#'
#' @details Different parameters can be defined, via the method.params named character list, to specify the behaviour of the such kind of comparisons:\cr
#' * the success.th parameter indicates the similarity or inclusion threshold\cr
#' * the success.p parameter indicates the expected proportion of marker successes\cr
#' * the nbcells.th parameter indicates the number of cells per cluster below which the marker expression density of a cell cluster profile will be approximated by a normal distribution\cr
#' * the cluster.quantiles parameter indicates the quantiles that will define the marker expression ranges for the cell cluster profile 
#' 
#' In the case of comparisons between cell profiles, the marker similarity measures are calculated based on the Euclidean distance. The parameter success.th is set to 0.20 by default and the parameter success.p is set to 0.75 default. 
#'
#' In the case of comparisons between cell cluster profiles, the marker similarity measures are calculated based on the Kolmogorov-Smirnov distance. The parameter success.th is set to 0.30 by default and the parameter success.p is set to 0.75 default. The nbcells.th parameter indicates the number of cells per cluster below which the density will be approximated by a normal distribution (set to 50 by default)
#'
#' In the case of comparisons between gate profiles, gates are modeled by uniform distributions, and the marker similarity measures are calculated based on the Kolmogorov-Smirnov distance. The parameter success.th is set to 0.30 by default and the parameter success.p is set to 0.75 default.
#'
#' In the case of comparisons between cell profiles and gate profiles, the marker inclusion measures are defined as the minimal distances between the cell expression and the gate ranges. The marker inclusion measures are equals to zero if the cell marker expressions are included in the gate marker ranges. The parameter success.th is set to 0.10 by default and the parameter success.p is set to 0.75 default.
#'
#' In the case of comparisons between cell cluster profiles and cluster profiles, the marker inclusion measures are defined as the minimal distances between the cell expressions and the cell cluster ranges. The marker inclusion measures are equals to zero if the cell marker expressions are included in the cell cluster marker ranges. The parameter success.th is set to 0.10 by default and the parameter success.p is set to 0.75 default. The cluster.quantiles parameter indicates the quantiles that will define the marker expression ranges for the cell cluster profile (set to 0.10 and 0.90 by default).
#'
#' In the case of comparisons between cell cluster profiles and gate profiles, the inclusion measures used are equals to the maximal absolute distances between the cluster and gate lower bounds and the cluster and gate upper bounds. The inclusion measures are equals to zero if the cell cluster marker expressions are included in the gate marker ranges. The parameter success.th is set to 0.10 by default and the parameter success.p is set to 0.75 default. The cluster.quantiles parameter indicates the quantiles that will define the marker expression ranges for the cell cluster profile (set to 0.10 and 0.90 by default).
#'
#' @param object1 a CELL, CLUSTER or GATE object
#' @param object2 a CELL, CLUSTER or GATE object
#' @param mweights a MWEIGHTS object specifying the markers to use in the comparison procedure with theirs associated weights
#' @param method a function or a character specifying the name of a function to use when performing the statistical comparisons between the cytometry profiles
#' @param method.params a named character list used to parametrize the comparison function (please see the details section)
#' @param ... other parameters
#'
#' @return a S4 object of class RES
#'
#' @name compare
#' @rdname compare-methods
#'
#' @export
setGeneric("compare", function(object1,object2,...){ standardGeneric("compare")})

#' @rdname compare-methods
#' @export
setMethod("compare",c("CELL","missing"),
    function(object1,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object1)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.marker.measures  <- NULL
        comb.marker.successes <- NULL
        
        message(paste0("Comparing CELL:",object1@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object1@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object1@profiles.nb*100)
                message(paste0(percent,"% CELL:",object1@name,":",object1@profiles[i]," vs. ","CELL:",object1@name,":",object1@profiles[j]))
                
                params <- list(profile1.type = "CELL",
                        profile2.type        = "CELL",
                        profile1.intensities = object1@intensities[i,object1@markers %in% markers],
                        profile2.intensities = object1@intensities[j,object1@markers %in% markers],
                        weights              = mweights@weights)
                        
                res.sub               <- do.call(method,c(params,method.params))
                comparison            <- data.frame(profile1=paste("CELL",object1@name,object1@profiles[i],sep=":"),profile2=paste("CELL",object1@name,object1@profiles[j],sep=":"),check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.marker.measures  <- rbind(comb.marker.measures,res.sub$marker.measures)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count  <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.marker.measures,comb.marker.successes)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("CLUSTER","missing"),
    function(object1,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object1)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.marker.measures  <- NULL
        comb.marker.successes <- NULL
        object.density_tmp    <- prod(dim(object1@densities) == c(object1@profiles.nb,length(object1@markers)))
        
        message(paste0("Comparing CLUSTER:",object1@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object1@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object1@profiles.nb*100)
                message(paste0(percent,"% CLUSTER:",object1@name,":",object1@profiles[i]," vs. ","CLUSTER:",object1@name,":",object1@profiles[j]))
                
                params <- list(profile1.type = "CLUSTER",
                        profile2.type       = "CLUSTER",
                        profile1.mean       = object1@means[i,object1@markers %in% markers],
                        profile1.sd         = object1@sd[i,object1@markers %in% markers],
                        profile1.nbcells    = object1@profiles.sizes[i],
                        profile2.mean       = object1@means[j,object1@markers %in% markers],
                        profile2.sd         = object1@sd[j,object1@markers %in% markers],
                        profile2.nbcells    = object1@profiles.sizes[i],
                        weights             = mweights@weights)
                        
                if(object.density_tmp)
                    params <- c(params,list(profile1.density=object1@densities[i,object1@markers %in% markers],profile2.density=object1@densities[j,object1@markers %in% markers]))
                
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("CLUSTER",object1@name,object1@profiles[i],sep=":"),profile2=paste("CLUSTER",object1@name,object1@profiles[j],sep=":"),check.names = FALSE)

                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.marker.measures  <- rbind(comb.marker.measures,res.sub$marker.measures)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.marker.measures,comb.marker.successes)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("GATE","missing"),
    function(object1,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object1)
        comparisons           <- NULL
        markers               <- mweights@markers
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.marker.measures  <- NULL
        comb.marker.successes <- NULL
        
        message(paste0("Comparing GATE:",object1@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object1@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object1@profiles.nb*100)
                message(paste0(percent,"% GATE:",object1@name,":",object1@profiles[i]," vs. ","GATE:",object1@name,":",object1@profiles[j]))
                
                params <- list(profile1.type = "GATE",
                        profile2.type        = "GATE",
                        profile1.range       = object1@ranges[i,object1@markers %in% markers,],
                        profile2.range       = object1@ranges[j,object1@markers %in% markers,],
                        weights              = mweights@weights)
                        
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("GATE",object1@name,object1@profiles[i],sep=":"),profile2=paste("GATE",object1@name,object1@profiles[j],sep=":"),check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.marker.measures  <- rbind(comb.marker.measures,res.sub$marker.measures)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.marker.measures,comb.marker.successes)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("CELL","CELL"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.marker.measures  <- NULL
        comb.marker.successes <- NULL
        
        message(paste0("Comparing CELL:",object1@name," vs. ","CELL:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% CELL:",object1@name,":",object1@profiles[i]," vs. ","CELL:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "CELL",
                        profile2.type        = "CELL",
                        profile1.intensities = object1@intensities[i,object1@markers %in% markers],
                        profile2.intensities = object2@intensities[j,object2@markers %in% markers],
                        weights              = mweights@weights)

                res.sub               <- do.call(method,c(params,method.params))
                comparison            <- data.frame(profile1=paste("CELL",object1@name,object1@profiles[i],sep=":"),profile2=paste("CELL",object2@name,object2@profiles[j],sep=":"),check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.marker.measures  <- rbind(comb.marker.measures,res.sub$marker.measures)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.marker.measures,comb.marker.successes)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("CLUSTER","CLUSTER"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.marker.measures  <- NULL
        comb.marker.successes <- NULL
        profile1.density_tmp  <- prod(dim(object1@densities) == c(object1@profiles.nb,length(object1@markers)))
        profile2.density_tmp  <- prod(dim(object2@densities) == c(object2@profiles.nb,length(object2@markers)))
        
        message(paste0("Comparing CLUSTER:",object1@name," vs. ","CLUSTER:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% CLUSTER:",object1@name,":",object1@profiles[i]," vs. ","CLUSTER:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "CLUSTER",
                        profile2.type        = "CLUSTER",
                        profile1.mean        = object1@means[i,object1@markers %in% markers],
                        profile1.sd          = object1@sd[i,object1@markers %in% markers],
                        profile1.nbcells     = object1@profiles.sizes[i],
                        profile2.mean        = object2@means[j,object2@markers %in% markers],
                        profile2.sd          = object2@sd[j,object2@markers %in% markers],
                        profile2.nbcells     = object2@profiles.sizes[j],
                        weights              = mweights@weights)
                
                if(profile1.density_tmp)
                    params <- c(params,list(profile1.density=object1@densities[i,object1@markers %in% markers]))
                if(profile2.density_tmp)
                    params <- c(params,list(profile2.density=object2@densities[j,object2@markers %in% markers]))
                    
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("CLUSTER",object1@name,object1@profiles[i],sep=":"),profile2=paste("CLUSTER",object2@name,object2@profiles[j],sep=":"),check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.marker.measures  <- rbind(comb.marker.measures,res.sub$marker.measures)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.marker.measures,comb.marker.successes)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("GATE","GATE"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.marker.measures  <- NULL
        comb.marker.successes <- NULL
        
        message(paste0("Comparing GATE:",object1@name," vs. ","GATE:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% GATE:",object1@name,":",object1@profiles[i]," vs. ","GATE:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "GATE",
                        profile2.type        = "GATE",
                        profile1.range       = object1@ranges[i,object1@markers %in% markers,],
                        profile2.range       = object2@ranges[j,object2@markers %in% markers,],
                        weights              = mweights@weights)
                
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("GATE",object1@name,object1@profiles[i],sep=":"),profile2=paste("GATE",object2@name,object2@profiles[j],sep=":"),check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.marker.measures  <- rbind(comb.marker.measures,res.sub$marker.measures)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.marker.measures,comb.marker.successes)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("CELL","CLUSTER"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.marker.measures  <- NULL
        comb.marker.successes <- NULL
        profile2.density_tmp  <- prod(dim(object2@densities) == c(object2@profiles.nb,length(object2@markers)))
        
        message(paste0("Comparing CELL:",object1@name," vs. ","CLUSTER:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% CELL:",object1@name,":",object1@profiles[i]," vs. ","CLUSTER:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "CELL",
                        profile2.type        = "CLUSTER",
                        profile1.intensities = object1@intensities[i,object1@markers %in% markers],
                        profile2.mean        = object2@means[j,object2@markers %in% markers],
                        profile2.sd          = object2@sd[j,object2@markers %in% markers],
                        profile2.nbcells     = object2@profiles.sizes[j],
                        weights              = mweights@weights)
                
                if(profile2.density_tmp)
                    params <- c(params,list(profile2.density=object2@densities[j,object2@markers %in% markers]))
                    
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("CELL",object1@name,object1@profiles[i],sep=":"),profile2=paste("CLUSTER",object2@name,object2@profiles[j],sep=":"),check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.marker.measures  <- rbind(comb.marker.measures,res.sub$marker.measures)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.marker.measures,comb.marker.successes)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("CLUSTER","CELL"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        comparisons           <- NULL
        markers               <- mweights@markers
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.marker.measures  <- NULL
        comb.marker.successes <- NULL
        profile1.density_tmp  <- prod(dim(object1@densities) == c(object1@profiles.nb,length(object1@markers)))
       
       message(paste0("Comparing CLUSTER:",object1@name," vs. ","CELL:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% CLUSTER:",object1@name,":",object1@profiles[i]," vs. ","CELL:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "CLUSTER",
                        profile2.type        = "CELL",
                        profile1.mean        = object1@means[i,object1@markers %in% markers],
                        profile1.sd          = object1@sd[i,object1@markers %in% markers],
                        profile1.nbcells     = object1@profiles.sizes[i],
                        profile2.intensities = object2@intensities[j,object2@markers %in% markers],
                        weights              = mweights@weights)
                
                if(profile1.density_tmp)
                    params <- c(params,list(profile1.density=object1@densities[i,object1@markers %in% markers]))
                    
                res.sub              <- do.call(method,c(params, method.params))
                comparison           <- data.frame(profile1=paste("CLUSTER",object1@name,object1@profiles[i],sep=":"),profile2=paste("CELL",object2@name,object2@profiles[j],sep=":"),check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.marker.measures  <- rbind(comb.marker.measures,res.sub$marker.measures)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.marker.measures,comb.marker.successes)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("CELL","GATE"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.marker.measures  <- NULL
        comb.marker.successes <- NULL
        
        message(paste0("Comparing CELL:",object1@name," vs. ","GATE:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% CELL:",object1@name,":",object1@profiles[i]," vs. ","GATE:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "CELL",
                        profile2.type        = "GATE",
                        profile1.intensities = object1@intensities[i,object1@markers %in% markers],
                        profile2.range       = object2@ranges[j,object2@markers %in% markers,],
                        weights              = mweights@weights)
                
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("CELL",object1@name,object1@profiles[i],sep=":"),profile2=paste("GATE",object2@name,object2@profiles[j],sep=":"),check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.marker.measures  <- rbind(comb.marker.measures,res.sub$marker.measures)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.marker.measures,comb.marker.successes)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("GATE","CELL"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.marker.measures  <- NULL
        comb.marker.successes <- NULL
        
        message(paste0("Comparing GATE:",object1@name," vs. ","CELL:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% GATE:",object1@name,":",object1@profiles[i]," vs. ","CELL:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "GATE",
                        profile2.type        = "CELL",
                        profile1.range       = object1@ranges[i,object1@markers %in% markers,],
                        profile2.intensities = object2@intensities[j,object2@markers %in% markers],
                        weights              = mweights@weights)
                        
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("GATE",object1@name,object1@profiles[i],sep=":"),profile2=paste("CELL",object2@name,object2@profiles[j],sep=":"),check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.marker.measures  <- rbind(comb.marker.measures,res.sub$marker.measures)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.marker.measures,comb.marker.successes)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("CLUSTER","GATE"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.marker.measures  <- NULL
        comb.marker.successes <- NULL
        profile1.density_tmp  <- prod(dim(object1@densities) == c(object1@profiles.nb,length(object1@markers)))
        
        message(paste0("Comparing CLUSTER:",object1@name," vs. ","GATE:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% CLUSTER:",object1@name,":",object1@profiles[i]," vs. ","GATE:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "CLUSTER",
                        profile2.type        = "GATE",
                        profile1.mean        = object1@means[i,object1@markers %in% markers],
                        profile1.sd          = object1@sd[i,object1@markers %in% markers],
                        profile1.nbcells     = object1@profiles.sizes[i],
                        profile2.range       = object2@ranges[j,object2@markers %in% markers,],
                        weights              = mweights@weights)
                        
                if(profile1.density_tmp)
                    params <- c(params,list(profile1.density=object1@densities[i,object1@markers %in% markers]))
                    
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("CLUSTER",object1@name,object1@profiles[i],sep=":"),profile2=paste("GATE",object2@name,object2@profiles[j],sep=":"),check.names = FALSE)

                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.marker.measures  <- rbind(comb.marker.measures,res.sub$marker.measures)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.marker.measures,comb.marker.successes)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("GATE","CLUSTER"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.marker.measures  <- NULL
        comb.marker.successes <- NULL
        profile2.density_tmp  <- prod(dim(object2@densities) == c(object2@profiles.nb,length(object2@markers)))
        
        message(paste0("Comparing GATE:",object1@name," vs. ","CLUSTER:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% GATE:",object1@name,":",object1@profiles[i]," vs. ","CLUSTER:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "GATE",
                        profile2.type        = "CLUSTER",
                        profile1.range       = object1@ranges[i,object1@markers %in% markers,],
                        profile2.mean        = object2@means[j,object2@markers %in% markers],
                        profile2.sd          = object2@sd[j,object2@markers %in% markers],
                        profile2.nbcells     = object2@profiles.sizes[j],
                        weights              = mweights@weights)
                        
                if(profile2.density_tmp)
                    params <- c(params,list(profile2.density=object2@densities[j,object2@markers %in% markers]))
                    
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("GATE",object1@name,object1@profiles[i],sep=":"),profile2=paste("CLUSTER",object2@name,object2@profiles[j],sep=":"),check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.marker.measures  <- rbind(comb.marker.measures,res.sub$marker.measures)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.marker.measures,comb.marker.successes)
        return(res)
    }
)


# title Internal - Creation of a RES object based the comparison results 
#
# @description This function is used internally to create a RES object based on calculated marker measures, marker p-values, aggregated measure and aggregated p-value. 
#
# @param comparisons a data.frame with two-columns containing the names of the profiles to compare
# @param markers a character vector containing the marker names
# @param comb.measure a numeric value containing the comparison measure
# @param comb.pvalue a numeric value containing the comparison p-value
# @param comb.marker.measures a numeric vector containing the marker comparison measures
# @param comb.marker.successes a numeric vector containing the marker comparison success
#
# @return a RES object containing
RES.format <- function(comparisons,markers,comb.measure,comb.pvalue,comb.marker.measures,comb.marker.successes){
    
    comparisons[,1]    <- as.character(comparisons[,1])
    comparisons[,2]    <- as.character(comparisons[,2])
    comparisons$measure <- comb.measure
    comparisons$pvalue <- comb.pvalue
    
    rownames(comparisons)           <- paste(comparisons[,1],comparisons[,2],sep="/vs/")
    comb.marker.measures            <- data.frame(comb.marker.measures)
    comb.marker.successes           <- data.frame(comb.marker.successes)
    colnames(comb.marker.measures)  <- markers
    rownames(comb.marker.measures)  <- rownames(comparisons) 
    colnames(comb.marker.successes) <- markers
    rownames(comb.marker.successes) <- rownames(comparisons) 
    
    res <- RES(comparisons = comparisons,
        comparisons.nb     = nrow(comparisons),
        markers            = markers,
        marker.measures    = comb.marker.measures,
        marker.successes   = comb.marker.successes
        )
    return(res)
}
