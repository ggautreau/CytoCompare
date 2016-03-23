#' @title CELL class definition
#'
#' @description CELL is a S4 object containing one or several cell profiles.
#'
#' @details This object mainly stores for each cell profile, the intensities of each marker.
#'
#' @slot name a character indicating the internal name of the CELL object
#' @slot profiles a character vector containing the names of the cell profiles 
#' @slot profiles.nb an integer value indicating the number of cell profiles
#' @slot markers a character vector containing the marker names
#' @slot markers.nb an integer value indicating the number of markers
#' @slot intensities a numeric matrix containing the intensities of each marker for each cell profile
#' @slot overview.function a character specifying the name of a function to call when plotting the CELL object overview (please refer to the documentation of the 'plot()' function)
#' @slot layout a numeric matrix that can be used to store the positions of cells in a 2-dimensional space (e.g. tSNE1 and tSNE2 dimensions provided by viSNE)
#'
#' @import methods
#'
#' @name CELL-class
#' @rdname CELL-class
#' @exportClass CELL
CELL <- setClass("CELL",
    slots=c(name          = "character",
        profiles          = "character",
        profiles.nb       = "integer",
        markers           = "character",
        markers.nb        = "integer",
        intensities       = "matrix",
        overview.function = "character",
        layout            = "matrix"),
    validity=function(object){
        if(length(object@profiles)!=length(unique(object@profiles)))
            stop("Error in profiles slot: profile names are not unique")
        if(length(object@markers)!=length(unique(object@markers)))
            stop("Error in markers slot: marker names are not unique")
        if(object@profiles.nb!=length(object@profiles))
            stop("Error in profiles.nb slot: profiles.nb do not correspond to the number of profile names")
        if(object@markers.nb!=length(object@markers))
            stop("Error in markers.nb slot: markers.nb do not correspond to the number of marker names")
        return(TRUE)
    }
)
setMethod("initialize",c("CELL"),
    function(.Object,
            name              = "",
            profiles          = "",
            profiles.nb       = 0,
            markers           = "",
            markers.nb        = 0,
            intensities       = as.matrix(0),
            overview.function = "",
            layout            = as.matrix(0)){
        if(profiles.nb==0)
            stop("Error can not create a CELL object with no profile")
        .Object@name              = name     
        .Object@profiles          = profiles
        .Object@profiles.nb       = profiles.nb
        .Object@markers           = markers
        .Object@markers.nb        = markers.nb
        .Object@intensities       = intensities
        .Object@overview.function = overview.function
        .Object@layout            = layout
        return(.Object)
    }
)


#' @title CLUSTER class definition
#'
#' @description CLUSTER is a S4 object containing one or several cell cluster profiles.
#'
#' @details This object mainly stores for each cell cluster profile, the means, the standard deviations and the densities of each marker.
#'
#' @slot name a character indicating the internal name of the CLUSTER object
#' @slot profiles a character vector containing the names of the cell cluster profiles
#' @slot profiles.nb an integer value indicating the number of cell cluster profiles
#' @slot profiles.sizes an integer vector indicating the number of cells associated to each cluster profile
#' @slot markers a character vector containing the marker names
#' @slot markers.nb an integer value indicating the number of markers
#' @slot markers.clustering a logical vector specifying the makers used as clustering markers
#' @slot means a numeric matrix containing the means of each maker for each cluster profile
#' @slot sd a numeric matrix containing the standard deviations of each maker for each cluster profile
#' @slot densities a matrix of DENSITY objects containing the densities of each marker for each cluster profile
#' @slot overview.function a character specifying the name of a function to call when plotting the CLUSTER object overview (please refer to the documentation of the 'plot()' function)
#' @slot graph an object that can be used to store a visual representation of the cell clusters (e.g. a SPADE tree)
#' @slot graph.layout a numeric matrix that can be used to store the positions of cell clusters in a 2-dimensional space (e.g. a SPADE tree layout)
#' 
#' @import methods
#'
#' @name CLUSTER-class
#' @rdname CLUSTER-class
#' @exportClass CLUSTER
#' @export CLUSTER
CLUSTER <- setClass("CLUSTER",
    slots=c(name           = "character",
        profiles           = "character",
        profiles.nb        = "integer",
        profiles.sizes     = "integer",
        markers            = "character",
        markers.nb         = "integer",
        markers.clustering = "logical",
        means              = "matrix",
        sd                 = "matrix",
        overview.function  = "character",
        graph              = "ANY",
        graph.layout       = "matrix",
        densities          = "matrix"),
    validity=function(object){
        if(length(object@profiles)!=length(unique(object@profiles)))
            stop("Error in profiles slot: profile names are not unique")
        if(length(object@markers)!=length(unique(object@markers)))
            stop("Error in markers slot: marker names are not unique")
        if(object@profiles.nb!=length(object@profiles))
            stop("Error in profiles.nb slot: profiles.nb do not correspond to the number of profile names")
        if(object@markers.nb!=length(object@markers))
            stop("Error in markers.nb slot: markers.nb do not correspond to the number of marker names")
        if(object@profiles.nb==0)
            stop("Error can not create a CLUSTER object with no profile")
        return(TRUE)
    }
)
setMethod("initialize",c("CLUSTER"),
    function(.Object,
            name               = "",
            profiles           = "",
            profiles.nb        = 0,
            profiles.sizes     = 0,
            markers            = "",
            markers.nb         = 0,
            markers.clustering = FALSE,
            means              = as.matrix(0),
            sd                 = as.matrix(0),
            overview.function  = "",
            graph              = NULL,
            graph.layout       = as.matrix(0),
            densities          = as.matrix(0)){
        if(profiles.nb==0)
            stop("Error can not create a CLUSTER object with no profile")
        .Object@name               <- name     
        .Object@profiles           <- profiles
        .Object@profiles.nb        <- profiles.nb
        .Object@profiles.sizes     <- profiles.sizes
        .Object@markers            <- markers
        .Object@markers.nb         <- markers.nb
        .Object@markers.clustering <- markers.clustering
        .Object@means              <- means
        .Object@sd                 <- sd
        .Object@overview.function  <- overview.function
        .Object@graph              <- graph
        .Object@graph.layout       <- graph.layout
        .Object@densities          <- densities
        return(.Object)
    }
)

#' @title GATE class definition
#'
#' @description GATE is a S4 object containing one or several gate profiles.
#'
#' @details This object mainly stores for each gate profile, the intensity ranges of each marker. 
#'
#' @slot name a character indicating the internal name of the GATE object
#' @slot profiles a character vector containing the names of the gate profiles
#' @slot profiles.nb an integer value indicating the number of cell gate profiles
#' @slot markers a character vector containing the marker names
#' @slot markers.nb an integer value indicating the number of markers
#' @slot ranges a 3-dimensional numeric array containing the intensity ranges of each marker for each gate profile
#'
#' @import methods
#'
#' @name GATE-class
#' @rdname GATE-class
#' @exportClass GATE
#' @export GATE
GATE <- setClass("GATE",
    slots=c(name    = "character",
        profiles    = "character",
        profiles.nb = "integer",
        markers     = "character",
        markers.nb  = "integer",
        ranges      = "array"),
    validity=function(object){
        if(length(object@profiles)!=length(unique(object@profiles)))
            stop("Error in profiles slot: profile names are not unique")
        if(length(object@markers)!=length(unique(object@markers)))
            stop("Error in markers slot: marker names are not unique")
        if(object@profiles.nb!=length(object@profiles))
            stop("Error in profiles.nb slot: profiles.nb do not correspond to the number of profile names")
        if(object@markers.nb!=length(object@markers))
            stop("Error in markers.nb slot: markers.nb do not correspond to the number of marker names")
        if(object@profiles.nb==0)
            stop("Error can not create a GATE object with no profile")
        return(TRUE)
    }
)
setMethod("initialize",c("GATE"),
    function(.Object,
            name        = "",
            profiles    = "",
            profiles.nb = 0,
            markers     = "",
            markers.nb  = 0,
            ranges      = as.matrix(0)){
        if(profiles.nb==0)
            stop("Error can not create a GATE object with no profile")
        .Object@name        <- name     
        .Object@profiles    <- profiles
        .Object@profiles.nb <- profiles.nb
        .Object@markers     <- markers
        .Object@markers.nb  <- markers.nb
        .Object@ranges      <- ranges
        return(.Object)
    }
)


#' @title DENSITY class definition
#'
#' @description DENSITY is a S4 object used to stores a marker expression density.
#'
#' @details This object mainly stores for each marker: the bin characteristics, the negative and positive marker densities values, and the number of cells used in the density estimation. Densities are stored using two numeric vectors: values.neg for the negative densities and values.pos for the positive densities. This strategy allows to compute and store densities without defining an absolute minimal value.
#'
#' @slot name a character indicating the internal name of the CELL object
#' @slot bin.interval a numeric vector of two values specifying the density boundaries
#' @slot bin.nb a numeric vector of two values specifying the numbers of negative and positive bins
#' @slot values.pos a numeric vector containing the positive density bins
#' @slot values.neg a numeric vector containing the negative density bins
#' @slot point.nb a numeric value indicating the number of point used to compute the expression density
#' @slot bin.width a numeric value indicating the width of the bins used in the density estimation
#'
#' @import methods
#'
#' @name DENSITY-class
#' @rdname DENSITY-class
#' @exportClass DENSITY
#' @export DENSITY
DENSITY <- setClass("DENSITY",
    slots=c(name         = "character",
            bin.interval = "numeric",
            bin.nb       = "numeric",
            values.pos   = "numeric",
            values.neg   = "numeric",
            point.nb     = "numeric",
            bin.width    = "numeric")
)
setMethod("initialize",c("DENSITY"),
    function(.Object,bin.width=0.05,bin.interval,point.nb,values.neg,values.pos,values=NULL,name=""){
        if(is.null(values)){
            .Object@bin.interval <- bin.interval
            .Object@bin.nb       <- c(length(values.neg),length(values.pos))
            .Object@point.nb     <- point.nb
            .Object@values.pos   <- values.pos
            .Object@values.neg   <- values.neg
        }else{
            density              <- compute.density(values,bin.width)
            .Object@bin.interval <- density$bin.interval
            .Object@bin.nb       <- density$bin.nb
            .Object@values.pos   <- density$values.pos
            .Object@values.neg   <- density$values.neg 
        }
        .Object@bin.width <- bin.width
        .Object@name      <- name
        return(.Object)
    }
)


#' @title MWEIGHTS class definition
#'
#' @description MWEIGHTS is a S4 object containing the marker weights to use in the comparison computations. 
#'
#' @details This object mainly stores for each marker: the markers names and marker weights. 
#'
#' @slot markers a character vector containing the marker names
#' @slot weights a numeric vector containing the marker weights
#'
#' @import methods
#'
#' @name MWEIGHTS-class
#' @rdname MWEIGHTS-class
#' @exportClass MWEIGHTS
#' @export MWEIGHTS
MWEIGHTS <- setClass("MWEIGHTS",
    slots = c(markers = "character",
              weights = "numeric")
)


#' @title RES class definition
#'
#' @description RES is a S4 object containing one or several comparison results.
#'
#' @details This object mainly stores for each comparison result: the similarity or inclusion measure, the associated similarity or inclusion p-value, the markers similarity or inclusion measures and the marker similarity or inclusion successes (measures below a specific threshold).
#'
#' @slot comparisons a data.frame containing for each comparison: the profile names, the similarity or inclusion measure, and the associated p-value
#' @slot comparisons.nb is an integer indicating the number of comparisons  
#' @slot markers a character vector containing the marker names used in the comparisons
#' @slot marker.measures a data.frame containing the marker similarity or inclusion measures for each comparison
#' @slot marker.successes a data.frame containing the marker successes for each comparison
#' 
#' @import methods
#'
#' @name RES-class
#' @rdname RES-class
#' @exportClass RES
#' @export RES
RES <- setClass("RES",
    slots=c(comparisons  = "data.frame",
        comparisons.nb   = "numeric",
        markers          = "character",
        marker.measures  = "data.frame",
        marker.successes = "data.frame"),
    validity=function(object){ 
        if(length(object@comparisons.nb)!=nrow(object@comparisons))
            stop("Error in comparisons.nb slot: comparisons.nb do not correspond to the number of comparisons")
         if(nrow(object@marker.measures)!=nrow(object@comparisons))
            stop("Error in marker.measures slot: marker.measures do not correspond to the number of comparisons")
         if(nrow(object@marker.successes)!=nrow(object@comparisons))
            stop("Error in marker.successes slot: marker.successes do not correspond to the number of comparisons")
        return(TRUE)
    }
)
setMethod("initialize",c("RES"),
          function(.Object,
                   comparisons      = data.frame(),
                   comparisons.nb   = 0,
                   markers          = "",
                   marker.measures  = data.frame(),
                   marker.successes = data.frame()){
            .Object@comparisons      <- comparisons     
            .Object@comparisons.nb   <- comparisons.nb
            .Object@markers          <- markers
            .Object@marker.measures  <- marker.measures
            .Object@marker.successes <- marker.successes
            return(.Object)
          }
)


# title Internal - Common markers between two cytometry objects
#
# @description Identifies the common markers between two cytometry objects (CELL, CLUSTER or GATE objects) and store the results in a MWEIGHTS object.
#
# @details The weight of each marker is set to '1' but can be changed afterwards using the 'set()' function.
#
# @param object1 a CELL, CLUSTER or GATE object
# @param object2 another CELL, CLUSTER or GATE object
#
# @return a S4 object of class MWEIGHTS.
inner.intersect <- function(object1,object2){
    cmarkers <- intersect(object1@markers,object2@markers)
    weights  <- rep(1,length(cmarkers))
    mweights <- MWEIGHTS(markers=cmarkers,weights=weights)
    return(mweights)
}


#' @title Identification of common markers between two cytometry objects
#'
#' @description Identifies the common markers between two cytometry objects (CELL, CLUSTER or GATE objects) and store the results in a MWEIGHTS object.
#'
#' @details The weight of each marker is set to 1 but can be changed afterwards using the function set.
#'
#' @param x a CELL, CLUSTER or GATE object
#' @param y a CELL, CLUSTER or GATE object
#'
#' @return a S4 object of class MWEIGHTS
#' @name intersect
#' @rdname intersect-methods
NULL

#' @rdname intersect-methods
#' @export
setMethod("intersect",c("CELL","CELL"),
    function(x,y){
        inner.intersect(x,y)
    }
)

#' @rdname intersect-methods
#' @export
setMethod("intersect",c("CLUSTER","CLUSTER"),
    function(x,y){
        inner.intersect(x,y)
    }
)

#' @rdname intersect-methods
#' @export
setMethod("intersect",c("GATE","GATE"),
    function(x,y){
        inner.intersect(x,y)
    }
)          

#' @rdname intersect-methods
#' @export
setMethod("intersect",c("CELL","CLUSTER"),
    function(x,y){
        inner.intersect(x,y)
    }
)

#' @rdname intersect-methods
#' @export
setMethod("intersect",c("CELL","GATE"),
    function(x,y){
        inner.intersect(x,y)
    }
)

#' @rdname intersect-methods
#' @export
setMethod("intersect",c("CLUSTER","GATE"),
    function(x,y){
        inner.intersect(x,y)
    }
)

#' @rdname intersect-methods
#' @export
setMethod("intersect",c("CLUSTER","CELL"),
    function(x,y){
        inner.intersect(x,y)
    }
)

#' @rdname intersect-methods
#' @export
setMethod("intersect",c("GATE","CELL"),
    function(x,y){
        inner.intersect(x,y)
    }
)

#' @rdname intersect-methods
#' @export
setMethod("intersect",c("GATE","CLUSTER"),
    function(x,y){
        inner.intersect(x,y)
    }
)


#' @title Extraction of subsets of data from CytoCompare objects
#'
#' @description Extracts subsets of CELL, CLUSTER, GATE, MWEIGHTS or RES object.
#'
#' @details For cytometry objects (CELL, CLUSTER, or GATE objects), the parameter i represents a vector of profiles to extract and the parameter j represents a vector of markers to extract.
#'
#' For MWEIGHTS objects, the parameter i represents a vector of markers to extract.
#'
#' For RES objects, the parameter i represents a vector of comparisons to extract.
#'
#' @param x a CELL, CLUSTER, GATE, MWEIGHTS or RES object
#' @param i a numeric, logical or character vector
#' @param j a numeric, logical or character vector
#'
#' @return a S4 object of class CELL, CLUSTER, GATE, MWEIGHTS or RES
#'
#' @name extract
#' @rdname extract-methods
NULL

#' @rdname extract-methods
#' @export
setMethod("[",c("CELL","ANY","ANY"),
    function(x,i,j){

        if(!missing(i) && length(i)>x@profiles.nb)
            stop("Too many cluster profiles to extract")
        if(!missing(j) && length(j)>length(x@markers))
            stop("Too many markers to extract")
        if(!missing(i) && !(is.logical(i) || is.numeric(i) || is.character(i)))
            stop(paste0("Wrong types in i: ",typeof(i)))
        if(!missing(j) && !(is.logical(j) || is.numeric(j) || is.character(j)))
            stop(paste0("Wrong types in j: ",typeof(j)))
                
        if(missing(i)){
            k <- 1:x@profiles.nb
        }else{
            if(is.character(i)){
                k <- which(x@profiles %in% i)
            }else{
                k <- i
            }
        }
        
        if(missing(j)){
            l <- 1:length(x@markers)
        }else{
            if(is.character(j)){
                l <- unlist(sapply(paste0("^",j,"$"),grep,x=x@markers))
            }else{
                l <- j
            }    
        }

        if(missing(i)){
            layout            <- x@layout
            overview.function <- x@overview.function 
        }else{
            layout            <- as.matrix(0)
            overview.function <- ""
        }
        
        cell <- CELL(name      = x@name,
            profiles           = x@profiles[k],
            profiles.nb        = length(x@profiles[k]),
            markers            = x@markers[l],
            markers.nb         = length(x@markers[l]),
            intensities        = x@intensities[k,l,drop=FALSE],
            layout             = layout,
            overview.function  = overview.function )
            
        return(cell)
    }
)

#' @rdname extract-methods
#' @export
setMethod("[",c("CLUSTER","ANY","ANY"),
    function(x,i,j){
    
        if(!missing(i) && length(i)>x@profiles.nb)
            stop("Too many cluster profiles to extract")
        if(!missing(j) && length(j)>length(x@markers))
            stop("Too many markers to extract")
        if(!missing(i) && !(is.logical(i) || is.numeric(i) || is.character(i)))
            stop(paste0("Wrong types in i: ",typeof(i)))
        if(!missing(j) && !(is.logical(j) || is.numeric(j) || is.character(j)))
            stop(paste0("Wrong types in j: ",typeof(j)))
            
        if(missing(i)){
            k <- 1:x@profiles.nb
        }else{
            if(is.character(i)){
                k <- which(x@profiles %in% i)
            }else{
                k <- i
            }
        }
        
        if(missing(j)){
            l <- 1:length(x@markers)
        }else{
            if(is.character(j)){
                l <- unlist(sapply(paste0("^",j,"$"),grep,x=x@markers))
            }else{
                l <- j
            }
        }
        
        if(missing(i)){
            overview.function <- x@overview.function    
            graph             <- x@graph
            graph.layout      <- x@graph.layout
        }else{
            overview.function <- ""
            graph             <- igraph::graph.empty(0,directed=FALSE)
            graph.layout      <- igraph::layout.auto(graph)
        }
            
        cluster <- CLUSTER(name = x@name,
            profiles            = x@profiles[k],
            profiles.nb         = length(x@profiles[k]),
            profiles.sizes      = as.integer(x@profiles.sizes[k]),
            markers             = x@markers[l],
            markers.nb          = length(x@markers[l]),
            markers.clustering  = x@markers.clustering[l],
            means               = x@means[k,l,drop=FALSE],
            sd                  = x@sd[k,l,drop=FALSE],
            densities           = x@densities[k,l,drop=FALSE],
            graph               = graph,
            graph.layout        = graph.layout,
            overview.function   = overview.function    )
        return(cluster)
    }
)

#' @rdname extract-methods
#' @export
setMethod("[",c("GATE","ANY","ANY"),
    function(x,i,j){
        
        if(!missing(i) && length(i)>x@profiles.nb)
            stop("Too many cluster profiles to extract")
        if(!missing(j) && length(j)>length(x@markers))
            stop("Too many markers to extract")
        if(!missing(i) && !(is.logical(i) || is.numeric(i) || is.character(i)))
            stop(paste0("Wrong types in i: ",typeof(i)))
        if(!missing(j) && !(is.logical(j) || is.numeric(j) || is.character(j)))
            stop(paste0("Wrong types in j: ",typeof(j)))
                
        if(missing(i)){
            k  <- 1:x@profiles.nb
        }else{
            if(is.character(i)){
                k <- which(x@profiles %in% i)
            }else{
                k <- i
            }
        }
        
        if(missing(j)){
            l <- 1:length(x@markers)
        }else{
            if(is.character(j)){
                l <- unlist(sapply(paste0("^",j,"$"),grep,x=x@markers))
            }else{
                l <- j
            }
        }
        
        gate <- GATE(name = x@name,
            profiles      = x@profiles[k],
            profiles.nb   = length(x@profiles[k]),
            markers       = x@markers[l],
            markers.nb    = length(x@markers[l]),
            ranges        = x@ranges[k,l,,drop=FALSE])
        return(gate)
    }
)

#' @rdname extract-methods
#' @export
setMethod("[",c("MWEIGHTS","ANY"),
    function(x,i){
    
        if(!(is.character(i) || is.numeric(i) || is.logical(i)))
            stop("Wrong types in parameter i")
    
        if(is.character(i)){
            k <- unlist(sapply(paste0("^",i,"$"),grep,x=x@markers))
        }else{
            k <- i
        }
        
        return(MWEIGHTS(markers=x@markers[k],weights=x@weights[k]))    
    }
)

#' @rdname extract-methods
#' @export
setMethod("[",c("RES","ANY"),
    function(x,i){
        
        if(!(is.logical(i) || is.numeric(i) || is.character(i)))
            stop(paste0("Wrong types in i: ",typeof(i)))
        if(is.logical(i) && length(i)!=nrow(x@comparisons))
            stop("logical vector i must have the same length as the number of comparisons in the RES object")
            
        if(is.character(i)){
            idx <- rownames(x@comparisons) %in% i
        }else{
            idx <- i
        }
        
        comparisons <- x@comparisons[idx,]
        
        res <- RES(comparisons = comparisons,
            comparisons.nb     = nrow(comparisons),
            markers            = x@markers,
            marker.measures    = x@marker.measures[idx,],
            marker.successes   = x@marker.successes[idx,])
        return(res)
    }
)


#' @title Change marker weights in a MWEIGHTS object
#'
#' @description Sets the weights of the markers in a MWEIGHTS object.
#'
#' @details Marker weights can be set based on their indexes or names in the MWEIGHTS object.
#'
#' @param x a MWEIGHTS object
#' @param i a numeric, logical or character vector
#' @param value a numeric vector containing the new marker weight values
#' @return a S4 object of class MWEIGHTS
#'
#' @name set
#' @rdname set-methods
NULL

#' @rdname set-methods
#' @export
setMethod("[<-",c("MWEIGHTS","numeric","missing","numeric"),
    function(x,i,value){
    
        if(min(i)<1 || max(i)>length(x@weights))
            stop("Indices out of range")
        if(length(i)>length(x@weights))
            stop("Too many markers to modify")
            
        x@weights[i] <- value
        return(x)
    }
)

#' @rdname set-methods
#' @export
setMethod("[<-",c("MWEIGHTS","character","missing","numeric"),
    function(x,i,value){
    
        if(sum(!i %in% x@markers)>0)
            stop("Unknown marker in i")
        if(length(i)>length(x@weights))
            stop("Too many markers to modify")
            
        idx            <- match(i,x@markers)
        x@weights[idx] <- value
        return(x)
    }
)

#' @rdname set-methods
#' @export
setMethod("[<-",c("MWEIGHTS","logical","missing","numeric"),
    function(x,i,value){
    
        if(length(i)!=length(x@weights))
            stop("logical vector i must have the same length as the number of weights in the MWEIGHTS object")
        
        x@weights[i] <- value
        return(x)
    }
)


#' Combination of CytoCompare objects
#'
#' @description Combines two or several CELL, CLUSTER, GATE or RES objects.
#'
#' @details All the different objects to combine must be of the same type.
#'
#' This function is especially useful when combining comparison results from different RES objects into a single RES object. RES objects can be combined to an empty RES object (i.e. RES()).
#' 
#' This function is also especially useful when combining cell profiles obtained from different FCS files into one single CELL object. 
#'
#' @param x a first CELL, CLUSTER, GATE or RES object
#' @param ... further objects of the same class as x to be combined
#' @param recursive a logical value indicating if the function recursively descends through lists combining all their elements into a vector. Not implemented and should be set to FALSE
#'
#' @return a S4 object of class CELL, CLUSTER, GATE or RES
#'
#' @name c
#' @rdname c-methods
NULL

#' @rdname c-methods
#' @export
setMethod("c",c("CELL"),
    function(x,...){
    
        other.CELL  <- list(x,...)
        name        <- c()
        markers     <- x@markers
        profiles.nb <- 0
        
        i <- 1
        for(cell in other.CELL){
            if(class(cell) != "CELL")
                stop(paste0("Cannot combine objects of different classes (element at position ",i," is of type ",class(cell)))
            name        <- c(name,cell@name)
            markers     <- union(markers,cell@markers)
            profiles.nb <- profiles.nb+cell@profiles.nb
            i <- i+1
        }
        name        <- paste0(name,collapse=";")
        profiles    <- as.character(1:profiles.nb)
        
        intensities <- matrix(NA,ncol=length(markers),nrow=profiles.nb,dimnames=list(profiles,markers))
        
        i <- 1
        for(cell in other.CELL){
            nb                                   <- cell@profiles.nb
            intensities[i:(i+nb-1),cell@markers] <- cell@intensities
            i <- i+nb
        }
        dimnames(intensities) <- NULL
        
        cell <- CELL(name = name,
            profiles      = profiles,
            profiles.nb   = length(profiles),
            markers       = markers,
            markers.nb    = length(markers),
            intensities   = intensities)
        return(cell)
    }
)

#' @rdname c-methods
#' @export
setMethod("c",c("CLUSTER"),       
    function(x,...){
    
        other.CLUSTER  <- list(x,...)
        name           <- c()
        markers        <- c()
        profiles       <- c()
        profiles.sizes <- c()
        profiles.nb    <- 0
        
        i <- 1 
        for(cluster in other.CLUSTER){
            if(class(cluster) != "CLUSTER")
                stop(paste0("Cannot combine objects of different classes (element at position ",i," is of type ",class(cluster)))
            name           <- c(name,cluster@name)
            markers        <- union(markers,cluster@markers)
            profiles.sizes <- c(profiles.sizes,cluster@profiles.sizes)
            profiles.nb    <- profiles.nb+cluster@profiles.nb
            profiles       <- c(profiles,cluster@profiles)
            i <- i+1
        }
        
        if(length(profiles)!=length(unique(profiles))){
            stop("Error: profile names are not unique")
        }
        
        name <- paste0(name,collapse=";")
        
        means          <- matrix(0,ncol=length(markers),nrow=profiles.nb,dimnames=list(profiles,markers))
        sd             <- matrix(0,ncol=length(markers),nrow=profiles.nb,dimnames=list(profiles,markers))
        densities      <- matrix(list(),ncol=length(markers),nrow=profiles.nb,dimnames=list(profiles,markers))
        
        i <- 1
        for(cluster in other.CLUSTER){
            nb                                    <- cluster@profiles.nb
            means[i:(i+nb-1),cluster@markers]     <- cluster@means
            sd[i:(i+nb-1),cluster@markers]        <- cluster@sd
            densities[i:(i+nb-1),cluster@markers] <- cluster@densities
            i <- i+nb
        }
        dimnames(means)     <- NULL
        dimnames(sd)        <- NULL
        dimnames(densities) <- NULL
        
        markers.clustering  <- rep(FALSE,length(markers))
        graph               <- igraph::graph.empty(0,directed=FALSE)
        graph.layout        <- igraph::layout.auto(graph)
        
        cluster <- CLUSTER(name = name,
            profiles            = profiles,
            profiles.nb         = length(profiles),
            profiles.sizes      = profiles.sizes,
            markers             = markers,
            markers.nb          = length(markers),
            markers.clustering  = markers.clustering,
            means               = means,
            sd                  = sd,
            densities           = densities,
            graph               = graph,
            graph.layout        = graph.layout)
        return(cluster)
    }
)

#' @rdname c-methods
#' @export
setMethod("c",c("GATE"),       
    function(x,...){

        other.GATE     <- list(x,...)
        name           <- c()
        markers        <- c()
        profiles       <- c()
        profiles.nb    <- 0
        
        i <- 1
        for(gate in other.GATE){
            if(class(gate) != "GATE")
                stop(paste0("Cannot combine objects of different classes (element at position ",i," is of type ",class(gate)))
            name           <- c(name,gate@name)
            profiles       <- c(profiles,gate@profiles)
            markers        <- union(markers,gate@markers)
            profiles.nb    <- profiles.nb+gate@profiles.nb
            i <- i+1
        }
        if(length(profiles)!=length(unique(profiles))){
            stop("Error: profile names are not unique")
        }
        
        name <- paste0(name,collapse=";")
        
        ranges           <- numeric(length(profiles)*length(markers)*2)*NA
        dim(ranges)      <- c(length(profiles),length(markers),2)
        dimnames(ranges) <- list(profiles,markers,c("inf","sup"))
        i  <- 1
        for(gate in other.GATE){ 
            nb                               <- gate@profiles.nb
            ranges[i:(i+nb-1),gate@markers,] <- gate@ranges
            i <- i+nb
        }

        gate <- GATE(name = name,
            profiles      = profiles,
            profiles.nb   = length(profiles),
            markers       = markers,
            markers.nb    = length(markers),
            ranges        = ranges)
        return(gate)
    }
)

#' @rdname c-methods
#' @export
setMethod("c",c("RES"),   
function(x,...){

    other.RES        <- list(x,...)
    comparisons      <- NULL
    markers          <- NULL
    
    i <- 1
    for(res in other.RES){
        if(class(res) != "RES")
            stop(paste0("Cannot combine objects of different classes (element at position ",i," is of type ",class(res)))
        if(any(rownames(res@comparisons) %in% rownames(comparisons)))
            stop(paste0("Duplicated comparison results found"))
        comparisons <- rbind(comparisons,res@comparisons)
        if(res@comparisons.nb>0)
            markers <- union(markers,res@markers)
        i <- i+1
    }
    rownames(comparisons) <- paste(comparisons[,1],comparisons[,2],sep="/vs/")
    
    marker.measures  <- matrix(NA,nrow=nrow(comparisons),ncol=length(markers),dimnames=list(rownames(comparisons),markers))
    marker.successes <- matrix(NA,nrow=nrow(comparisons),ncol=length(markers),dimnames=list(rownames(comparisons),markers))
    i <- 1
    for(res in other.RES){ 
        if(nrow(res@marker.measures)>1 && nrow(res@marker.successes)>1){
            nb                                       <- res@comparisons.nb
            marker.measures[i:(i+nb-1),res@markers]  <- as.matrix(res@marker.measures)
            marker.successes[i:(i+nb-1),res@markers] <- as.matrix(res@marker.successes)
            i <- i+nb
        }
    }
    marker.measures  <- data.frame(marker.measures)
    marker.successes <- data.frame(marker.successes)
    comparisons[,1]  <- as.character(comparisons[,1])
    comparisons[,2]  <- as.character(comparisons[,2])
    
    res <- RES(comparisons = comparisons,
        comparisons.nb     = nrow(comparisons),
        markers            = markers,
        marker.measures    = marker.measures,
        marker.successes   = marker.successes
        )
    return(res)
    }
)


#' @title Coercion to a CELL object
#'
#' @description Coerces a numeric matrix into a CELL object.
#'
#' This function transforms a numeric matrix into one or several cell profiles.
#'
#' @details The matrix must have its column names corresponding to the cell markers.
#'
#' @param object a numeric matrix
#' @param name a character specifying the internal name of the CELL object to create
#'
#' @return a S4 object of class CELL
#'
#' @name as.CELL
#' @rdname as.CELL-methods
#'
#' @export
setGeneric("as.CELL", function(object,name="cell") { standardGeneric("as.CELL") })

#' @rdname as.CELL-methods
#' @export
setMethod("as.CELL",c("matrix"),
    function(object,name){  
        data           <- object
        dimnames(data) <- NULL
        cell <- CELL(name = name,
            profiles      = as.character(1:nrow(object)),
            profiles.nb   = nrow(object),
            markers       = colnames(object),
            markers.nb    = ncol(object),
            intensities   = data)
        return(cell)
    }
)

#' @title Coercion to a CLUSTER object
#'
#' @description Coerces a CELL object or a numeric matrix into a CLUSTER object.
#'
#' This function transforms the cell profiles from a CELL object or from a numeric matrix into one or several cell cluster profiles by computing the means, the standard deviations, and the densities of each marker.
#'
#' @details The 'cluster' parameter is especially useful when importing FCS files containing the SPADE clustering results (where an additional channel is used to indicate the associations between cells and cell clusters) or FCS files from any other automatic gating algorithm.
#'
#' In the context of a numeric matrix coercion, the matrix must have its column names corresponding to the cell markers.
#'
#' @importFrom igraph graph.empty layout.auto
#'
#' @param object a CELL object or a numeric matrix
#' @param name a character specifying the internal name of the CLUSTER object to create
#' @param cluster a character indicating a channel name that can be used to gather the cell profiles into several cluster profiles. If a channel named is specified then the created CLUSTER object will contain as many profiles as different values present in this channel. If this parameter is NULL then a CLUSTER object with only one profile will be created
#' @param bin.width a numeric value indicating the width of the bins in the density estimation computations (default=0.05)
#'
#' @return a S4 object of class CLUSTER
#'
#' @name as.CLUSTER
#' @rdname as.CLUSTER-methods
#'
#' @export
setGeneric("as.CLUSTER",function(object,name=object@name,cluster=NULL,bin.width=0.05){standardGeneric("as.CLUSTER")})

#' @rdname as.CLUSTER-methods
#' @export
setMethod("as.CLUSTER",c("CELL"),
    function(object,name=object@name,cluster=NULL,bin.width=0.05){
        colnames(object@intensities) <- object@markers
        cluster                      <- as.CLUSTER(object=object@intensities,name=object@name,cluster=cluster,bin.width=bin.width)
        return(cluster)
    }
)

#' @rdname as.CLUSTER-methods
#' @export
setMethod("as.CLUSTER",c("matrix"),
    function(object,name="cell_cluster",cluster=NULL,bin.width=0.05){
        markers <- colnames(object)
        if(is.null(cluster)){
            profiles       <- "1"
            profiles.nb    <- as.integer(1)
            profiles.sizes <- integer(profiles.nb)
            message("Computation of marker means, standard deviations and densities...")
            means          <- matrix(apply(object,2,mean),nrow=profiles.nb,ncol=length(markers))
            sd             <- matrix(apply(object,2,function(x){sqrt(var(x))}),nrow=profiles.nb,ncol=length(markers))
            densities      <- matrix(list(),nrow=1,ncol=length(markers))
            for(i in 1:length(markers)){
                densities[[i]] <- DENSITY(bin.width=bin.width,values=object[,i],name=markers[i])
                message(paste0(round(i/length(markers)*100),"%"))
            }
            message("done")
            profiles.sizes <- nrow(object)
        }else{
            val.cluster    <- sort(unique(object[,cluster]))
            profiles       <- as.character(val.cluster)
            profiles.nb    <- length(val.cluster)
            profiles.sizes <- integer(profiles.nb)
            markers        <- markers[markers!=cluster]
            message("Computation of marker means, standard deviations and densities...")
            means          <- matrix(0,nrow=profiles.nb,ncol=length(markers))
            sd             <- matrix(0,nrow=profiles.nb,ncol=length(markers))
            densities      <- matrix(list(),nrow=profiles.nb,ncol=length(markers))
            for(i in 1:profiles.nb){
                current_profile   <- object[object[,cluster]==val.cluster[i],markers,drop=FALSE]
                profiles.sizes[i] <- nrow(current_profile)
                means[i,]         <- apply(current_profile,2,mean)
                sd[i,]            <- apply(current_profile,2,sd)
                densities[i,]     <- apply(current_profile,2,function(x){ return(DENSITY(bin.width=bin.width,values=x))})
                for(j in 1:length(markers)){
                    densities[i,j][[1]]@name <- markers[j]
                }
                message(paste0(round(i/profiles.nb*100),"%"))
            }
            message("done")
        }
        
        markers.clustering     <- rep(FALSE,length(markers))
        graph                  <- igraph::graph.empty(0,directed=FALSE)
        graph.layout           <- igraph::layout.auto(graph)
        colnames(graph.layout) <- NULL
        cluster <- CLUSTER(name = name,
            profiles            = profiles,
            profiles.nb         = profiles.nb,
            profiles.sizes      = profiles.sizes,
            markers             = markers,
            markers.nb          = length(markers),
            markers.clustering  = markers.clustering,
            means               = means,
            sd                  = sd,
            densities           = densities,
            graph               = graph,
            graph.layout        = graph.layout)
        return(cluster)
    }
)


#' @title Coercion to a GATE object
#'
#' @description Coerces a CELL or CLUSTER object into a GATE object.
#'
#' This function transforms cell or cell cluster profiles from a CELL or CLUSTER object into one or several gate profiles by computing the range of each marker. 
#'
#' @details By default the expression ranges are computed based on the 0.01 and 0.99 quantiles of the marker expression densities, but can be specified by the user. If quantiles is set to 0 and 1, then the expression ranges will correspond to the minimal and maximal values of the expression markers.
#'
#' If several cell cluster profiles are present in a CLUSTER object, then the resulting GATE object will contain as many gate profiles.
#'
#' @param object a CELL or CLUSTER object
#' @param name a character specifying the internal name of the GATE object to create
#' @param quantiles a numeric vector of two values specifying the quantiles to use for the computations of marker expression ranges
#'
#' @return a S4 object of class GATE
#'
#' @name as.GATE
#' @rdname as.GATE-methods
#'
#' @export
setGeneric("as.GATE", function(object,name=object@name,quantiles=c(0.01,0.99)) { standardGeneric("as.GATE") })

#' @rdname as.GATE-methods
#' @export
setMethod("as.GATE",c("CELL"),
    function(object,name,quantiles){
        ranges           <- apply(object@intensities,2,FUN=quantile,probs=quantiles)
        ranges           <- c(ranges[1,],ranges[2,])
        dim(ranges)      <- c(1,length(object@markers),2)
        dimnames(ranges) <- list(NULL,NULL,c("inf","sup"))
        gate <- GATE(name = name,
            profiles      = object@name,
            profiles.nb   = length(object@name),
            markers       = object@markers,
            markers.nb    = length(object@markers),
            ranges        = ranges)
        return(gate)
    }
)

#' @rdname as.GATE-methods
#' @export
setMethod("as.GATE",c("CLUSTER"),
    function(object,name,quantiles){
        ranges_inf       <- qnorm(quantiles[1],object@means,object@sd)
        ranges_sup       <- qnorm(quantiles[2],object@means,object@sd)
        ranges           <- cbind(ranges_inf,ranges_sup)
        dim(ranges)      <- c(object@profiles.nb,object@markers.nb,2)
        dimnames(ranges) <- list(NULL,NULL,c("inf","sup"))
        gate <- GATE(name = name,
            profiles      = object@profiles,
            profiles.nb   = length(object@profiles),
            markers       = object@markers,
            markers.nb    = length(object@markers),
            ranges        = ranges)
        return(gate)
    }
)
