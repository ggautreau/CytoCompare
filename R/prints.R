#' @title Textual preview for all S4 CytoCompare objects
#'
#' @description Prints a preview for a CELL, CLUSTER, GATE, RES, MWEIGHTS or DENSITY object.
#'
#' @param x a CELL, CLUSTER, GATE, RES, MWEIGHTS or DENSITY object
#' 
#' @return none
#' 
#' @name print
#' @rdname print-methods
NULL

#' @rdname print-methods
#' @export
setMethod("print","CELL",
    function(x){
        cat("Object class: CELL\n")
        cat(paste0("Object name: ",x@name,"\n"))
        cat(paste0("Number of cell profiles: ",x@profiles.nb,"\n"))
        cat(paste0("Number of markers: ",length(x@markers),"\n"))
        cat("Markers: ")
        cat(x@markers,sep="; ")
        cat("\n")
    }
)

#' @rdname print-methods
#' @export
setMethod("print","CLUSTER",
    function(x){
        cat("Object class: CLUSTER\n")
        cat(paste0("Object name: ",x@name,"\n"))
        cat(paste0("Number of clusters profiles: ",x@profiles.nb,"\n"))
        cat("Markers: ")
        cat(x@markers,sep="; ")
        cat("\n")
        cat(paste0("Number of markers: ",length(x@markers),"\n"))
        cat("Clustering markers: ")
        if(length(x@markers[x@markers.clustering])>1)
            cat(x@markers[x@markers.clustering],sep="; ")
        else
            cat("none")
        cat("\n")
        cat(paste0("Number of clustering markers: ",sum(x@markers.clustering),"\n"))
        if(prod(dim(x@densities) == c(x@profiles.nb,length(x@markers))))
            cat(paste0("Density bin width: ",x@densities[[1]]@bin.width,"\n"))
        else
            cat(paste0("No Densities associated to marker expressions\n"))
        cat("Cluster profile names and number of associated cells:\n")
        cat(paste0(" ",x@profiles,": ",ifelse(!is.na(x@profiles.sizes),x@profiles.sizes,"no associated number of")," cells",collapse="\n"))
        cat("\n")
    }
)

#' @rdname print-methods
#' @export
setMethod("print","GATE",
    function(x){
        cat("Object class: GATE\n")
        cat(paste0("Object name: ",x@name,"\n"))
        cat(paste0("Number of gate profiles: ",x@profiles.nb,"\n"))
        cat("Markers: ")
        cat(x@markers,sep="; ")
        cat("\n")
        cat(paste0("Number of markers: ",length(x@markers),"\n"))
        cat("Gate profile names: \n")
        cat(paste0(" ",x@profiles,collapse="\n"))
        cat("\n")
    }
)

#' @rdname print-methods
#' @export
setMethod("print","RES",
    function(x){
        cat("Object class: RES\n")
        cat(paste0("Number of comparisons: ",x@comparisons.nb,"\n"))
        cat("Markers: ")
        cat(x@markers,sep="; ")
        cat("\n")
        cat(paste0("Number of markers: ",length(x@markers),"\n"))
        cat("Profiles present in the comparisons:\n")
        cat(unique(c(as.character(x@comparisons[,1]),as.character(x@comparisons[,2]))),sep="; ")
        cat("\n")
        cat("Comparisons available:\n")
        cat(rownames(x@comparisons))
        cat("\n")
    }
)

#' @rdname print-methods
#' @export
setMethod("print","MWEIGHTS",
    function(x){
        cat("Object class: MWEIGHTS\n")
        cat(paste0("Number of markers: ",length(x@weights),"\n"))
        cat("Markers with associated weights:\n")
        cat(paste0(" ",x@markers,": ",x@weights,collapse="\n"))
        cat("\n")
    }
)

#' @rdname print-methods
#' @export
setMethod("print","DENSITY",
    function(x){
        cat("Object class: DENSITY\n")
        cat(paste0("Bin widths: ",x@bin.width,"\n"))
        cat(paste0("Number of bins: ",(x@bin.nb[1]+x@bin.nb[2]),"\n"))
        cat(paste0("Minimal bin: ",x@bin.interval[1],"\n"))
        cat(paste0("Maximal bin: ",x@bin.interval[2],"\n"))
        cat(paste0("Number of negative bins: ",x@bin.nb[1],"\n"))
        cat(paste0("Number of positive bins: ",x@bin.nb[2],"\n"))
        cat(paste0("Number of values used in the density estimation: ",x@point.nb,"\n"))
    }
)


#' @title Textual preview for all S4 CytoCompare objects
#'
#' @description Shows a preview for a CELL, CLUSTER, GATE, DENSITY, MWEIGHTS or RES object.
#'
#' @param object a CELL, CLUSTER, GATE, DENSITY, MWEIGHTS or RES object
#' 
#' @return none
#' 
#' @name show
#' @rdname show-methods
NULL

#' @rdname show-methods
#' @export
setMethod("show","CELL", function(object) print(object))

#' @rdname show-methods
#' @export
setMethod("show","CLUSTER", function(object) print(object))

#' @rdname show-methods
#' @export
setMethod("show","GATE", function(object) print(object))

#' @rdname show-methods
#' @export
setMethod("show","RES", function(object) print(object))

#' @rdname show-methods
#' @export
setMethod("show","MWEIGHTS", function(object) print(object))

#' @rdname show-methods
#' @export
setMethod("show","DENSITY", function(object) print(object))
