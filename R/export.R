#' @title Exportation of CELL or GATE objects
#'
#' @description Exports a CELL object into a FCS file or export a GATE object into a GatingML-XML file.
#'
#' @param object a CELL or GATE object
#' @param filename a character indicating the location of the FCS or GatingML-XML file to save
#'
#' @return none
#' 
#' @importFrom flowCore write.FCS rectangleGate
#' @import flowUtils
#'
#' @name export
#' @rdname export-methods
setGeneric("export",function(object,filename) { standardGeneric("export") })

#' @rdname export-methods
#' @export
setMethod("export",c("CELL"),
    function(object,filename="cells.fcs"){
        data                      <- sinh_coeff(object@intensities)
        flowFrame                 <- createFlowFrame(data,object@markers)
        suppressWarnings(filename <- flowCore::write.FCS(flowFrame,filename))
    }
)

#' @rdname export-methods
#' @export
setMethod("export",c("GATE"),
    function(object,filename="gates.xml"){
        flowEnv <- new.env()
        for(i in 1:object@profiles.nb){
            gate.name              <- paste(object@name,object@profiles[i], sep="")
            gate.ranges            <- t(object@ranges[i,,])
            colnames(gate.ranges)  <- object@markers
            row.names(gate.ranges) <- c("min","max")
            flowEnv[[gate.name]]   <- flowCore::rectangleGate(filterId=gate.name, .gate=gate.ranges)
        }
        flowUtils::write.gatingML(flowEnv,filename)
    }
)
