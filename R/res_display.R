#' @title Circular graph representation of a RES object
#'
#' @description Creates a circular graph representation of the comparison results. In such graph representation, each cytometry profile is represented by a node and links between the nodes represent significant similarities or inclusions between the profiles. Nodes are positioned on a circular layout and organized based on their object names.
#'
#' @details The representation is provided as an interactive HTML file via a Scalable Vector Graphics (SVG) element created with the D3.js library. Users can interact with the representation by modifying the link tensions and the p-value cutoff for the significant similarities or inclusions.
#'
#' @importFrom RJSONIO toJSON
#'
#' @param res a RES object
#' @param filename a character specifying a file location where to save the HTML file of the representation
#' @param pvalue.th a numeric value specifying a p-value cutoff. Only the associations below this specific value will be returned
#' @param svgsize a numeric value specifying the size of the SVG representation in pixels
#'
#' @return none
#'
#' @export
res.graph <- function(res,filename="res.html",svgsize=1000,pvalue.th=0.2){
    
    profiles1 <- res@comparisons$profile1
    profiles2 <- res@comparisons$profile2
    profiles  <- unique(profiles1,profiles2)
    maxlength <- max(nchar(profiles))
    
    graph_res <- create.graph(res,pvalue.th)
    
    html <- ""
    html <- paste0(html,'<!doctype html>\n')
    html <- paste0(html,'<html>\n')
    html <- paste0(html,'<head>\n<meta charset="utf-8">\n</head>\n')
    html <- paste0(html,'<body>\n')
    html <- paste0(html,'<script type="text/javascript">\n')
    html <- paste0(html,d3.js,'\n')
    html <- paste0(html,'</script>\n')
    html <- paste0(html,'tension of links: \n')
    html <- paste0(html,'<input id="tension" type="range" min="0" max="1" value="0.80" step=".01"/>\n')
    html <- paste0(html,'<output id="label_tension" readonly>0.80</output>')
    html <- paste0(html,'<br>\n')
    html <- paste0(html,'p-value cutoff: \n')
    html <- paste0(html,'<input id="pvaluecutoff" type="range" min="0" max="1" value="0.05" step=".01">\n')
    html <- paste0(html,'<output id="label_pvaluecutoff" readonly>0.05</output>')
    html <- paste0(html,'</div>\n')
    html <- paste0(html,'<div id="graph">\n')
    html <- paste0(html,'</div>\n')
    html <- paste0(html,'<script type="text/javascript">\n')
    html <- paste0(html,'var graph = ',RJSONIO::toJSON(graph_res),';\n')
    html <- paste0(html,'var svgsize = ',svgsize,';\n')
    html <- paste0(html,'var maxlength = ',maxlength,';\n')
    html <- paste0(html,res_graph.js,'\n')
    html <- paste0(html,'drawGraph();\n')
    html <- paste0(html,'</script>\n')
    
    cat(iconv(html,to="UTF-8"),file = filename)
    
}

# title Internal - Graph representation of the comparison results
#
# @description Create a graph structure of the comparison results. In such graph structure each cytometry profile is represented by a node and links between the nodes represents significant similarities or inclusions between the profiles.
#
# @details This function is used internally when generating the D3.js circular graph representation of the comparison results. 
#
# @param res a RES object
# @param pvalue.th a numeric value specifying a p-value cutoff. Only the associations below this specific value will be returned
#
# @return a named list containing the nodes and links of the graph
create.graph <- function(res,pvalue.th){

    nodes.init <- function(res){
        return(list(name=as.character(res),children=list()))
    }
    
    nodes          <- nodes.init("root")
    nodes$children <- lapply("profiles", nodes.init)
    
    profile        <- unique(c(as.character(res@comparisons[,1]), as.character(res@comparisons[,2])))
    
    profile.type   <- do.call("rbind", strsplit(profile,":"))[,1]
    profile.name   <- do.call("rbind", strsplit(profile,":"))[,2]
    types          <- unique(paste0(profile.type,":",profile.name))
    
    nodes$children[[1]]$children <- lapply(types, nodes.init)
    
    for(i in 1:length(types)){
        objects <- profile[grep(types[i],profile)]
        nodes$children[[1]]$children[[i]]$children <- lapply(objects, nodes.init)
    }
    
    res <- res[res@comparisons[,1]!=res@comparisons[,2]]
    res <- res[res@comparisons$pvalue<pvalue.th] 
    
    links <- mapply(function(source,target,pvalue){
            return(list(source=source,target=target,pvalue=pvalue))
            },as.character(res@comparisons[,1]),as.character(res@comparisons[,2]),res@comparisons[,4],SIMPLIFY=FALSE)
    names(links) <- NULL

    return(list(node = nodes, link = links))
}

#' @title Create a Multidimensional scaling (MDS) representation of a RES object
#'
#' @description Creates a Multidimensional scaling (MDS) representation of the comparison results. In such MDS representation each cytometry profile is represented by a dot in a two-dimensional space and the distances between the nodes are proportional to the similarity measures between the profiles. The Kruskal Stress displayed at the left bottom of the representation quantifies the quality of the representation as the percentage of information lost in the dimensionality reduction process.
#'
#' @details The representation is provided as a HTML file via a Scalable Vector Graphics (SVG) element created with the D3.js library.
#'
#' @importFrom RJSONIO toJSON
#'
#' @param res a RES object
#' @param filename a file location where to save the objects
#' @param cols a character vector specifying the colours of the node in the SVG representation
#' @param sizes a numeric vector specifying the sizes of the nodes in pixels in the SVG representation
#' @param svgsize a numeric value specifying the size of the SVG representation in pixels
#'
#' @return none
#'
#' @export
res.mds <- function(res,filename="res.html",cols=NULL,sizes=NULL,svgsize=1000){

    mds_res  <- create.mds(res,cols,sizes)
    
    html <- ""
    html <- paste0(html,'<!doctype html>\n')
    html <- paste0(html,'<html>\n')
    html <- paste0(html,'<head>\n<meta charset="utf-8">\n</head>\n')
    html <- paste0(html,'<body>\n')
    html <- paste0(html,'<script type="text/javascript">\n')
    html <- paste0(html,d3.js,'\n')
    html <- paste0(html,'</script>\n')
    html <- paste0(html,'<div id="mds" align="center">\n')
    html <- paste0(html,'</div>\n')
    html <- paste0(html,'<script type="text/javascript">\n')
    html <- paste0(html,'var profiles = ',RJSONIO::toJSON(mds_res$profiles),';\n')
    html <- paste0(html,'var pos = ',RJSONIO::toJSON(mds_res$positions),';\n')
    html <- paste0(html,'var nodes_cols = ',RJSONIO::toJSON(mds_res$nodes_cols),';\n')
    html <- paste0(html,'var nodes_sizes = ',RJSONIO::toJSON(mds_res$nodes_sizes),';\n')
    html <- paste0(html,'var stress = ',mds_res$stress,';\n')
    html <- paste0(html,'var svgsize = ',svgsize,';\n')
    html <- paste0(html,res_mds.js,'\n')
    html <- paste0(html,'drawMDS();\n')
    html <- paste0(html,'</script>\n')
    html <- paste0(html,'</body>\n</html>')
    
    cat(iconv(html,to="UTF-8"),file = filename)
}

# title Internal - Multidimensional scaling representation of the comparisons results
#
# @description Computes a Multidimensional scaling (MDS) representation of the comparison results. In such MDS representation each cytometry profile is represented by a dot in a two-dimensional space and the distances between the nodes are proportional to the similarity measures between the profiles. The Kruskal Stress  quantifies the quality of the representation as the percentage of information lost in the dimensionality reduction process.
#
# @details This function is used internally when generating the D3.js MDS representation of the comparison results. The 'cols' and 'sizes' parameters can be used to specify the colors and sizes of the dots in the MDS representation.
#
# @param res a RES object
# @param cols a named character vector specifying the colours of the node in the representation
# @param sizes a named numeric vector specifying the sizes of the node in the representation
#
# @return a named list containing the positions of the profiles in the 2-dimensional MDS space and the Kruskal Stress of the MDS representation
# 
#' @import MASS 
create.mds <- function(res,cols=NULL,sizes=NULL){
    
    profiles_source <- as.character(res@comparisons[,1])
    profiles_target <- as.character(res@comparisons[,2])
    profiles        <- unique(c(profiles_source,profiles_target))
    profiles.types  <- do.call("rbind",strsplit(profiles,":"))[,1]
    
    if(length(unique(profiles.types))>1)
        stop("MDS can only be performed on profiles of same types")
    
    if(!is.null(cols) && is.null(names(cols)))
        stop("character vector cols is not named")
    if(!is.null(sizes) && is.null(names(sizes)))
        stop("character vector sizes is not named")
    
    if(!is.null(cols) && length(profiles)!=length(cols))
        stop("character vector cols do not have the same length as the number of profiles in the comparisons")
    if(!is.null(sizes) && length(profiles)!=length(sizes))
        stop("character vector sizes do not have the same length as the number of profiles in the comparisons")
        
    dist   <- matrix(0,nrow=length(profiles),ncol=length(profiles),dimnames=list(profiles,profiles))
    measure <- res@comparisons[,3]
    for(i in 1:length(measure)){
        profile1                <- as.character(res@comparisons[i,1])
        profile2                <- as.character(res@comparisons[i,2])
        dist[profile1,profile2] <- measure[i]
        dist[profile2,profile1] <- measure[i]
    }
    capture.output(mds <- MASS::isoMDS(dist))
    
    if(is.null(cols)){
        cols         <- rep("blue",length(profiles))
        names(cols)  <- profiles
    }else{
        cols  <- cols[profiles]
    }
    
    if(is.null(sizes)){
        sizes        <- rep(6,length(profiles))
        names(sizes) <- profiles
    }else{
        sizes <- sizes[profiles]
    }
    
    return(list(profiles = profiles, positions = mds$point, stress = mds$stress, nodes_cols=unname(cols), nodes_sizes=unname(sizes)))
}

