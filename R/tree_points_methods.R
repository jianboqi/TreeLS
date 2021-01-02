#' Tree points algorithm: voronoi polygons.
#' @description This function is meant to be used inside \code{\link{treePoints}}. Assign **all** points to a \code{TreeID} based on their closest \code{\link{treeMap}} coordinate.
#' @importFrom dismo voronoi
#' @importFrom sp over SpatialPoints
#' @export
trp.voronoi = function(){

  func = function(las, xymap){
    xt = extent(las) + c(-1,1,-1,1)
    v_poly = voronoi(xymap[,2:3], xt)
    v_poly$id = xymap$TreeID
    names(v_poly) = 'TreeID'
    crs(v_poly) = crs(las)

    xysp = las@data[,.(X,Y)] %>% SpatialPoints
    crs(xysp) = crs(las)
    xysp %<>% over(y = v_poly)

    xysp = xysp[,1] %>% as.double
    xysp[is.na(xysp)] = 0

    las@data$TreeID = xysp
    las@data$TreeID[las@data$Classification == 2] = 0

    las %<>% setAttribute('tree_points')

    return(las)
  }

  func %<>% setAttribute('tpt_mtd')
  return(func)
}


#' Tree points algorithm: fixed size patches.
#' @description This function is meant to be used inside \code{\link{treePoints}}. Assign points to a \code{TreeID} inside circles/squares of fixed area around \code{\link{treeMap}} coordinates.
#' @param l \code{numeric} - circle radius or square side length.
#' @param circle \code{logical} - assign \code{TreeID}s to circular (\code{TRUE}) or squared (\code{FALSE}) patches.
#' @export
trp.crop = function(l = 1, circle=TRUE){
  func = function(las, xymap){
    las@data$TreeID = treeIdsFromMap(las@data[,.(X,Y)] %>% as.matrix, xymap[,.(X,Y)] %>% as.matrix, xymap$TreeID %>% as.integer, l, circle)
    las %<>% setAttribute('tree_points')
    return(las)
  }
  func %<>% setAttribute('tpt_mtd')
  return(func)
}

#' Tree points algorithm: auto segmentation.
#' @description This function is meant to be used inside \code{\link{treePoints}}. Assign points to a nearest \code{TreeID} in \code{\link{treeMap}} coordinates.
#' @export
trp.autoseg = function(){
  minfun <- function(x){
    if(any(x < 1.797693e+308)){
      return(which.min(x));
    }else{
      return(0);
    }
  }

  func = function(las, xymap){
    xymap[["Z"]]=0;
    numoftrees = nrow(xymap);

    xyz <- las@data[,.(X,Y,Z)];
    xyz <- rbind(xymap[,c(2,3,4)], xyz);
    nnself <- knn(data=xyz, k=5, radius=0.3);
    edges <- data.frame(from_vertex=rep(1:nrow(nnself$nn.idx),4), to_vertex=as.vector(nnself$nn.idx[,2:5]),cost=as.vector(nnself$nn.dists[,2:5]))
    non_directed<-makegraph(edges,directed=FALSE);

    # simple = cpp_simplify(non_directed);
    # origin <- unique(c(edges$from_vertex,edges$to_vertex));
    # origin <- origin[-c(1:numoftrees)];
    origin <- (numoftrees+1):nrow(nnself$nn.idx);
    nodes <- 1:numoftrees;
    non_dir_dist<-get_distance_matrix(Graph=non_directed, from=origin, to=nodes, allcores=TRUE);
    IDs = apply(non_dir_dist, 1, minfun);
    las@data$TreeID = IDs
    las@data$TreeID[las@data$Classification == 2] = 0

    las %<>% setAttribute('tree_points')
    return(las)
  }
  func %<>% setAttribute('tpt_mtd')
  return(func)
}
