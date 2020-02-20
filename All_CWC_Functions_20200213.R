######## CWC Simulation Functions ########

##### functions directly from gsl-functions-6-8-2018.R #####
gsl.cluster <- function(X, mdsfun, pruning.criterion = "size",
                        pruning.threshold = 0, gsl.cluster.out = NULL,
                        assign.fluff = T,
                        obs.density.lower.confidence.bounds = NULL) {
  n <- attr(mdsfun, "n")
  if (is.null(gsl.cluster.out)) {
    mast <- gsl.mast(mdsfun)
    mast.dendogram <- gsl.mast.dendogram(mast)
    obs.density <- rep(0, n)
    for (i in 1:n) obs.density[i] <- mdsfun(i, i)
    pc.values <- gsl.pruning.criteria(mast.dendogram, obs.density,
                                      obs.density.lower.confidence.bounds =
                                        obs.density.lower.confidence.bounds)
  }
  else {
    mast.dendogram <- gsl.cluster.out$mast.dendogram
    pc.values <- gsl.cluster.out$pc.values
    obs.density <- gsl.cluster.out$obs.density
  }
  cluster.tree <- gsl.compute.cluster.tree(mast.dendogram, pc.values,
                                           pruning.criterion = pruning.criterion,
                                           pruning.threshold = pruning.threshold)
  afo.out <- NULL
  if (assign.fluff) afo.out <- gsl.assign.fluff.oneshot(X, cluster.tree,
                                                        mast.dendogram, obs.density)
  
  return(list(cluster.tree = cluster.tree, leaf.code = afo.out$leaf.code,
              cluster.core = afo.out$cluster.core, mast.dendogram = mast.dendogram,
              pc.values = pc.values, obs.density = obs.density))
}
gsl.runt.size <- function(gsl.cluster.out) {
  return(rev(sort((gsl.cluster.out$cluster.tree[, "runt.size"]))))
}
gsl.runt.excess.mass <- function(gsl.cluster.out) {
  return(rev(sort((gsl.cluster.out$cluster.tree[, "runt.excess.mass"]))))
}
gsl.runt.rise <- function(gsl.cluster.out) {
  return(rev(sort((gsl.cluster.out$cluster.tree[, "runt.rise"]))))
}
gsl.runt.pruning.criterion <- function(gsl.cluster.out) {
  return(rev(sort((gsl.cluster.out$cluster.tree[, "runt.pruning.crit"]))))
}
gsl.observation.labels <- function(gsl.cluster.out) {
  return(gsl.cluster.out$leaf.code)
}
gsl.observation.in.core <- function(gsl.cluster.out) {
  return(gsl.cluster.out$cluster.core %% 2 == 1)
}
gsl.cluster.tree <- function(gsl.cluster.out) {
  cluster.tree <- gsl.cluster.out$cluster.tree
  columns <- c("leftson", "rightson", "runt.size", "runt.excess.mass",
               "leaf.code", "level")
  return(cluster.tree[,columns])
}
gsl.mdsmat <- function(X, density, min.density.ut = NULL, kmin = 1,
                       kmax = nrow(X) - 1, ngrid = 10, interactive = F,
                       include.mst = T) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  n <- nrow(X)
  p <- ncol(X)
  all <- (kmin == 1) & (kmax == (n - 1))
  if(is.null (min.density.ut)) min.density <- matrix(-1, nrow = n, ncol = n)
  else {
    min.density <- matrix(0, nrow = n, ncol = n)
    min.density[upper.tri(min.density, diag = T)] <- min.density.ut
    min.density <- min.density + t(min.density)
    diag(min.density) <- diag(min.density) / 2
    min.density[min.density == min(min.density)] <- -1
  }
  if (!all | include.mst) Dist <- gsl.interpoint.distance.matrix(X, X)
  if (!all) {
    ## Dist <- gsl.interpoint.distance.matrix(X, X)
    NN <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
      NN[i,] <- order(Dist[i,])
    }
  }
  if (!all & include.mst) {
    mst.mat.out <- gsl.mst.mat(Dist)
    mst.edges <- mst.mat.out$edges
    for (i in 1:(n-1)) {
      v1 <- mst.edges[i, 1]
      v2 <- mst.edges[i, 2]
      Xeval <- gsl.find.evaluation.points(X[v1,], matrix(X[v2,], nrow = 1),
                                          ngrid)
      phat.eval <- density(Xeval)
      min.density[v1, v2] <- min(phat.eval)
      min.density[v2, v1] <- min(phat.eval)
    }
  }
  for (i in 1:n) {
    if (all) closest <- (1:n)[-i] else closest <- NN[i, (kmin + 1):(kmax + 1)]
    n.closest <- kmax - kmin + 1
    Xeval <- gsl.find.evaluation.points(X[i,], matrix(X[closest,], nrow = n.closest),
                                        ngrid)
    phat.eval <- matrix(density(Xeval), ncol = ngrid, byrow = T)
    mins <- apply(phat.eval, 1, min) 
    min.density[i, closest] <- mins
    min.density[closest, i] <- mins
    if (interactive) cat(" ", i)
  }
  diag(min.density) <- density(X)
  actual.min <- min(min.density[min.density > 0])
  min.density[min.density < 0] <- 0.5 * actual.min
  up <- as.vector(min.density[upper.tri(min.density, diag = T)])
  return(up)
}
gsl.repeat.rows <- function(X, nrep) {
  n <- nrow(X)
  p <- ncol(X)
  vX <- as.vector(X)
  vX.rep <- rep(vX, rep(nrep, n*p))
  R <- matrix(vX.rep, ncol = p)
  return(R)
}
gsl.find.evaluation.points <- function(Xi, X.closest, ngrid) {
  n.closest <- nrow(X.closest)
  p <- length(Xi)
  tgrid <- seq(0, 1, length = ngrid)
  Xi.mat <- gsl.repeat.rows(matrix(Xi, nrow = 1), ngrid * n.closest)
  X.closest.mat <- gsl.repeat.rows(X.closest, ngrid)
  t <- rep(seq(0, 1, length = ngrid), n.closest)
  E <- (1 - t) * Xi.mat + t * X.closest.mat
  return(E)
}
gsl.make.nn.mdsfun <- function(X, on.the.fly = T) {
  if(on.the.fly) {
    X <- as.matrix(X)
    m <- ncol(X)
    mdsfun <- function(i, js) {
      Y <- matrix(X[js,], ncol = m)
      Z <- matrix(X[i,], nrow = length(js), ncol = m, byrow = T)
      d <- 1/apply((Y - Z)^2, 1, sum)
      d[d == Inf] <- 1.e50
      return(d)
    }
  }
  else {
    D <- 1/as.matrix(dist(X))
    D[D == Inf] <- 1.e50
    mdsfun <- function(i, js) {
      return(D[i, js])
    }
  }
  attr(mdsfun, "n") <- nrow(X)
  return(mdsfun)
}
gsl.make.mdsfun.from.mdsmat <- function(ut.mdsmat) {
  k <- length(ut.mdsmat)
  n <- (sqrt(1 + 8*k) - 1) / 2
  mdsmat <- matrix(0, nrow = n, ncol = n)
  mdsmat[upper.tri(mdsmat, diag = T)] <- ut.mdsmat
  mdsmat <- mdsmat + t(mdsmat)
  diag(mdsmat) <- diag(mdsmat) / 2
  mdsfun <- function(i, j) {
    return(mdsmat[i,j])
  }
  attr(mdsfun, "n") <- n
  return(mdsfun)
}
gsl.make.kernel.mdsfun <- function(X, bandwidth, kmax = 20, ngrid = 10,
                                   interactive = F) {
  density <- make.gaussian.kernel.density.estimate(X, bandwidth)
  mdsmat <- gsl.mdsmat(X, density, kmin = 1, kmax = kmax, ngrid = ngrid,
                       interactive = interactive)
  mdsfun <- gsl.make.mdsfun.from.mdsmat(mdsmat)
  return(mdsfun)
}
gsl.mst.mat <- function(D) {
  n <- nrow(D)
  edges <- matrix(0, nrow = n-1, ncol = 2)
  edge.weights <- rep(0, n-1)
  ##D[diag(D)] <- max(D) + 1
  out.points <- 2:n
  closest.in.points <- rep(1, n-1)
  dist.to.closest.in.point <- D[1, (2:n)]
  for (i in (n-1):2) {
    o <- order(dist.to.closest.in.point)
    new.in.point <- out.points[o[1]]
    edges[i, 1] <- closest.in.points[o[1]]
    edges[i, 2] <- new.in.point
    edge.weights[i] <- dist.to.closest.in.point[o[1]]
    out.points[o[1]] <- out.points[i]
    closest.in.points[o[1]] <- closest.in.points[i]
    dist.to.closest.in.point[o[1]] <- dist.to.closest.in.point[i]
    out.points <- out.points[1:(i-1)]
    closest.in.points <- closest.in.points[1:(i-1)]
    dist.to.closest.in.point <- dist.to.closest.in.point[1:(i-1)]
    dist.to.new.in.point <- D[new.in.point, out.points]
    new.in.point.closer <- (dist.to.new.in.point < dist.to.closest.in.point)
    closest.in.points[new.in.point.closer] <- new.in.point
    dist.to.closest.in.point[new.in.point.closer] <-
      dist.to.new.in.point[new.in.point.closer]
  }
  edges[1, 1] <- out.points[1]
  edges[1, 2] <- closest.in.points[1]
  edge.weights[1] <- dist.to.closest.in.point[1]
  return(list(edges = edges, edge.weights = edge.weights))
}
gsl.mst <- function(Dfun) {
  n <- attr(Dfun, "n")
  edges <- matrix(0, nrow = n-1, ncol = 2)
  edge.weights <- rep(0, n-1)
  out.points <- 2:n
  closest.in.points <- rep(1, n-1)
  dist.to.closest.in.point <- Dfun(1, (2:n))
  for (i in (n-1):2) {
    o <- order(dist.to.closest.in.point)
    new.in.point <- out.points[o[1]]
    edges[i, 1] <- closest.in.points[o[1]]
    edges[i, 2] <- new.in.point
    edge.weights[i] <- dist.to.closest.in.point[o[1]]
    out.points[o[1]] <- out.points[i]
    closest.in.points[o[1]] <- closest.in.points[i]
    dist.to.closest.in.point[o[1]] <- dist.to.closest.in.point[i]
    out.points <- out.points[1:(i-1)]
    closest.in.points <- closest.in.points[1:(i-1)]
    dist.to.closest.in.point <- dist.to.closest.in.point[1:(i-1)]
    dist.to.new.in.point <- Dfun(new.in.point, out.points)
    new.in.point.closer <- (dist.to.new.in.point < dist.to.closest.in.point)
    closest.in.points[new.in.point.closer] <- new.in.point
    dist.to.closest.in.point[new.in.point.closer] <-
      dist.to.new.in.point[new.in.point.closer]
    ## cat(" ", i)
  }
  edges[1, 1] <- out.points[1]
  edges[1, 2] <- closest.in.points[1]
  edge.weights[1] <- dist.to.closest.in.point[1]
  return(list(edges = edges, edge.weights = edge.weights))
}
gsl.mast <- function(Sfun) {
  Dfun <- function(...) {
    return(- Sfun(...))
  }
  attr(Dfun, "n") <- attr(Sfun, "n")
  mst <- gsl.mst(Dfun)
  mast <- mst
  mast$edge.weights <- -mst$edge.weights
  return(mast)
}
gsl.mst.dendogram <- function(gsl.mst.out) {
  mst.edges <- gsl.mst.out$edges
  mst.edge.weights <- gsl.mst.out$edge.weights
  n <- nrow(mst.edges) + 1
  weight.order <- order(mst.edge.weights)
  group.members <- vector("list", n)
  for (i in 1:n) group.members[[i]] <- i
  group.size <- rep(1, n)
  group.index <- 1:n
  merge <- matrix(0, nrow = n-1, ncol = 2)
  edges <- matrix(0, nrow = n-1, ncol = 2)
  height <- rep(0, n-1)
  ## runt.size <- rep(0, n-1)
  group.row <- rep(0, n)
  for (i in 1:(n-1)) {
    e <- mst.edges[weight.order[i],]
    v1 <- e[1]
    v2 <- e[2]
    edges[i, 1] <- v1
    edges[i, 2] <- v2
    height[i] <- mst.edge.weights[weight.order[i]]
    gi1 <- group.index[v1]
    gi2 <- group.index[v2]
    s1 <- group.size[gi1]
    s2 <- group.size[gi2]
    ## runt.size[i] <- min(s1, s2)
    if (s1 == 1) merge[i, 1] <- -group.members[[gi1]]
    else merge[i, 1] <- group.row[gi1]
    if (s2 == 1) merge[i, 2] <- -group.members[[gi2]]
    else merge[i, 2] <- group.row[gi2]
    group.members[[gi1]] <- c(group.members[[gi1]], group.members[[gi2]])
    group.size[gi1] <- s1 + s2
    group.index[group.index == gi2] <- gi1
    group.row[gi1] <- i
  }
  return(list(merge = merge, height = height, edges = edges))
}
gsl.mast.dendogram <- function(gsl.mast.out) {
  gsl.mast.out$edge.weights <- - gsl.mast.out$edge.weights
  dend <- gsl.mst.dendogram(gsl.mast.out)
  dend$height <- - dend$height
  return(dend)
}
gsl.pruning.criteria <- function(mast.dendogram, obs.density,
                                 obs.density.lower.confidence.bounds = NULL) {
  merge <- mast.dendogram$merge
  height <- mast.dendogram$height
  left.size <- rep(0, nrow(merge))
  right.size <- rep(0, nrow(merge))
  left.excess.mass <- rep(0, nrow(merge))
  right.excess.mass <- rep(0, nrow(merge))
  left.rise <- rep(0, nrow(merge))
  right.rise <- rep(0, nrow(merge))
  left.bootstrap.rise <- rep(0, nrow(merge))
  right.bootstrap.rise <- rep(0, nrow(merge))
  pc.int <- function(inode) {
    ileft <- merge[inode, 1]
    iright <- merge[inode, 2]
    if (ileft < 0) left.leaves <- -ileft
    if (ileft > 0) left.leaves <- pc.int(ileft)
    if (iright < 0) right.leaves <- -iright
    if (iright > 0) right.leaves <- pc.int(iright)
    left.high.density.cluster <- left.leaves[obs.density[left.leaves] >=
                                               height[inode]]
    right.high.density.cluster <- right.leaves[obs.density[right.leaves] >=
                                                 height[inode]]
    ## if (!is.null(left.high.density.cluster)) {
    if (length(left.high.density.cluster) == 0) browser()
    if (length(left.high.density.cluster) > 0) {
      left.size[inode] <<- length(left.high.density.cluster)
      left.excess.mass[inode] <<- sum(1 - height[inode] /
                                        obs.density[left.high.density.cluster])/
        (nrow(merge) + 1)
      left.max <- max(obs.density[left.high.density.cluster])
      if (left.max == Inf) left.rise[inode] <<- Inf
      else left.rise[inode] <<- left.max - height[inode]
      ##left.rise[inode] <<- max(obs.density[left.high.density.cluster]) -
      ##  height[inode]
      ## if (left.rise[inode] == -Inf) browser()
      ##if (left.rise[inode] == Inf) browser()
      if (!is.null(obs.density.lower.confidence.bounds)) 
        left.bootstrap.rise[inode] <<-
        max(obs.density.lower.confidence.bounds[left.high.density.cluster]) -
        height[inode]
    }
    ## if (!is.null(right.high.density.cluster)) {    ## Changed 6-6-2018
    if (length(right.high.density.cluster) == 0) browser()
    if (length(right.high.density.cluster) > 0) {
      right.size[inode] <<- length(right.high.density.cluster)
      right.excess.mass[inode] <<- sum(1 - height[inode] /
                                         obs.density[right.high.density.cluster])/
        (nrow(merge) + 1)
      right.max <- max(obs.density[right.high.density.cluster])
      if (right.max == -Inf) browser()  ## 6-6-2018
      if (right.max == Inf) right.rise[inode] <<- Inf
      else right.rise[inode] <<- right.max - height[inode]
      ##right.rise[inode] <<- max(obs.density[right.high.density.cluster]) -
      ##  height[inode]
      ## if (right.rise[inode] == -Inf) browser()  
      ## if (right.rise[inode] == Inf) browser()
      if (!is.null(obs.density.lower.confidence.bounds))
        right.bootstrap.rise[inode] <<-
        max(obs.density.lower.confidence.bounds[right.high.density.cluster]) -
        height[inode]
    }
    return(c(left.leaves, right.leaves))
  }
  pc.int(nrow(merge))
  res <- cbind(left.size, right.size, left.excess.mass, right.excess.mass,
               left.rise, right.rise, left.bootstrap.rise, right.bootstrap.rise)
  colnames(res) <- c("left.size", "right.size", "left.excess.mass", "right.excess.mass",
                     "left.rise", "right.rise", "left.bootstrap.rise",
                     "right.bootstrap.rise")
  return(res)
}
gsl.compute.cluster.tree <- function(dendogram, gsl.pc.out, pruning.criterion = "size",
                                     pruning.threshold = 0) {
  runt.size <- pmin(gsl.pc.out[,"left.size"], gsl.pc.out[,"right.size"])
  runt.excess.mass <- pmin(gsl.pc.out[,"left.excess.mass"],
                           gsl.pc.out[,"right.excess.mass"])
  runt.rise <- pmin(gsl.pc.out[,"left.rise"], gsl.pc.out[,"right.rise"])
  runt.bootstrap.rise <- pmin(gsl.pc.out[,"left.bootstrap.rise"],
                              gsl.pc.out[,"right.bootstrap.rise"])
  if (pruning.criterion == "size") runt.crit <- runt.size
  if (pruning.criterion == "excess.mass") runt.crit <- runt.excess.mass
  if (pruning.criterion == "rise") runt.crit <- runt.rise
  if (pruning.criterion == "bootstrap.rise") runt.crit <- runt.bootstrap.rise
  ##browser()
  merge <- dendogram$merge
  level <- dendogram$height
  nobs <- nrow(merge) + 1
  
  cluster.tree <- matrix(0, nrow = 2 * nobs - 1, ncol = 14)
  colnames(cluster.tree) <- c("leftson", "rightson", "runt.size", "runt.excess.mass",
                              "leaf.code", "level", "dendogram.node", "size",
                              "excess.mass", "rise", "bootstrap.rise",
                              "runt.rise", "runt.bootstrap.rise", "runt.pruning.crit")
  nnodes <- 1
  cluster.tree[1, "size"] <- nobs
  cluster.tree[1, "excess.mass"] <- 1
  cluster.tree[1, "level"] <- 0
  cluster.tree[1, "leaf.code"] <- 1
  cluster.tree[1, "dendogram.node"] <- nrow(merge)
  
  large.runt <- (runt.crit >= pruning.threshold) & (runt.rise > 0)
  ## large.runt <- (runt.crit >= pruning.threshold)
  
  find.large.runt.node.in.subdend <- function(root) {
    ## Note: Because of the monotonicity of the pruning criteria either
    ## the root has large runt or at most one of the two subdendograms.
    if (root < 0) return(0)
    if (large.runt[root]) return(root)
    large.runt.node <- 0
    if (merge[root, 1] > 0)
      large.runt.node <- find.large.runt.node.in.subdend(merge[root, 1])
    if ((large.runt.node == 0) & (merge[root, 2] > 0))
      large.runt.node <-
      find.large.runt.node.in.subdend(merge[root, 2])
    return(large.runt.node)
  }
  make.cluster.tree.sons <- function(parent) {
    ## browser()
    if (cluster.tree[parent, "dendogram.node"] < 0) return()
    large.runt.node <-
      find.large.runt.node.in.subdend(cluster.tree[parent, "dendogram.node"])
    if (large.runt.node > 0) {
      cluster.tree[parent, "dendogram.node"] <<- large.runt.node
      cluster.tree[parent, "runt.size"] <<- runt.size[large.runt.node]
      cluster.tree[parent, "runt.excess.mass"] <<- runt.excess.mass[large.runt.node]
      cluster.tree[parent, "runt.rise"] <<- runt.rise[large.runt.node]
      cluster.tree[parent, "runt.bootstrap.rise"] <<- runt.bootstrap.rise[large.runt.node]
      cluster.tree[parent, "runt.pruning.crit"] <<- runt.crit[large.runt.node]
      nnodes <<- nnodes + 1
      cluster.tree[parent, "leftson"] <<- nnodes
      cluster.tree[nnodes, "leaf.code"] <<- 2 * cluster.tree[parent, "leaf.code"]
      cluster.tree[nnodes, "size"] <<- gsl.pc.out[large.runt.node, "left.size"]
      cluster.tree[nnodes, "excess.mass"] <<-
        gsl.pc.out[large.runt.node, "left.excess.mass"]
      cluster.tree[nnodes, "rise"] <<-
        gsl.pc.out[large.runt.node, "left.rise"]
      cluster.tree[nnodes, "bootstrap.rise"] <<-
        gsl.pc.out[large.runt.node, "left.bootstrap.rise"]
      cluster.tree[nnodes, "level"] <<- level[large.runt.node]
      cluster.tree[nnodes, "dendogram.node"] <<- merge[large.runt.node, 1]
      left.parent <- nnodes
      make.cluster.tree.sons(left.parent)
      nnodes <<- nnodes + 1
      cluster.tree[parent, "rightson"] <<- nnodes
      cluster.tree[nnodes, "leaf.code"] <<- 2 * cluster.tree[parent, "leaf.code"] + 1
      cluster.tree[nnodes, "size"] <<- gsl.pc.out[large.runt.node, "right.size"]
      cluster.tree[nnodes, "excess.mass"] <<-
        gsl.pc.out[large.runt.node, "right.excess.mass"]
      cluster.tree[nnodes, "rise"] <<-
        gsl.pc.out[large.runt.node, "right.rise"]
      cluster.tree[nnodes, "bootstrap.rise"] <<-
        gsl.pc.out[large.runt.node, "right.bootstrap.rise"]
      cluster.tree[nnodes, "level"] <<- level[large.runt.node]
      cluster.tree[nnodes, "dendogram.node"] <<- merge[large.runt.node, 2]
      right.parent <- nnodes
      make.cluster.tree.sons(right.parent)
    }
    return()
  }
  make.cluster.tree.sons(1)
  if (nnodes == 1) {
    ct.out <- matrix(cluster.tree[1,], nrow = 1, byrow = T)
    colnames(ct.out) <- c("leftson", "rightson", "runt.size", "runt.excess.mass",
                          "leaf.code", "level", "dendogram.node", "size",
                          "excess.mass", "rise", "bootstrap.rise",
                          "runt.rise", "runt.bootstrap.rise", "runt.pruning.crit")
    return(ct.out)
  }
  return(cluster.tree[1:nnodes,])
}
gsl.assign.fluff.oneshot <- function(X, cluster.tree, dendogram, obs.density) {
  ## browser()
  merge <- dendogram$merge
  n <- nrow(merge) + 1
  m <- ncol(X)
  obs.leaf.code <- rep(1, n)
  obs.in.core <- rep(1, n)
  collect.obs.in.core <- function(dendogram.node, level) {
    is.in.core <- rep(F, n)
    collect.descendents.of.dendogram.node <- function(dendogram.node) {
      ## cat(" ", dendogram.node)
      if (dendogram.node < 0) {
        if (obs.density[-dendogram.node] > level)
          is.in.core[-dendogram.node] <<- T
        return()
      }
      leftson <- merge[dendogram.node, 1]
      rightson <- merge[dendogram.node, 2]
      collect.descendents.of.dendogram.node(leftson)
      collect.descendents.of.dendogram.node(rightson)
    }
    collect.descendents.of.dendogram.node(dendogram.node)
    return(is.in.core)
  }
  preserve.ncol <- function(x, m) {
    if (is.matrix(x)) return(x)
    else return(matrix(x, ncol = m, byrow = T))
  }
  afo.internal <- function(node.index) {
    ##browser()
    if (cluster.tree[node.index, "leftson"] == 0) return()
    dendogram.node <- cluster.tree[node.index, "dendogram.node"]
    node.leaf.code <- cluster.tree[node.index, "leaf.code"]
    leftson <- cluster.tree[node.index, "leftson"]
    rightson <- cluster.tree[node.index, "rightson"]
    ## cat("\n", node.index, leftson, rightson, "\n")
    leftson.leaf.code <- cluster.tree[leftson, "leaf.code"]
    rightson.leaf.code <- cluster.tree[rightson, "leaf.code"]
    in.node <- obs.leaf.code == node.leaf.code
    ## browser()
    in.left.core <- collect.obs.in.core(dendogram$merge[dendogram.node, 1],
                                        dendogram$height[dendogram.node])
    in.right.core <- collect.obs.in.core(dendogram$merge[dendogram.node, 2],
                                         dendogram$height[dendogram.node])
    obs.leaf.code[in.left.core] <<- leftson.leaf.code
    obs.in.core[in.left.core] <<- 2 * obs.in.core[in.left.core] + 1
    obs.leaf.code[in.right.core] <<- rightson.leaf.code
    obs.in.core[in.right.core] <<- 2 * obs.in.core[in.right.core] + 1
    ## browser()
    if (sum(in.left.core) + sum(in.right.core) < sum(in.node)) {
      X.train <- rbind(preserve.ncol(X[in.left.core,], m),
                       preserve.ncol(X[in.right.core,], m))
      y.train <- c(rep(1, sum(in.left.core)), rep(2, sum(in.right.core)))
      in.fluff <- in.node & !(in.left.core | in.right.core)
      obs.in.core[in.fluff] <<- 2 * obs.in.core[in.fluff]
      X.test <- preserve.ncol(X[in.fluff,], m)
      y.test <- as.vector(gsl.knn.classifier(X.train, y.train, X.test, kk = 1)$ypred)
      obs.leaf.code[in.fluff][y.test == 1] <<- leftson.leaf.code
      obs.leaf.code[in.fluff][y.test == 2] <<- rightson.leaf.code
    }
    afo.internal(leftson)
    afo.internal(rightson)
  }
  afo.internal(1)
  return(list(leaf.code = obs.leaf.code, cluster.core = obs.in.core))
}
gsl.cluster.dendogram <- function(gsl.cluster.out, leaf.level.rule = "max", hang = 0.1) {
  cluster.tree <- gsl.cluster.out$cluster.tree
  mast.dendogram <- gsl.cluster.out$mast.dendogram
  obs.density <- gsl.cluster.out$obs.density
  if (sum(obs.density) == Inf) leaf.level.rule <- "constant"
  highest.level <- max(cluster.tree[, "level"])
  offset <- hang * highest.level
  obs.leaf.code <- gsl.cluster.out$leaf.code
  obs.in.core <- gsl.cluster.out$cluster.core
  cluster.dendogram <- cluster.tree
  convert.internal <- function(inode) {
    leftson <- cluster.tree[inode, "leftson"]
    rightson <- cluster.tree[inode, "rightson"]
    if (leftson == 0) {
      node.leaf.code <- cluster.tree[inode, "leaf.code"]
      obs.in.node <- gsl.select.obs.in.node(obs.leaf.code, obs.in.core, node.leaf.code)
      cluster.core.density <- obs.density[obs.in.node$in.core]
      if (leaf.level.rule == "mean") leaf.level <- mean(cluster.core.density)
      if (leaf.level.rule == "max") leaf.level <- max(cluster.core.density)
      if (leaf.level.rule == "constant")
        leaf.level <- (1 + hang) * highest.level
      cluster.dendogram[inode, "level"] <<- leaf.level
    }
    else {
      cluster.dendogram[inode, "level"] <<- cluster.tree[leftson, "level"]
      convert.internal(leftson)
      convert.internal(rightson)
    }
  }
  convert.internal(1)
  return(cluster.dendogram)
}
gsl.select.obs.in.node <- function(obs.leaf.code, obs.in.core, node.leaf.code) {
  n <- length(obs.leaf.code)
  node.leaf.code <- rep(node.leaf.code, n)
  delta.floor <- floor(log(obs.leaf.code, 2)) - floor(log(node.leaf.code, 2))
  in.node <- floor(obs.leaf.code / 2^delta.floor) == node.leaf.code
  in.core <- in.node & (floor(obs.in.core / 2^delta.floor) %% 2 == 1)
  return(list(in.node = in.node, in.core = in.core))
}
gsl.safe.lsfit <- function(X, y, tol = 1.e-6) {
  n <- nrow(X)
  p <- ncol(X)
  X1 <- cbind(rep(1,n), X)
  X1.svd <- svd(X1)
  U <- X1.svd$u
  V <- X1.svd$v
  d <- X1.svd$d
  dpos <- d[d/d[1] > tol]
  rank <- length(dpos)
  ytilde <- t(U[,1:rank]) %*% y
  z <- ytilde / dpos
  b <- V[,1:rank] %*% z
  return(list(coef = b, residuals = y - X1 %*% b))
}
gsl.knn.classifier <- function(X.train, y.train, X.test, kk = 1,
                               pi = rep(1/K, K), CV = F, block.size = 200) {
  ## browser()
  n.train <- nrow(X.train)
  n.test <- nrow(X.test)
  n.block <- n.test %/% block.size
  m <- ncol(X.train)
  if (n.block * block.size < n.test) n.block <- n.block + 1
  K <- max(y.train)
  ypred <- matrix(0, nrow = n.test, ncol = length(kk))
  weight <- vector("numeric", K)
  for (i in 1:K) {
    ni <- sum(y.train==i)
    if (ni > 0) {weight[i] <- pi[i] * n.train / ni}
    else {weight[i] <- 1}
  } 
  for (block in 1:n.block) {
    i.start <- (block - 1) * block.size + 1
    i.end <- min(block * block.size, n.test)
    ## cat("\ni.start = ", i.start)
    ## if (i.end > i.start) X.test.block <- X.test[i.start:i.end,]
    ## Changed 5-20-09
    if (i.end > i.start) {
      X.test.block <- matrix(X.test[i.start:i.end,], ncol = m)
    } else {
      X.test.block <- matrix(X.test[i.start,], ncol = m)
    }
    ## if (ncol(X.train) != ncol(X.test.block)) browser()
    dist <- gsl.interpoint.distance.matrix(X.train, X.test.block)
    if (CV) {dist[row(dist) == (col(dist) + (i.start - 1))] <- max(dist) + 1}
    # cat("Computed distance matrix \n")
    permut <- apply(dist, 2, order)
    ordered.class.ids <- matrix(y.train[permut], nrow = n.train)
    counts <- matrix(0, nrow = K, ncol = i.end - i.start + 1)
    ypred.block <- matrix(0, nrow = i.end - i.start + 1, ncol = length(kk))
    # cat("Sorted distances \n")
    # browser()
    for(j in 1:length(kk)) {
      k <- kk[j]
      for(i in 1:K){
        #browser()
        if(k==1) {counts[i,] <- weight[i] * (ordered.class.ids[1,] == i)}
        else {counts[i,] <- weight[i] * 
          apply((ordered.class.ids[1:k,] == i), 2, sum)}
      }
      #browser()
      ypred.block[,j] <- apply(counts, 2, which.max)
      counts <- sweep(counts, 2, apply(counts, 2, sum), "/")
    }
    ##browser()
    ypred[i.start:i.end,] <- ypred.block
  }
  return(list(ypred = ypred))
}
gsl.interpoint.distance.matrix <- function(X, Y) {
  XY <- rbind(X, Y)
  XY.dist <- as.matrix(dist(XY))
  row.ind <- 1:nrow(X)
  col.ind <- (nrow(X) + 1):(nrow(X) + nrow(Y))
  Dist <- XY.dist[row.ind, col.ind]
  if (nrow(X) == 1) Dist <- matrix(Dist, nrow = 1)
  if (nrow(Y) == 1) Dist <- matrix(Dist, ncol = 1)
  return(Dist^2)
}
gsl.draw.cluster.tree <- function(gsl.cluster.out, draw.cluster.numbers = F,
                                  draw.cluster.leaf.codes = F,
                                  draw.runt.excess.mass = F,
                                  draw.runt.size = F) {
  ## browser()
  dct.out <- gsl.draw.tree(gsl.cluster.out$cluster.tree,
                           draw.cluster.numbers = draw.cluster.numbers,
                           draw.cluster.leaf.codes = draw.cluster.leaf.codes,
                           draw.runt.excess.mass = draw.runt.excess.mass,
                           draw.runt.size = draw.runt.size)
  return(dct.out)
}
gsl.draw.cluster.dendogram <- function(gsl.cluster.out, leaf.level.rule = "max",
                                       hang = 0.1, ...) {
  cd.out <- gsl.cluster.dendogram(gsl.cluster.out,
                                  leaf.level.rule = leaf.level.rule,
                                  hang = hang)
  dct.out <- gsl.draw.tree(cd.out, ...)
  return(dct.out)
}
gsl.gkl.diagnostic.plot <- function(X, gsl.cluster.out, tree.type = "cluster.tree",
                                    n.breaks = 20, ...){
  par(mfrow = c(3,1))
  cluster.tree <- gsl.cluster.out$cluster.tree
  obs.leaf.code <- gsl.cluster.out$leaf.code
  obs.cluster.core <- gsl.cluster.out$cluster.core
  if (tree.type == "cluster.tree") {
    node.xy <- gsl.draw.cluster.tree(gsl.cluster.out, ...)
    title("Cluster tree")
  }
  else {
    node.xy <- gsl.draw.cluster.dendogram(gsl.cluster.out, ...)
    title("Cluster dendogram")
  }
  #title(sub = "(c)")
  selected.node <- identify(node.xy$x, node.xy$y, labels = rep(" ", length(node.xy$x)),
                            n = 1)
  points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
         col = "blue1")
  left.son <- cluster.tree[selected.node, "leftson"]
  right.son <- cluster.tree[selected.node, "rightson"]
  if (left.son == 0) {
    cat("\nYou should select an interior node\n")
    return()
  }
  ## browser()
  m <- ncol(X)
  selected.leaf.code <- cluster.tree[selected.node, "leaf.code"]
  left.son.leaf.code <- cluster.tree[left.son, "leaf.code"]
  right.son.leaf.code <- cluster.tree[right.son, "leaf.code"]
  selected <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, selected.leaf.code)
  left <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, left.son.leaf.code)
  right <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, right.son.leaf.code)
  is.left.in.core <- left$in.core
  is.right.in.core <- right$in.core
  is.selected <- selected$in.node
  is.in.core <- selected$in.core
  if (sum(is.left.in.core) == 1) left.X <- matrix(X[is.left.in.core,], nrow = 1) else
    ## left.X <- X[is.left.in.core,]
    left.X <- matrix(X[is.left.in.core,], ncol = m)
  if (sum(is.right.in.core) == 1) right.X <- matrix(X[is.right.in.core,], nrow = 1) else
    ## right.X <- X[is.right.in.core,]
    right.X <- matrix(X[is.right.in.core,], ncol = m)
  X.train <- rbind(left.X, right.X)
  y.train <- c(rep(0, nrow(left.X)), rep(1, nrow(right.X)))
  lsmod <- gsl.safe.lsfit(X.train, y.train)
  yhat.selected <- matrix(X[is.selected,], ncol = m) %*% lsmod$coef[2:(ncol(X)+1)]
  yhat.in.core <- matrix(X[is.in.core,], ncol = m) %*% lsmod$coef[2:(ncol(X)+1)]
  breaks <- seq(from = min(yhat.selected), to = max(yhat.selected),
                length = n.breaks)
  hist(yhat.selected, breaks = breaks,  xaxt = "n",
       xlab = "", col = "blue1", main = "Observations in node")
  hist(yhat.in.core, breaks = breaks,  xaxt = "n",
       xlab = "", col = "blue1", main = "Observations in core")
  par(mfrow = c(1,1))
  return()
}
gsl.all.gkl.diagnostic.plots <- function(X, gsl.cluster.out, tree.type = "cluster.tree",
                                         n.breaks = 20,...){
  par(mfrow = c(3,1))
  cluster.tree <- gsl.cluster.out$cluster.tree
  obs.leaf.code <- gsl.cluster.out$leaf.code
  obs.cluster.core <- gsl.cluster.out$cluster.core
  agdp.internal <- function(selected.node) {
    left.son <- cluster.tree[selected.node, "leftson"]
    right.son <- cluster.tree[selected.node, "rightson"]
    if (left.son == 0) return()
    ##browser()
    if (tree.type == "cluster.tree") {
      node.xy <- gsl.draw.cluster.tree(gsl.cluster.out, ...)
      title("Cluster tree")
    }
    else {
      node.xy <- gsl.draw.cluster.dendogram(gsl.cluster.out, ...)
      title("Cluster dendogram")
    }
    points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
           col = "blue1")
    selected.leaf.code <- cluster.tree[selected.node, "leaf.code"]
    left.son.leaf.code <- cluster.tree[left.son, "leaf.code"]
    right.son.leaf.code <- cluster.tree[right.son, "leaf.code"]
    selected <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, selected.leaf.code)
    left <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, left.son.leaf.code)
    right <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, right.son.leaf.code)
    is.left.in.core <- left$in.core
    is.right.in.core <- right$in.core
    is.selected <- selected$in.node
    is.in.core <- selected$in.core
    if (sum(is.left.in.core) == 1) left.X <- matrix(X[is.left.in.core,], nrow = 1)
    else left.X <- X[is.left.in.core,]
    if (sum(is.right.in.core) == 1) right.X <- matrix(X[is.right.in.core,], nrow = 1)
    else right.X <- X[is.right.in.core,]
    X.train <- rbind(left.X, right.X)
    y.train <- c(rep(0, nrow(left.X)), rep(1, nrow(right.X)))
    lsmod <- gsl.safe.lsfit(X.train, y.train)
    yhat.selected <- X[is.selected,] %*% lsmod$coef[2:(ncol(X)+1)]
    yhat.in.core <- X[is.in.core,] %*% lsmod$coef[2:(ncol(X)+1)]
    breaks <- seq(from = min(yhat.selected), to = max(yhat.selected),
                  length = n.breaks)
    hist(yhat.selected, breaks = breaks,  xaxt = "n",
         xlab = "", col = "blue1", main = "Observations in node")
    hist(yhat.in.core, breaks = breaks,  xaxt = "n",
         xlab = "", col = "blue1", main = "Observations in core")
    bla <- readline("\nHit return to see next plot or any character and return to exit ")
    if (bla != "") stop
    agdp.internal(left.son)
    agdp.internal(right.son)
  }
  agdp.internal(1)
  par(mfrow = c(1,1))
}
gsl.two.dim.diagnostic.plot <- function(X, gsl.cluster.out, tree.type = "cluster.tree",
                                        ...){
  cluster.tree <- gsl.cluster.out$cluster.tree
  obs.leaf.code <- gsl.cluster.out$leaf.code
  obs.cluster.core <- gsl.cluster.out$cluster.core
  if (tree.type == "cluster.tree") {
    node.xy <- gsl.draw.cluster.tree(gsl.cluster.out, ...)
    title("Cluster tree")
  }
  else {
    node.xy <- gsl.draw.cluster.dendogram(gsl.cluster.out, ...)
    title("Cluster dendogram")
  }
  selected.node <- identify(node.xy$x, node.xy$y, labels = rep(" ", length(node.xy$x)),
                            n = 1)[1]
  points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
         col = "black")
  left.son <- cluster.tree[selected.node, "leftson"]
  right.son <- cluster.tree[selected.node, "rightson"]
  if (left.son == 0) {
    cat("\nYou should select an interior node\n")
    return()
  }
  points(node.xy$x[left.son], node.xy$y[left.son], pch = 18, cex = 2, 
         col = "red1")
  points(node.xy$x[right.son], node.xy$y[right.son], pch = 18, cex = 2, 
         col = "blue1")
  left.son.leaf.code <- cluster.tree[left.son, "leaf.code"]
  right.son.leaf.code <- cluster.tree[right.son, "leaf.code"]
  left <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, 
                                 left.son.leaf.code)
  right <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, 
                                  right.son.leaf.code)
  is.left.not.core <- left$in.node & !left$in.core
  is.left.core <- left$in.node & left$in.core
  is.right.not.core <- right$in.node & !right$in.core
  is.right.core <- right$in.node & right$in.core
  plot(X[,1], X[,2], type = "n", xlab = "", ylab = "")
  points(X[left$in.node, 1], X[left$in.node, 2], pch = 20, col = "red1")
  points(X[left$in.core, 1], X[left$in.core, 2], pch = 19, col = "red1")
  points(X[right$in.node, 1], X[right$in.node, 2], pch = 20, col = "blue1")
  points(X[right$in.core, 1], X[right$in.core, 2], pch = 19, col = "blue1")
  points(X[!left$in.node & !right$in.node, 1],
         X[!left$in.node & !right$in.node, 2], pch = 18, col = "black")
  return(selected.node)
}
gsl.visually.select.observations.in.node <- function(gsl.cluster.out) {
  cluster.tree <- gsl.cluster.out$cluster.tree
  obs.leaf.code <- gsl.cluster.out$leaf.code
  obs.cluster.core <- gsl.cluster.out$cluster.core
  node.xy <- gsl.draw.cluster.tree(gsl.cluster.out, 1)
  selected.node <- identify(node.xy$x, node.xy$y,
                            labels = rep(" ", length(node.xy$x)), n = 1)
  points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
         col = "black")
  selected.leaf.code <- cluster.tree[selected.node, "leaf.code"]
  selected <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, selected.leaf.code)
  return(list(in.node = selected$in.node, in.core = selected$in.core))
}
gsl.subtree.gkl.diagnostic.plot <- function(X, gsl.cluster.out, ...){
  n <- nrow(X)
  m <- ncol(X)
  cluster.tree <- gsl.cluster.out$cluster.tree
  obs.leaf.code <- gsl.cluster.out$leaf.code
  obs.cluster.core <- gsl.cluster.out$cluster.core
  par(mfrow = c(2, 1))
  node.xy <- gsl.draw.cluster.tree(gsl.cluster.out, ...)
  title("Cluster tree")
  selected.node <- identify(node.xy$x, node.xy$y, labels = rep(" ", length(node.xy$x)),
                            n = 1)[1]
  points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
         col = "blue1")
  left.daughter <- cluster.tree[selected.node, 1]
  if (left.daughter == 0) {
    cat("\nYou should select an interior node\n")
    return()
  }
  find.leaf.inodes.for.subtree <- function(inode, ifs = NULL) {
    leftson <- cluster.tree[inode, 1]
    rightson <- cluster.tree[inode, 2]
    if (leftson == 0) ifs <- c(ifs, inode)
    else {
      ifs <- find.leaf.inodes.for.subtree(leftson, ifs)
      ifs <- find.leaf.inodes.for.subtree(rightson, ifs)
    }
    return(ifs)
  }
  leaf.inodes <- find.leaf.inodes.for.subtree(selected.node) 
  leaf.codes <- cluster.tree[leaf.inodes, "leaf.code"]
  n.leaves <- length(leaf.codes)
  X.in.leaves <- matrix(0, nrow = n, ncol = m)
  y.in.leaves <- rep(0, n)
  core.in.leaves <- rep(0, n)
  n.in.leaves <- 0
  for (i in 1:n.leaves) {
    lc <- leaf.codes[i]
    selected <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, lc)
    in.leaf <- selected$in.node
    in.leaf.core <- selected$in.core
    n.in.leaf <- sum(in.leaf)
    ind.range <- (n.in.leaves + 1):(n.in.leaves + n.in.leaf)
    X.in.leaves[ind.range,] <- X[in.leaf,]
    y.in.leaves[ind.range] <- i
    core.in.leaves[ind.range] <- in.leaf.core[in.leaf]
    n.in.leaves <- n.in.leaves + n.in.leaf
    ## browser()
  }
  X.in.leaves <- X.in.leaves[1:n.in.leaves,]
  y.in.leaves <- y.in.leaves[1:n.in.leaves]
  core.in.leaves <- core.in.leaves[1:n.in.leaves]
  ## browser()
  lda.out <- lda(X.in.leaves[(core.in.leaves == 1), ],
                 y.in.leaves[core.in.leaves == 1])
  if (n.leaves == 2) {
    par(mfrow = c(3, 1))
    node.xy <- gsl.draw.cluster.tree(gsl.cluster.out, ...)
    title("Cluster tree")
    points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
           col = "blue1")
    X.proj <- X.in.leaves %*% lda.out$scaling
    hist(X.proj, xlim = range(X.proj), nclass = 20,  xaxt = "n", xlab = "",
         col = "blue1", main = "Observations in node")
    X.proj.core <- X.proj[core.in.leaves == 1]
    hist(X.proj.core, xlim = range(X.proj), nclass = 20,  xaxt = "n",
         xlab = "", col = "blue1", main = "Observations in core")
    par(mfrow = c(1, 1))
    return()
  }
  par(mfrow = c(2, 1))
  node.xy <- gsl.draw.cluster.tree(gsl.cluster.out, ...)
  title("Cluster tree")
  points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
         col = "blue1")
  X.proj <- X.in.leaves %*% lda.out$scaling[,1:2]
  plot(X.proj, col = "grey", pch = 20, xaxt = "n", yaxt = "n")
  points(X.proj[(core.in.leaves == 1),], col = "black", pch = 20)
  readline("\nHit return to see colored plot ")
  
  node.xy <- gsl.draw.cluster.tree(gsl.cluster.out, ...)
  title("Cluster tree")
  points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
         col = "black")
  points(node.xy$x[leaf.inodes], node.xy$y[leaf.inodes], pch = 19,
         col = rainbow(n.leaves))
  col <- rainbow(n.leaves)[y.in.leaves]     
  plot(X.proj, pch = ".", col = col, xaxt = "n", yaxt = "n")
  points(X.proj[(core.in.leaves == 1),], col = col[core.in.leaves == 1],
         pch = 20)
  return()
}
gsl.draw.tree <- function(cluster.tree, draw.cluster.numbers = F,
                          draw.cluster.leaf.codes = F,
                          draw.runt.excess.mass = F,
                          draw.runt.size = F) {
  spacing <- 1
  nnodes <- nrow(cluster.tree)
  left.width <- rep(0, nnodes)
  right.width <- rep(0, nnodes)
  node.x <- rep(0, nnodes)
  node.y <- cluster.tree[, "level"]
  runt.excess.mass <- cluster.tree[, "runt.excess.mass"]
  runt.size <- cluster.tree[, "runt.size"]
  n <- cluster.tree[1, "size"]
  plot.symbol <- rep(0, nnodes)
  current <- 1
  for (i in 1:nnodes) {
    if (cluster.tree[i, 1] == 0) {
      if (draw.cluster.numbers) {
        plot.symbol[i] <- toString(current)
        current <- current + 1
      }
      if (draw.cluster.leaf.codes) 
        plot.symbol[i] <- toString(cluster.tree[i, "leaf.code"])
    }
  }
  is.leaf <- cluster.tree[, 1] == 0
  compute.width <- function(i, spacing) {
    left.width[i] <<- 0
    right.width[i] <<- 0
    if(cluster.tree[i, 1] > 0) {
      left.width[i] <<- compute.width(cluster.tree[i, 1], spacing)
      right.width[i] <<- compute.width(cluster.tree[i, 2], spacing)
    }
    return(left.width[i] + right.width[i] + spacing)
  }
  compute.node.xy <- function(i, lr, mother.x, mother.y, spacing) {
    if(lr == "left") {
      node.x[i] <<- mother.x - (right.width[i] + 0.5 * spacing)
    }
    else {
      node.x[i] <<- mother.x + (left.width[i] + 0.5 * spacing)
    }
    if (cluster.tree[i, 1] > 0) {
      compute.node.xy(cluster.tree[i, 1], "left", node.x[i], node.y[i], spacing)
      compute.node.xy(cluster.tree[i, 2], "right", node.x[i], node.y[i], spacing)
    }
  }
  draw.edges <- function(i) {
    if(cluster.tree[i, 1] > 0) {
      leftson <- cluster.tree[i, 1]
      rightson <- cluster.tree[i, 2]
      lines(c(node.x[i], node.x[leftson]), c(node.y[i], node.y[leftson]))
      lines(c(node.x[i], node.x[rightson]), c(node.y[i], node.y[rightson]))
      draw.edges(leftson)
      draw.edges(rightson)
    }
  }
  compute.width(1, spacing)
  compute.node.xy(1, "left", 0, 0, spacing)
  ylim = c(-0.05 * max(node.y), 1.05 * max(node.y))
  if (draw.cluster.numbers | draw.cluster.leaf.codes) {
    plot(node.x, node.y, pch = 19, xaxt = "n", xlab = "", ylab = "level",
         ylim = ylim, cex.axis = 1.5, cex.lab = 1.5)
    text(node.x[is.leaf], node.y[is.leaf] + 0.05 * max(node.y),
         labels = plot.symbol[is.leaf], cex = 1.5)
    points(node.x[!is.leaf], node.y[!is.leaf], pch = 19)
  }
  else {
    plot(node.x, node.y, pch = 19, xaxt = "n", xlab = "", ylab = "level", ylim = ylim,
         cex.axis = 1.5, cex.lab = 1.5)
  }
  if(draw.runt.excess.mass) {
    ##    text(node.x[!is.leaf] + 0.2, node.y[!is.leaf],
    ##       round(n * runt.excess.mass[!is.leaf], 0))
    text(node.x[!is.leaf], node.y[!is.leaf] - 0.05 * max(node.y),
         round(n * runt.excess.mass[!is.leaf], 0), cex = 1.5)
  }
  if(draw.runt.size) {
    text(node.x[!is.leaf], node.y[!is.leaf] - 0.05 * max(node.y),
         runt.size[!is.leaf], cex = 1.5)
  }
  
  draw.edges(1)
  return(list(x = node.x, y = node.y))
}
sphere <- function(X) {
  Centered <- sweep(X, 2, apply(X, 2, mean))
  Sigma  <- var(Centered)
  R <- chol(Sigma)
  Sphered <- t(solve(t(R), t(Centered)))
  return(Sphered)
}
standardize <- function(X) {
  centered <- sweep(X, 2, apply(X, 2, mean))
  sds <- sqrt(apply(centered, 2, var))
  scaled <- sweep(centered, 2, sds, FUN = "/")
  return(scaled)
}
make.gaussian.kernel.density.estimate <- function(X.train, h) {
  density.estimate <- function(X.eval) {
    phat <- gaussian.kernel.density.estimate(X.train, X.eval, h, cv = F)
    return(phat)
  }
  return(density.estimate)
}
gaussian.kernel.density.estimate <- function(X.obs, X.eval, h, cv = F) {
  ## K <- function(tsquared, h) {
  ##     return(exp(-tsquared / (2 * h^2)) / h^p)
  ##   }
  K <- function(tsquared, h) {
    return(exp(-tsquared / (2 * h^2)) / (((2 * pi)^(p/2)) * h^p))
  }
  if (is.vector(X.obs)) {
    X.obs <- matrix(X.obs, ncol = 1)
    X.eval <- matrix(X.eval, ncol = 1)
  }
  n.eval <- nrow(X.eval)
  n.obs <- nrow(X.obs)
  p <- ncol(X.obs)
  dens <- rep(0, n.eval)
  squared.dist <- gsl.interpoint.distance.matrix(X.eval, X.obs) 
  for (i in 1:n.eval) {
    kernel.values <- K(squared.dist[i,], h)
    dens[i] <- sum(kernel.values) / n.obs
    if (cv) {
      dens[i] <- (sum(kernel.values) - kernel.values[i]) / (n.obs - 1)
    }
  }
  return(dens)
}
gaussian.kernel.least.squares.cv <- function(X, h) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  n <- nrow(X)
  p <- ncol(X)
  nh <- length(h)
  cv <- rep(0, nh)
  K <- function(tsquared, h) {
    return(exp(-tsquared / (2 * h^2)) / h^p)
  }
  Kstar <- function(tsquared, h) {
    return(K(tsquared, h*sqrt(2)))
  }
  D <- gsl.interpoint.distance.matrix(X, X)
  for (i in 1:nh) {
    cv[i] <- sum(Kstar(D, h[i]))/n^2 - 2 * sum(K(D, h[i]))/(n * (n-1)) +
      2 * K(0, h[i])/(n-1)
  }
  return(cv)
}
golden.section <- function(fun, xleft = 0.05, xright = 1, ngrid = 10, maxit = 8){
  xgrid <- seq(xleft, xright, length = ngrid)
  fgrid <- rep(0, ngrid)
  for (i in 1:ngrid) fgrid[i] <- fun(xgrid[i])
  iopt <- which.min(fgrid)
  if ((iopt == 1)|(iopt == ngrid)) return(-1)
  xleft <- xgrid[iopt-1]; fleft = fgrid[iopt-1]
  xmiddle <- xgrid[iopt]; fmiddle = fgrid[iopt]
  xright <- xgrid[iopt+1]; fright = fgrid[iopt+1]
  ## We have initial bracket
  for (i in 1:maxit) {
    if ((xmiddle - xleft) > (xright - xmiddle)) {
      xnew <- xmiddle - 0.38 * (xmiddle - xleft)
      fnew <- fun(xnew)
      if (fnew < fmiddle) {
        xright <- xmiddle; fright <- fmiddle; xmiddle <- xnew; fmiddle <- fnew
      }
      else {
        xleft <- xnew; fleft = fnew
      }
    }
    else {
      xnew <- xright - 0.38 * (xright - xmiddle)
      fnew <- fun(xnew)
      if (fnew < fmiddle) {
        xleft <- xmiddle; fleft <- fmiddle; xmiddle <- xnew; fmiddle <- fnew
      }
      else{
        xright <- xnew; fright <- fnew
      }
    }
  }
  return(xmiddle)
}
cv.search <- function(X, trial.par = c(0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4,
                                       0.6, 0.8, 1.0, 1.5, 2.0, 4, 8)) {
  cv.function <- gaussian.kernel.least.squares.cv
  search.for.h <- T
  cv <- cv.function(X, trial.par)
  finite <- is.finite(cv)
  opt.smopar <- NA
  finite.par <- NA
  finite.cv <- NA
  if (sum(finite) > 3) {
    finite.cv <- cv[finite]
    finite.par <- trial.par[finite]
    imin <- which.min(finite.cv)
    if((imin != 1) & (imin != length(finite.cv))) {
      opt.smopar <- finite.par[imin]
      if (search.for.h) {
        fun <- function(h) {
          return(cv.function(X, h))
        }
        opt.smopar <- golden.section(fun, xleft = finite.par[imin - 1],
                                     xright = finite.par[imin + 1], maxit = 8)
      }
    }
  }
  return(list(opt.smopar = opt.smopar, smopars = finite.par, cv = finite.cv))
}
rand.index <- function(true.labels, estimated.labels) {
  tab <- table(true.labels, estimated.labels)
  nidot <- apply(tab, 1, sum)
  ndotj <- apply(tab, 2, sum)
  n <- sum(tab)
  num <- sum(choose(tab, 2)) - 
    sum(choose(nidot, 2)) * sum(choose(ndotj, 2)) / choose(n, 2)
  denom <- 0.5 * (sum(choose(nidot, 2)) + sum(choose(ndotj, 2))) -
    sum(choose(nidot, 2)) * sum(choose(ndotj, 2)) / choose(n, 2)
  return(num/denom)
}

##### michael's original functions #####
get.empirical.conf.band <- function(hs,target.cov.prob){
  data <- rbind(hs$resample.edge.weights,hs$resample.vertex.weights)
  nc <- dim(data)[2]
  nr <- dim(data)[1]
  prop <- 0
  out <- nc*target.cov.prob/2
  while(prop < target.cov.prob){
    out <- out+1
    lower <- floor(nc/2+1-out)
    upper <- ceiling(nc/2+out)
    lows <- apply(data,1,FUN=function(x){sort(x)[lower]})
    upps <- apply(data,1,FUN=function(x){sort(x)[upper]})
    prop <- mean(apply(data,2,FUN=function(x){all(x>=lows & x<=upps)}))
  }
  lower <- max(2,lower)
  upper <- min(nc-1,upper)
  lows <- apply(data,1,FUN=function(x){sort(x)[lower]})
  upps <- apply(data,1,FUN=function(x){sort(x)[upper]})
  prop <- mean(apply(data,2,FUN=function(x){all(x>=lows & x<=upps)}))
  #print(prop)
  #print(c(lower,upper))
  lows <- apply(data,1,FUN=function(x){sort(x)[lower]})
  upps <- apply(data,1,FUN=function(x){sort(x)[upper]})
  list(upper.vw = upps[(nr+1-nrow(hs$resample.vertex.weights)):nr],
       lower.vw = lows[(nr+1-nrow(hs$resample.vertex.weights)):nr],
       upper.ew = upps[1:nrow(hs$resample.edge.weights)],
       lower.ew = lows[1:nrow(hs$resample.edge.weights)],
       cov.parameter = c(lower,upper),
       prop = prop)
}
get.bc.conf.band <- function(hs,target.cov.prob){
  data <- rbind(hs$resample.edge.weights,hs$resample.vertex.weights)
  
  #find necessary cov.prob for desired simultaneous coverage
  get.simul.cov <- function(cov.prob){
    alpha <- (1-cov.prob)/2
    vals <- c()
    alpha1<-pnorm(2*apply(data,1,function(x){qnorm(mean(x<mean(x)))}) + qnorm(alpha))
    alpha2<-pnorm(2*apply(data,1,function(x){qnorm(mean(x<mean(x)))}) + qnorm(1-alpha))
    vals <- unlist(lapply(1:nrow(data),function(x){ci <- quantile(data[x,],c(alpha1[x],alpha2[x]))
    which(data[x,]>=ci[1] & data[x,]<=ci[2])}))
    prop <- length(which(table(vals) == nrow(data)))/nrow(data)
    return(prop-target.cov.prob)
  }
  cov.prob <- uniroot(f=get.simul.cov,interval=c(0,1))$root
  
  #get intervals
  alpha <- (1-cov.prob)/2
  alpha1<-pnorm(2*apply(data,1,function(x){qnorm(mean(x<mean(x)))}) + qnorm(alpha))
  alpha2<-pnorm(2*apply(data,1,function(x){qnorm(mean(x<mean(x)))}) + qnorm(1-alpha))
  CIs <- unlist(lapply(1:nrow(data),function(x){as.vector(quantile(data[x,],c(alpha1[x],alpha2[x])))}))
  intervals <- matrix(CIs,nrow=nrow(data),ncol=2,byrow = TRUE)
  
  list(upper.vw = intervals[(nrow(hs$resample.edge.weights)+1):nrow(data),2],
       lower.vw = intervals[(nrow(hs$resample.edge.weights)+1):nrow(data),1],
       upper.ew = intervals[1:nrow(hs$resample.edge.weights),2],
       lower.ew = intervals[1:nrow(hs$resample.edge.weights),1],
       prop = cov.prob)
}
generate_simplex_data <- function(n,M,N,d,distribution,scale){
  #### Code to Generate Data in a Simplex
  # n = sample size, i.e. number of points to be generated
  # M = number of signal dimensions = number of modes
  # N = number of noise dimensions, meaning that M+N = total number of dimensions
  # d = distance between modes (equal distances assumed)
  # distribution = {"uniform", "normal"}
  # scale = scale parameter for each distribution. For normal distributions, it is the standard deviation.
  # For uniform distributions, it is the length of the support of the distribution. If the vector
  # has length 1, the the same value is used for all distributions. Else, a different value is used for each.
  library(MASS)
  
  if(length(scale) == 1){scale <- rep(scale,N+M)}
  if(length(scale) != (N+M)){stop("Scale vector of incorrect length")}
  
  if(distribution == "normal"){
    if(M>0){
      vals <- lapply(1:M,function(i){
        lower <- round(n*(i-1)/M)+1
        upper <- round(n*i/M)
        mu <- rep(0,M)
        mu[i]<-d/sqrt(2)
        mvrnorm(n=upper-lower+1,mu=mu,Sigma = diag(scale[i]^2,nrow=M))
      })
      signal_data<-matrix(NA,nrow=0,ncol=M)
      for(i in 1:M){signal_data<- rbind(signal_data,vals[[i]])}
    }else{signal_data <- c()}
    if(N>0){noise_data  <- mvrnorm(n=n,mu=rep(0,N),Sigma = diag(scale[(M+1):(M+N)]^2))
    }else{noise_data <- c()}
    data <- c(unlist(signal_data),unlist(noise_data))
  }
  if(distribution == "uniform"){
    if(M>0){
      vals <- lapply(1:M,function(i){
        lower <- round(n*(i-1)/M)+1
        upper <- round(n*i/M)
        data <- matrix(NA,nrow=upper-lower+1,ncol=M)
        data[,i] <- runif(n=upper-lower+1,min = d/sqrt(2)-0.5*scale[i],max = d/sqrt(2) + 0.5*scale[i])
        data[,-i] <- runif(n=length(data[,-i]),min = -0.5*scale[i],max =0.5*scale[i])
        data
      })
      signal_data<-matrix(NA,nrow=0,ncol=M)
      for(i in 1:M){signal_data<- rbind(signal_data,vals[[i]])}
    }else{signal_data <- c()}
    if(N>0){
      noise_data  <- matrix(NA,nrow=n,ncol=N)
      for(i in 1:N){noise_data[,i]<-runif(n=n,min=-0.5*scale[M+i],max=0.5*scale[M+i])}
    }else{noise_data <- c()}
    data <- c(unlist(signal_data),unlist(noise_data))
  }
  
  X <- matrix(data,nrow = n, ncol = M + N)
  
  return(X)
}

cluster_results <- function(hs,X,group.id,confidence.bounds,pruning.criterion){
  uc.mdsfun <- cwc.make.upper.confidence.mdsfun(hs, confidence.bounds)
  obs.lower <- confidence.bounds$lower.vw
  gsl.cluster.out <- gsl.cluster(X, uc.mdsfun, pruning.criterion = pruning.criterion,
                                 pruning.threshold = 0, gsl.cluster.out = NULL,
                                 assign.fluff = T, obs.density.lower.confidence.bounds = obs.lower)
  numClusters <- length(unique(gsl.cluster.out$leaf.code))
  list(table = table(group.id, gsl.cluster.out$leaf.code),
       numClusters = numClusters,
       lambda = unique(gsl.cluster.out$cluster.tree[,"level"]))
}
get_results <- function(hs,X,confidence.bounds,pruning.criterion){
  uc.mdsfun <- cwc.make.upper.confidence.mdsfun(hs, confidence.bounds)
  obs.lower <- confidence.bounds$lower.vw
  gsl.cluster.out <- gsl.cluster(X, uc.mdsfun, pruning.criterion = pruning.criterion,
                                 pruning.threshold = 0, gsl.cluster.out = NULL,
                                 assign.fluff = T, obs.density.lower.confidence.bounds = obs.lower)
  length(unique(gsl.cluster.out$leaf.code))
}
full_cluster_results <- function(data,nres,target.cov.prob,
                                 kmax,ngrid,pruning.criterion){
  #### perform half-sampling
  X<-as.matrix(unlist(data[1]),ncol=1)
  group.id<-unlist(data[2])
  hs <- cwc.resample.kernel.ve.weights(X,kmax=kmax, # number of nearest neighbors in spanning tree
                                       nres=nres,  #number of resamples
                                       ngrid = ngrid, #number of tested values in grid search?
                                       bandwidth = "cv", #get cv bandwidth for each resample
                                       resample.type = "half.sampling")
  vw.hs <- hs$resample.vertex.weights
  ew.hs <- hs$resample.edge.weights
  
  #### confidence bounds
  empirical.confidence.bounds <- get.empirical.conf.band(hs,target.cov.prob)
  gaussian.confidence.bounds <-cwc.confidence.bounds(hs,target.cov.prob,simultaneous = T,band.type = "gaussian")
  poisson.confidence.bounds <-cwc.confidence.bounds(hs,target.cov.prob,simultaneous = T,band.type = "poisson")
  bc.confidence.bounds <- get.bc.conf.band(hs,target.cov.prob)
  
  #### results
  e_res <- cluster_results(hs,X,group.id,empirical.confidence.bounds,pruning.criterion)
  g_res <- cluster_results(hs,X,group.id,gaussian.confidence.bounds,pruning.criterion)
  p_res <- cluster_results(hs,X,group.id,poisson.confidence.bounds,pruning.criterion)
  b_res <- cluster_results(hs,X,group.id,bc.confidence.bounds,pruning.criterion)
  
  return(list(e_table = e_res[1],e_nc = e_res[2],
              g_table = g_res[1],g_nc = g_res[2],
              p_table = p_res[1],p_nc = p_res[2],
              b_table = b_res[1],b_nc = b_res[2]))
}



##### updated werner functions #####
cwc.resample.kernel.ve.weights <- function(X, tg.vertices = NULL, tg.edges = NULL,kmax = 20,
                                           include.mst = T,nres = 0,resamples = NULL, 
                                           trial.h = NULL,ngrid = 10, bandwidth = "cv",
                                           resample.type = "half.sampling",debug.filename = NULL){
  debug <- !is.null(debug.filename)
  if (is.null(tg.vertices)){
    tg.edges <- cwc.determine.edges(X, kmax = kmax, include.mst = include.mst)
    tg.vertices <- X
  }
  if (is.null(trial.h))
    trial.h <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80,
                 1.00, 1.50, 2.00)
  nvertices <- nrow(tg.vertices)
  nedges <- nrow(tg.edges)
  n <- nrow(X)
  ##
  ## Compute vertex and edge weights for original sample
  ##
  if (debug) cat("\n\n\nComputing edge weights for original sample",
                 file = debug.filename, append = T)
  cvs.orig <- NULL
  if (is.numeric(bandwidth)) h.orig <- bandwidth
  if (!is.numeric(bandwidth)) {
    orig.cvs <- cv.search(X, trial.h)
    orig.h <- orig.cvs$opt.smopar
  }
  density <- make.gaussian.kernel.density.estimate(X, orig.h)
  cwcvw.out <- cwc.ve.weights(tg.vertices, tg.edges, density, ngrid)
  vertex.weights <- cwcvw.out$vertex.weights
  edge.weights <- cwcvw.out$edge.weights
  
  resample.vertex.weights <- NULL
  resample.edge.weights <- NULL
  resample.h <- NULL
  resample.cvs <- NULL
  
  if ((nres > 0) | (!is.null(resamples))) {
    if (is.null(resamples)) {
      if (resample.type == "half.sampling") rs.size <- round(nrow(X) / 2)
      if (resample.type == "bootstrap") rs.size <- nrow(X)
      resamples <- matrix(0, nrow = rs.size, ncol = nres)
      for (i in 1:nres) {
        if (resample.type == "half.sampling") rs <- sample(1:n, rs.size, replace = F)
        if (resample.type == "bootstrap") rs <- sample(1:n, rs.size, replace = T)
        resamples[,i] <- rs
      }
    }
    nres <- ncol(resamples)
    resample.type <- "half.sampling"
    if (nrow(resamples) == n) resample.type <- "bootstrap"
    resample.vertex.weights <- matrix(0, nrow = nvertices, ncol = nres)
    resample.edge.weights <- matrix(0, nrow = nedges, ncol = nres)
    resample.cvs <- vector("list", nres)
    resample.h <- rep(0, nres)
    
    if (debug) cat("\n\nFinding bandwidths for resamples",
                   file = debug.filename, append = T)
    if ((resample.type == "bootstrap") & !is.numeric(bandwidth)) {
      cvso <- cv.search(X, trial.h)
      resample.cvs[[1]] <- cvso
      resample.h[1:nres] <- cvso$opt.smopar
      if (debug) cat("\n", resample.h[1], file = debug.filename, append = T)
    } 
    if ((resample.type == "half.sampling") & !is.numeric(bandwidth)) {
      for (i in 1:nres) {
        rs <- resamples[,i]
        X.rs <- as.matrix(X[rs, ], nrow = rs.size)
        cvso <- cv.search(X.rs, trial.h)
        resample.cvs[[i]] <- cvso
        resample.h[i] <- cvso$opt.smopar
        if (debug) cat("\n", i, resample.h[i], file = debug.filename, append = T)
      }
      mean.h <- mean(resample.h[!is.na(resample.h)])
    }
    if (debug) cat("\n\nComputing resample vertex and edge weights\n",
                   file = debug.filename, append = T)
    for (i in 1:nres) {
      if (is.numeric(bandwidth)) h <- bandwidth
      if (!is.numeric(bandwidth)) {
        if (resample.type == "bootstrap") h <- resample.h[1]
        if (resample.type == "half.sampling") {
          h <- mean.h
          if ((bandwidth == "cv") & (!is.na(resample.h[i]))) h <- resample.h[i]
        }
      }
      rs <- resamples[,i]
      X.rs <- as.matrix(X[rs, ], nrow = rs.size)
      density <- make.gaussian.kernel.density.estimate(X.rs, h)
      cwcvw.out <- cwc.ve.weights(tg.vertices, tg.edges, density, ngrid)
      ##  browser()
      resample.vertex.weights[,i] <- cwcvw.out$vertex.weights
      resample.edge.weights[,i] <- cwcvw.out$edge.weights
      if(debug) cat(" ", i, file = debug.filename, append = T)
    }
  }
  return(list(tg.vertices = tg.vertices,
              tg.edges = tg.edges,
              vertex.weights = vertex.weights,
              edge.weights = edge.weights,
              orig.cvs = orig.cvs,
              resample.vertex.weights = resample.vertex.weights,
              resample.edge.weights = resample.edge.weights,
              resamples = resamples,
              resample.h = resample.h,
              resample.cvs = resample.cvs))
}
cwc.determine.edges <- function(X, kmax = 20, include.mst = T,debug = F) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  n <- nrow(X)
  p <- ncol(X)
  already.computed <- matrix(F, nrow = n, ncol = n)
  edges <- matrix(0, nrow = n * kmax + (n-1), ncol = 2)
  nedges <- 0
  Dist <- gsl.interpoint.distance.matrix(X, X)
  if (include.mst) {
    mst.mat.out <- gsl.mst.mat(Dist)
    mst.edges <- mst.mat.out$edges
    for (i in 1:(n-1)) {
      nedges <- nedges + 1
      edges[nedges,] <- mst.edges[i,]
      v1 <- mst.edges[i, 1]
      v2 <- mst.edges[i, 2]
      already.computed[v1, v2] <- T
      already.computed[v2, v1] <- T
    }
  }
  if(kmax >= 1) {
    for (i in 1:n) {
      v1 <- i
      NN <- order(Dist[i,])
      for (j in 2:(kmax + 1)) {
        v2 <- NN[j]
        if (!already.computed[v1, v2]) {
          nedges <- nedges + 1
          edges[nedges, 1] <- v1
          edges[nedges, 2] <- v2
          already.computed[v1, v2] <- T
          already.computed[v2, v1] <- T
        }
      }
    }
  }
  return(edges[1:nedges, ])
}
cwc.ve.weights <- function(tg.vertices, tg.edges, density, ngrid) {
  nedges <- nrow(tg.edges)
  m <- ncol(tg.vertices)
  edge.weights <- rep(0, nedges)
  vertex.weights <- density(tg.vertices)
  t <- seq(0, 1, length = ngrid)[2:(ngrid-1)]
  Tmat <- matrix(t, nrow = ngrid-2, ncol = m)
  for (i in 1:nedges) {
    v1.index <- tg.edges[i, 1]
    v2.index <- tg.edges[i, 2]
    v1 <- tg.vertices[v1.index,]
    v2 <- tg.vertices[v2.index,]
    V1 <- matrix(v1, nrow = ngrid-2, ncol = m, byrow = T)
    V2 <- matrix(v2, nrow = ngrid-2, ncol = m, byrow = T)
    E <- (1-Tmat) * V1 + Tmat * V2
    phat <- c(vertex.weights[v1.index], density(E), vertex.weights[v2.index])
    edge.weights[i] <- min(phat)
  }
  return(list(vertex.weights = vertex.weights, edge.weights = edge.weights))
}
cwc.confidence.bounds <- function(resample.weights, target.cov.prob,band.type = "empirical",
                                  simultaneous = F){
  if(simultaneous)
    confidence.bounds <- cwc.simul.confidence.bounds(resample.weights, target.cov.prob,
                                                     band.type = band.type)
  if(!simultaneous)
    confidence.bounds <- cwc.nonsimul.confidence.bounds(resample.weights, target.cov.prob,
                                                        band.type = band.type)
  return(confidence.bounds)
}
cwc.simul.confidence.bounds <- function(resample.weights, target.simul.cov.prob,
                                        band.type = "empirical") {
  make.fun <- function(rvw, rew, band.type, target.sim.cov.prob) {
    fun <- function(cov.parameter) {
      return(cwc.simul.cov.prob(resample.weights, cov.parameter,
                                band.type)$simul.cov.prob
             - target.sim.cov.prob)
    }
    return(fun)
  }
  rvw <- resample.weights$resample.vertex.weights
  rew <- resample.weights$resample.edge.weights
  ne <- nrow(rew)
  nv <- nrow(rvw)
  nres <- ncol(rvw)
  
  if ((band.type == "gaussian") | (band.type == "poisson")) {
    xmin <- 0
    ## xmax <- 4
    xmax <- 10  ## changed 11-02-2018
    fun <- make.fun(rvw, rew,  band.type, target.simul.cov.prob)
    gamma <- binary.root.search(fun, xmin, xmax, maxit = 30) #maxit set to 30 by michael
    #print(gamma) #added by michael
    ew.mean <- apply(rew, 1, mean)
    ew.sd <- apply(rew, 1, sd)
    vw.mean <- apply(rvw, 1, mean)
    vw.sd <- apply(rvw, 1, sd)
  }
  if (band.type == "gaussian") {
    upper.ew <- ew.mean + gamma * ew.sd
    upper.vw <- vw.mean + gamma * vw.sd
    lower.vw <- vw.mean - gamma * vw.sd
  }
  if (band.type == "poisson") {
    vw.sd <- sqrt(vw.mean)
    upper.ew <- ew.mean + gamma * sqrt(ew.mean)
    upper.vw <- vw.mean + gamma * sqrt(vw.mean)
    lower.vw <- vw.mean - gamma * sqrt(vw.mean)
  }
  # if (band.type == "empirical") {
  #   xmin <- 0
  #   xmax <- nres / 2
  #   fun <- make.fun(rvw, rew, band.type, target.simul.cov.prob)
  #   nout <- round(binary.root.search(fun, xmin, xmax, maxit = 8))
  #   print(cwc.simul.cov.prob(resample.weights,nout,band.type)$simul.cov.prob) #added by michael
  #   gamma <- nout
  #   lower <- nout + 1
  #   upper <- nres - nout
  #   if (lower >= upper) {
  #     lower <- lower - 1
  #     upper <- upper + 1
  #   }
  #   upper.ew <- rep(0, ne)
  #   upper.vw <- rep(0, nv)
  #   lower.vw <- rep(0, nv)
  #   for (i in 1:ne) {
  #     sew <- sort(rew[i,])
  #     upper.ew[i] <- sew[upper]
  #   }
  #   for (i in 1:nv) {
  #     svw <- sort(rvw[i,])
  #     upper.vw[i] <- svw[upper]
  #     lower.vw[i] <- svw[lower]
  #   }
  #  print(c(lower,upper)) #added by michael
  #}
  return(list(upper.ew = upper.ew, upper.vw = upper.vw, lower.vw =
                lower.vw, cov.parameter = gamma))
}
binary.root.search <- function(fun, xmin, xmax, maxit = 8) {
  nit = 0
  f.xmin <- fun(xmin)
  f.xmax <- fun(xmax)
  if (sign(f.xmin * f.xmax) > 0)
    stop("binary.root.search: initial bounds not a bracket")
  if (sign(f.xmin * f.xmax) == 0) {
    if (f.xmin == 0) return(xmin) else return(xmax)
  }
  for (i in 1:maxit) { 
    xmiddle <- (xmin + xmax) / 2
    f.xmiddle <- fun(xmiddle)
    if (f.xmiddle == 0) return(xmiddle)
    if (sign(f.xmiddle * f.xmin) < 0) {
      xmax <- xmiddle
      f.xmax <- f.xmiddle
      next
    }
    if (sign(f.xmiddle * f.xmax) < 0) {
      xmin <- xmiddle
      f.xmin <- f.xmiddle
      next
    }
  }
  return(xmiddle)
}
cwc.simul.cov.prob <- function(resample.weights, cov.parameter,
                               band.type = "empirical") {
  rvw <- resample.weights$resample.vertex.weights
  rew <- resample.weights$resample.edge.weights
  nres <- ncol(rew)
  ne <- nrow(rew)
  nv <- nrow(rvw)
  
  if (band.type == "empirical") {
    nout <- cov.parameter
    lower <- nout + 1
    upper <- nres - nout
    if (lower >= upper) {
      lower <- lower - 1
      upper <- upper + 1
    }
    upper.ew <- rep(0, ne)
    upper.vw <- rep(0, nv)
    lower.vw <- rep(0, nv)
    for (i in 1:ne) {
      sew <- sort(rew[i,])
      upper.ew[i] <- sew[upper]
    }
    for (i in 1:nv) {
      svw <- sort(rvw[i,])
      upper.vw[i] <- svw[upper]
      lower.vw[i] <- svw[lower]
    }
  }
  if (band.type != "empirical") {
    gamma <- cov.parameter
    ew.mean <- apply(rew, 1, mean)
    ew.sd <- apply(rew, 1, sd)
    vw.mean <- apply(rvw, 1, mean)
    vw.sd <- apply(rvw, 1, sd)
  }
  if (band.type == "gaussian") {
    upper.ew <- ew.mean + gamma * ew.sd
    upper.vw <- vw.mean + gamma * vw.sd
    lower.vw <- vw.mean - gamma * vw.sd
  }
  if (band.type == "poisson") {
    vw.sd <- sqrt(vw.mean)
    upper.ew <- ew.mean + gamma * sqrt(ew.mean)
    upper.vw <- vw.mean + gamma * sqrt(vw.mean)
    lower.vw <- vw.mean - gamma * sqrt(vw.mean)
  }
  nout <- 0
  out.resamples <- rep(0, nres)
  for (i in 1:nres) {
    out <- sum(rew[,i] > upper.ew) + sum(rvw[,i] > upper.vw) +
      sum(rvw[,i] < lower.vw)
    if (out > 0) {
      nout <- nout + 1
      out.resamples[nout] <- i
    }
  }
  return(list(simul.cov.prob = 1 - nout/nres, out.resamples =
                out.resamples[1:nout]))
  
}
cwc.make.upper.confidence.mdsfun <- function(resample.weights, confidence.bounds) {
  tg.vertices <- resample.weights$tg.vertices
  tg.edges <- resample.weights$tg.edges
  n <- nrow(tg.vertices)
  nedges <- nrow(tg.edges)
  uc.mdsmat <- matrix(0, nrow = n, ncol = n)
  diag(uc.mdsmat) <- confidence.bounds$upper.vw
  for (i in 1:nedges) {
    v1 <- tg.edges[i, 1]
    v2 <- tg.edges[i, 2]
    uc.mdsmat[v1, v2] <- confidence.bounds$upper.ew[i]
    uc.mdsmat[v2, v1] <- confidence.bounds$upper.ew[i]
  }
  up <- as.vector(uc.mdsmat[upper.tri(uc.mdsmat, diag = T)])
  upper.confidence.mdsfun <- gsl.make.mdsfun.from.mdsmat(up)
  return(upper.confidence.mdsfun)
}
gsl.pruning.criteria <- function(mast.dendogram, obs.density, ## modified to remove "browser()"
                                 obs.density.lower.confidence.bounds = NULL) {
  merge <- mast.dendogram$merge
  height <- mast.dendogram$height
  left.size <- rep(0, nrow(merge))
  right.size <- rep(0, nrow(merge))
  left.excess.mass <- rep(0, nrow(merge))
  right.excess.mass <- rep(0, nrow(merge))
  left.rise <- rep(0, nrow(merge))
  right.rise <- rep(0, nrow(merge))
  left.bootstrap.rise <- rep(0, nrow(merge))
  right.bootstrap.rise <- rep(0, nrow(merge))
  pc.int <- function(inode) {
    ileft <- merge[inode, 1]
    iright <- merge[inode, 2]
    if (ileft < 0) left.leaves <- -ileft
    if (ileft > 0) left.leaves <- pc.int(ileft)
    if (iright < 0) right.leaves <- -iright
    if (iright > 0) right.leaves <- pc.int(iright)
    left.high.density.cluster <- left.leaves[obs.density[left.leaves] >=
                                               height[inode]]
    right.high.density.cluster <- right.leaves[obs.density[right.leaves] >=
                                                 height[inode]]
    ## if (!is.null(left.high.density.cluster)) {
    #if (length(left.high.density.cluster) == 0) browser()
    if (length(left.high.density.cluster) > 0) {
      left.size[inode] <<- length(left.high.density.cluster)
      left.excess.mass[inode] <<- sum(1 - height[inode] /
                                        obs.density[left.high.density.cluster])/
        (nrow(merge) + 1)
      left.max <- max(obs.density[left.high.density.cluster])
      if (left.max == Inf) left.rise[inode] <<- Inf
      else left.rise[inode] <<- left.max - height[inode]
      ##left.rise[inode] <<- max(obs.density[left.high.density.cluster]) -
      ##  height[inode]
      ## if (left.rise[inode] == -Inf) browser()
      ##if (left.rise[inode] == Inf) browser()
      if (!is.null(obs.density.lower.confidence.bounds)) 
        left.bootstrap.rise[inode] <<-
        max(obs.density.lower.confidence.bounds[left.high.density.cluster]) -
        height[inode]
    }
    ## if (!is.null(right.high.density.cluster)) {    ## Changed 6-6-2018
    #if (length(right.high.density.cluster) == 0) browser()
    if (length(right.high.density.cluster) > 0) {
      right.size[inode] <<- length(right.high.density.cluster)
      right.excess.mass[inode] <<- sum(1 - height[inode] /
                                         obs.density[right.high.density.cluster])/
        (nrow(merge) + 1)
      right.max <- max(obs.density[right.high.density.cluster])
      #if (right.max == -Inf) browser()  ## 6-6-2018
      if (right.max == Inf) right.rise[inode] <<- Inf
      else right.rise[inode] <<- right.max - height[inode]
      ##right.rise[inode] <<- max(obs.density[right.high.density.cluster]) -
      ##  height[inode]
      ## if (right.rise[inode] == -Inf) browser()  
      ## if (right.rise[inode] == Inf) browser()
      if (!is.null(obs.density.lower.confidence.bounds))
        right.bootstrap.rise[inode] <<-
        max(obs.density.lower.confidence.bounds[right.high.density.cluster]) -
        height[inode]
    }
    return(c(left.leaves, right.leaves))
  }
  pc.int(nrow(merge))
  res <- cbind(left.size, right.size, left.excess.mass, right.excess.mass,
               left.rise, right.rise, left.bootstrap.rise, right.bootstrap.rise)
  colnames(res) <- c("left.size", "right.size", "left.excess.mass", "right.excess.mass",
                     "left.rise", "right.rise", "left.bootstrap.rise",
                     "right.bootstrap.rise")
  return(res)
}