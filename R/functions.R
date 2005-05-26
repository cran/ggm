"adjMatrix" <-
function (A) 
{
### From the edge matrix to the adjacency matrix
  E <- t(A)
  diag(E) <- 0
  E
}

"allEdges" <-
function(amat){
### Finds all the edges of a graph with edge matrix amat.
    nn <- 1:nrow(amat)
    E <- c()
    if(all(amat == t(amat))) { 
      amat[lower.tri(amat)] <- 0
    }
    for(i in nn) {
      e <- nn[amat[i,] == 1]
      if(length(e) == 0) next
      li <- cbind(i,  e)
      dimnames(li) <- list(rep("", length(e)), rep("", 2))
      E <- rbind(E, li) 
    }
    E
  }

"basiSet" <-
function(amat){
### Basis set of a DAG with adjacency matrix amat.
    amat <- topSort(amat)
    nod <- vertices(amat)
    dv <- length(nod)
    ind <- NULL
    ## NOTE. This is correct if the adj mat is upper triangular.
    for(r in 1:dv){
      for(s in r:dv) {
        if((amat[r,s] != 0) | (s==r))
          next
        else{
          ed <- nod[c(r,s)]
          pa.r <- nod[amat[,r] == 1]
          pa.s <- nod[amat[,s] == 1] 
          dsep <- union(pa.r, pa.s) 
          dsep <- setdiff(dsep, ed)
          b <- list(c(ed, dsep))
          ind <- c(ind, b)
        }
      }
    }
    ##      ind <- lapply(ind, function(x) nn[x])
    ind
  }

"bd" <-
function (nn, amat) 
{
### Boundary of the nodes nn for a graph with adjacency matrix amat.
  nod <- vertices(amat)
  if(is.null(nod)) stop("The edge matrix must have dimnames!")
  if(!all(is.element(nn, nod))) stop("Some of the nodes are not among the vertices.")
  b <- vector(length(nn), mode="list")
  diag(amat) <- 0  # As you do not want the node itself in the list
  k <- length(nn)
  for(i in 1:k) {
    b[[i]] <- c(  nod[amat[nn[i], ]==1 ],nod[amat[,nn[i]]==1 ] )
  }
  b <- unique(unlist(b))
  setdiff(b, nn)
}

"bfs" <-
function(amat, v=1) {
### Breadth-first search of a connected UG with adjacency matrix amat.
    n <- nrow(amat)
    indices <- 1:n
    if(n==1) return(NULL)
    visited <- rep(0, n)
    Q <- c()
    tree <- matrix(0, n,n)
    dimnames(tree) <- dimnames(amat)
    E <- c()
    visited[v] <- 1
    Q <- c(v, Q)
    while(!all(visited)){
      x <- Q[1]
      Q <- Q[-1]
     ## b <- bd(x, amat)
      b <- indices[amat[x,]==1] # Boundary of x. Assumes that amat is symmetric
      for(y in b){
        if(visited[y] == 0){
          visited[y] <- 1
          Q <- c(Q, y)
          tree[x,y]<- 1 ; tree[y,x] <- 1
          E <- rbind(E, c(x,y))
        }
      }
    }
    cross <- amat - tree
    V <- allEdges(cross)
    dimnames(E) <- list(rep("", nrow(E)), rep("", 2))
    list(tree = tree, branches = E,chords = V ) 
  }

"ch" <-
function (nn, amat) 
{
### List of the children of nodes nn for a given with adjacency matrix amat.
  nod <- vertices(amat)
  if(is.null(nod)) stop("The adjacency matrix must have dimnames!")
  if(!all(is.element(nn, nod))) stop("Some of the nodes are not among the vertices.")
  k <- length(nn)
  p <- vector(k, mode="list")
  A <- 0 + ((amat != t(amat)) & (amat == 1)) # Select the directed edges
  for(i in 1:k) {
    p[[i]] <- nod[A[nn[i],]==1 ]
  }
  setdiff(unique(unlist(p)), nn)
}

"checkIdent" <-
function(amat, latent) {
### Checks SW sufficient conditions for identifiability of a DAG
### with adjacency matrix edge amat and one latent variable.
   "allSubsets" <-
     function (n) 
       {
         ## Returns all subsets of n
         p <- length(n)
         H <- data.matrix(expand.grid(rep(list(1:2), p))) - 1
         H <- split(H==1, row(H))
         lapply(H, function(i) n[i])
       }
   
    nod <- vertices(amat)
    if(is.null(nod)) stop("The adjacency matrix must have dimnames.")
    gcov <- inducedCovGraph(amat, sel=nod, cond=NULL)
    L <- latent
    if(length(L) > 1)
      stop("I have an answer only for one latent variable.")
    O <- setdiff(nod, L)
    m <-  bd(L, gcov)
    ## Theorem 1
    if(length(m) > 2){ 
      G <- inducedCovGraph(amat, sel=O, cond=L)
      cond.i <- isGident(G[m,m,drop=FALSE])
    }
    else
      cond.i <- FALSE
    gcon <- inducedConGraph(amat, sel=nod, cond=NULL) 
    cc <- bd(L, gcon)
    if(length(cc) > 2) {
      cond.ii <- isGident(gcon[cc, cc, drop=FALSE])
    } 
    else  
      cond.ii <- FALSE
    ## Theorem 2 (revised)
    a <- union(pa(L, amat), ch(L, amat))
    if(length(a) > 2){
      Oa <- setdiff(O, a)   # O \ a
      S.a <- inducedCovGraph(amat, sel=union(Oa, L), cond=a)
      first <- S.a[Oa, L]
      first <- Oa[first==0]  # variables that satisfy condition (i) of thm. 2.
      if(length(first)==0){
        cond.iii <- FALSE
      }
      else {
        cond.iii <- FALSE
        H <- allSubsets(first) # in all subsets of these variables
        for(h in H){           # look for a G-identifiable cov. or con. graph
          isgid <- isGident(inducedCovGraph(amat, sel=a, cond=union(L, h)))
          if(isgid){
            cond.iii <- TRUE
            break
          }
          else{
            isgid <- isGident(inducedConGraph(amat, sel=a, cond=union(L, h)))
            if(isgid){
              cond.iii <- TRUE
              break
            }
          }
        }
      }
      second <- setdiff(O,m) # variables that satisfy condition (ii) of thm. 2.
      if(length(second)==0){
        cond.iv <- FALSE
      }
      else {
        cond.iv <- FALSE
        H <- allSubsets(second)  # in all subsets of these variables
        for(h in H){             # look for a G-identifiable cov. or con. graph
          isgid <- isGident(inducedCovGraph(amat, sel=a, cond=union(L, h)))
          if(isgid){
            cond.iv <- TRUE
            break
          }
          else{
            isgid <- isGident(inducedConGraph(amat, sel=a, cond=union(L, h)))
            if(isgid){
              cond.iv <- TRUE
              break
            }
          }
        }
      }  
    }
    else{
      cond.iii <- FALSE
      cond.iv <- FALSE
    }
   c(T1.i = cond.i, T1.ii = cond.ii,
     T2.i = cond.iii, T2.ii = cond.iv)
 }

"cliques" <-
function(amat){
### Finds the cliques of an UG. Mathias Drton, 2003.
    is.subset <- function(x,y)
      all(is.element(x,y))
    find.clique <- function(al, i, j){
      ## Local function
      ali <- c(i, al[[i]])
      if(length(ali)==1){
        return(c(i))
      }
      else{
        clique.lower <- c(i,j)
        clique.upper <-
          intersect(ali,c(j, al[[j]]))
        nb.remain <- setdiff(clique.upper, clique.lower)
        if(length(nb.remain)==0){
          return(clique.lower)
        }
        else{
          while(length(nb.remain)>=1){
            new.clique.lower <- c(clique.lower, nb.remain[1])
            new.clique.upper <-
              intersect(clique.upper,
                        c(nb.remain[1], al[[nb.remain[1]]]))
            if(!is.subset(new.clique.lower, new.clique.upper)){
              nb.remain <- setdiff(nb.remain, nb.remain[1])
            }
            else{
              clique.lower <- new.clique.lower
              clique.upper <- new.clique.upper
              nb.remain <- setdiff(clique.upper, clique.lower)
            }
          }
          return(sort(clique.lower))
        }
      }
    }
    ## Finds the adjacency list
    alist <- split(amat, row(amat))
    nod <- 1:nrow(amat) #  nod <- vertices(amat)
    alist <- lapply(alist, function(x) nod[x==1])

    ## Start
    cli <- c()
    index <- 1:length(alist)
    al <- alist
    while(length(index)>=1){
      i <- index[1]
      nbi <- al[[i]]
      if(length(nbi)==0){
        cli <- c(cli, c(i))
        index <- index[-1]
      }
      else{
        while(length(nbi)>=1){
          new.clique <- find.clique(alist,i, nbi[1])
          cli <- c(cli, list(new.clique))
          for(j in setdiff(new.clique, i)){
            al[[j]] <- setdiff(al[[j]], new.clique)
            if(length(al[[j]])==0){
              index <- setdiff(index, j)
            }
          }
          nbi <- setdiff(nbi, new.clique)
        }
        index <- setdiff(index, i)
      }
    }
    nod <- vertices(amat)
    lapply(cli, function(x) nod[x])
  }

"cmpGraph" <-
function(amat){
### Adjacency matrix of the complementary graph
    g <- 1*!amat
    diag(g) <- 0
    g
  }

"conComp" <-
function (amat) 
{
### Finds the connected components of an UG graph from its adjacency matrix amat.
  if(!all(amat == t(amat)))
    stop("Not an undirected graph.")
  A  <- transClos(amat)
  diag(A) <- 1
  n <- nrow(A)
  A <- sign(A + t(A))
  u <- A %*% 2^((n - 1):0)
  return(match(u, unique(u)))
}

"correlations" <-
function (x)
{
### Marginal correlations (lower half) and
### partial correlations given all remaining variables (upper half).
  
  if(is.data.frame(x))
    r <- cor(x)
  else  { # Recomputes the corr matrix
    Dg <- 1/sqrt(diag(x))
    r <- x * outer(Dg, Dg)
  }
  rp <- parcor(r)
  r[upper.tri(r)] <- rp[upper.tri(rp)]
  r
}

"cycleMatrix" <-
function(amat){
### Fundamental Cycle matrix of the UG amat.
    fc <- fundCycles(amat)  # Edges of the fundamental cycles
    E <- allEdges(amat)     # All the edges of the graph
    n <- nrow(E)            # Number of edges
    k <- length(fc)         # Number of FC
    if(k == 0) return(NULL)
    cmat <- matrix(0, k, n)
    for(cy in 1:k) {
      M <- fc[[cy]]         # Edges in cycle cy
      for(j in 1:nrow(M)) {
        e <- sort(M[j,])   
        for(i in 1:n){          
          cmat[cy, i] <- cmat[cy, i] | all(E[i,] == e)
        }
      }
    }
    dimnames(cmat) <- list(1:k, paste(E[,1], E[,2]))
    cmat       
  }

"DAG" <-
function (...,order=FALSE) 
{
### Defines a DAG from a set of equations (defined with model formulae).
  f <- list(...)
  nb <- length(f)  # nb is the number of model formulae (of blocks)
  nod <- c()       # Counts the number of nodes
  for (k in 1:nb) {
    tt <- terms(f[[k]])
    vars <- dimnames(attr(tt, "factors"))[[1]]
    nod <- c(nod, vars)
  }
  N <- unique(nod) # set of nodes
  dN <- length(N)  # number of nodes
  amat <- matrix(0,dN,dN)
  for (k in 1:nb) {
    tt <- terms(f[[k]])
    vars <- dimnames(attr(tt, "factors"))[[1]]
    if (attr(tt, "response") == 1) {
      j <- match(vars[1], N)
      i <- match(vars[-1], N)
      amat[i, j] <- 1
    }
    else if (attr(tt, "response") == 0) 
      stop("Some equations have no response")
  }
    if(!isAcyclic(amat))
      warning("The graph contains directed cycles!")
  dimnames(amat) <- list(N, N)
  if(order){
    amat <- topSort(amat)
  }
  amat
}

"drawGraph" <-
function(amat, coor=NULL, adjust=TRUE, alpha = 1.2, beta = 3, lwd=1, ecol = "blue", left=FALSE){
### A small function to plot mixed graphs with adjancency matrix amat.

    ## Some local functions
    
    plot.dots <- function(xy, v, dottype, n, beta) {
      for (i in 1:n) {
        if (dottype[i] == 1) {
          points(xy[i, 1], xy[i, 2], pch = 1, cex = 1.2, lwd=lwd)
        }
        else if (dottype[i] == 2) {
          points(xy[i, 1], xy[i, 2], pch = 16, cex = 1.2)
        }
      }
      text(xy[, 1] - beta, xy[, 2] + beta, v, cex = 1.2)
    }
    double.edges <- function(x1, x2, y1, y2, lwd, ecol){
      ## Double edges  for summary graphs
      arrows(x1, x2, y1, y2, lty=1, code=1, angle=20, length=.1, lwd=lwd, col=ecol)
      m1 <- (x1+y1)/2; m2 <- (x2+y2)/2
      a1 <- c(x1, m1+alpha, y1);  a2 <- c(x2, m2+alpha, y2)
      lines(a1, a2,  lwd=lwd, lty=2, col=ecol)
    }
    draw.edges <-  function(coor, u, alpha, type, lwd, ecol) {
      ## type = 0, 1, 2, 3 for undirected, directed, bidirected, double edges
      for (k in 1:nrow(u)) {
        a <- coor[u[k, 1], ]
        b <- coor[u[k, 2], ]
        ba <- b-a
        ba <- ba/sqrt(sum(ba*ba))
        x <- a + ba*alpha
        y <- b - ba*alpha
        switch(type+1, 
               segments(x[1], x[2], y[1], y[2], lty = 1, lwd=lwd, col=ecol),
               arrows(x[1], x[2], y[1], y[2], code=2, angle=20, length=.1, lty = 1, lwd=lwd, col=ecol),
##               {arrows(x[1], x[2], y[1], y[2], code=3,length = 0.1, lty = 1, lwd=lwd);
##                segments(x[1], x[2], y[1], y[2], lty=1, lwd=lwd, col="white");
##                segments(x[1], x[2], y[1], y[2], lty=2, lwd=lwd)},
               arrows(x[1], x[2], y[1], y[2], code=3, angle=20, length = 0.1, lty = 5, lwd=lwd, col=ecol),
               double.edges(x[1], x[2], y[1], y[2], lwd=lwd, ecol)
               )
      }
    }
    def.coor <- function(ce, k, h, w) {
      ## Default coords of the dots for an UG.
      if (k == 1) 
        return(ce)
      else if (k == 2) {
        r1 <- c(ce[1], ce[1])
        r2 <- c(ce[2] + h * 0.3, ce[2] - h * 0.3)
      }
      else if (k == 3) {
        r1 <- c(ce[1], ce[1], ce[1])
        r2 <- c(ce[2] + h * 0.25, ce[2], ce[2] - h * 0.25)
      }
      else if (k == 4) {
        r1 <- c(ce[1] - w * 0.3, ce[1] + w * 0.3, ce[1] + 
                w * 0.3, ce[1] - w * 0.3)
        r2 <- c(ce[2] - h * 0.3, ce[2] - h * 0.3, ce[2] + 
                h * 0.3, ce[2] + h * 0.3)
      }
      else {
                                        # a <- runif(1, 0, pi)
        a <- 1
        z <- seq(a, a + 2 * pi, len = k + 1)
        z <- z[-1]
        r1 <- ce[1] + w/2.5 * cos(z)
        r2 <- ce[2] + h/2.5 * sin(z)
      }
      cbind(r1, r2)
    }
    def.coor.dag <- function(amat, w, h, left){
      nod <- rownames(amat)
      o <- topOrder(amat)
      if(left)
        o <- rev(o)
      k <- length(nod)
      x <- seq(0, 100, len=k)
      y <- rep(c(20, 40, 60, 80),ceiling(k/4))[1:k] 
      xy <- cbind(x, y)
      rownames(xy) <- nod[o]
      xy[nod,]
    }    
    
    ## Basic definitions
    
    v <- rownames(amat)
    n <- length(v)

    ## The dot type is fixed to 1.
    
    dottype <- rep(1,n)

    ## Initialize plot
    old <- par(mar=c(0,0,0,0))
    on.exit(par(old))
    
    plot(c(0, 100), c(0, 100), type = "n", axes = FALSE, xlab = "", 
         ylab = "")

    center <- matrix(c(50, 50), ncol = 2)
    
    ## Default coordinates of the dots    
    if (is.null(coor)) {
      isdag <- isAcyclic(amat)
      if(isdag)
        coor <- def.coor.dag(amat, 100, 100, left=left)
      else
        coor <- def.coor(center, n, 100, 100)
    }
    ## Find the undirected edges
    g0 <- amat * ((amat == 1) & (t(amat) == 1))
    g0[lower.tri(g0)] <- 0

    ## Find the directed edges                                    
    g1 <- amat * ((amat == 1)  & !((amat > 0) & (t(amat) > 0)))
    
    ## Find the bidirected edges (for ancestral graphs)
    g2 <- amat * ((amat == 2) & (t(amat) == 2))
    g2[lower.tri(g2)] <- 0    

    ## Find the double edges (for summary graphs)    
    g3 <-  amat *( ((amat == 2)  & (t(amat) == 3)))
    
    i <- expand.grid(1:n, 1:n)
    u0 <- i[g0 > 0, ]
    u1 <- i[g1 > 0, ]
    u2 <- i[g2 > 0, ]
    u3 <- i[g3 > 0, ]
    
    ## Plot the dots and the edges
    if(nrow(coor) !=  length(v))
      stop("Wrong coordinates of the vertices.")  
    plot.dots(coor, v, dottype, n, beta)  
    draw.edges(coor, u0, alpha, type=0, lwd=lwd, ecol) 
    draw.edges(coor, u1, alpha, type=1, lwd=lwd, ecol) 
    draw.edges(coor, u2, alpha, type=2, lwd=lwd, ecol)
    draw.edges(coor, u3, alpha, type=3, lwd=lwd, ecol) 

    ## Adjust the plot
    if(adjust){
      repeat {
        xnew <- unlist(locator(1))
        if (length(xnew) == 0) {
          break
        }
        d2 <- (xnew[1] - coor[, 1])^2 + (xnew[2] - coor[, 2])^2
        i <- (1:n)[d2 == min(d2)]
        coor[i, 1] <- xnew[1]
        coor[i, 2] <- xnew[2]
        plot(c(0, 100), c(0, 100), type = "n", axes = FALSE,
             xlab = "", ylab = "")
        plot.dots(coor, v, dottype, n, beta)
        draw.edges(coor, u0, alpha, type=0, lwd=lwd, ecol)
        draw.edges(coor, u1, alpha, type=1, lwd=lwd, ecol) 
        draw.edges(coor, u2, alpha, type=2, lwd=lwd, ecol)
        draw.edges(coor, u3, alpha, type=3, lwd=lwd, ecol) 
      }
    }
    colnames(coor) <- c("x", "y")
    return(invisible(coor))
  }

"dSep" <-
function(amat, first, second, cond) {
### Are first and second d-Separated by cond in a DAG? 
    e <- inducedCovGraph(amat, sel=c(first,second), cond=cond)
    all(e[first,second] == 0)
  }

"edgeMatrix" <-
function (E, inv=FALSE) 
{
### From the adjacency matrix to the edge matrix
  if(inv){
    ord <- topOrder(E)
    ord <- rev(ord) # Inverse topological order: Nanny ordering.
    E <- E[ord, ord]
  }
  A <- t(E)
  diag(A) <- 1
  A
}

"essentialGraph" <-
function(dagx){
### Converts a DAG into Essential Graph. 
### Is implemented by the algorithm by D.M.Chickering (1995).
### A transformational characterization of equivalent Bayesian network
### structures. Proceedings of Eleventh Conference on Uncertainty in
### Artificial Intelligence, Montreal, QU, pages 87-98. Morgan Kaufmann 
### http://research.microsoft.com/~dmax/publications/uai95.pdf 
### Implemented in Matlab by Tomas Kocka, AAU.
### Translated in R by Giovanni Marchetti, University of Florence.
  
  ord <- topOrder(dagx);      # get the topological order of nodes 
  n <- nrow(dagx)             # gets the number of nodes
  i <- expand.grid(1:n, 1:n)  # finds all nonzero elements in the adj matrix
  IJ <- i[dagx==1,]           # sort the arcs from lowest possible y
  I <- IJ[, 1]; J <- IJ[, 2]  # and highest possible x, arcs are x->y
  e <- 1
  for(y in 1:n){
    for(x in n:1){
      if(dagx[ord[x], ord[y]] == 1) { 
        I[e] <- ord[x]
        J[e] <- ord[y]
        e <- e + 1
      }
    }
  }
  ## Now we have to decide which arcs are part of the essential graph and
  ## which are undirected edges in the essential graph.
  ## Undecided arc in the DAG are 1, directed in EG are 2 and undirected in EG are 3.
  
  for(e in 1:length(I)){
    if(dagx[I[e],J[e]] == 1){
      cont <- TRUE
      for(w in 1:n){ 
        if(dagx[w,I[e]] == 2){
          if(dagx[w,J[e]] != 0)
            dagx[w,J[e]] <- 2
          else {
            for(ww in 1:n){
              if(dagx[ww,J[e]] != 0)
                      dagx[ww,J[e]] <- 2
            } # skip the rest and start with another arc from the list
            w <- n
            cont <- FALSE
          }
        }
      }
      if(cont){
        exists <- FALSE
        for(z in 1:n){
          if((dagx[z,J[e]] != 0) & (z != I[e]) & (dagx[z,I[e]] == 0)){
            exists <- TRUE
            for(ww in 1:n){
              if(dagx[ww,J[e]] == 1){
                dagx[ww,J[e]] <- 2
              }
            }
          }
        }
        if(!exists){
          for(ww in 1:n){
            if(dagx[ww,J[e]] == 1){
              dagx[ww,J[e]] <- 3
            }
          } 
        }
      }
    }          
  }
  (dagx==2) + (dagx==3) + t(dagx==3)
}

"findPath" <-
function (amat, st, en, path = c()) 
{
### Find a path between nodes st and en in a UG with adjacency mat. amat.
  indices <- 1:nrow(amat)
  if(st == en) # st is 'node' in recursive calls
    return(c(path, st))
  if(sum(amat[st,]) == 0 ) 
    return(NULL)
  ## ne <- bd(st,amat)
  ne <- indices[amat[st,]==1] # Boundary of x. Assumes that amat is symmetric
  for(node in ne){
    if(!is.element(node, c(path, st))){
      newpath <- findPath(amat, node, en, c(path, st))
      if(!is.null(newpath))
        return(newpath)
    }
  }
}

"fitAncestralGraph" <-
function (amat, S, n, tol = 1e-06){
### Fit Ancestral Graph. Mathias Drton, 2003.
    nam <- rownames(S)
    nod <- vertices(amat)
    ## permute graph to have same layout as S
    if(is.null(nod)){
      stop("The adjacency matrix has no labels!")
    }
    if(!all(is.element(nod, nam)))
      stop("The nodes of the graph do not match the names of the variables")
    else
      sek <- intersect(nam, nod)
    S <- S[sek,sek, drop=FALSE]              # Resizes eventually S
    amat <- amat[sek,sek, drop=FALSE]        # and reorders amat
    
    temp <- icfmag(amat, S, tol)
    lik <-  function(K, S, n, k){
      ## Computes the deviance
      trace <- function(a) sum(diag(a)) 
      SK <- S %*% K 
      (trace(SK) - log(det(SK)) - k) * n
    }
    p <- ncol(S)
    df <- p*(p-1)/2 - sum(In(amat+t(amat)))/2   # Degrees of freedom 
    dev <- lik(solve(temp$Sigmahat), S, n, p)
    if(is.null(temp$Bhat)){
      Beta <- NULL
    }
    else{
      ## Beta <- diag(1,p)-temp$Bhat
      Beta <- temp$Bhat
    }
    return(list(Shat=temp$Sigmahat, Lhat=temp$Lambdahat, Bhat=Beta,
                Ohat=temp$Omegahat, dev = dev, df = df, it=temp$iterations))
  }

"fitConGraph" <-
function (amat, S, n, pri=FALSE, alg=2, tol = 1e-06)
{
### Fits a concentration graph G.
  nam <- rownames(S)
  nod <- vertices(amat)
  if(is.null(nod)){
    stop("The adjacency matrix has no labels.")
  }
  if(!all(is.element(nod, nam)))
    stop("The nodes of the graph do not match the names of the variables.")
  else
    sek <- intersect(nam, nod)
  S <- S[sek,sek, drop=FALSE]              # Resizes eventually S
  amat <- amat[sek,sek, drop=FALSE]        # and reorders amat
  nod <- vertices(amat)        
  cli <- cliques(amat)              # Finds the cliques

  lik <-  function(K, S, n, k){
    ## Computes the deviance
    trace <- function(a) sum(diag(a)) 
    SK <- S %*% K 
    (trace(SK) - log(det(SK)) - k) * n
  }
  
  nc <- length(cli)
  if(nc==1)
    return(list(Shat = S, dev = 0, df = 0, it=1))
  k <- ncol(S)
  it <- 0
  if(alg == 1){
    V <-   diag(diag(S))  # Starting value
    dimnames(V) <- dimnames(S)
    repeat {
      V.old <- V
      it <- it+1
      for (i in 1:nc) {
        a <- cli[[i]]
        b <- setdiff(nod, a)
        Saa <- S[a, a]
        Vaa <- V[a, a]
        Vba <- V[b, a]
        Vbb <- V[b, b]
        B <- Vba %*% solve(Vaa)
        Spar <- Vbb - B %*% Vaa %*% t(B)
        BV <- B %*% Saa
        V[b, a] <- BV
        V[a, b] <- t(BV)
        V[a, a] <- Saa
        V[b, b] <- Spar + BV %*% t(B)
        if(pri) {
          dev <- lik(solve(V), S, n, k)
          cat(dev, "\n")
        }
      }
      if(sum(abs(V-V.old)) < tol) break
    }
  }
  else if(alg==2) {
     K <-   solve(diag(diag(S)))  # Starting value
    dimnames(K) <- dimnames(S)
    repeat {
      K.old <- K
      it <- it+1
      for (i in 1:nc) {
        a <- cli[[i]]
        b <- setdiff(nod, a)
        K[a,a] <- solve(S[a,a]) + K[a,b] %*% solve(K[b,b]) %*% K[b,a] 
        if(pri) {
          dev <- lik(K, S, n, k)
          cat(dev, "\n")
        }
      }
      if(sum(abs(K-K.old)) < tol) break
    }
     V <- solve(K)
  }
  df <- sum(amat[upper.tri(amat)] == 0) # Degrees of freedom
  dev <- lik(solve(V), S, n, k) 
  list(Shat = V, dev = dev, df = df, it=it)
}

"fitCovGraph" <-
function (amat, S, n, alg="icf", dual.alg=2, start.icf=NULL, tol = 1e-06){
### Fits a Covariance Graph. Mathias Drton, 2003
### amat: adjacency matrix; S: covariance matrix; n: sample size.
    amat <- In(amat) # Forces the ones in a bidirected graph defined with makeAG
    nam <- rownames(S)
    nod <- vertices(amat)
    if(is.null(nod)){
      stop("The adjacency matrix has no labels.")
    }
    if(!all(is.element(nod, nam)))
      stop("The nodes of the graph do not match the names of the variables.")
    else
      sek <- intersect(nam, nod)
    S <- S[sek,sek, drop=FALSE]              # Resizes eventually S
    amat <- amat[sek,sek, drop=FALSE]        # and reorders amat
    
    if(alg=="icf"){
      temp <- icf(amat, S, start.icf, tol)
    }
    else{
      if(alg == "dual"){
        Sinv <- solve(S)
        temp <- fitConGraph(amat, Sinv, n, pri=FALSE, alg=dual.alg, tol = 1e-06)
        temp <- list(Sigmahat=zapsmall(solve(temp$Shat)), iterations=temp$it)
      }
      else{
        stop("Algorithm misspecified!")
      }
    }
    lik <-  function(K, S, n, k){
      ## Computes the deviance
      trace <- function(a) sum(diag(a)) 
      SK <- S %*% K 
      (trace(SK) - log(det(SK)) - k) * n
    }
    df <- sum(amat[upper.tri(amat)] == 0) # Degrees of freedom
    k <- ncol(S)
    dev <- lik(solve(temp$Sigmahat), S, n, k) 
    return(list(Shat=temp$Sigmahat, dev = dev, df = df, it=temp$iterations))
  }

"fitDag" <-
function (amat, S, n)
{
### Fits linear recursive regressions with independent residuals. 
### amat: the adjacency matrix of the DAG. S: cov matrix. n: sample size.
  if(missing(amat)){ # saturated model
    amat <-  lower.tri(diag(ncol(S)), diag=FALSE) * 1
    dimnames(amat) <- dimnames(S)
  }
  nam <- rownames(S)
  nod <- vertices(amat)
  if(is.null(nod))
    stop("The adjacency matrix has no labels.")
  if(!all(is.element(nod, nam)))
    stop("The nodes of the graph do not match the names of the variables.")
  else
    sek <- intersect(nam, nod) 
  S <- S[sek,sek]              # Resizes eventually S 
  amat <- amat[sek,sek]        # and reorders amat
  Delta <- rep(length(sek),0)
  
  emat <- edgeMatrix(amat) 
  A <- emat
  p <- ncol(S)
  ip <- 1:p
  for(i in 1:p) {
    u <- emat[i,]
    v <- ip[u == 1 & ip != i]
    M <- swp(S, v)[i,]
    A[i, ] <- - A[i, ] * M
    A[i,i] <- 1
    k <- sum(A[i,])
    Delta[i] <- M[i]
  }
  names(Delta) <- sek
  B <- solve(A)
  Shat <- B %*% diag(Delta) %*% t(B)
  dimnames(Shat) <- dimnames(S)
  Khat <- solve(Shat)
  H <- S %*% Khat
  trace <- function(A) sum(diag(A))
  dev <- (trace(H) - log(det(H)) - p) * n
  df <- p*(p+1)/2 - sum(amat==1) - p
  list(Shat=Shat, Ahat = A, Dhat=Delta, dev=dev, df=df)
}

"fitDagLatent" <-
function (amat, Syy, n, latent, norm = 1,  seed=144, maxit=9000, tol=1e-6, pri=FALSE) 
{
### Fits linear recursive regressions with independent residuals and one latent variable.
### Syy: covariance matrix, n: sample size, amat: adjacency matrix.
### NOTE: both amat and Syy must have rownames.
### latent is the "name" of the latent in the rownames of Syy. norm = normalisation type.

  ## Local functions
  setvar1 <- function (V, z, paz, norm)
    ## paz are the parents of the latent (needed if norm=2)
    {
      ## Normalizes V
      if(norm == 1){
        ## Rescales V forcing V[z,z] = 1
        a <- 1 / sqrt(V[z,z])
      }
      else if(norm== 2) {
        ## Rescales V forcing Delta[z,z] = 1
        if(sum(paz) > 0)
          sig <- V[z,z] - V[z, paz] %*% solve(V[paz,paz]) %*% V[paz,z]
        else
          sig <- V[z,z]
        a <- 1/sqrt(sig)
      }
      V[z,] <- V[z,] * a
      V[,z] <- V[,z] * a
      V
    }

  cmqi <- function (Syy, Sigma, z) 
    {
      ## Computes the matrix C(M | Q) by Kiiveri (1987), Psychometrika.
      ## It is a slight generalization in which Z is not the last element.
      ## z is a Boolean vector indicating the position of the latent variable in X.
      y <- ! z
      Q <- solve(Sigma)
      Qzz <- Q[z,z]
      Qzy <- Q[z,y]
      B <- - solve(Qzz) %*% Qzy
      BSyy <- B %*% Syy
      E <- Sigma*0 
      E[y,y] <- Syy
      E[y,z] <- t(BSyy)
      E[z,y] <- BSyy
      E[z,z] <- BSyy %*% t(B) + solve(Qzz)
      dimnames(E) <- dimnames(Sigma) 
      E
    }
  fitdag <- function (amat, S,n, constr=NULL)
    {
      ## Fits linear recursive regressions with independent residuals (fast version).
      ## NOTE. amat and S must have the same size and variables. constr is a matrix
      ## indicating the edges that must be constrained to 1.
      emat <- edgeMatrix(amat)
      A <- emat
      p <- ncol(S)
      Delta <- rep(p,0)
      ip <- 1:p
      for(i in ip) {
        u <- emat[i,]
        v <- ip[(u == 1) & (ip != i)]           # Parents 
        if(length(v) == 0){                     # If pa is empty
          Delta[i] <- S[i,i]
          next
        }
        M <- lmfit(S, y=i, x=v, z=(constr[i, ] == 1))
        A[i, ] <- - A[i, ] * M
        A[i,i] <- 1
        k <- sum(A[i,])
        Delta[i] <- M[i]
      }
      Khat <- t(A) %*% diag(1/Delta) %*% A
      Shat <- solve(Khat)
      list(A = A, Delta=Delta, Shat=Shat, Khat=Khat)
    }
  lmfit <- function(S, y, x, z=NULL){
    ## Regression coefficients of y given x eventually with z constrained to 1.
    ## Residual variance in position y.
    Sxy <- S[x, y, drop=FALSE] - apply(S[x, z, drop=FALSE], 1, sum) 
    Sxx <- S[x, x, drop=FALSE]
    bxy <- solve(Sxx,Sxy)
    out<- rep(0, nrow(S))
    out[x] <- bxy
    out[z] <- 1
    names(out) <- rownames(S)
    xz <- c(x,z)
    b <- out[xz, drop=FALSE]
    res <- S[y,y] + t(b) %*% (S[xz, xz, drop=FALSE] %*% b - 2 * S[xz,y, drop=FALSE])
    out[y] <- res
    out
  }
  lik <-  function(K, S, n, p){
    ## Computes the deviance
    trace <- function(a) sum(diag(a)) 
    SK <- S %*% K 
    (trace(SK) - log(det(SK)) - p) * n
  }
  ## Beginning  of the main function
  nam <- rownames(Syy)        # Names of the variables (they can be more than the nodes)
  nod <- vertices(amat)       # Names of the nodes of the DAG (that contains the latent)
 
  if(is.null(nod))
    stop("The adjacency matrix has no labels.")
  if(!all(is.element(setdiff(nod, latent), nam)))
    stop("The observed nodes of the graph do not match the names of the variables.")
  else
    sek <- intersect(nam, nod)
  Syy <- Syy[sek,sek, drop=FALSE]          # Resizes eventually S
  sek <- c(sek, latent)
  amat <- amat[sek,sek, drop=FALSE]        # and reorders amat
  nod <- vertices(amat)
  paz <- pa(latent, amat)                  # Parents of the latent
  paz <- is.element(nod, paz)
  dn <- list(nod,nod)
  wherez <- is.element(nod, latent)        # Index of the latent
  if(is.null(wherez))
    stop("Wrong name of the latent variable!")
  wherey <- ! wherez
  p <- ncol(Syy)
  df <- p*(p+1)/2  - sum(amat==1) - p      # Degrees of freedom 
  
  if(df <= 0)
    warning(paste("The degrees of freedom are ", df))
  
  set.seed(seed)     # For random starting value
  Sigma.old <- rcorr(p+1)
  Sigma.old <- setvar1(Sigma.old, wherez, paz, norm=norm)
  dimnames(Sigma.old) <- dn    
  it <- 0
  repeat{
    it <- it+1
    if(it > maxit){
      warning("Maximum number of iterations reached.")
      break
    }
    Q <- cmqi(Syy, Sigma.old, wherez) # E-step. See Kiiveri, 1987                                        
    fit <- fitdag(amat, Q, n)         # M-step
    Sigma.new <- fit$Shat
                                       
    if(pri) {
      dev <-  lik(fit$Khat, Q, n, (p+1)) # Monitoring progress of iterations
      cat(dev, "\n")
    }
    else{
      if(0==(it %% 80))
        cat("\n")
      else
        cat(".")
    }                                     
    if(sum(abs(Sigma.new - Sigma.old)) < tol) break  # Test convergence
    Sigma.old <- Sigma.new
  }
  cat("\n")
  dev <-  lik(fit$Khat, Q, n, (p+1))
  Shat <- setvar1(Sigma.new, wherez, paz, norm=norm) # Normalize Shat
  dimnames(Shat) <- dn
  Khat <- solve(Shat)
  fit <- fitDag(amat, Shat, n)
  list(Shat=Shat, Ahat=fit$Ahat, Dhat=fit$Dhat, dev=dev, df=df, it=it) 
}

"fundCycles" <-
function(amat){
### Finds a set of fundamental cycles for an UG with adj. matrix amat.
    fc <- c()
    tr <- bfs(amat) # Spanning tree
    if(is.null(tr)) return(NULL)
    if(is.null(tr$chords)) return(NULL)
    co <- tr$chords # edges of the cospanning tree
    for(i in 1:nrow(co)) {
      e <- co[i,] 
      g <- tr$tree # edge matrix of the spanning tree
      cy <- findPath(g, st=e[1], en=e[2])
       splitCycle <- function(v){
         ## Splits a cycle v into a matrix of edges.
         v <- c(v, v[1])
         cbind(v[-length(v)], v[-1])
       }
      cy <- splitCycle(cy)
      fc <- c(fc, list(cy))
    }
    fc 
  }

"icf" <-
function(bi.graph, S, start=NULL, tol = 1e-06){
### Iterative conditional fitting for bidirected graphs. Mathias Drton, 2003
    if(!is.matrix(S)){
      stop("Second argument is not a matrix!")
    }
    if(dim(S)[1]!=dim(S)[2]){
      stop("Second argument is not a square matrix!")
    }
    if(min(eigen(S)[[1]])<=0){
      stop("Second argument is not a positive definite matrix!")
    }

    p <- nrow(S)
    i <- 0

    ## prep spouses and non-spouses
    
    pa.each.node <-function(amat){
      ## List of the parents of each node.
      ## If amat is symmetric it returns the boundaries.
      p <- nrow(amat)
      b <- vector(p, mode="list")
      ip <- 1:p
      for(i in 1:p)
        b[[i]] <- ip[amat[,i]==1]
      b 
    }
    spo <- pa.each.node(bi.graph)
    nsp <- pa.each.node(cmpGraph(bi.graph))
    number.spouses <- unlist(lapply(spo, length))
    no.spouses <- (1:p)[number.spouses==0]
    all.spouses <- (1:p)[number.spouses==(p-1)]
    nontrivial.vertices <- setdiff((1:p), no.spouses)

    if(length(nontrivial.vertices)==0){
      if(p==1){
        Sigma <- S
      }
      else{
        Sigma <- diag(diag(S))
        dimnames(Sigma) <- dimnames(S)
      }
      return(list(Sigmahat=Sigma, iterations=1))
    }

    if(is.null(start)){
      Sigma <- as.matrix(diag(diag(S))) # starting value
    }
    else{
      temp <- diag(start)
      start[bi.graph==0] <- 0
      diag(start) <- temp
      diag(start)[no.spouses] <- diag(S)[no.spouses]
      Sigma <- as.matrix(start)
      if(min(eigen(Sigma)$values) <= 0){
        stop("Starting value is not feasible!")
      }
    }
    
    repeat{
      i <- i+1
      Sigma.old <- Sigma
      for(v in nontrivial.vertices){
        if(is.element(v, all.spouses)){
          B <- S[v,-v]%*%solve(S[-v,-v])
          lambda <- S[v,v]-B%*%S[-v,v]
        }
        else{
          B.spo.nsp <-
            Sigma[spo[[v]],nsp[[v]]]%*%solve(Sigma[nsp[[v]],nsp[[v]]])
          YZ <- S[v,spo[[v]]]-S[v,nsp[[v]]]%*%t(B.spo.nsp) 
          B.spo <- 
            YZ %*%
              solve( S[spo[[v]],spo[[v]]]
                    -S[spo[[v]],nsp[[v]]]%*%t(B.spo.nsp)
                    -B.spo.nsp%*%S[nsp[[v]],spo[[v]]]
                    +B.spo.nsp%*%S[nsp[[v]],nsp[[v]]]%*%t(B.spo.nsp) )
          lambda <- S[v,v]-B.spo%*%t(YZ)
          B.nsp <- -B.spo%*%B.spo.nsp
          B <- rep(0, p)
          B[spo[[v]]] <- B.spo
          B[nsp[[v]]] <- B.nsp
          B <- B[-v]
        }
        ## here I can improve by only using B[spo[[v]]]!
        Sigma[v,-v] <- B%*%Sigma[-v,-v]
        Sigma[v,nsp[[v]]] <- 0
        Sigma[-v,v] <- t(Sigma[v,-v])
        Sigma[v,v] <- lambda + B%*%Sigma[-v,v]
      }
      if(sum(abs(Sigma.old-Sigma)) < tol) break
    }
    dimnames(Sigma) <- dimnames(S)
    return(list(Sigmahat=Sigma, iterations=i))
  }

"icfmag" <-
function(mag, S, tol = 1e-06){
    ## Iterative conditional fitting for ancestral graphs. Mathias Drton, 2003.
    if(!is.matrix(S)){
      stop("Second argument is not a matrix!")
    }
    if(dim(S)[1]!=dim(S)[2]){
      stop("Second argument is not a square matrix!")
    }
    if(min(eigen(S)[[1]])<=0){
      stop("Second argument is not a positive definite matrix!")
    }
    p <- nrow(S)  # Dimensionality
    temp <- unmakeAG(mag)
    mag.ug <- temp$ug
    mag.dag <- temp$dag
    mag.bg <- temp$bg
    
    ## Catch trivial case
    if(p==1){
      return(list(Sigmahat=S, Omegahat=S, Bhat=NULL, Lambdahat=NULL, iterations=1))
    }
    
    ## Starting value
    Omega <- diag(diag(S))
    dimnames(Omega) <- dimnames(S)
    B <- diag(p)
    dimnames(B) <- dimnames(S)

    ## IPS for UG
    
    UG.part <- (1:p)[0==apply(mag.dag + mag.bg,2,sum)] 
    if(length(UG.part)> 0){
      Lambda.inv <-
        fitConGraph(mag.ug[UG.part,UG.part, drop=FALSE],
                    S[UG.part,UG.part, drop=FALSE], p+1,tol = tol)$Shat
      Omega[UG.part,UG.part] <- Lambda.inv
    }

    ## Prepare list of spouses, parents, etc.
    pa.each.node <-function(amat){
      ## List of the parents of each node.
      ## If the adjacency matrix is symmetric it gives the boundary.
      p <- nrow(amat)
      b <- vector(p, mode="list")
      ip <- 1:p
      for(i in 1:p)
        b[[i]] <- ip[amat[,i]==1]
      b 
    }
    spo <- pa.each.node(mag.bg)
    nsp <- pa.each.node(cmpGraph(mag.bg))
    pars <- pa.each.node(mag.dag) 
    
    i <- 0
    repeat{
      i <- i+1
      Omega.old <- Omega
      B.old <- B
      for(v in setdiff(1:p, UG.part)){
        parv <- pars[[v]]
        spov <- spo[[v]]
        if(length(spov)==0){
          if(length(parv)!=0){
            if(i == 1){ # do it only once
              ## attention: B = - beta
              B[v,parv] <- -S[v,parv]%*%solve(S[parv,parv])
              Omega[v,v] <- S[v,v]+B[v,parv]%*%S[parv,v]
            }
          }
        }
        else{
          if(length(parv)!=0){
            O.inv <- matrix(0, p,p)
            O.inv[-v,-v] <- solve(Omega[-v,-v])
            Z <- O.inv[spov,-v] %*%B[-v,]
            lpa <- length(parv)
            lspo <- length(spov)
            XX <- matrix(0, lpa+lspo, lpa+lspo)
            XX[1:lpa, 1:lpa] <- S[parv,parv]
            XX[1:lpa,(lpa+1):(lpa+lspo)] <- S[parv,]%*%t(Z)
            XX[(lpa+1):(lpa+lspo),1:lpa] <- t(XX[1:lpa,(lpa+1):(lpa+lspo)])
            XX[(lpa+1):(lpa+lspo),(lpa+1):(lpa+lspo)] <- Z%*%S%*%t(Z)
            YX <- c(S[v,parv], S[v,]%*%t(Z))
            temp <- YX %*% solve(XX)
            B[v,parv] <- -temp[1:lpa]
            Omega[v,spov] <- temp[(lpa+1):(lpa+lspo)]
            Omega[spov,v] <- Omega[v,spov]
            
            temp.var <- S[v,v] - temp %*% YX
            Omega[v,v] <- temp.var +
              Omega[v,spov] %*% O.inv[spov,spov] %*% Omega[spov,v]
          }
          else{
            O.inv <- matrix(0, p,p)
            O.inv[-v,-v] <- solve(Omega[-v,-v])
            Z <- O.inv[spov,-v] %*%B[-v,]
            XX <- Z%*%S%*%t(Z)
            YX <- c(S[v,]%*%t(Z))
            Omega[v,spov] <- YX %*% solve(XX)
            Omega[spov,v] <- Omega[v,spov]
            
            temp.var <- S[v,v] -  Omega[v,spov] %*% YX
            Omega[v,v] <- temp.var +
              Omega[v,spov] %*% O.inv[spov,spov] %*%
                Omega[spov,v]
          }
        }
      }
      if(sum(abs(Omega.old-Omega)) + sum(abs(B.old-B)) < tol) break
    }
    Sigma <- solve(B)%*%Omega%*%solve(t(B))
    Lambda <- Omega
    Lambda[-UG.part,-UG.part] <- 0
    Omega[UG.part,UG.part] <- 0
    return(list(Sigmahat=Sigma, Bhat=B, Omegahat=Omega, Lambdahat=Lambda,
                iterations=i))
  }

"In" <-
function (A) 
{
### Indicator matrix of structural zeros.
  abs(sign(A)) 
}

"inducedChainGraph" <-
function(amat, cc=vertices(amat), cond=NULL, type="LWF"){
### Induced chain graph with chain components cc.
    inducedBlockGraph <- function(amat, sel, cond){
      S <- inducedConGraph(amat, sel=union(cond, sel), cond=NULL)
      S[cond, sel, drop=FALSE]
    }
    nod <- vertices(amat)
    nam <- c(list(cond), cc)
    if(!all(unlist(nam) %in% nod))
      stop("The chain components or the conditioning set are wrong.")
    for(h in nam){
      for(k in nam){
        h <- unlist(h)
        k <- unlist(k)
        if (setequal(h,k))
          next
        if(length(intersect(h,k) > 0))
          stop("The sets are not disjoint!")
      }
    }
    nam <- unlist(nam)
    cg <- matrix(0, length(nam),length(nam))
    dimnames(cg) <- list(nam,nam)
    kc <- length(cc)
    if(type=="AMP"){
      for(i in 1:kc){
        past <-  unlist(cc[0:(i-1)]) 
        Past <- union(cond,past)
        g <- cc[[i]]
        
        Sgg.r <- inducedConGraph(amat, sel=g, cond=Past)
        cg[g, g] <- Sgg.r
        if(length(past) !=0){
          Pgr <- inducedRegGraph(amat, sel=g, cond=Past)
          cg[Past, g] <- Pgr
        }
      }
    }
    else if(type=="LWF"){
      for(i in 1:kc){
        past <-  unlist(cc[0:(i-1)]) 
        Past <- union(cond,past)
        g <- cc[[i]]
        
        Sgg.r <- inducedConGraph(amat,sel=g, cond=Past)
        cg[g, g] <- Sgg.r
        if(length(past) !=0){
          Cgr <- inducedBlockGraph(amat, sel=g, cond=Past)
          cg[Past, g] <- Cgr
        }
      }
    }
    else if(type=="MRG"){
      for(i in 1:kc){
        past <-  unlist(cc[0:(i-1)]) 
        Past <- union(cond,past)
        g <- cc[[i]]
        
        Sgg.r <- inducedCovGraph(amat, sel=g, cond=Past)
        cg[g, g] <- Sgg.r
        if(length(past) != 0){
          Pgr <- inducedRegGraph(amat, sel=g, cond=Past)
          cg[Past, g] <- Pgr
        }
      }
    }
    else
      stop("Wrong type.")
    n <- unlist(cc)
    cg[n,n, drop=FALSE]
  }

"inducedConGraph" <-
function(amat, sel=vertices(amat), cond=NULL){
### Induced concentration graph for a set of nodes given a conditioning set.
    ancGraph <- function(A) {
      ## Edge matrix of the overall ancestor graph.
      if(sum(dim(A)) == 0)
        return(A)
      else
        return(In(solve(2*diag(nrow(A)) - A)))
    }
   
   trclos <- function(M) {
      ## Transitive closure of an UG with edge matrix M. See Wermuth and Cox (2004). 
     edgeMatrix(transClos(adjMatrix(M)))
    }

    A <- edgeMatrix(amat) # From the adjacency matrix to edge matrix    
    nod <- vertices(A)
      if(!all(cond %in% nod))
        stop("Conditioning nodes are not among the vertices.")
    if(!all(sel %in% nod))
      stop("Selected nodes are not among the vertices.")
    
    if(length(intersect(sel,cond) > 0))
      stop("The sets are not disjoint!")

    l <- setdiff(nod, union(sel, cond))  # Marginal nodes
    g <- sel
    r <- cond
    L <- union(l,g)
    R <- union(g,r)
    
    Al <- ancGraph(A[l,l,drop=FALSE])
    ARR.l <- In(A[R,R, drop=FALSE] +
                A[R,l, drop=FALSE]%*% Al %*% A[l,R, drop=FALSE])
    TRl <- In(A[R,l,drop=FALSE] %*% Al)
    DRl <- In(diag(length(R)) + TRl %*% t(TRl))
    out <- In(t(ARR.l) %*% trclos(DRl) %*% ARR.l)
    out <- out[g,g, drop=FALSE]
    adjMatrix(out)
  }

"inducedCovGraph" <-
function(amat, sel=vertices(amat), cond=NULL){
### Induced covariance graph for a set of nodes given a conditioning set.
    ancGraph <- function(A) {
      ## Edge matrix of the overall ancestor graph.
      if(sum(dim(A)) == 0)
        return(A)
      else
        return(In(solve(2*diag(nrow(A)) - A)))
    }
     trclos <- function(M) {
      ## Transitive closure of an UG with edge matrix M. See Wermuth and Cox (2004). 
     edgeMatrix(transClos(adjMatrix(M)))
    }
    A <- edgeMatrix(amat) # From the adjacency matrix to edge matrix
    nod <- vertices(A)
    if(!all(cond %in% nod))
      stop("Conditioning nodes are not among the vertices.")
    if(!all(sel %in% nod))
      stop("Selected nodes are not among the vertices.")
    if(length(intersect(sel,cond) > 0))
      stop("The sets are not disjoint!")
    l <- setdiff(nod, union(sel, cond))  # Marginalized over nodes
    g <- sel
    r <- cond
    L <- union(l,g)
    R <- union(g,r)
    
    AL <-  ancGraph(A[L,L,drop=FALSE]) # In(solve(2*diag(length(L)) - A[L,L]))
    TrL <- In(A[r,L,drop=FALSE] %*% AL)
    DLr <- In(diag(length(L)) + t(TrL) %*% TrL)
    cl <- trclos(DLr)
    out <- In(AL %*% cl %*% t(AL))
    out <- out[g,g, drop=FALSE]
    adjMatrix(out)
  }

"inducedDAG" <-
function(amat, order, cond=NULL){
### Induced DAG in a new ordering.
    cc=as.list(order)
    inducedChainGraph(amat, cc=cc, cond=cond)
  }

"inducedRegGraph" <-
function(amat, sel=vertices(amat), cond=NULL){
### Induced regression graph for a set of nodes given a conditioning set.
    ancGraph <- function(A) {
      ## Edge matrix of the overall ancestor graph.
      if(sum(dim(A)) == 0)
        return(A)
      else
        return(In(solve(2*diag(nrow(A)) - A)))
    }
    trclos <- function(M) {
      ## Transitive closure of an UG with edge matrix M. See Wermuth and Cox (2004). 
     edgeMatrix(transClos(adjMatrix(M)))
    }
    A <- edgeMatrix(amat) # From the adjacency matrix to edge matrix
    nod <- vertices(A)
    if(!all(cond %in% nod))
      stop("Conditioning nodes are not among the vertices.")
    if(!all(sel %in% nod))
      stop("Selected nodes are not among the vertices.")
    if(length(intersect(sel,cond) > 0))
      stop("The sets are not disjoint!")
    l <- setdiff(nod, union(sel, cond))  # Nodes marginalized over
    g <- sel
    r <- cond
    L <- union(l,g)
    R <- union(g,r)
    
    AL <-  ancGraph(A[L,L,drop=FALSE])          # A^{LL} 
    TrL <- In(A[r,L,drop=FALSE] %*% AL)         # T_{rL}
    DrL <- In(diag(length(r)) + TrL %*% t(TrL)) # D_{rr-L}
    Arr.L <-  In(A[r,r, drop=FALSE] + A[r,L, drop=FALSE] %*% AL 
                 %*% A[L,r, drop=FALSE])        # A_{rr.L}
    FLr <- In(AL %*% A[L, r, drop=FALSE])       # F_{Lr} 
    out <- In(AL %*% t(TrL) %*% trclos(DrL) %*% Arr.L + FLr)
    t(out[g,r, drop=FALSE])
  }

"isAcyclic" <-
function (amat) 
{
### Tests if the graph is acyclic.
  B <- transClos(amat)
  l <- B[lower.tri(B)]
  u <- t(B)[lower.tri(t(B))]
  com <- (l&u)
  all(!com)
}

"isGident" <-
function(amat){
### Is the UG with adjacency matrix amat G-identifiable?
    is.odd <- function(x) (x %% 2) == 1
    cmpGraph <- function(amat){
      ## Adjacency matrix of the complementary graph.
      g <- 0+ (!amat)
      diag(g) <- 0
      g
    }
    cgr <- cmpGraph(amat) 
    cc <- conComp(cgr) 
    l <- unique(cc)
    k <- length(l)
    g <- rep(k, 0)
    for(i in 1:k){
      subg <- cgr[cc==i, cc==i, drop=FALSE]
      m <- cycleMatrix(subg)
      if(is.null(m))
        rt <- 0
      else        
        rt <- apply(m, 1, sum)
      g[i] <- any(is.odd(rt))
    }
    all(g)
  }

"makeAG" <-
function(dag = NULL, ug = NULL, bg = NULL) {
### Builds an ancestral graph from the adj.mat of: dag, ug, bg
    dag.nodes <- vertices(dag)
    ug.nodes <- vertices(ug) # ; print(ug.nodes)
    bg.nodes <- vertices(bg) # ; print(bg.nodes)
    if(all((!is.null(ug)) & (!is.null(bg)))){
      if(length(intersect(ug.nodes, bg.nodes)) > 0)
        warning("Undirected edges meeting an arrowhead.")
    }
    anteriorGraph <- function (amat) 
      {
        ## Adjacency matrix of the anterior graph from an AG.
        A <- 0 + (amat == 1) # Select the directed and undirected edges
        transClos(A)
      }
    vertices <- unique(c(bg.nodes, ug.nodes, dag.nodes))
    d <- length(vertices)
    amat <- matrix(0, d, d)
    dimnames(amat) <- list(vertices, vertices)
    amat.dag <- amat
    amat.ug <- amat
    amat.bg <- amat
    if(!is.null(dag))
      amat.dag[dag.nodes, dag.nodes] <- dag # directed edges
    if(!is.null(ug)){
      if(any(amat.dag[ug.nodes, ug.nodes] & ug))
        warning("Overlapping directed and undirected edges.")
    }
    if(!is.null(bg)){
      if(any(amat.dag[bg.nodes, bg.nodes] & bg))
        warning("Overlapping directed and bi-directed edges.")
    }
    if(!is.null(ug))
      amat.ug[ug.nodes, ug.nodes] <- ug      # undirected edges
    if(!is.null(bg))
      amat.bg[bg.nodes, bg.nodes] <- bg*2    # bi-directed edges
    if(any(apply(amat.dag, 2, sum) & apply(amat.ug, 2, sum)))
      warning("Undirected edges meeting an arrowhead...")
    amat <- amat.dag + amat.ug + amat.bg
    H <- anteriorGraph(amat)[bg.nodes, bg.nodes]
    K <- amat[bg.nodes, bg.nodes]
    if(any(H & (K == 2)))
      warning("Spouses cannot be ancestors.")
    amat
  }

"pa" <-
function (nn, amat) 
{
### List of the parents of nodes nn for a given with adjacency matrix amat.
  nod <- rownames(amat)
  if(is.null(nod)) stop("The adjacency matrix must have dimnames!")
  if(!all(is.element(nn, nod))) stop("Some of the nodes are not among the vertices.")
  k <- length(nn)
  p <- vector(k, mode="list")
  A <- 0 + ((amat != t(amat)) & (amat == 1)) # Select the directed edges
  for(i in 1:k) {
    p[[i]] <- nod[A[,nn[i]]==1 ]
  }
  setdiff(unique(unlist(p)), nn)
}

"parcor" <-
function (S)
{
### Finds the partial correlation matrix of the variables given the rest.
### S is the covariance matrix.  
  p <- ncol(S)
  K <- solve(S)
  a <- 1/sqrt(diag(K))
  K <- K * outer(a, a)
  out <- 2 * diag(p) - K
  dimnames(out) <- dimnames(S)
  out
}

"pcor" <-
function (u, S) 
{
### Partial correlation between u[1:2], given th rest of u. S: cov matrix.
  k <- solve(S[u,u])
  -k[1,2]/sqrt(k[1,1]*k[2,2])
}

"pcor.test" <-
function(r, q, n){
                df = n - 2 - q
                tval <- r * sqrt(df)/sqrt(1-r*r)
                pv <- 2 * pt(-abs(tval), df)
  list(tval = tval, df = df, pvalue = pv)

}

"rcorr" <-
function(d)
{
# Generates a random correlation matrix of dimension d
# with the method of Marsaglia and Olkin (1984).
 h<-rsphere(d,d)
 h %*% t(h)
}

"rnormDag" <-
function (n, A, Delta) 
{
### Generates n observations from a multivariate normal with mean 0
### and a covariance matrix A^-1 Delta (A^-1)'.
  p <- length(Delta)
  E <- matrix(0, n, p)
  for(j in 1:p) { 
    E[,j] <- rnorm(n, 0, sqrt(Delta[j]))
  }
  B <- solve(A)
  Y <- E %*% t(B) 
  colnames(Y) <- colnames(A)
  Y
}

"rsphere" <-
function(n, d)
{
## Generates n random vectors uniformly dist. on the
## surface of a sphere, in d dimensions.
  X <- matrix(rnorm(n*d),n,d)
  d <- apply(X, 1, function(x) sqrt(sum(x*x)))
  sweep(X, 1, d, "/")
}

"shipley.test" <-
function (amat, S, n) 
{
### Overall d-separation test. See Shipley (2000).
### amat: adjacency matrix; S: covariance matrix;  n: observations.
  pval <- function(r, q, n){
    ## See pcor
    df = n - 2 - q
    tval <- r * sqrt(df)/sqrt(1-r*r)
    2 * pt(-abs(tval), df)
  }
  l <- basiSet(amat)
  k <- length(l)
  p <- rep(0, k)
  for(i in 1:k){
    r <- pcor(l[[i]], S)
    q <- length(l[[i]]) - 2
    p[i] <- pval(r, q, n)
  }
  ctest <- -2 * sum(log(p))
  df <- 2*k
  pv <- 1 - pchisq(ctest, df)
  list(ctest=ctest, df=df, pvalue=pv)
}

"swp" <-
function (V, b) 
{
### SWP operator. V is the covariance matrix, b  is a  subset of indices.
  p <- ncol(V)
  u <- is.na(match(1:p, b))
  a <- (1:p)[u]
  out <- 0 * V
  dimnames(out) <- dimnames(V)
  if (length(a) == 0) 
    return(-solve(V))
  else if (length(a) == p) 
    return(V)
  else{
    Saa <- V[a, a, drop = FALSE]
    Sab <- V[a, b, drop = FALSE]
    Sbb <- V[b, b, drop = FALSE]
    B <- Sab %*% solve(Sbb)
    out[a, a] <- Saa - B %*% t(Sab)
    out[a, b] <- B
    out[b, a] <- t(B)
    out[b, b] <- -solve(Sbb)
    return(out)
  }
  ## list(swept = out, coef = out[a, b], rss = out[a, a, drop = F])
}

"topOrder" <-
function (amat) 
{
### Return the nodes in topological order (parents before children).
### Translated from: Kevin Murphy's BNT.
  if(!isAcyclic(amat)) stop("The graph is not acyclic!")
  n <- nrow(amat)
  nod <- 1:n
  indeg <- rep(0, n)
  up <- !amat[lower.tri(amat)]
  if(all(up))
    return(nod)
  zero.indeg <- c() #  a stack of nodes with no parents
  for(i in nod) {
    indeg[i] <- sum(amat[,i])
    if(indeg[i] == 0)
      zero.indeg <- c(i,  zero.indeg)
  }
  s <- 1
  ord <- rep(0, n)
  while(length(zero.indeg) > 0){
    v <- zero.indeg[1]  #  pop v
    zero.indeg <- zero.indeg[-1]
    ord[s] <- v
    s <- s + 1
    cs <- nod[amat[v,]==1]
    if(length(cs) == 0) next
    for(j in 1:length(cs)){
      k <- cs[j]
      indeg[k] <- indeg[k] - 1
      if(indeg[k] == 0)
        zero.indeg <- c(k,  zero.indeg) # push k    
    }
  }
  ord
}

"topSort" <-
function (amat) {
### Topological sort of the DAG with adjacency matrix amat.
    ord <- topOrder(amat)
    amat[ord, ord]
}

"transClos" <-
function (amat) 
{
### Transitive closure of the relation with adjacency matrix amat.
  if (nrow(amat) == 1) 
    return(amat)
  A <- amat
  diag(A) <- 1
  repeat {
    B <- sign(A %*% A)
    if (all(B == A))
      break
    else A <- B
  }
  diag(A) <- 0
  A
}

"triDec" <-
function(Sigma){
### Triangular decomposition of covariance matrix Sigma.  
  R = chol(solve(Sigma))
  dimnames(R) = dimnames(Sigma)
  D = diag(R)
  A = diag(1/D) %*% R
  dimnames(A) <- dimnames(Sigma)
  B = solve(A) 
  list(A = A, B = B, Delta = 1/(D^2))
}

"UG" <-
function (f) 
{
### Defines an UG from a model formula. Returns the adj. matrix.  
  tt <- terms(f)
  if (attr(tt, "response") == 1)
    stop("You should not specify a response!")
  nod <- dimnames(attr(tt, "factors"))[[1]]
  
  N <- unique(nod) # set of nodes
  dN <- length(N)  # number of nodes
  amat <- matrix(0, dN, dN)
  o <- attr(tt, "order") <= 2
  v <- attr(tt, "factors")[, o, drop = FALSE]
  m <- match(dimnames(v)[[1]], N)
  for (i in 1:sum(o)) {
    ij <- m[v[, i] == 1]
    amat[ij[1], ij[2]] <- 1
    amat[ij[2], ij[1]] <- 1
  }
  dimnames(amat) <- list(N, N)
  amat
}

"unmakeAG" <-
function(mag){
    ### Returns a list with the three components of an Ancestral Graph.
    temp <- mag+t(mag)
    ug <- dag <- bg <- matrix(0, nrow(mag), ncol(mag))
    dimnames(ug) <- dimnames(dag) <- dimnames(bg) <- dimnames(mag)
    ug[temp==2] <- 1
    dag[(temp==1) & (mag==1)] <- 1
    bg[temp==4] <- 1
    return(list(dag = dag, ug = ug, bg = bg))
  }

"vertices" <-
function(amat) rownames(amat)

