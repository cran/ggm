"ancGraph" <-
function(A) {
### Edge matrix of the overall ancestor graph.
    In(solve(2*diag(nrow(A)) - A))
  }
"anGraph" <-
function (A) 
{
### Edge matrix of the ancestor graph from edge matrix A.
  if (length(A) == 1) 
    return(1)
  repeat {
    B <- sign(A %*% A)
    if (all(B == A)) 
      return(B)
    else A <- B
  }
}
"basiSet" <-
function(A){
### Basis set of a DAG. A is the edge matrix.
    A <- topSort(A)$triu
    nam <- rownames(A)
    dv <- ncol(A)
    ind <- NULL
    for(r in 1:dv){
      for(s in r:dv) {
        if(A[r,s] != 0)
          next
        else{
          ed <- nam[c(r,s)]
          pa.r <- nam[A[r, ] == 1]
          pa.s <- nam[A[s, ] == 1] 
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
function (nn, A) 
{
### Boundary of the nodes nn for a graph with edge matrix A.
  if(!is.character(nn)) # Names of the nodes
    nod <- 1:nrow(A)
  else {
    nod <- rownames(A)
    if(is.null(nod)) stop("The edge matrix must have dimnames!")
  }
  b <- vector(length(nn), mode="list")
  diag(A) <- 0  # As you do not want the node itself in the list
  k <- length(nn)
  for(i in 1:k) {
    b[[i]] <- c(  nod[A[nn[i], ]==1 ],nod[A[,nn[i]]==1 ] )
  }
  b <- unique(unlist(b))
  setdiff(b, nn)
}
"bfs" <-
function(gmat, v=1) {
### Breadth-first search of a connected UG with edge matrix gmat.
    n <- nrow(gmat)
    if(n==1) return(NULL)
    visited <- rep(0, n)
    Q <- c()
    tree <- diag(n)
    dimnames(tree) <- dimnames(gmat)
    E <- c()
    visited[v] <- 1
    Q <- c(v, Q)
    while(!all(visited)){
      x <- Q[1]
      Q <- Q[-1]      
      b <- bd(x, gmat)
      for(y in b){
        if(visited[y] == 0){
          visited[y] <- 1
          Q <- c(Q, y)
          tree[x,y]<- 1 ; tree[y,x] <- 1
          E <- rbind(E, c(x,y))
        }
      }
    }
    cross <- gmat - tree
    diag(cross) <- 1
    V = edges(cross)
    dimnames(E) <- list(rep("", nrow(E)), rep("", 2))
    list(tree = tree, branches = E,chords = V ) 
  }
"ch" <-
function (nn, A) 
{
### Children of nodes nn for a given with edge matrix A.
  if(!is.character(nn)) # Names of the nodes
    nod <- 1:nrow(A)
  else {
    nod <- rownames(A)
    if(is.null(nod)) stop("The edge matrix must have dimnames!")
  }
  p <- vector(length(nn), mode="list")
  diag(A) <- 0    # As you do not want the node itself in the list
  k <- length(nn)
  for(i in 1:k) {
    p[[i]] <- nod[A[,nn[i] ]==1 ]
  }
  setdiff(unique(unlist(p)), nn)
}
"checkIdent" <-
function(gmat, hidden) {
### Checks SW sufficient conditions for identifiability of a DAG
### with edge matrix gmat and one hidden variable L.
    if(!is.character(hidden)) # Names of the nodes
      nod <- 1:nrow(gmat)
    else {
      nod <- rownames(gmat)
      if(is.null(nod)) stop("The edge matrix must have dimnames!")
    }
    gcov <- inducedCovGraph(gmat, sel=1:nrow(gmat), cond=NULL)
    L <- hidden
    if(length(L) > 1)
      stop("Only one hidden variable, please!")
    O <- setdiff(nod, hidden)
    m <-  bd(L, gcov)
    
    ## Theorem 1
    if(length(m) > 2){ 
      G <- inducedCovGraph(gmat, sel=O, cond=L)
      cond.i <- is.Gident(G[m,m,drop=FALSE])
    }
    else
      cond.i <- FALSE
    gcon <- inducedConGraph(gmat, sel=1:nrow(gmat), cond=NULL) 
    cc <- bd(L, gcon)
    if(length(cc) >2) { 
      cond.ii <- is.Gident(gcon[cc, cc, drop=FALSE])
    } 
    else  
      cond.ii <- FALSE
    ## Theorem 2
    a <- union(pa(L, gmat), ch(L, gmat))
    ma <- setdiff(m,a)
    ca <- setdiff(cc,a)
    if( 0 != length(intersect(ma,ca))) {
      cond.iii <- FALSE
      cond.iv <- FALSE
    }
    else {
      if(length(a) > 2){
        c1 <- is.Gident(inducedCovGraph(gmat, sel=a, cond=union(L, ma)))
        c2 <- is.Gident(inducedConGraph(gmat, sel=a, cond=union(L, ma)))
        cond.iii <- c1 | c2
        
        c3 <- is.Gident(inducedCovGraph(gmat, sel=a, cond=union(L, ca)))
        c4 <- is.Gident(inducedConGraph(gmat, sel=a, cond=union(L, ca)))
        cond.iv <- c3 | c4      
      }
      else{
        cond.iii <- FALSE
        cond.iv <- FALSE
      }
    }
    c(cond.i = cond.i, cond.ii = cond.ii,
      cond.iii=cond.iii, cond.iv=cond.iv)  
  }
"cliques" <-
function (gmat) 
{
### Finds the cliques of an UG.
  is.complete <- function (model, cli, j) 
    {
      isOK <- TRUE
      for (i in cli){
        if (!model[i, j]) 
          isOK <- FALSE
      }
      isOK
    }
  addEdge <- function (gmat, cliqs, cliq, n) 
    {
      maximal <- TRUE
      if (max(cliq) < n) 
        for (j in (max(cliq) + 1):n) if (is.complete(gmat,cliq, j)) {
          maximal <- FALSE
          cliqs <- addEdge(gmat, cliqs, append(cliq, j), n)
        }
      if (maximal) {
        isOK <- TRUE 
        if (length(cliq) != n) 
          for (j in (1:n)[-cliq]) if (is.complete(gmat, cliq, j)) 
            isOK <- FALSE 
        if (isOK) {
          cliqs <- append(cliqs, 0)
          cliqs[[length(cliqs)]] <- cliq
        }   
      }
      cliqs
    } 
  cliqs <- list()
  n <- nrow(gmat)
  for (i in 1:n) {
    cliqs <- addEdge(gmat, cliqs, i, n)
  }
  nc <- length(cliqs)
  nod <- rownames(gmat)
  if(is.null(length(nod)))
    stop("The nodes should have a name.")
  for(i in 1:length(cliqs)){
    cliqs[[i]] <- nod[cliqs[[i]]]
  }
  cliqs
}
"clos" <-
function(M) {
### See Wermuth and Cox (2003).
    triu <- function(A) {
      ## Upper tri of A including the diagonal.
      A[lower.tri(A)] <- 0
      A
    }
    M1 <- In(triu(M))
    M2 <- In(crossprod(M1))
    p <- nrow(M)
    M3 <- ancGraph(triu(M2)) #  In(solve(2*diag(p) - triu(M2)))
    In(M3 %*% t(M3))
  }
"cmpGraph" <-
function(gmat){
### Edge matrix of the complementary graph
    g <- 1*!gmat
    diag(g) <- 1
    g
  }
"conComp" <-
function (gmat) 
{
### Finds the connected components of an UG graph from its edge matrix gmat.
  if(!all(gmat == t(gmat)))
    gmat <- 0 + (gmat | t(gmat))
  gmat <- anGraph(gmat)
  n <- nrow(gmat)
  gmat <- sign(gmat + t(gmat))
  u <- gmat %*% 2^((n - 1):0)
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
function(gmat){
### Fundamental Cycle matrix of the UG gmat.
    fc <- fundCycles(gmat)  # Edges of the fundamental cycles
    E <- edges(gmat)        # All the edges of the graph
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
function (...,order=FALSE, test=TRUE) 
{
### Defines a DAG from a set of equations(defined with model formulae).
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
  gmat <- diag(dN)
  for (k in 1:nb) {
    tt <- terms(f[[k]])
    vars <- dimnames(attr(tt, "factors"))[[1]]
    if (attr(tt, "response") == 1) {
      i <- match(vars[1], N)
      j <- match(vars[-1], N)
      gmat[i, j] <- 1
    }
    else if (attr(tt, "response") == 0) 
      stop("Error! Some equations have no response")
  }
  if(test) {
    if(!is.acyclic(gmat))
      stop("The graph contains directed cycles!")
  }
  dimnames(gmat) <- list(N, N)
  if(order){
    gmat <- topSort(gmat)$triu
  }
  gmat
}
"dSep" <-
function(gmat, first, second, cond) {
### Are first and second d-Separated by cond in a DAG with edge matrix gmat?. 
    e <- inducedCovGraph(gmat, sel=c(first,second), cond=cond)
    all(e[first,second] == 0)
  }
"edges" <-
function(gmat){
### Finds all the edges of a graph with edge matrix gmat.
    nn <- 1:nrow(gmat)
    E <- c()
    diag(gmat) <- 0
    if(all(gmat == t(gmat))) { 
      gmat[lower.tri(gmat)] <- 0
    }
    for(i in nn) {
      e <- nn[gmat[i,] == 1]
      if(length(e) == 0) next
      li <- cbind(i,  e)
      dimnames(li) <- list(rep("", length(e)), rep("", 2))
      E <- rbind(E, li) 
    }
    E
  }
"findPath" <-
function (gmat, st, en, path = c()) 
{
### Find a path between nodes st and en in a graph gmat.
  if(st == en) # st is 'node' in recursive calls
    return(c(path, st))
  if(sum(gmat[st,]) == 1 ) 
    return(NULL)
  ne <- bd(st,gmat)
  
  for(node in ne){
    if(!is.element(node, c(path, st))){
      newpath <- findPath(gmat, node, en, c(path, st))
      if(!is.null(newpath))
        return(newpath)
    }
  }
}
"fitDag" <-
function (gmat, Syy,n)
{
### Fits linear recursive regressions with independent residuals. 
### gmat: the edge matrix of the DAG. Syy: cov matrix. n: sample size.
  if(missing(gmat)) # saturated model
    gmat <-  upper.tri(diag(ncol(Syy)), diag=TRUE) * 1
  nam <- rownames(Syy)
  nod <- rownames(gmat)
  if(is.null(nod))
    stop("Edge matrices must have rownames!")
  if(!all(is.element(nod, nam)))
    stop("The nodes of the graph do not match the names of the variables")
  else
    sek <- intersect(nam, nod) # What if nod is not a subset of nam? FIXME
  Syy <- Syy[sek,sek]   # Resizes eventually Syy 
  gmat <- gmat[sek,sek] # and reorders gmat
  Delta <- rep(length(sek),0)
  A <- gmat
  p <- ncol(Syy)
  ip <- 1:p
  for(i in 1:p) {
    u <- gmat[i,]
    v <- ip[u == 1 & ip != i]
    M <- swp(Syy, v)[i,]
    A[i, ] <- - A[i, ] * M
    A[i,i] <- 1
    k <- sum(A[i,])
    Delta[i] <- M[i]
  }
  names(Delta) <- sek
  B <- solve(A)
  Shat <- B %*% diag(Delta) %*% t(B)
  dimnames(Shat) <- dimnames(Syy)
  Khat <- solve(Shat)
  H <- Syy %*% Khat
  Trace <- function(A) sum(diag(A))
  l2 <- (Trace(H) - log(det(H)) - p) * n
  gmat <- topSort(gmat)$triu
  df <- sum(!gmat[upper.tri(gmat)])
  list(A = A, B=B, Delta=Delta, Shat=Shat, Khat=Khat, dev=l2, df=df)
}
"fitDagLatent" <-
function (gmat, Syy, n, latent, norm = 1,  seed=144, maxit=9000, tol=1e-6, pri=FALSE) 
{
### Fits linear recursive regressions with independent residuals and one latent variable.
### Syy: covariance matrix, n: sample size, gmat: edge matrix.
### NOTE: both gmat and Syy must have rownames.
### latent is the "name" of the latent in the rownames of Syy. norm = normalisation type.

  ## Local functions
  setvar1 <- function (V, z, norm) 
    {
      ## Normalizes V
      if(norm == 1){
        ## Rescales V forcing V[z,z] = 1
        a <- 1 / sqrt(V[z,z])
      }
      else if(norm== 2) {
        ## Rescales V forcing Delta[z,z] = 1
        a <- 1/sqrt((triDec(V)$Delta)[z]) # Note: this needs an upper trianguar edge matrix
      }
      V[z,] <- V[z,] * a
      V[,z] <- V[,z] * a
      V
    }

  cmqi <- function (Syy, Sigma, z) 
    {
      ## Computes the matrix C(M | Q) by Kiiveri (1987), Psychometrika.
      ## It is a slight generalization in which Z is not the last element.
      ## z is a Boolean vector indicating the position of the latent
      ## variable in X.
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
  fitdag <- function (gmat, Syy,n, constr=NULL)
    {
      ## Fits linear recursive regressions with independent residuals (fast version).
      ## NOTE. gmat and Syy must have the same size and variables. constr is a matrix
      ## indicating the edges that must be constrained to 1.
      A <- gmat
      p <- ncol(Syy)
      Delta <- rep(p,0)
      ip <- 1:p
      for(i in 1:p) {
        u <- gmat[i,]
        v <- ip[u == 1 & ip != i] # Parents without response
        if(length(v) == 0){  # if pa is empty
          Delta[i] <- Syy[i,i]
          next
        }
        M <- lmfit(Syy, y=i, x=v, z=(constr[i, ] == 1))
        A[i, ] <- - A[i, ] * M
        A[i,i] <- 1
        k <- sum(A[i,])
        Delta[i] <- M[i]
      }
      B <- solve(A)
      Shat <- B %*% diag(Delta) %*% t(B)
      Khat <- solve(Shat)
      H <- Syy %*% Khat
      Trace <- function(A) sum(diag(A))
      l2 <- (Trace(H) - log(det(H)) - p) * n
      list(A = A, B=B, Delta=Delta, Shat=Shat, Khat=Khat, dev=l2)
    }
  lmfit <- function(V, y, x, z=NULL){
    ## Regression coefficients of y given x eventually with z constrained to 1.
    ## Residual variance in position y.
    Vxy <- V[x, y, drop=FALSE] - apply(V[x, z, drop=FALSE], 1, sum) 
    Vxx <- V[x, x, drop=FALSE]
    bxy <- solve(Vxx,Vxy)
    out<- rep(0, nrow(V))
    out[x] <- bxy
    out[z] <- 1
    names(out) <- rownames(V)
    xz <- c(x,z)
    b <- out[xz, drop=FALSE]
    res <- V[y,y] + t(b) %*% (V[xz, xz, drop=FALSE] %*% b - 2 * V[xz,y, drop=FALSE])
    out[y] <- res
    out
  }
 
  o <- topSort(gmat)$perm # It is safer to have an upper triangular edge matrix.
  gmat <- gmat[o,o]
  nam <- rownames(Syy)    # Names of the variables (they can be more than the nodes)
  nod <- rownames(gmat)   # Names of the nodes of the DAG (that contains the latent)

  if(is.null(nod))
    stop("Edge matrices must have rownames!")
  if(!setequal(latent, setdiff(nod, nam)))
    stop("The nodes of the graph do not match the names of the variables")
  else
    sek <- intersect(nod, nam) 
  Syy <- Syy[sek,sek]   # Reorders and eventually resizes Syy. Leaves untouched gmat. 
  a <- rownames(gmat)                          # Ordering of the dag
  dn <- list(a,a)
  wherez <- is.element(a, latent)              # Index of the latent
  if(is.null(wherez))
    stop("Wrong name of the latent variable!")
  wherey <- ! wherez
                                       
  df <- sum(!gmat[upper.tri(gmat)])  - sum(wherey)  # Degrees of freedom
  if(df <= 0)
    warning(paste("The degrees of freedom are ", df))
  p <- ncol(Syy)
  dev.old <- Inf
  converge <- FALSE
                                        # Starting values (random)
  set.seed(seed)
  Sigma.old <- rcorr(p+1)
  Sigma.old <- setvar1(Sigma.old, wherez, norm=norm)
  
  dimnames(Sigma.old) <- dn

  it <- 0
  while(!converge){
    it <- it+1
    if(it > maxit){
      warning("Maximum number of iterations reached!")
      break
    }
                                        # E-step
    Q <- cmqi(Syy, Sigma.old, wherez) # See Kiiveri, 1987 
    
                                        # M-step
    fit <- fitdag(gmat, Q, n)
    
    Sigma.new <- fit$Shat
    dev.new <- fit$dev
    
                                        # Monitoring progress of iterations
    if(pri) cat(dev.new, "\n")
    else{
      if(0==(it %% 80))
        cat("\n")
      else
        cat(".")
    }
                                        # Test convergence
    converge <- (abs(dev.old - dev.new) < tol)
    Sigma.old <- Sigma.new
    dev.old <- dev.new
  }
  cat("\n")
  Shat <- setvar1(Sigma.new, wherez, norm=norm) # Normalize Shat
  dimnames(Shat) <- dn
  Khat <- solve(Shat)
  dimnames(Shat) <- dn
  fit <- fitDag(gmat, Shat, n)
  Delta <- fit$Delta
  names(Delta) <- a
  list(Shat=Shat, Khat=Khat, A=fit$A, Delta=Delta, dev=dev.new,it=it, df=df) 
}
"fitUg" <-
function (gmat, V, n, cov=FALSE, pri=FALSE, tol = 1e-06)
{
### Fits a concentration graph or a covariance graph.  G. Marchetti, 1994, 2002, 2003.
### gmat is the edge matrix; cov = TRUE fits a covariance graph.
  nam <- rownames(V)
  nod <- rownames(gmat)
  if(!all(is.element(nod, nam)))
    stop("The nodes of the graph do not match the names of the variables")
  else
    sek <- intersect(nam, nod)
  V <- V[sek,sek]              # Resizes eventually Syy
  gmat <- gmat[sek,sek]        # and reorders gmat
                                        # Finds the cliques
  cli <- cliques(gmat)
                                        # Inverts matrix V when cov = TRUE 
  type <- "concentration graph"
  if(cov){
    V <- solve(V) # NOTE !
    type <- "covariance graph"
  }
                                        # Finds the indicators of the cliques
  maxnodes <- unique(unlist(cli))       
  l <- lapply(cli , function(x) maxnodes %in% x)
  cliq <- matrix(unlist(l), nrow=length(l), byrow=TRUE)
  trace <- function(a) sum(diag(a))
  lold <- +Inf
  converge <- FALSE
  nc <- nrow(cliq)
  k <- ncol(cliq)
                                        # Starting value
  Vg <- diag(k)
  while (!converge) {
    for (i in 1:nc) {
      a <- cliq[i, ]
      b <- !a
      Vaa <- V[a, a]
      Vgaa <- Vg[a, a]
      Vgba <- Vg[b, a]
      Vgbb <- Vg[b, b]
      B <- Vgba %*% solve(Vgaa)
      Vpar <- Vgbb - B %*% Vgaa %*% t(B)
      BV <- B %*% Vaa
      Vg[b, a] <- BV
      Vg[a, b] <- t(BV)
      Vg[a, a] <- Vaa
      Vg[b, b] <- Vpar + BV %*% t(B)
      sd <- V %*% solve(Vg)
      l2 <- (trace(sd) - log(det(sd)) - k) * n
      if (tol > abs(l2 - lold)) {
        converge <- TRUE
        break
      }
      else {
        if(pri) cat(l2, "\n")
        lold <- l2
      }
    }
  }
                                        # Degrees of freedom
  df <-  !crossprod(cliq)
  df <- sum(df[row(df) > col(df)])
  if(cov) {
    Vg <- solve(Vg)
    V <- solve(V)
    sd <- V %*% solve(Vg)
    l2 <- (trace(sd) - log(det(sd)) - k) * n	
  }
  dimnames(Vg) <- dimnames(V)
  sigma <- sqrt(diag(Vg))
  Rhat <-  Vg / outer(sigma,sigma)
  Phat <- parcor(Vg)
  out <- list(type, cliques=cli, Vhat = Vg, Rhat = Rhat,
              Phat = Phat, dev = l2, df = df)
  out
}
"fundCycles" <-
function(gmat){
### Finds a set of fundamental cycles for an UG gmat.
    fc <- c()
    tr <- bfs(gmat) # Spanning tree
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
"In" <-
function (A) 
{
### Indicator matrix of structural zeros.
   (A != 0) + 0 
}
"inducedConGraph" <-
function(A, sel=1:nrow(A), cond=NULL){
### Edge matrix of a con graph for sel|cond deduced by a DAG.
    nod <- rownames(A)
    if(is.null(nod))
      stop("The edge matrix must have dimnames!") 
    if(!is.character(sel))
      sel <- nod[sel]
    if(!is.character(cond))
      cond <- nod[cond]
    if(!all(cond %in% nod))
      stop("Wrong names!")
    if(!all(sel %in% nod))
      stop("Wrong names!")
    if(length(intersect(sel,cond) > 0))
      stop("The sets are not disjoint!")
    l <- setdiff(nod, union(sel, cond))  # Marginal nodes
    g <- sel
    r <- cond
    L <- union(l,g)
    R <- union(g,r)
    zeroLen <- function(a) length(a) == 0
    ## Compute B
    p <- nrow(A)
    B <- ancGraph(A)
    ## Case cond is empty
    if(zeroLen(l)  & zeroLen(cond)){
      out <- In(crossprod(A))
    }
    else if((length(l)>0) & zeroLen(cond)){
      Al <- ancGraph(A[l,l,drop=FALSE])
      All.g <- In(A[g,g,drop=FALSE] +
                  A[g,l,drop=FALSE] %*% Al %*% A[l,g,drop=FALSE])
      Tgl <- In(A[g,l,drop=FALSE] %*% Al)
      Dgl <- In(diag(length(g)) + Tgl %*% t(Tgl))
      out <- In(t(All.g) %*% clos(Dgl) %*% All.g)
    }
    ## Case cond not NULL
    else if(zeroLen(l) & (! zeroLen(cond))) { # if l is empty
      out <- In(crossprod(A[g,g,drop=FALSE]) + crossprod(A[r,g,drop=FALSE]))
    }
    else { # if l is not empty
      Al <- ancGraph(A[l,l,drop=FALSE])
      ARR.l <- In(A[R,R, drop=FALSE] +
                  A[R,l, drop=FALSE]%*% Al %*% A[l,R, drop=FALSE])
      TRl <- In(A[R,l,drop=FALSE] %*% Al)
      DRl <- In(diag(length(R)) + TRl %*% t(TRl))
      out <- In(t(ARR.l) %*% clos(DRl) %*% ARR.l)
      out <- out[g,g, drop=FALSE]
    }
    out
  }
"inducedCovGraph" <-
function(A, sel=1:nrow(A), cond=NULL){
### Induced edge matrix of a cov graph for sel|cond from a DAG.
     
    nod <- rownames(A)
    if(is.null(nod))
      stop("The edge matrix must have dimnames!") 
    if(!is.character(sel))
      sel <- nod[sel]
    if(!is.character(cond))
      cond <- nod[cond]
    if(!all(cond %in% nod))
      stop("Wrong names!")
    if(!all(sel %in% nod))
      stop("Wrong names!")
    if(length(intersect(sel,cond) > 0))
      stop("The sets are not disjoint!")
    l <- setdiff(nod, union(sel, cond))  # Marginal nodes
    g <- sel
    r <- cond
    L <- union(l,g)
    R <- union(g,r); 
    ## Compute B
    B <- ancGraph(A) 
    ## Case cond is empty
    if(length(cond)==0){
      Bgl <- B[g,l, drop=FALSE]
      Bgg <- B[g,g, drop=FALSE]
      out <- In(Bgg %*% t(Bgg) + Bgl %*% t(Bgl))
    }
    ## Case cond not NULL
    else{
      AL <-  ancGraph(A[L,L,drop=FALSE]) # In(solve(2*diag(length(L)) - A[L,L]))
      TrL <- In(A[r,L,drop=FALSE] %*% AL)
      DLr <- In(diag(length(L)) + crossprod(TrL))
      out <- In(AL %*% clos(DLr) %*% t(AL))
      out <- out[g,g, drop=FALSE]
    }
    out
  }
"is.acyclic" <-
function (A) 
{
### Tests if the graph is acyclic.
  B <- anGraph(A)
  l <- B[lower.tri(B)]
  u <- t(B)[lower.tri(t(B))]
  com <- (l&u)
  all(!com)
}
"is.Gident" <-
function(gmat){
### Is the UG with edge matrix gmat G-identifiable?
    is.odd <- function(x) (x %% 2) == 1
    cgr <- cmpGraph(gmat)
    cc <- conComp(cgr)
    l <-unique(cc)
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
"pa" <-
function (nn, A) 
{
### List of the parents of nodes nn for a given with edge matrix A.
  if(!is.character(nn)) # Names of the nodes
    nod <- 1:nrow(A)
  else {
    nod <- rownames(A)
    if(is.null(nod)) stop("The edge matrix must have dimnames!")
  }
  p <- vector(length(nn), mode="list")
  diag(A) <- 0    # As you do not want the node itself in the list
  k <- length(nn)
  for(i in 1:k) {
    p[[i]] <- nod[A[nn[i], ]==1 ]
  }
  setdiff(unique(unlist(p)), nn)
}
"parcor" <-
function (Smat)
{
### Finds the partial correlation matrix of the variables given the rest.
### Smat is the covariance matrix.  
  p <- ncol(Smat)
  K <- solve(Smat)
  a <- 1/sqrt(diag(K))
  K <- K * outer(a, a)
  out <- 2 * diag(p) - K
  dimnames(out) <- dimnames(Smat)
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
function (S, n, A) 
{
### Overall d-separation test. See Shipley (2000).
### S = covariance matrix; A = edge matrix; n = observations.
  pval <- function(r, q, n){
    ## See pcor
    df = n - 2 - q
    tval <- r * sqrt(df)/sqrt(1-r*r)
    2 * pt(-abs(tval), df)
  }
  l <- basiSet(A)
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
"topSort" <-
function (gmat) 
{
### Topological sort of the DAG with edge matrix gmat.
  if(!is.acyclic(gmat)) stop("The graph is not acyclic!")
  A <- gmat
  n <- ncol(A)
  if(all(!A[lower.tri(A)]))
    return(list(triu=A, perm=1:n))
  v <- 1:ncol(A)
  o <- rep(0, n)
  for(i in n:1) { 
    id <- apply(A, 1, sum)
    o[i] <- v[id == 1][1]
    A[o[i],]  <- 0
    A[, o[i]] <- 0
  }
  list(triu = gmat[o,o], perm=o)
}
"triDec" <-
function(Sigma){
### Triangular decomposition of covariance matrix Sigma.  
  R = chol(solve(Sigma))
  dimnames(R) = dimnames(Sigma)
  D = diag(R)
  A = diag(1/D) %*% R
  B = solve(A) 
  list(A = A, B = B, Delta = 1/(D^2))
}
"UG" <-
function (f) 
{
### Defines an UG from a model formula.  
  tt <- terms(f)
  if (attr(tt, "response") == 1)
    stop("You should not specify a response!")
  nod <- dimnames(attr(tt, "factors"))[[1]]
  
  N <- unique(nod) # the set of notes
  dN <- length(N) # the number of nodes
  gmat <- diag(dN)
  o <- attr(tt, "order") <= 2
  v <- attr(tt, "factors")[, o, drop = FALSE]
  m <- match(dimnames(v)[[1]], N)
  for (i in 1:sum(o)) {
    ij <- m[v[, i] == 1]
    gmat[ij[1], ij[2]] <- 1
    gmat[ij[2], ij[1]] <- 1
  }
  dimnames(gmat) <- list(N, N)
  gmat
}
