pmin2 <- function(xvec, yvec){
    #res <- pmin3(xvec, yvec)
    #print(c(length(xvec), length(yvec), length(res)))
    #return(res)
    return(((xvec + yvec) - (abs(xvec - yvec))) / 2)
}

websearch_NODF_fast <- function(mtx){
    k = matrix(0,3,3)
    k[1,2] = 1
    k[2,1] = 1
    k[2,3] = 1
    k[3,2] = 1

    NodesA = nrow(mtx)
    NodesB = ncol(mtx)
    mtx2 <- 0.0*mtx
    mtx2[1:NodesA - 1,] <- mtx2[1:NodesA -1,] + mtx[2:NodesA,]
    mtx2[2:NodesA,] <- mtx2[2:NodesA,] + mtx[1:NodesA-1,]
    mtx2[,1:NodesB-1] <- mtx2[,1:NodesB-1] + mtx[,2:NodesB]
    mtx2[,2:NodesB] <- mtx2[,2:NodesB] + mtx[,1:NodesB-1]
    mtx3 <- matrix(as.numeric((mtx2 <= 4.0) * (mtx >= 0.5)), nrow(mtx), ncol(mtx))
    mtx4 <- 0.0*mtx3
    mtx4[2:NodesA,] <- mtx4[2:NodesA,] + mtx3[1:NodesA-1,]
    mtx4[,2:NodesB] <- mtx4[,2:NodesB] + mtx3[,1:NodesB-1]
    mtx5 = (mtx4 >= 0.5)
    mtx6 = matrix(as.numeric(mtx5 * (mtx == 0.0) * (mtx4 == 2.0)), nrow(mtx), ncol(mtx))
    posList = which(mtx6 == 1, arr.ind = T)
    return(posList)
}

# calculates the combined NODF statistic
# inputs: web = mutualistic network, raw NODF, maximum raw NODF
# output: the combined NODF statistic
comb_nest <- function(web,NODF,max_NODF){
    C <- sum(web)/(ncol(web)*nrow(web))
    S <- sqrt(ncol(web) * nrow(web) )
    out <- NODF / (max_NODF * C * log10(S))
    return(out)
}

# Computes the row marginal totals
marginal_total0 <- function(mtx){
    mt_0 <- rowSums(mtx, na.rm = FALSE, dims = 1)
    return(mt_0)
}

# Computes the column marginal totals
marginal_totalt <- function(mtx){
    mt_t <- colSums(mtx, na.rm = FALSE, dims = 1)
    return(mt_t)
}

# Assembles the list of marginal totals
compute_marginal_totals <- function(mtx){
    mt_0  <- marginal_total0(mtx)
    mt_t  <- marginal_totalt(mtx)
    return(list(mt_0, mt_t))
}

# Computes a list containing both Fill matrices
compute_fill_factors <- function(mtx){
    F0 <- mtx %*% t(mtx)
    Ft <- t(mtx) %*% mtx
    return(list(F0, Ft))
}

# Compute a list containing both degree matrices
compute_deg_mtx <- function (mtx) {
    NodesA <- nrow(mtx)
    NodesB <- ncol(mtx)
    mt_0 <- marginal_total0(mtx)
    mt_t <- marginal_totalt(mtx)
    deg_mtx0 <- matrix(mt_0, nrow=length(mt_0),ncol=length(mt_0),byrow=TRUE)
    deg_mtxt <- matrix(mt_t, nrow=length(mt_t),ncol=length(mt_t),byrow=TRUE)
    return(list(deg_mtx0, deg_mtxt))
}

# Compute a list containing both degree minima matrices
compute_deg_minima <- function(mtx){
    my_ans <- compute_deg_mtx(mtx)
    DM0 <- my_ans[[1]]
    DMt <- my_ans[[2]]
    deg_min0 <- pmin2(DM0, t(DM0))
    deg_mint <- pmin2(DMt, t(DMt))
    return(list(deg_min0, deg_mint))
}

# Compute a list containing both negative delta matrices
compute_neg_deltas <- function(mtx){
    mt_0 <- marginal_total0(mtx)
    mt_t <- marginal_totalt(mtx)
    neg_delta0 = outer(mt_0, mt_0, FUN = ">")
    neg_deltat = outer(mt_t, mt_t, FUN = ">")
    return(list(neg_delta0, neg_deltat))
}

# Compute a list containing both the row and column sum values
compute_sums <- function(mtx){
    my_res <- compute_fill_factors(mtx)
    F0 <- my_res[[1]]
    Ft <- my_res[[2]]
    my_res <- compute_deg_minima(mtx)
    DM0 <- my_res[[1]]
    DMt <- my_res[[2]]
    my_res <- compute_neg_deltas(mtx)
    ND0 <- my_res[[1]]
    NDt <- my_res[[2]]
    n_paris0 = F0[ND0] / (DM0[ND0])
    n_parist = Ft[NDt] / (DMt[NDt])
    sum0 = sum(n_paris0)
    sumt = sum(n_parist)
    return(c(sum0, sumt))
}

# Assembles the list containing all the additional data required
# for fast_nodf computations
init_nodf <- function(mtx){
    MT <- compute_marginal_totals(mtx)
    Fill <- compute_fill_factors(mtx)
    DM <- compute_deg_minima(mtx)
    ND <- compute_neg_deltas(mtx)
    S <- compute_sums(mtx)
    return(list(MT, Fill, DM, ND, S))
}

get_contributions <- function(Fill, ND, DM, idx){
    A1 <- Fill[idx,][ND[idx,]] / (DM[idx,][ND[idx,]])
    A2 <- Fill[,idx][ND[,idx]] / (DM[,idx][ND[,idx]])
    return(sum(A1) + sum(A2))
}

# Efficient way to compute the nodf value of a matrix where the
# link at position pos = list(xpos, ypos) is removed.
nodf_one_link_removed <- function(mtx, pos, support_data){
    # Unpack and update all the support data:
    xpos <- pos[[1]]
    ypos <- pos[[2]]
    MT <- support_data[[1]]
    Fill <- support_data[[2]]
    DM <- support_data[[3]]
    ND <- support_data[[4]]
    S <- support_data[[5]]
    # print("B")
    # Unpack even futher to get a set a variables that we can work with:
    mt_0 <- MT[[1]]
    mt_t <- MT[[2]]
    # print("B1")
    F0 <- Fill[[1]]
    Ft <- Fill[[2]]
    # print("B2")
    DM0 <- DM[[1]]
    DMt <- DM[[2]]
    # print("B3")
    ND0 <- ND[[1]]
    NDt <- ND[[2]]
    # print("B4")
    S0 <- S[[1]]
    St <- S[[2]]
    # print("B5")
    NodesA <- nrow(mtx)
    NodesB <- ncol(mtx)
    my_norm <- 0.5*((NodesA *(NodesA -1)) + (NodesB*(NodesB-1)))
    # print("B6")
    # Modify the matrix appropriately (remove a link):
    put_val(mtx, xpos, ypos, 0.0)
    # print("C")

    # Compute old contributions
    old_contrib_0 <- get_contributions(F0, ND0, DM0, xpos)
    old_contrib_t <- get_contributions(Ft, NDt, DMt, ypos)

    # modify the marginal totals
    mt_0[xpos] <- mt_0[xpos] - 1
    mt_t[ypos] <- mt_t[ypos] - 1

    # Update the degree matrix:
    m0 <- rep(mt_0[xpos], length(mt_0))
    mt <- rep(mt_t[ypos], length(mt_t))

    # Update the degree minima:
    DM0[xpos,] <- pmin3(m0, mt_0)
    DM0[,xpos] <- pmin3(m0, mt_0)
    DMt[ypos,] <- pmin3(mt, mt_t)
    DMt[,ypos] <- pmin3(mt, mt_t)

    # Update negative deltas:
    ND0[xpos, ] <- (m0 > mt_0)
    ND0[, xpos] <- (m0 < mt_0)
    NDt[ypos, ] <- (mt > mt_t)
    NDt[, ypos] <- (mt < mt_t)

    # Update the fill factors
    F0[,xpos] <- F0[,xpos] - mtx[,ypos]
    F0[xpos,] <- F0[xpos,] - t(mtx[,ypos])
    F0[xpos,xpos] <- F0[xpos,xpos] - 1
    #
    Ft[, ypos] <- Ft[,ypos] - mtx[xpos,]
    Ft[ypos,] <- Ft[ypos,] - mtx[xpos,]
    Ft[ypos,ypos] <- Ft[ypos,ypos] - 1

    # Compute the new contributions:
    new_contrib_0 <- get_contributions(F0, ND0, DM0, xpos)
    new_contrib_t <- get_contributions(Ft, NDt, DMt, ypos)

    # Update the sums
    S0 <- S0 - old_contrib_0 + new_contrib_0
    St <- St - old_contrib_t + new_contrib_t

    # Compute the new nodf:
    nodf <- (S0+St) / my_norm

    # Pack up the results
    # Unpack even futher to get a set a variables that we can work with:
    MT <- list(mt_0, mt_t)
    Fill <- list(F0, Ft)
    DM <- list(DM0, DMt)
    ND <- list(ND0, NDt)
    S <- c(S0, St)
    # Unpack all the support data:
    support_data <- list(MT, Fill, DM, ND, S)
    return(list(nodf, mtx, support_data))
}

# Efficient way to compute the nodf value of a matrix where the
# link at position pos = list(xpos, ypos) is added:
nodf_one_link_added <- function(mtx, pos, support_data){
    # Unpack and update all the support data:
    xpos <- pos[[1]]
    ypos <- pos[[2]]
    MT <- support_data[[1]]
    Fill <- support_data[[2]]
    DM <- support_data[[3]]
    ND <- support_data[[4]]
    S <- support_data[[5]]
    # print("B")
    # Unpack even futher to get a set a variables that we can work with:
    mt_0 <- MT[[1]]
    mt_t <- MT[[2]]
    # print("B1")
    F0 <- Fill[[1]]
    Ft <- Fill[[2]]
    # print("B2")
    DM0 <- DM[[1]]
    DMt <- DM[[2]]
    # print("B3")
    ND0 <- ND[[1]]
    NDt <- ND[[2]]
    # print("B4")
    S0 <- S[[1]]
    St <- S[[2]]
    # print("B5")
    NodesA <- nrow(mtx)
    NodesB <- ncol(mtx)
    my_norm <- 0.5*((NodesA *(NodesA -1)) + (NodesB*(NodesB-1)))
    # print("B6")
    # Modify the matrix appropriately (remove a link):
    mtx[xpos, ypos] <- 1.0
    # print("C")

    # Compute old contributions
    old_contrib_0 <- get_contributions(F0, ND0, DM0, xpos)
    old_contrib_t <- get_contributions(Ft, NDt, DMt, ypos)

    # modify the marginal totals
    mt_0[xpos] <- mt_0[xpos] + 1
    mt_t[ypos] <- mt_t[ypos] + 1

    # Update the degree matrix:
    # m0 <- matrix(mt_0[xpos], nrow=length(mt_0), ncol=1)
    # mt <- matrix(mt_t[ypos], nrow=length(mt_t), ncol=1)
    m0 <- rep(mt_0[xpos], length(mt_0))
    mt <- rep(mt_t[ypos], length(mt_t))

    # Update the degree minima:
    DM0[xpos,] <- pmin3(m0, mt_0)
    DM0[,xpos] <- pmin3(m0, mt_0)
    DMt[ypos,] <- pmin3(mt, mt_t)
    DMt[,ypos] <- pmin3(mt, mt_t)

    # Update negative deltas:
    ND0[xpos, ] <- (m0 > mt_0)
    ND0[, xpos] <- (m0 < mt_0)
    NDt[ypos, ] <- (mt > mt_t)
    NDt[, ypos] <- (mt < mt_t)

    # Update the fill factors
    F0[,xpos] <- F0[,xpos] + mtx[,ypos]
    F0[xpos,] <- F0[xpos,] + t(mtx[,ypos])
    F0[xpos,xpos] <- F0[xpos,xpos] - 1
    #
    Ft[,ypos] <- Ft[,ypos] + t(mtx[xpos,])
    Ft[ypos,] <- Ft[ypos,] + mtx[xpos,]
    Ft[ypos,ypos] <- Ft[ypos,ypos] - 1

    # Compute the new contributions:
    new_contrib_0 <- get_contributions(F0, ND0, DM0, xpos)
    new_contrib_t <- get_contributions(Ft, NDt, DMt, ypos)
    #print(c(S0, old_contrib_0, new_contrib_0))
    #print(c(St, old_contrib_t, new_contrib_t))

    # Update the sums
    S0 <- S0 - old_contrib_0 + new_contrib_0
    St <- St - old_contrib_t + new_contrib_t
    #print(c(S0, St))

    # Compute the new nodf:
    nodf <- (S0+St) / my_norm

    # Pack up the results
    # Unpack even futher to get a set a variables that we can work with:
    MT <- list(mt_0, mt_t)
    Fill <- list(F0, Ft)
    DM <- list(DM0, DMt)
    ND <- list(ND0, NDt)
    S <- c(S0, St)
    # Unpack all the support data:
    support_data <- list(MT, Fill, DM, ND, S)
    return(list(nodf, mtx, support_data))
}

# Efficient way to compute the nodf value of a neighbor graph
nodf_neighbor <- function(mtx, oPos, zPos, support_data){
    my_res <- nodf_one_link_removed(mtx, oPos, support_data)
    mtx <- my_res[[2]]
    support_data <- my_res[[3]]
    my_res  <- nodf_one_link_added(mtx, zPos, support_data)
    # my_res = list(nodf, mtx, support_data)
    return(my_res)
}

# Efficient way to compute the nodf value of a neighbor graph
nodf_neighbor2 <- function(mtx, oPos, zPos, mt_0, mt_t, F0, Ft, DM0, DMt, ND0, NDt, S){
    xpos <- zPos[[1]]
    ypos <- zPos[[2]]
    my_nodf0<- nodf_one_link_added_cpp(mtx, xpos, ypos, mt_0, mt_t, F0, Ft,
                                        DM0, DMt, ND0, NDt, S)
    xpos <- oPos[[1]]
    ypos <- oPos[[2]]
    my_nodf <- nodf_one_link_removed_cpp(mtx, xpos, ypos, mt_0, mt_t, F0, Ft,
                                        DM0, DMt, ND0, NDt, S)
    return(my_nodf)
}

get_valid_ones <- function(mtx){
    NodesA <- nrow(mtx)
    NodesB <- ncol(mtx)
    sub_mtx <- mtx[2:NodesA,2:NodesB]
    one_pos <- which(sub_mtx == 1, arr.ind = T)
    shift_mtx <- matrix(1, nrow = nrow(one_pos), ncol = 2)
    return(one_pos + shift_mtx)
}

get_zeros <- function(mtx){
    zero_pos <- which(mtx == 0, arr.ind = T)
    return(zero_pos)
}

accept_probability <- function(new_cost, old_cost, temp){
    if(new_cost < old_cost){
        result <- 1.0
    }
    else{
        a <- -1.0*(new_cost- old_cost) / temp
        result <- exp(a)
    }
    return(result)
}
