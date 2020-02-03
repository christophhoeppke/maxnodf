greedy_add_link <- function(mtx, support_data){
    # Find the most promising zeros:
    zPosList <- websearch_NODF_fast(mtx)
    opt_nodf <- -100.0
    opt_pos <- c(-1,-1)
    opt_res  <- 0
    for(idx in 1:nrow(zPosList)){
        zPos <- zPosList[idx, ]
        my_res <- nodf_one_link_added(mtx, zPos, support_data)
        new_nodf <- my_res[[1]]
        # mtx <- my_res[[2]]
        # support_data <- my_res[[3]]
        if(new_nodf > opt_nodf){
            opt_nodf <- new_nodf
            opt_pos <- zPos
            opt_res <- my_res
        }
        # Reverting is not necessary as the changes are never accepted!
        # Revert the change:
        # my_res <- nodf_one_link_removed(mtx, zPos, support_data)
        # new_nodf <- my_res[[1]]
        # mtx <- my_res[[2]]
        # support_data <- my_res[[3]]
    }
    return(opt_res)
}

greedy_add_link_smart <- function(mtx, support_data){
    # Find the most promising zeros:
    NodesA = nrow(mtx)
    NodesB = ncol(mtx)
    zPosList <- websearch_NODF_fast(mtx)
    opt_nodf <- -100.0
    opt_pos <- c(-1,-1)
    for(idx in 1:nrow(zPosList)){
        zPos <- zPosList[idx, ]
        my_res <- nodf_one_link_added(mtx, zPos, support_data)
        new_nodf <- my_res[[1]]
        mtx <- my_res[[2]]
        support_data <- my_res[[3]]
        if(new_nodf > opt_nodf){
            opt_nodf <- new_nodf
            opt_pos <- zPos
        }
        else if(new_nodf == opt_nodf){
            score_old <- (opt_pos[1] / NodesA) + (opt_pos[2] / NodesB)
            score_new <- (zPos[1] / NodesA) + (zPos[2] / NodesB)
            if(score_new < score_old){
                opt_nodf <- new_nodf
                opt_pos <- zPos
            }
        }
        # Revert the change:
        my_res <- nodf_one_link_removed(mtx, zPos, support_data)
        new_nodf <- my_res[[1]]
        mtx <- my_res[[2]]
        support_data <- my_res[[3]]
    }
    # actually perform the update
    my_res <- nodf_one_link_added(mtx, opt_pos, support_data)
    nodf <- my_res[[1]]
    mtx <- my_res[[2]]

    return(my_res)
}

greedy_solve <- function(NodesA, NodesB, Edges){
    mtx <- matrix(0.0, nrow=NodesA, ncol=NodesB)
    mtx[1,] = 1.0
    mtx[,1] = 1.0
    mtx[2,2] = 1.0
    support_data <- init_nodf(mtx)
    nodf <- nodf_cpp(mtx)
    tp <- utils::txtProgressBar(min = sum(mtx), max = Edges, style = 3)
    while(sum(mtx) < Edges){
        my_res <- greedy_add_link_smart(mtx, support_data)
        nodf <- my_res[[1]]
        mtx <- my_res[[2]]
        support_data <- my_res[[3]]
        utils::setTxtProgressBar(tp, sum(mtx))
    }
    return(mtx)
}

greedy_solve2 <- function(NodesA, NodesB, Edges){
    mtx <- matrix(0.0, nrow=NodesA, ncol=NodesB)
    mtx[1,] = 1.0
    mtx[,1] = 1.0
    mtx[2,2] = 1.0
    support_data <- init_nodf(mtx)
    nodf <- nodf_cpp(mtx)
    Edges_Left <- Edges - NodesA - NodesB
    if(Edges_Left >= 1){
        tp <- utils::txtProgressBar(min = 0, max = Edges_Left, style = 3)
        # Unpack the support data:
        MT <- support_data[[1]]
        Fill <- support_data[[2]]
        DM <- support_data[[3]]
        ND <- support_data[[4]]
        S <- support_data[[5]]
        # Unpack even futher to get a set a variables that we can work with:
        mt_0 <- as.vector(MT[[1]])
        mt_t <- as.vector(MT[[2]])
        # print("B1")
        F0 <- Fill[[1]]
        Ft <- Fill[[2]]
        # print("B2")
        DM0 <- DM[[1]]
        DMt <- DM[[2]]
        # print("B3")
        ND0 <- ND[[1]]*1
        NDt <- ND[[2]]*1

        for(i in 1:Edges_Left){
            zPosList <- websearch_cpp(mtx)
            # zPosList <- websearch_NODF_fast(mtx)
            opt_nodf <- -100.0
            opt_pos <- c(-1,-1)
            for(idx in 1:nrow(zPosList)){
                zPos <- zPosList[idx, ]
                xpos <- zPos[[1]]
                ypos <- zPos[[2]]
                new_nodf <- nodf_one_link_added_cpp(mtx,xpos,ypos,mt_0,mt_t,F0,Ft,DM0,DMt,ND0,NDt,S);
                if(new_nodf > opt_nodf){
                    opt_nodf <- new_nodf
                    opt_pos <- zPos
                }
                else if(new_nodf == opt_nodf){
                    score_old <- (opt_pos[1] / NodesA) + (opt_pos[2] / NodesB)
                    score_new <- (zPos[1] / NodesA) + (zPos[2] / NodesB)
                    if(score_new < score_old){
                        opt_nodf <- new_nodf
                        opt_pos <- zPos
                    }
                }
                # Revert the change:
                nodf_2 <- nodf_one_link_removed_cpp(mtx, zPos[1], zPos[2], mt_0, mt_t, F0, Ft, DM0, DMt, ND0, NDt, S);
            }
            # actually perform the update
            xpos <- opt_pos[[1]]
            ypos <- opt_pos[[2]]
            new_nodf <- nodf_one_link_added_cpp(mtx,xpos,ypos,mt_0,mt_t,F0,Ft,DM0,DMt,ND0,NDt,S);
            utils::setTxtProgressBar(tp, i)
        }
    }
    return(mtx)
}

