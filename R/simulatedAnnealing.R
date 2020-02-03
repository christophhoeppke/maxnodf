
# This file contains the final simulated annealing algorithm!
# at least it will once I am done implementing it :)

cost_function <- function(nodf_value){
    return(1.0 - nodf_value)
}

sim_anneal_step <- function(mtx, temp, iters, support_data){
    opt_mtx <- mtx
    opt_cost <- cost_function(nodf_cpp(mtx))
    old_cost <- opt_cost
    oPosList <- get_valid_ones(mtx)
    zPosList <- get_zeros(mtx)

    for(i in 1:iters){
        oidx <- sample(1:nrow(oPosList), 1)
        zidx <- sample(1:nrow(zPosList), 1)

        oPos <- oPosList[oidx, ]
        zPos <- zPosList[zidx, ]
        my_res <- nodf_neighbor(mtx, oPos, zPos, support_data)
        new_nodf <- my_res[[1]]
        new_cost <- cost_function(new_nodf)
        new_mtx <- my_res[[2]]
        new_support_data <- my_res[[3]]

        # Test if the new solution is optimal:
        if(new_cost < opt_cost){
            opt_cost <- new_cost
            opt_mtx <- new_mtx
        }

        # Test if we should accept the new solution
        new_cost <- cost_function(new_nodf)
        acc_prob <- accept_probability(new_cost, old_cost, temp)

        if(stats::runif(1,0,1) <= acc_prob){
            # Accept the new solution:
            mtx <- new_mtx
            support_data <- new_support_data
            old_cost <- new_cost
            oPosList <- get_valid_ones(mtx)
            zPosList[zidx, ] <- oPos
        }
        # Matrix is rejected by not updating mtx and support_data
    }
    # End of the iteration:
    return(list(opt_mtx, opt_cost, mtx, old_cost, support_data))
}

sim_anneal_optimisation <- function(mtx, alpha = 0.998, iters = 6, init_temp = 0.25, min_temp =1e-4){

    NodesA <- nrow(mtx)
    NodesB <- ncol(mtx)
    cool_steps <- round((log(min_temp / init_temp) / log(alpha))+ 0.5)

    # Allocate space for the optimal variables:
    opt_mtx <- mtx
    opt_nodf <- nodf_cpp(mtx)
    opt_cost <- cost_function(opt_nodf)

    temp <- init_temp
    eps <- 1e-12

    # Initialise support data
    support_data <- init_nodf(mtx)
    opt_support_data <- init_nodf(opt_mtx)

    pb <- utils::txtProgressBar(min = 0, max = cool_steps, style = 3)
    for(i in 1:cool_steps){
        utils::setTxtProgressBar(pb, i)
        my_res <- sim_anneal_step(mtx, temp, iters, support_data)
        new_opt_mtx <- my_res[[1]]
        new_opt_cost <- my_res[[2]]
        mtx <- my_res[[3]]
        new_cost <- my_res[[4]]
        support_data <- my_res[[5]]
        # Test if a new optimal solution was found:
        if(new_opt_cost + eps < opt_cost){
            opt_mtx <- new_opt_mtx
            # Hill climb at this point once I have implemented it!
            print("Hill Climbing!")
            opt_mtx <- full_hill_climb(opt_mtx)
            opt_nodf <- nodf_cpp(opt_mtx)
            opt_cost <- cost_function(opt_nodf)
            opt_support_data <- init_nodf(opt_mtx)
            # print(opt_nodf)
        }
        # Test to see if we should go back to the optimal solution
        acc_prob <- accept_probability(opt_cost, new_cost, temp)
        if(stats::runif(1,0,1) > acc_prob){
            mtx <- opt_mtx
            support_data <- opt_support_data
        }
        temp <- temp * alpha
    }
    return(opt_mtx)
}

sim_anneal_optimisation2 <- function(mtx, alpha= 0.998, iters=6, init_temp = 0.25, min_temp = 1e-4){
    NodesA <- nrow(mtx)
    NodesB <- ncol(mtx)
    cool_steps <- round((log(min_temp / init_temp) / log(alpha))+ 0.5)

    # Allocate space for the optimal variables:
    opt_mtx <- mtx
    opt_nodf <- nodf_cpp(mtx)
    opt_cost <- cost_function(opt_nodf)

    temp <- init_temp
    eps <- 1e-12

    # Initialise support data
    support_data <- init_nodf(mtx)
    opt_support_data <- init_nodf(opt_mtx)

    # Unpack the support data:
    MT <- support_data[[1]]
    Fill <- support_data[[2]]
    DM <- support_data[[3]]
    ND <- support_data[[4]]
    S <- support_data[[5]]
    # Unpack even futher to get a set a variables that we can work with:
    mt_0 <- as.vector(MT[[1]])
    mt_t <- as.vector(MT[[2]])
    F0 <- Fill[[1]]
    Ft <- Fill[[2]]
    DM0 <- DM[[1]]
    DMt <- DM[[2]]
    ND0 <- ND[[1]]*1
    NDt <- ND[[2]]*1

    pb <- utils::txtProgressBar(min = 0, max = cool_steps, style = 3)
    for(i in 1:cool_steps){
        utils::setTxtProgressBar(pb, i)
        # The inner sim_anneal_step:

        inner_opt_mtx <- mtx
        inner_opt_cost <- cost_function(nodf_cpp(mtx))
        inner_old_cost <- inner_opt_cost
        oPosList <- get_valid_ones_cpp(mtx)
        zPosList <- get_zeros(mtx)

        for(i in 1:iters){
            oidx <- sample(1:nrow(oPosList), 1)
            zidx <- sample(1:nrow(zPosList), 1)

            oPos <- oPosList[oidx, ]
            zPos <- zPosList[zidx, ]
            new_nodf <- nodf_neighbor2(mtx,oPos,zPos,mt_0,mt_t,F0,Ft,DM0,DMt,ND0,NDt,S)
            new_cost <- cost_function(new_nodf)

            # Test if the new solution is optimal:
            if(new_cost < inner_opt_cost){
                inner_opt_cost <- new_cost
                inner_opt_mtx <- mtx
            }

            # Test if we should accept the new solution
            new_cost <- cost_function(new_nodf)
            acc_prob <- accept_probability(new_cost, inner_old_cost, temp)

            if(stats::runif(1,0,1) <= acc_prob){
                # Accept the new solution:
                # Accept by not reverting back and updating the old_cost
                inner_old_cost <- new_cost
                oPosList <- get_valid_ones_cpp(mtx)
                zPosList[zidx,] <- oPos
            }else{
                # Reject by reverting back and not update the old_cost
                new_nodf <- nodf_neighbor2(mtx,zPos,oPos,mt_0,mt_t,F0,Ft,DM0,DMt,ND0,NDt,S)
            }
        }
        # End of the inner iteration.

        # Test if a new optimal solution was found:
        if(inner_opt_cost + eps < opt_cost){
            opt_mtx <- inner_opt_mtx
            # Hill climb at this point once I have implemented it!
            print("Hill Climbing!")
            opt_mtx <- full_hill_climb(opt_mtx)
            opt_nodf <- nodf_cpp(opt_mtx)
            opt_cost <- cost_function(opt_nodf)
            opt_support_data <- init_nodf(opt_mtx)
            print(opt_nodf)
        }
        # Test to see if we should go back to the optimal solution
        acc_prob <- accept_probability(new_cost, opt_cost, temp)
        print(c("acc_prob", acc_prob, 1.0 - opt_cost, 1.0 - new_cost))
        if(stats::runif(1,0,1) > acc_prob){
            print("Back to the optimum")
            mtx <- opt_mtx
            print(c(nodf_cpp(mtx), 1.0 - opt_cost))
            new_cost <- opt_cost
            # support_data <- init_nodf(mtx)
            # Re-unpack the support data:
            # Unpack the support data:
            MT <- opt_support_data[[1]]
            Fill <- opt_support_data[[2]]
            DM <- opt_support_data[[3]]
            ND <- opt_support_data[[4]]
            S <- opt_support_data[[5]]
            # Unpack even futher to get a set a variables that we can work with:
            mt_0 <- as.vector(MT[[1]])
            mt_t <- as.vector(MT[[2]])
            F0 <- Fill[[1]]
            Ft <- Fill[[2]]
            DM0 <- DM[[1]]
            DMt <- DM[[2]]
            ND0 <- ND[[1]]*1
            NDt <- ND[[2]]*1
        }
        temp <- temp * alpha
    }
    return(opt_mtx)
}
