test_that("check FDR control (fdrControl)",{
    prob <- seq(0,1,0.001)
    expect_equal(sum(fdrControl(prob,0.05)),101)
})

test_that("test loglikDifferentialHMM (nb)",{

    # set seed
    set.seed(210427)
    N <- 1000
    M <- 2
    betas <- c(0,1)
    lambdas <- c(0,1)
    dt <- data.table::data.table(Window = rep(seq_len(N),times = M),
                                 ChIP = 1,
                                 offsets = -1,
                                 Dsg.Mix1 = rep(1:0,each = N),
                                 Dsg.Mix2 = rep(0:1,each = N),
                                 Group = 1)
    ll <- loglikDifferentialHMM(dt = dt,
                          theta = list('psi' = list(betas,lambdas),
                                       'delta' = c(0.1,0.9)),
                          N = 2,
                          M = 1000,
                          minZero = 1e-32,
                          dist = 'nb')
    
    unique_ll <- unique(ll)
    
    expect_true(length(unique_ll) == 3)
    
    # Consensus NB
    expect_equal(unique_ll[1],
                 dnbinom(1,mu = exp(-1),size = exp(1),log = TRUE)*2)
    # Mixture NB
    expect_equal(unique_ll[2],
                 dnbinom(1,mu = exp(-1),size = exp(1),log = TRUE) +
                     dnbinom(1,mu = exp(-1),size = exp(2),log = TRUE))
    # Consensus NB
    expect_equal(unique_ll[3],
                 dnbinom(1,mu = exp(-1),size = exp(2),log = TRUE)*2)
})

test_that("test loglikDifferentialHMM (zinb)",{
    
    # set seed
    set.seed(210427)
    N <- 1000
    M <- 2
    zip = 0
    betas <- c(0,1)
    lambdas <- c(0,1)
    dt <- data.table::data.table(Window = rep(seq_len(N),times = M),
                                 ChIP = 1,
                                 offsets = -1,
                                 Dsg.Mix1 = rep(1:0,each = N),
                                 Dsg.Mix2 = rep(0:1,each = N),
                                 Group = 1)
    ll <- loglikDifferentialHMM(dt = dt,
                                theta = list('delta' = c(0.5,0.5),
                                             'psi' = list(c(zip,betas),lambdas)),
                                N = 2,
                                M = 1000,
                                minZero = 1e-32,
                                dist = 'zinb')
    
    unique_ll <- unique(ll)
    
    expect_true(length(unique_ll) == 3)
    
    # Consensus ZINB
    expect_equal(unique_ll[1],
                 log((1/(1+exp(zip)))*(stats::dnbinom(1,mu = exp(-1),size = exp(1))))*2)
    # Differential mixture ZINB/NB
    expect_equal(unique_ll[2],
                 log((1/(1+exp(zip)))*(stats::dnbinom(1,mu = exp(-1),size = exp(1)))) +
                     dnbinom(1,mu = exp(-1),size = exp(2),log = TRUE))
    # Consensus NB
    expect_equal(unique_ll[3],
                 dnbinom(1,mu = exp(-1),size = exp(2),log = TRUE)*2)
})

test_that("test optimizers",{
    
    # set seed
    set.seed(210427)
    rcWeights <- data.table::data.table(id = rep(1:4,each = 250),
                                        weights = 1,
                                        Group = sample(25,1000,replace = TRUE),
                                        ChIP = 1,
                                        offsets = -1,
                                        Dsg.Mix1 = rep(1:0,each = 500),
                                        Dsg.Mix2 = rep(0:1,each = 500),
                                        Intercept = 1)
    
    maxDisp <- log(controlEM()$maxDisp)
    
    opt.nb <- optimDifferential(par = c(1,2,3,4),rcWeights = rcWeights,control = controlEM(),dist = 'nb')
    
    expect_equal(opt.nb$par[1],1,tolerance = 1e-5)
    expect_equal(opt.nb$par[2],maxDisp,tolerance = 1e-5)
    expect_equal(opt.nb$par[3],0,tolerance = 1e-5)
    expect_true(opt.nb$par[4]>0)

    opt.zinb <- optimDifferential(par = c(0,1,2,3,4),rcWeights = rcWeights,control = controlEM(),dist = 'zinb')
    
    expect_equal(1/(1+exp(-opt.zinb$par[1])),0,tolerance = 1e-5)
    expect_equal(opt.zinb$par[2],1,tolerance = 1e-5)
    expect_equal(opt.zinb$par[3],maxDisp,tolerance = 1e-5)
    expect_equal(opt.zinb$par[4],0,tolerance = 1e-5)
    expect_equal(opt.zinb$par[5],0,tolerance = 1e-5)
})
