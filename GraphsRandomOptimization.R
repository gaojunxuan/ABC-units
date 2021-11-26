KLdiv <- function(p, q) {
  # p and q are two pmfs of discrete probability distributions
  # with the same outcomes, which are nowhere 0.
  # Value:  Kullback-Leibler divergence  sum(p * log( p / q))).

  if (length(p) != length(q)) {
    print(length(p))
    print(length(q))
    stop("PANIC: input vector lengths differ!")
  }
  if (any(c((p == 0), (q == 0)))) {
    stop("PANIC: 0's found in input vectors!")
  }

  return(sum(p * log( p / q )))
}

graphKLDiv <- function(x) {
  scNet <- readRDS("./data/scCCnet.rds")
  scGraph <- igraph::graph_from_data_frame(scNet, directed = FALSE)

  # calculate the degree of each node
  scDeg <- igraph::degree(scGraph)
  xDeg <- igraph::degree(x)

  # calculate the maximum degree
  maxDeg <- max(max(scDeg), max(xDeg))

  # we want to count everything from 0 to the maximum degree
  # so we use "factor()" to set the explicitly specify the categories
  # to include, so that degrees with 0 count also occur in the final
  # count when passed to "table()"
  scCount <- table(factor(scDeg, levels=0:maxDeg))
  xCount <- table(factor(xDeg, levels=0:maxDeg))

  # add pseudocount of 0.5 in case of 0
  scCount <- scCount + 0.5
  xCount <- xCount + 0.5

  scFreq <- scCount / sum(scCount)
  xFreq <- xCount / sum(xCount)

  return(KLdiv(scFreq, xFreq))
}

simForestFire <- function() {
  fw.prob = seq(from = 0,
                to = 0.9,
                by = 0.01)
  bw.factor = seq(from = 0,
                  to = 1,
                  by = 0.1)
  palette(rainbow(11))

  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,3))
  axis(1)
  axis(2)
  title(xlab="fw.prob", ylab="KL divergence")

  minKL <- Inf

  bestFwProb <- 0
  bestBwFactor <- 0

  for (ii in seq_along(bw.factor)) {
    klDiv <- numeric()
    for (j in fw.prob) {
      i = bw.factor[ii]

      set.seed(112358)
      ff <- igraph::sample_forestfire(283, fw.prob=j, bw.factor=i, directed=F)
      set.seed(NULL)

      ffKL <- graphKLDiv(ff)
      if (ffKL < minKL) {
        minKL <- ffKL
        bestFwProb <- j
        bestBwFactor <- i
      }
      klDiv <- rbind(klDiv, ffKL)
    }
    points(fw.prob, klDiv, col=ii, pch=16)
  }

  fwProbLegend = paste("bw.factor =", bw.factor)
  legend("topleft",
         inset=c(0.05,0),
         cex = 0.5,
         legend=fwProbLegend,
         col=1:length(bw.factor),
         pch=16)

  print(bestFwProb)
  print(bestBwFactor)
  print(minKL)
}
simForestFire()

simER <- function() {
  prob <- seq(0, 0.04, 0.001)

  minKL <- Inf

  bestProb <- 0

  klDiv <- numeric()

  for (p in prob) {
    set.seed(112358)
    erGraph <- igraph::sample_gnp(283, p=p)
    set.seed(NULL)

    kl <- graphKLDiv(erGraph)
    klDiv <- rbind(klDiv, kl)

    if (kl < minKL) {
      minKL <- kl
      bestProb <- p
    }
  }
  print(bestProb)
  print(minKL)
  plot(prob, klDiv, xlab="p", ylab="KL divergence")
}

simER()

simSmallWorld <- function() {
  p = seq(from = 0,
          to = 0.9,
          by = 0.01)
  nei = seq(from = 1,
            to = 10,
            by = 1)
  palette(rainbow(11))

  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,3))
  axis(1)
  axis(2)
  title(xlab="p", ylab="KL divergence")

  minKL <- Inf

  bestP <- 0
  bestNei <- 0

  for (ii in seq_along(nei)) {
    klDiv <- numeric()
    for (j in p) {
      i = nei[ii]

      set.seed(112358)
      sw <- igraph::sample_smallworld(dim=1, size=283,
                                      nei=i, p=j)
      set.seed(NULL)

      swKL <- graphKLDiv(sw)
      if (swKL < minKL) {
        minKL <- swKL
        bestP <- j
        bestNei <- i
      }
      klDiv <- rbind(klDiv, swKL)
    }
    points(p, klDiv, col=ii, pch=16)
  }
  neiLegend = paste("nei =", nei)
  legend("bottomleft",
         inset=c(0.05,0.05),
         cex = 0.5,
         legend=neiLegend,
         col=1:length(nei),
         pch=16)

  print(bestP)
  print(bestNei)
  print(minKL)
}
simSmallWorld()


