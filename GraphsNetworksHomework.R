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

scCCnet <- readRDS("./data/scCCnet.rds")
scCCGraph <- igraph::graph_from_data_frame(scCCnet, directed = FALSE)



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

scCCnet <- readRDS("./data/scCCnet.rds")
scCCGraph <- igraph::graph_from_data_frame(scCCnet, directed = FALSE)
d <- igraph::degree(scCCGraph)
ddistro <- igraph::degree.distribution(scCCGraph)
d2 <- 4 * d

sim <- igraph::sample_degseq(d2)
dSim <- table(factor(igraph::degree(sim)))

graphKLDiv(x = sim)

