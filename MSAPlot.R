#' Flatten a MsaAA/DNA/RNAMultipleAlignment object into an AA/DNA/RNAStringSet
#'
#' @param aln The msa object to flatten
#' @return The flattened BioStrings StringSet
flattenAlignment <- function(aln) {
  # check if aln is generated from msa
  if (class(aln) == "MsaAAMultipleAlignment" |
      class(aln) == "MsaDNAMultipleAlignment" |
      class(aln) == "MsaRNAMultipleAlignment") {
    # if so, unmask it
    aln_coerced <- Biostrings::unmasked(msaM)
    # split each sequence into a sequence of individual strings
    alnmat <- lapply(seq_along(aln_coerced), function(i) {
      strsplit(as.character(aln_coerced[[i]]), '')[[1]]
    })
    return(alnmat)
  }
  else {
    alnmat <- lapply(seq_along(aln), function(i) {
      strsplit(as.character(aln[[i]]), '')[[1]]
    })
    return(alnmat)
  }
}

#' Calculate the the step size required to reach the end of the current
#' alignment block
#'
#' @param j The starting index
#' @param gap A boolean indicating if the block is a sequence of gaps
#' @param alphabet A set of letters that could be in the block
#' @return The step size
calcStepSize <- function(seq, j, gap = FALSE, alphabet = c(letters, LETTERS)) {
  step <- 1
  lookAheadChar <- ifelse(gap, '-', alphabet)
  # look ahead until we reach the end of the block (either a gap, or
  # a character that is not in the alphabet)
  while (j + step <= length(seq) && any(seq[j+step] == lookAheadChar)) {
    step <- step + 1
  }
  return(step)
}

#' Plot the multiple sequence alignment
#' @param aln The alignment to be plotted
#' @param step_size The default step size. This will only be used if scores
#' is NULL
#' @param col The color of each aligned block, or the color palette if scores
#' is not NULL
#' @param scores The numeric vector returned from msa::msaConservationScore
plotAlignment <- function(aln, step_size = 1, col = "skyblue", scores = NULL) {
  if (missing(aln))
    stop("Missing input alignment data.")
  alnmat <- flattenAlignment(aln)

  # initialize the plotting env
  # adjust margin to fit sequence names
  op <- par(las=1, mar=c(5.1, 6.1, 4.1, 2.1))
  # try to read the name of each sequence
  y_names <- tryCatch(expr = { names(Biostrings::unmasked(aln)) },
                      error = function(e) { return(NULL) }
  )
  # if failed, default to sequence numbers
  if (is.null(y_names)) {
    y_names <- seq(1, length(alnmat), 1)
    # change margins back, we don't need them anymore
    op <- par(las=1, mar=c(5.1, 4.1, 4.1, 2.1))
  }
  plot.new()
  plot.window(xlim=c(0,length(alnmat[[1]])+50), ylim=c(0,length(alnmat)))
  axis(1, cex.axis=0.6)
  axis(2,
       at=seq(1, length(alnmat), 1),
       labels = y_names,
       cex.axis=0.5)
  # iterate over the sequences in the alignment
  for (i in seq_along(alnmat)) {
    segments(x0=0, y0=i, x1=length(alnmat[[i]]), y1=i)
    j <- 1
    while (j <= length(alnmat[[i]])) {
      residue <- alnmat[[i]][j]
      if (residue != "-") {
        # define the color palette if scores are provided
        if (!is.null(scores)) {
          myPal <- colorRampPalette(col)
          myCol <- myPal(max(scores))
        }
        # calculate the step size or use the default step size
        step <- ifelse(is.null(scores),
                       calcStepSize(alnmat[[i]], j),
                       step_size)
        # draw a rectangle for each aligned block
        rect(xleft = j,
             ybottom = i-0.25,
             xright = j+step,
             ytop = i+0.25,
             col = ifelse(is.null(scores), col, myCol[scores[j]]),
             border = NA)
        j <- j + step

      }
      else {
        step <- ifelse(is.null(scores),
                       calcStepSize(alnmat[[i]], j, gap = T),
                       step_size)
        j <- j + step
      }
    }
  }
  par(op)
}

plotAlignment(msaM)
