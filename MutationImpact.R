myFA <-             readFASTA("data/RAB39B_HSa_coding.fa")
myFA <- rbind(myFA, readFASTA("data/PTPN5_HSa_coding.fa"))
myFA <- rbind(myFA, readFASTA("data/PTPN11_HSa_coding.fa"))
myFA <- rbind(myFA, readFASTA("data/KRAS_HSa_coding.fa"))

# Set row names to gene names
rownames(myFA)<- c('RAB39B', 'PTPN5', 'PTPN11', 'KRAS')

#' Get the type of mutation given the sequence before and after the mutation
#' @param prior: the sequence before mutation
#' @param after: the sequence after mutation
#'
#' @return c("nonsense", "missense", "silent"), the type of mutation
getMutationType <- function(prior, after) {
  i <- 1
  while (i < nchar(prior)) {
    # start from the beginning
    # go by step of 3
    # check if the current block is the same in both prior and after
    if (substr(prior, i, i + 2) != substr(after, i, i + 2)) {
      oldAA <- Biostrings::GENETIC_CODE[[substr(prior, i, i + 2)]]
      newAA <- Biostrings::GENETIC_CODE[[substr(after, i, i + 2)]]
      if (newAA != oldAA) {
        if (newAA == "*")
          return("nonsense")
        return("missense")
      }
      else {
        return("silent")
      }
    }
    i <- i + 3
  }
  return("silent")
}

#' Randomly mutate a sequence and count the number each type of mutations
#' @param seqIn: the sequence to be mutated
#' @param N: number of iterations
#'
#' @return c(silentCount, nonsenseCount, missenseCount), count of each type
#' in the given order
sampleMutations <- function(seqIn, N=100000) {
  nuc <- c('A', 'T', 'C', 'G')
  silentCount <- 0
  nonsenseCount <- 0
  missenseCount <- 0
  for (i in 1:N) {
    pos <- sample(1:nchar(seqIn), size = 1)
    originalNuc <- substr(seqIn, pos, pos)
    # filter out the nucleotide that is currently at that position
    newNuc <- sample(nuc[nuc != originalNuc], size = 1)
    newSeq <- seqIn
    # apply the mutation
    substr(newSeq, pos, pos) <- newNuc

    mutType <- getMutationType(seqIn, newSeq)
    # count each type of mutation
    if (mutType == "silent")
      silentCount <- silentCount + 1
    else if (mutType == "nonsense")
      nonsenseCount <- nonsenseCount + 1
    else
      missenseCount <- missenseCount + 1
  }
  return(c(silentCount, nonsenseCount, missenseCount))
}

print(sampleMutations("ATGATGATGATGATGATG"))
print(sampleMutations("CCCCCCCCCCCCCCCCCC"))
print(sampleMutations("TATTACTATTACTATTAC"))
print(sampleMutations("TGGTGGTGGTGGTGGTGGTGGTGG"))
print(sampleMutations("TGTTGTTGTTGTTGTTGTTGTTGT"))

for (i in rownames(myFA)) {
  geneName <- i
  originalSeq <- myFA[i, ]$seq
  print(sprintf("Gene name: %s", geneName))
  sampleResult <- sampleMutations(as.character(originalSeq))
  silentCount <- sampleResult[1]
  nonsenseCount <- sampleResult[2]
  missenseCount <- sampleResult[3]

  print(sprintf("silent: %d, nonsense: %d, missense: %d",
                silentCount, nonsenseCount, missenseCount))
}

originalSeq <- myFA["KRAS",]$seq
sampleResult <- sampleMutations(as.character(originalSeq), N=2374)
sampleResult
