rm(list = ls())
library(Biostrings)



seq_score <- readDNAStringSet("V:/PRG_BIN/exercise_07/seq_score.fasta", format = "fasta")
starting_indices <- c(8, 30, 1, 19, 44)
motiv_length <- 6


Score <- function(starting_indices, seq_score, motiv_length, l){
  motivs <- DNAStringSet()
  counter <- 1
  for (i in starting_indices){
    seq_cs <- seq_score[counter]
    motivs <- c(motivs, subseq(seq_cs, start = i, end = i+motiv_length-1))
  }
  

  freq_matrix <- consensusMatrix(motivs)
  
  max_values <- apply(freq_matrix[1:l, ], 2, max)
  
  
  score <- sum(max_values)
  return(score)
  
}

Score(starting_indices, seq_score, motiv_length)


NextLeaf <- function(array_start, L, k){
  for (i in L:1){
    if (array_start[i] < k){
      array_start[i] <- array_start[i] + 1
      return(array_start)
    }
    array_start[i] <- 1
  }
  return(array_start)
}


array_start <- c(8, 30, 1, 19, 44)
L <- length(seq_score)
k <- length(seq_score[[1]]) - length(array_start) + 1

NextLeaf(array_start, L, k)


BFMotifSearch <- function(DNA, t, n, l){
  s <- rep(1, l)
  bestScore <- Score(s, DNA, l)
  bestMotif <- c()
  while(TRUE){
    s <- NextLeaf(s, t, n - l + 1)
    if (Score(s, DNA, l) > bestScore){
      bestScore <- Score(s, DNA, l)
      bestMotif <- s
    }
    if (all(s == rep(1, l))){
      return(bestMotif)
    }
      
  }
}




DNA <- readDNAStringSet("V:/PRG_BIN/exercise_07/seq_motif.fasta", format = "fasta")
t <- length(DNA)
n <- length(DNA[[1]])
l <- 6


BFMotifSearch(DNA, t, n, l)






NextVertex <- function(array_start, i, L, k){
  if (i < L){
    array_start[i+1] = 1
    return(list(array_start, i+1))
  }
  else{
    for (j in L:1){
      if (array_start[j] < k){
        array_start[j] <- array_start[j] + 1
        return(list(array_start, j))
      }
    }
  }
  return(list(array_start, 0))
}




i <- 2
NextVertex(array_start, i, L, k)



ByPass <- function(array_start, i, L, k){
  for (j in i:1){
    if (array_start[j] < k){
      array_start[j] <- array_start[j] + 1
      return(list(array_start, j))
    }
  }
  return(list(array_start, 0))
  
}


ByPass(array_start, i, L, k)



BBMotifSearch <- function(array_start, t, L, k){
  s <- rep(1, l)
  bestScore <- 0
  i <- 1
  while(i > 0){
    if (i < t){
      optimisticScore <- Score(s, i, DNA, l) + (t - i)*l
      if (optimisticScore < bestScore){
        list(s, i) <- ByPass(s, i, t, n - l + 1)
      }
      else {
        list(s, i) <- NextVertex(s, i, t, n - l + 1)
      }
    }
    else{
      if (Score(s, t, DNA, l) > bestScore){
        bestScore <- Score(s, t, DNA, l)
        bestMotif <- s
      }
      list(s, i) <- NextVertex(s, i, t, n - l +1)
    }
  }
  return(bestMotif)
}

BBMotifSearch(array_start, t, L, k)
