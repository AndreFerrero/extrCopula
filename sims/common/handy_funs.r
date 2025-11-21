rCopFrechet <- function(alpha, cop) {
  # Given a copula, return Frechet(alpha) sample
  U <- as.vector(rCopula(1, cop))
  X <- qfrechet(U, shape = alpha)
  return(X)
}

block_max <- function(x, block_size) {
  # return block maxima given the block size
  sapply(split(x, ceiling(seq_along(x) / block_size)), max)
}