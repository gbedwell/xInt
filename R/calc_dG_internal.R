calc_dG_internal <- function( seqs, vals, len ){
  sw <- function( x, win.size ) {
    res <- matrix( nrow = length(x) - win.size + 1, ncol = 2 )

    for (i in 1:(length(x) - win.size + 1)) {
      res[i,1] <- x[i]
      res[i,2] <- x[i + win.size - 1]
    }
    return(res)
  }

  inds <- sw( x = seq(1, len), win.size = 3 )

  trimer.sw <- sapply( X = 1:nrow(inds),
                       FUN = function(x){
                         as.character( Biostrings::subseq(seqs, start = inds[x,1], end = inds[x,2] ) )
                       }
  )

  dG.sw <- apply( X = trimer.sw,
                  MARGIN = 1,
                  FUN = function(x){
                    vals[ x ]
                  }
  )

  return( rowMeans( dG.sw, na.rm = TRUE ) )
}
