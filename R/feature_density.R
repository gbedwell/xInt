feature_density <- function( site.list, features, win.size = 1E6, min.overlap = 1, average = TRUE ){
  
  gd <- lapply( X = site.list,
                FUN = function(x){
                  expanded <- resize( x, width = win.size, fix = "center", ignore.strand = TRUE )
                  
                  ol <- countOverlaps( query = expanded,
                                       subject = features,
                                       minoverlap = min.overlap,
                                       type = "any",
                                       ignore.strand = TRUE )
                
                if( isTRUE( average ) ){
                  
                  return( mean( ol ) )
                  
                  } else{
                    
                  return( ol )
                  
                  }
                }
              )
  
  if( isTRUE( average ) ){
    
    return( do.call( c, gd ) )
    
    } else{
      
      return( gd )
      
    }
}
