collapse_sites <- function( site.list, group.index.list, group.names ){
  
  if( !is.list( group.index.list ) ){
    stop( "Group index values must be given as a list.",
          call. = FALSE )
  }
  
  if( length( group.index.list ) != length( group.names ) ){
    stop( "The number of groups should equal the given group names.", 
          call. = FALSE ) 
  }
  
  ll <- lapply( X = 1:length( group.names ),
                FUN = function(x){
                  ind <- group.index.list[[x]]
                  out <- do.call( c, site.list[ ind ] )
                }
  )
  
  names( ll ) <- group.names
  
  return( ll )
}