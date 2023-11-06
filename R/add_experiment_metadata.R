add_experiment_metadata <- function( xint.obj = NULL, 
                                     virus, 
                                     cell.type, 
                                     num.cells, 
                                     dpi, 
                                     frag.method, 
                                     moi,
                                     seq.platform,
                                     seq.type,
                                     read.length,
                                     date ){
  
  if ( !missing( virus ) ){
    v <- paste( "Virus:", virus )
  } else{
    v <- paste( "Virus:", NULL )
  }
  
  if ( !missing( cell.type ) ){
    ct <- paste( "Cell Type:", cell.type )
  } else{
    ct <- paste( "Cell Type:", NULL )
  }
  
  if ( !missing( num.cells ) ){
    nc <- paste( "Cell Number:", num.cells )
  } else{
    nc <- paste( "Cell Number:", NULL )
  }  
  
  if ( !missing( dpi ) ){
    days <- paste( "DPI:", dpi )
  } else{
    days <- paste( "DPI:", NULL )
  }    
  
  if ( !missing( frag.method ) ){
    fm <- paste( "Fragmentation Method:", frag.method )
  } else{
    fm <- paste( "Fragmentation Method:", NULL )
  }
  
  if ( !missing( moi ) ){
    m <- paste( "MOI:", moi )
  } else{
    m <- paste( "MOI:", NULL )
  }
  
  if ( !missing( seq.platform ) ){
    sp <- paste( "Sequencing Platform:", seq.platform )
  } else{
    sp <- paste( "Sequencing Platform:", NULL )
  }  
  
  if ( !missing( seq.type ) ){
    st <- paste( "Sequencing Type:", seq.type )
  } else{
    st <- paste( "Sequencing Type:", NULL )
  } 
  
  if ( !missing( seq.type ) ){
    rl <- paste( "Read Length:", read.length )
  } else{
    rl <- paste( "Read Length:", NULL )
  }
  
  if ( !missing( date ) ){
    d <- paste( "Date:", date )
  } else{
    d <- paste( "Date:", NULL )
  }
  
  out.ls <- list( v, ct, nc, days, fm, m, sp, st, rl, d  )
  out <- do.call( rbind, out.ls )
  
  if ( is.null( xint.obj ) ){
    return( out )
  } else{
    metadata( xint.obj )$experiment.info <- out
    return( xint.obj  )
  }
  
}
