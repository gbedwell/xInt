#' Add Metadata
#'
#' Adds basic experiment metadata to an xInt object.
#' This is intended to be used to attach experiment information (e.g., virus type, cell type, genome fragmentation method, etc.)
#' directly to the related xInt object or to the R environment itself.
#' Any/all metadata fields can be left blank, if desired.
#'
#'@param xint.obj The xInt object of interest.
#'@param virus The type of virus used in the infection.
#'@param cell.type The type of cell infected.
#'@param poi The protein-of-interest in transfected cells.
#'@param num.cells The number of cells infected.
#'@param dpi The days post infection.
#'@param moi The multiplicity of infection.
#'@param frag.method The method of genome fragmentation.
#'@param seq.platform The sequencing platform used.
#'@param seq.type The type of sequencing (e.g., paired- or single-end).
#'@param read.length The maximum read length.
#'@param date The date of sequencing.
#'
#'@return The input xInt object with metadata attached.
#'
#'@examples
#'add_metadata(virus = "HIV-1",
#'             cell.type = "HEK293T",
#'             frag.method = "fragmentase",
#'             seq.type = "paired-end",
#'             read.length = "150 bp")
#'
#'@import S4Vectors
#'@import SummarizedExperiment
#'@import GenomicRanges
#'
#'@export
#'
add_metadata <- function( xint.obj = NULL,
                          virus,
                          cell.type,
                          poi,
                          num.cells,
                          dpi,
                          moi,
                          frag.method,
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

  if ( !missing( poi ) ){
    prot <- paste( "Protein-of-interest:", poi )
  } else{
    prot <- paste( "Protein-of-interest:", NULL )
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

  if ( !missing( moi ) ){
    m <- paste( "MOI:", moi )
  } else{
    m <- paste( "MOI:", NULL )
  }

  if ( !missing( frag.method ) ){
    fm <- paste( "Fragmentation Method:", frag.method )
  } else{
    fm <- paste( "Fragmentation Method:", NULL )
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

  out.ls <- list( v, ct, prot, nc, days, fm, m, sp, st, rl, d  )
  out <- do.call( rbind, out.ls )

  if ( is.null( xint.obj ) ){
    return( out )
  } else{
    if( !validObject( xint.obj ) ){
      stop( "xint.obj is not a valid xIntObject.",
            call. = FALSE )
    }
    metadata( xint.obj )$experiment.info <- out
    return( xint.obj  )
  }
}
