#'@import methods
#'@export
setClass( "xIntObject",
          contains = "RangedSummarizedExperiment" )

setValidity(
  "xIntObject",
  function( object ){

    if ( !( "counts" %in% assayNames( object ) ) ){
      return( "The assay slot must contain a matrix named 'counts'." )
    }

    if ( any( is.na( assays( object ) ) ) ){
      return( "Count data cannot contain NA." )
    }

    if ( !is.integer( assays( object )$counts ) ){
      return( "Count data must be integers." )
    }

    if ( any( assays( object )$counts < 0 ) ){
      return( "Count data must be positive." )
    }

    if( !all( c( "sample", "total.sites", "overlapping.sites", "fraction.overlap", "condition" ) %in% names( colData( object ) ) ) ){

      msg <- paste( "colData must contain",
                    paste( c( "sample", "total.sites", "overlapping.sites", "fraction.overlap", "condition" ), collapse = ", " ),
                    collapse = " " )
      return( msg )
    }

    if( !is.integer( colData( object )$total.sites ) ){
      return( "'total.sites' must be integers." )
    }

    if( !is.integer( colData( object )$overlapping.sites ) ){
      return( "'overlapping.sites' must be integers." )
    }

    if( !is.numeric( colData( object )$fraction.overlap ) ){
      return( "'fraction.overlap' must be numeric." )
    }

    if( !is.factor( colData( object )$condition ) ){
      return( "'condition' column must be a factor.")
    }

    if( !all( colData( object )$sample == colnames( object ) ) ){
      return( "count matrix column names and 'sample' column in colData must match." )
    }

    return( TRUE )

  }
)
