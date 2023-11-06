#' Build Reference Model
#'
#' Build feature-wise negative binomial models for the reference or control samples.
#'
#'@param xint.obj An xInt object. The output of make_xint_dataset().
#'
#'@import DESeq2
#'
#'@export
#'
build_reference_model <- function( xint.obj ){

  # Create DESeq object
  tmp.dds <- DESeqDataSetFromMatrix(countData = assay( xint.obj ),
                                    colData = colData( xint.obj ),
                                    design = ~ 1 )

  # Manually define size factors.
  # Normalizes for library size only
  # Size factors defined relative to the geometric mean of library sizes
  sizes <- xint.obj$total.sites
  gm <- exp( mean( log( sizes ) ) )
  sf <- sizes * gm

  sizeFactors( tmp.dds ) <- sf

  # Use DESeq2 machinery to define feature-wise dispersions.
  # Uses the parameteric fit described in the 2014 DESeq2 paper
  tmp.dds <- estimateDispersions( tmp.dds, fitType = "parametric" )

  # Add size factors to xint.obj
  colData( xint.obj )$size.factor <- sizeFactors( tmp.dds )

  # Extract relevant information
  # Get MAP dispersion estimates -- ignore the dispersion outliers highlighted by DESeq2
  df <- rowData( tmp.dds )[ names( rowData( tmp.dds ) ) %in% c( "baseMean", "baseVar", "dispGeneEst",
                                                                "dispFit", "dispMAP" ) ]
  colnames( df ) <- c( "mean", "var", "disp.gw", "disp.fit", "disp.MAP" )

  # Define relative mean and dispersion values.
  # Use MAP dispersion estimates for this
  df$rel.mean <- df$mean / mean( sizes * sf )
  df$rel.disp <- df$disp.MAP / mean( sizes * sf )

  # Append data to xint.obj
  rowData( xint.obj ) <- cbind( rowData( xint.obj ), df )

  return( xint.obj )
}
