#' Local Integration Targeting Comparisons
#'
#' Perform per-feature integration targeting comparisons within a given feature-set. This is analogous to differential expression analyses for e.g., RNA-seq and uses the same machinery in limma. The output object can be further analyzed using other limma functions, such as <code>decideTests()</code> and <code>topTable()</code>. See the limma documentation for more information.
#'
#'@param xint.obj The xInt object of interest.
#'@param min.count Used to filter lowly targeted features. If not provided, defaults to 1. Used as an input to edgeR's <code>filterByExpr()</code>.
#'@param min.total.count Used alongside min.count to filter lowly targeted features. If not provided, defaults to one-half of the number of samples. Used as an input to edgeR's <code>filterByExpr()</code>.
#'@param norm.method If NULL, only library size normalization is performed. Otherwise, must be one of the valid options in edgeR's <code>calcNormFactors()</code>.
#'@param plot Boolean. Whether or not to output diagnostic plots related to sample biases and mean-variance modeling. Defaults to TRUE.
#'@param return.contrasts Boolean. Whether or not to return the contrast fits. Defaults to TRUE. If FALSE, the lmFit() output is returned after empirical Bayes moderation.
#'
#'@return A MArrayLM object containing the results of feature-wise linear fits.
#'
#'@importFrom stats model.matrix
#'@import SummarizedExperiment
#'@import limma
#'@import edgeR
#'
#'@export
#'
local_comparisons <- function( xint.obj, min.count, min.total.count,
                               norm.method = NULL, plot = TRUE, return.contrasts = TRUE ){

  if( !validObject( xint.obj ) ){
    stop( "xint.obj is not a valid xIntObject.",
          call. = FALSE )
  }

  cs <- assay(xint.obj)
  dat <- colData( xint.obj )

  if( missing( min.count ) ){
    min.count <- 1
  }

  if( missing( min.total.count ) ){
    min.total.count <- floor( ncol( cs ) / 2 )
  }

  dge <- DGEList(
    counts = cs,
    lib.size = dat$total.sites,
    norm.factors = rep( 1, ncol( cs ) ),
    remove.zeros = TRUE,
    group = dat$condition
  )

  group <- dat$condition

  design <- model.matrix( ~0 + group )
  colnames( design ) <- gsub( pattern = "group", replacement = "", x = colnames( design ) )

  keep <- filterByExpr(
    y = dge,
    design = design,
    lib.size = dat$total.sites,
    min.count = min.count,
    min.total.count = min.total.count
  )

  dge <- dge[ keep, ]

  if( !is.null( norm.method ) ){
    dge <- calcNormFactors( dge, method = norm.method )
  }

  if( isTRUE( plot ) ){
    lcpm <- cpm( dge, log = TRUE )
    plotMDS( lcpm, labels = group )
  }

  design.cols <- colnames( design )
  contrast.list <- list()

  for (i in 1:( length( design.cols ) - 1 ) ) {
    for ( j in ( i + 1 ):length( design.cols ) ) {
      if ( design.cols[i] != design.cols[j] ){
        contrast.name <- paste( design.cols[i], "-", design.cols[j], sep = "" )
        contrast.expression <- substitute( x - y, list( x = as.name( design.cols[i] ),
                                                        y = as.name( design.cols[j] ) ) )
        contrast.list[[contrast.name]] <- contrast.expression
      }
    }
  }

  contrasts <- do.call( "makeContrasts", c( contrast.list, list( levels = design.cols ) ) )

  v <- voom( dge, design, plot = plot )

  fit <- lmFit( v, design )

  if( isTRUE( return.contrasts ) ){
    fit <- contrasts.fit( fit, contrasts = contrasts )
  }

  fit <- eBayes( fit )

  if( isTRUE( plot ) ){
    plotSA( fit )
  }

  return( fit )
}
