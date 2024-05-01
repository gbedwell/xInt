library(BSgenome.Hsapiens.UCSC.hs1)
library(xInt)

all.sites <- import_sites(
  path = "~/Documents/github/xInt_data/example_data/",
  pattern = ".bed",
  genome.obj = BSgenome.Hsapiens.UCSC.hs1
)

set.seed(1)

sites <- lapply(
  X = all.sites[ !names(all.sites) %in% "C1" ],
  FUN = function(x){
    nn <- sample(seq(10,30), size = 1)
    gr <- sample(x, size = nn)
    return( gr )
  }
)

uninf <- all.sites[ names(all.sites) %in% "C1" ][[1]]

set.seed(1)
uninf <- sample(uninf, size = 5)

set.seed(1)
sites <- lapply(
  X = sites,
  FUN = function(x){
    nn <- sample(seq(0,5), size = 1)
    un <- sample(uninf, size = nn)
    gr <- c( x, un )
    gr <- sortSeqlevels(gr)
    sort(gr, ignore.strand = TRUE)
  }
)

sites <- c( sites, list( C1 = uninf ) )

set.seed(1)
sites2 <- sample( all.sites[[1]], size = 1E4, replace = TRUE )
sites2 <- sortSeqlevels( sites2 )
sites2 <- sort( sites2, ignore.strand = TRUE )
sites2 <- list( A1 = sites2 )

feats <- rtracklayer::import(con="~/Documents/github/xInt_data/chm13v2.0_genes.bed.gz", format="BED")
feats$score <- NULL
feats$name <- paste0("Feature", seq_along(feats))

xobj <- make_xIntObject(
  site.list = all.sites[ !names(all.sites) %in% "C1" ],
  features = feats,
  conditions = c(rep("A", 4), rep("B", 5)),
  condition.levels = c("A", "B"),
  id.col = "name"
)

xobj <- xobj[-which( rowSums( assay(xobj) ) == 0 ),]

usethis::use_data(sites, sites2, xobj, internal = FALSE, overwrite = TRUE)
