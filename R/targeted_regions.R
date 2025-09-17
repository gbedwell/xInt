#' Unbiased Target Region Identification
#'
#' @param sites A SiteList object
#' @param conditions A character vector of conditions corresponding to the individual datasets.
#' Must be equal to length(sites) and must be in the appropriate order.
#' When not NULL (default), replicates from the same condition are pooled.
#' @param genome.obj The genome object being used.
#' @param bin.size The bin size resolution. Smaller values offer more granularity but increase noise.
#' Defaults to 1000.
#' @param minseglen The minimum number of consecutive bins required for a region to be called.
#' Larger sizes can help decrease spurious calls, but at the expense of identifying smaller regions.
#' Defaults to 10.
#' @param min.gapwidth The maximum distance between adjacent regions for them to be merged.
#' Defaults to 1.
#' @param p Savitzky–Golay polynomial order. Must be less-than n. Defaults to 2.
#' @param n Savitzky–Golay window size. Must be odd and greater-than p. Defaults to 3.
#' @param difference Boolean. Whether or not to return the unique regions for each sample/condition.
#' Defaults to FALSE.
#' @param BiocParallel parameters. Passed to bplapply().
#' Defaults to SerialParam() for cross-platform and cluster compatibility.
#' @param omit A character vector of chromosomes to omit from the analysis.
#' Defaults to c("chrM", "MT")
#' @param min.width The minimum allowed width of an identified region. Defaults to 100.
#'
#' @return A GRanges object containing or GRangesList object containing targeted regions for each sample/condition.
#'
#' @import GenomicRanges
#' @import BiocParallel
#' @importFrom signal sgolayfilt
#' @importFrom changepoint cpt.meanvar cpt.mean cpts param.est
#'
#' @export
#'
targeted_regions <- function(sites, genome.obj, conditions = NULL, bin.size = 1000, minseglen = 2, 
                             min.gapwidth = 1, p = 2, n = 3, difference = FALSE, BPPARAM = SerialParam(),
                             omit = c("chrM", "MT"), min.width = 100) {
  if (!validObject(sites)) {
    stop("sites is not a valid SiteList object.", call. = FALSE)
  }

  if (!is.null(conditions) && length(conditions) != length(sites)) {
    stop("conditions must have the same length as sites.", call. = FALSE)
  }

  if(!is.null(conditions)) {
    cond.list <- lapply(
      X = unique(conditions),
      FUN = function(x) {
        which(conditions == x)
      }
    )

    names(cond.list) <- unique(conditions)

    sites <- collapse_sites(sites, cond.list, sorted = TRUE)
    conditions <- unique(conditions)
  } else {
    conditions <- names(sites)
  }

  sl <- seqlengths(genome.obj)
  sl <- sl[!names(sl) %in% omit]
  tot.sl <- sum(sl)

  regions <- bplapply(
    X = sites@sites,
    FUN = function(s) {
      chroms <- as.character(seqnames(s)@values)
      chroms <- chroms[!chroms %in% omit]
      hs <- lapply(
        X = chroms,
        FUN = function(chrom) {
          tmp.sl <- sl[names(sl) == chrom]
          theor.rand <- tmp.sl / tot.sl
          tmp.sites <- s[seqnames(s) == chrom]
          
          if(length(tmp.sites) == 0) {
            return(GRanges())
          }

          int.pos <- start(tmp.sites)
          max.pos <- max(int.pos)
          bins <- seq(1, tmp.sl + bin.size, by = bin.size)
          bin.counts <- hist(int.pos, breaks = bins, plot = FALSE)$counts
          cpm <- (bin.counts / sum(bin.counts)) * 1e6
          log.cpm <- log2(cpm + 1)
          smoothed.counts <- sgolayfilt(log.cpm, p = p, n = n)

          segments <- changepoint::cpt.meanvar(
            data = smoothed.counts,
            test.stat = "Normal",
            method = "PELT",
            penalty = "MBIC",
            minseglen = minseglen
          )

          cpts.idx <- changepoint::cpts(segments)
          seg.means <- changepoint::param.est(segments)$mean
          cpt.pos <- bins[cpts.idx + 1]

          gr <- GRanges(
            seqnames = chrom,
            ranges = IRanges(start = c(1, cpt.pos), end = c(cpt.pos - 1, min((max(bins) - 1), tmp.sl)))
          )
          gr$mean <- seg.means
          gr <- gr[gr$mean > theor.rand]
          gr <- reduce(gr, min.gapwidth = min.gapwidth)

          gr.ov <- findOverlaps(query = gr, subject = tmp.sites)
          query.idx <- queryHits(gr.ov)
          subject.idx <- subjectHits(gr.ov)
          overlap.groups <- split(subject.idx, query.idx)

          adj.start <- vapply(
            X = overlap.groups, 
            FUN = function(indices) {
              min(int.pos[indices])
            },
            FUN.VALUE = numeric(1)
          )

          adj.end <- vapply(
            X = overlap.groups, 
            FUN = function(indices) {
              max(int.pos[indices])
            },
            FUN.VALUE = numeric(1)
          )

          start(gr)[as.numeric(names(overlap.groups))] <- adj.start
          end(gr)[as.numeric(names(overlap.groups))] <- adj.end
          
          gr <- gr[width(gr) >= min.width]

          return(gr)
        }
      )
      suppressWarnings(hs <- do.call(c, hs))
      return(hs)
    },
    BPPARAM = BPPARAM
  )
  if(isTRUE(difference)) {
    get_diff <- function(i, gr.list) {
      current <- gr.list[[i]]
      others <- Reduce(c, gr.list[-i])
      other.union <- reduce(others)
      result <- subtract(current, other.union) |> unlist()
      return(result)
    }

    diff.regions <- lapply(seq_along(regions), get_diff, gr.list = regions)
    diff.regions <- unlist(as(diff.regions, "GRangesList"))
    regions <- reduce(diff.regions, min.gapwidth = min.gapwidth)

  } else {
    regions <- unlist(as(regions, "GRangesList"))
    regions <- reduce(regions, min.gapwidth = min.gapwidth)
  }

  region.counts <- lapply(
    X = seq_along(sites@sites), 
    FUN = function(x) {
    dat <- sites@sites[[x]]
    overlap.counts <- countOverlaps(query = regions, subject = dat)
    lcpm <- calculate_cpm(
      y = as.matrix(overlap.counts), 
      lib.size = length(dat), 
      log = TRUE, 
      prior.count = 0.5, 
      lib.prior = 1)
    return(c(lcpm))
    }
  )

  names(region.counts) <- names(sites@sites)
  for (nm in names(region.counts)) {
    mcols(regions)[[paste0(nm, ".lcpm")]] <- region.counts[[nm]]
  }

  mcols.sum <- rowSums(as.data.frame(mcols(regions)))
  regions <- regions[mcols.sum > ncol(mcols(regions)) * min(mcols.sum)]

  regions$name <- paste0("region_", seq_along(regions))

  return(regions)
}