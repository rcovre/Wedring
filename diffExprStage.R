wedrDiffExprVersion <- "0.3.2"

wedrDiffExprHelp <- function() {
  writeLines(":: diffExprStage - The Wedring pipeline DE Stage ::
This script corresponds to the Differential Expression stage of the Wedring
pipeline. It calls several function from DESeq (a R's library from Bioconductor
project) to achieve differentially expressed genes. One table with the data
about the genes expression is generated. This table is a direct DESeq's output.
The script may also output up to three graphics, depending on the options passed
to it. The graphics are: the dispersion estimates of the samples, a MA plot
and/or a Volcano plot with the genes data.

Usage:
Rscript diffExprStage.R [options] <infile> <conditions> [<outfile>]

<infile>    Input table containing the counts per genome feature in each
treatment.
<conditions>    Comma-separated values containg the conditions of the
experiment. For example, 'treated,treated,untreated,untreated'. Whitespaces
are not allowed.
<outfile> -- Name of the file which the output table will be stored. If any
file is provided the name of the output file will be based on the name of the
input file. For example, if the input is 'counts.txt' the output will be
'counts_diffexpr.txt'.

Options:
Options are used in the format '-option=value'.
For example '-method=per-condition', '-img.size=large' etc.

 -help    Display help message.
 -version    Display version information.
 -method    Empirical dispersion computation method. Accepted values: 'pooled'
            (default), 'per-condition' or 'blind'.
 -sharing.mode    Wich data should be used, empirical and/or fitted values?
                  Accepted values: 'maximum' (default), 'fit-only' or
                  'gene-est-only'.
 -fit.type    Data fitting type. 'parametric' (default) or 'local'.
 -ma.plot    Generate MA plot. Default: off.
 -volcano.plot    Generate Volcano plot. Default: off.
 -dispest.plot    Generate dispersion estimates plot. Default: off.
 -alpha    Significance level. Defaults to 0.05.
 -img.size    Size in pixel of the plotted images. Accepted values: 'small'
              (256x256), 'medium' (516x516, default) or 'large' (1024x1024)")
}

parseOptions <- function(arguments, options) {
  parsed.opts <- list()
  parsed.opts[options] <- rep(NA, length(options))
  opts.idx <- c()
  opts.idx.psb <- grep('-.+', arguments) # possible options indexes
  for (idx in opts.idx.psb) {
    if (grepl('=', arguments[idx])) {
      opt <- gsub('-(.+)=.+', '\\1', arguments[idx])
      if (opt %in% options) {
        parsed.opts[[opt]] <- gsub('-.+=(.+)', '\\1', arguments[idx])
        opts.idx <- append(opts.idx, idx)
      }
    } else {
      opt <- substr(arguments[idx], 2, nchar(arguments[idx]))
      if (opt %in% options) {
        parsed.opts[[opt]] <- TRUE
        opts.idx <- append(opts.idx, idx)
      }
    }
  }
  arg.list = c()
  for (arg.idx in setdiff(seq(length(arguments)), opts.idx)) {
    arg.list <- append(arg.list, arguments[arg.idx])
  }
  parsed.opts$arguments <- arg.list
  return(parsed.opts)
}

plotDispEsts <- function(cds, imgsize="medium") {
  sizes <- c(small = 1, medium = 1, large = 1.5)
  plot(rowMeans(counts(cds, normalized = TRUE)), fitInfo(cds)$perGeneDispEsts,
       pch = '.',
       log = "xy",
       cex.lab = sizes[imgsize][[1]],
       xlab = "Per-gene mean expression",
       ylab = "Per-gene dispersion estimate")
  xg <- 10^seq(-.5, 5, length.out = 300)
  lines(xg, fitInfo(cds)$dispFun(xg), col = "red", lwd = 3)
}

plotMA <- function(res, alpha = .05, imgsize="medium") {
  sizes <- c(small = 1, medium = 1, large = 1.5)
  plot(res$baseMean, res$log2FoldChange,
       log = "x",
       pch = 20,
       cex.lab = sizes[imgsize][[1]],
       col = ifelse(res$padj < alpha, "red", "black"),
       xlab = "Base mean",
       ylab = "log2 Fold change")
}

plotVolcano <- function(res, alpha = .05, imgsize="medium") {
  sizes <- c(small = 1, medium = 1, large = 1.5)
  plot(res$log2FoldChange, -log10(res$pval),
       pch = 20,
       cex.lab = sizes[imgsize][[1]],
       col = ifelse(res$padj < alpha, "red", "black"),
       xlab = "log2 Fold change",
       ylab = "-log10 P-value (BH adjusted)")
}

prepareOutfile <- function(infile, outfile = NA) {
  if (is.na(outfile)) {
    of.base <- strsplit(infile, "\\.")[[1]]
    if (length(of.base) > 1) {
      of.ext <- tail(of.base, n=1)
      of.base <- paste(head(of.base, n = -1), collapse = ".")
    } else {
      of.ext = ""
      of.base <- paste(of.base, collapse = ".")
    }
    outfile <- ifelse(of.ext != "",
        paste(of.base, "_diffexpr.", of.ext, sep = ""),
        paste(of.base, "_diffexpr", sep = ""))
  } else {
    of.base <- strsplit(outfile, "\\.")[[1]]
    of.base <- ifelse(length(of.base) > 1,
        paste(head(of.base, n = -1), collapse = "."),
        paste(of.base, collapse = "."))
  }
  return(list(filename = outfile, base.name = of.base))
}

processDESeq <-function(infile, conditions, method = "pooled",
                         sharingmode = "maximum", fittype = "parametric",
                         dispestplot = FALSE, imgsize = "medium") {
  cnt.table <- read.table(infile, header = TRUE)
  conditions <- factor(conditions)
  conds.levels <- levels(conditions)
  cnt.data.set <- newCountDataSet(cnt.table, conditions)
  cnt.data.set <- estimateSizeFactors(cnt.data.set)
  cnt.data.set <- estimateDispersions(cnt.data.set,
                                      method = method,
                                      sharingMode = sharingmode,
                                      fitType = fittype)
  result <- nbinomTest(cnt.data.set, conds.levels[1], conds.levels[2])
  return(list(result = result, data.set = cnt.data.set))
}

processDiffExpr <- function(infile, conditions, outfile = NA, method = "pooled",
                             sharingmode = "maximum", fittype = "parametric",
                             maplot = FALSE, volcanoplot = FALSE, alpha = .05,
                             dispestplot = FALSE, imgsize = "medium") {
  out <- prepareOutfile(infile, outfile)
  deseq.out <- processDESeq(infile, conditions,
                            method = method,
                            fittype = fittype,
                            sharingmode = sharingmode)
  
  write.table(deseq.out$result, out$filename, sep = "\t", quote = FALSE,
              row.names = FALSE)
  
  if (any(maplot, volcanoplot, dispestplot)) {
    sizes <- c(small = 256, medium = 512, large = 1024)
    if (maplot) {
      jpeg(paste(out$base.name, "_maplot.jpg", sep = ""), quality = 100,
           height = sizes[imgsize][[1]], width = sizes[imgsize][[1]])
      plotMA(deseq.out$result, alpha = alpha, imgsize = imgsize)
      dev.off()
    }
    if (volcanoplot) {
      jpeg(paste(out$base.name, "_volcanoplot.jpg", sep = ""), quality = 100,
           height = sizes[imgsize][[1]], width = sizes[imgsize][[1]])
      plotVolcano(deseq.out$result, alpha = alpha, imgsize = imgsize)
      dev.off()
    }
    if (dispestplot) {
      jpeg(paste(out$base.name, "_dispersion_estimates.jpg", sep = ""),
           quality = 100, height = sizes[imgsize][[1]],
          width = sizes[imgsize][[1]])
      plotDispEsts(deseq.out$data.set, imgsize = imgsize)
      dev.off()
    }
  }
}

args <- commandArgs(trailingOnly = TRUE)
opts <- c("method", "sharing.mode", "fit.type", "ma.plot", "volcano.plot",
          "dispest.plot", "alpha", "img.size")
parsed.opts <- parseOptions(args, opts)

if ("-help" %in% args) {
  wedrDiffExprHelp()
  q(save = "no", status = 0)
} else if ("-version" %in% args) {
  writeLines(paste('diffExprStage.R version', wedrDiffExprVersion))
  q(save = "no", status = 0)
}

if (length(parsed.opts$arguments) < 2) {
  writeLines("Usage:
    Rscript wedrdep.R [options] <infile> <conditions> [<outfile>]")
  stop("Insufficient number of arguments.")
}

infile <- parsed.opts$arguments[1]
conditions <- unlist(lapply(strsplit(parsed.opts$arguments[2], ",")[[1]],
                            FUN = function(c) if (c != "") c))
outfile <- ifelse(is.na(parsed.opts$arguments[3]), NA, parsed.opts$arguments[3])

if (!is.na(parsed.opts$method) && parsed.opts$method != "pooled") {
  if (parsed.opts$method %in% c("per-condition", "blind")) {
    method <- parsed.opts$method
  } else {
    stop(paste("Invalid method name:", parsed.opts$method))
  }
} else {
  method <- "pooled"
}

if (!is.na(parsed.opts$sharing.mode) && parsed.opts$sharing.mode != "maximum") {
  if (parsed.opts$sharing.mode %in% c("fit-only", "gene-est-only")) {
    sharing.mode <- parsed.opts$sharing.mode
  } else {
    stop(paste("Invalid sharing mode:", parsed.opts$fit.type))
  }
} else {
  sharing.mode <- "maximum"
}

if (!is.na(parsed.opts$fit.type) && parsed.opts$fit.type != "parametric") {
  if (parsed.opts$fit.type == "local") {
    fit.type <- parsed.opts$fit.type
  } else {
    stop(paste("Invalid fit type:", parsed.opts$fit.type))
  }
} else {
  fit.type <- "parametric"
}

if (!is.na(parsed.opts$ma.plot)) {
  ma.plot <- TRUE
} else {
  ma.plot <- FALSE
}

if (!is.na(parsed.opts$volcano.plot)) {
  volcano.plot <- TRUE
} else {
  volcano.plot <- FALSE
}

if (!is.na(parsed.opts$dispest.plot)) {
  dispest.plot <- TRUE
} else {
  dispest.plot <- FALSE
}

if (!is.na(parsed.opts$alpha) && as.numeric(parsed.opts$alpha) != .05) {
  alpha <- as.numeric(parsed.opts$alpha)
  if (is.na(alpha))
    stop(paste("Invalid alpha value:", parsed.opts$alpha))
} else {
  alpha <- .05
}

if (!is.na(parsed.opts$img.size) && parsed.opts$img.size != "medium") {
  if (parsed.opts$img.size %in% c("small", "large")) {
    img.size <- parsed.opts$img.size
  } else {
    stop(paste("Invalid image size:", parsed.opts$img.size))
  }
} else {
  img.size <- "medium"
}

suppressMessages(library(DESeq))

processDiffExpr(infile, conditions, outfile,
                method = method,
                sharingmode = sharing.mode,
                fittype = fit.type,
                alpha = alpha,
                maplot = ma.plot,
                volcanoplot = volcano.plot,
                dispestplot = dispest.plot,
                imgsize = img.size)