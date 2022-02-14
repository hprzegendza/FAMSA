#' @export
famsa <- function(stringset = NULL, help = FALSE, verbose = FALSE,
                  tree = FALSE, distance_matrix = FALSE, advanced_settings = "") {

  args <- strsplit(advanced_settings, " ")[[1]]
  argv <- NULL
  if (help) {
      print("To get help - type help(\"famsa\").")
      return()
  } else {
    ss_class <- class(stringset)
    if (length(grep("StringSet", ss_class)) == 0) {
      stop("Input must be an object of class XStringSet: DNAStringSet, RNAStringSet, or AAStringSet!\nSee Biostrings: https://bioconductor.org/packages/release/bioc/html/Biostrings.html\n")
    }

    argv <- append(argv, c("famsa"))
    if (verbose) {
      argv <- append(argv, c("-v"))
    }

    sequences_names <- names(stringset)
    if (!tree && !distance_matrix) {
      names(stringset) <- 1:length(stringset)
      temp_out <- tempfile(fileext = ".fasta")
    } else if (tree && distance_matrix) {
        stop("Both tree and distance matrix have been selected as the output!\nUse help(\"famsa\"), to print instruction.")
    } else if (tree) {
      argv <- append(argv, c("-gt_export"))
      temp_out <- tempfile()
    } else {
      argv <- append(argv, c("-dist_export"))
      temp_out <- tempfile(fileext = ".csv")
    }

    if (length(args) > 0) {
      for (i in 1:length(args)) {
        if (is.logical(args[[i]])) {
          argv <- append(argv, as.character(names(args)[i]))
        } else {
          argv <- append(argv, c(as.character(names(args)[i]), as.character(args[[i]])))
        }
      }
    }

    temp_in <- tempfile(fileext = ".fasta")
    Biostrings::writeXStringSet(stringset, filepath = temp_in, format = "fasta")

    argv <- append(argv, c(temp_in, temp_out))
  }


  nargs <- as.integer(length(argv))
  .C(.famsaCPP, nargs, as.character(argv))

  if (!tree && !distance_matrix) {
    ret <- switch(
            class(stringset),
            AAStringSet = Biostrings::readAAStringSet(temp_out, format = "fasta"),
            DNAStringSet = Biostrings::readDNAStringSet(temp_out, format = "fasta"),
            RNAStringSet = Biostrings::readRNAStringSet(temp_out, format = "fasta")
        )

    ret <- ret[order(as.integer(names(ret)))]
    names(ret) <- sequences_names

    ret <- switch(
            class(stringset),
            AAStringSet = Biostrings::AAMultipleAlignment(ret),
            DNAStringSet = Biostrings::DNAMultipleAlignment(ret),
            RNAStringSet = Biostrings::RNAMultipleAlignment(ret)
        )

  } else if (tree) {
    ret <- phytools::read.newick(temp_out)
  } else {
    ret <- as.dist(
      data.matrix(
        read.csv(file=temp_out, row.names = 1, header = F, fill=T, col.names=append(sequences_names,NA,))
      )
    )
  }
  file.remove(temp_in, temp_out)

  return(ret)
}
