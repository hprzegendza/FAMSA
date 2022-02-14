\name{famsa}
\docType{package}
\title{
  MultiplaeSequence Alignment using FAMSA
}
\description{
  FAMSA is a Fast and Accurate Multiple Sequence algorithm for large-scale multiple sequence alignments (400k proteins in 2 hours and     8GB of RAM). Can output aligned sequences in XStringSet as well as a distance matrix and a phylogenetic  tree in Newick format.
}
\details{
    Pass a XStringSet (DNA, RNA or AA) from Biostrings package into the famsa() function to get an aligned sequence set.

    famsa(string_set)

    To get a distance matrix add a distance_matrix = TRUE parameter.

    famsa(stringSet, distance_matrix = TRUE)

    To get a phylogenetic tree add a tree = TRUE parameter.

    famsa(stringSet, true = TRUE)

    To specify the internal parameters of FAMSA program execution, use advanced_settings.
    Accepted options are described at https://github.com/refresh-bio/FAMSA.

}
\author{
    Original FAMSA program:
    Sebastian Deorowicz,
    Agnieszka Debudaj-Grabysz,
    Adam Gudyś,

    R package implementation:
    Hubert Przegendza <hp.przegendza@gmail.com>

    Maintainer: Hubert Przegendza <hp.przegendza@gmail.com>
}
\references{
    Repository with the base program: https://github.com/refresh-bio/FAMSA

    Citing:
    Deorowicz, S., Debudaj-Grabysz, A., Gudyś, A. (2016) FAMSA: Fast and accurate multiple sequence alignment of huge protein families. Scientific Reports, 6, 33964
    https://www.nature.com/articles/srep33964
}
\keyword{ package, msa,  }
\seealso{
  Biostrings, BH
}
\examples{
    famsa(exampleDNAStringSet)
    famsa(exampleDNAStringSet, distance_matrix = TRUE)
    plot(famsa(exampleDNAStringSet, tree = TRUE))
}
