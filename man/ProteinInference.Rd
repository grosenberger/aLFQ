\name{ProteinInference}
\alias{ProteinInference}
\alias{ProteinInference.default}

\title{Protein inference for aLFQ import data frame}
\description{Protein inference for aLFQ import data frame.}
\usage{
\method{ProteinInference}{default}(data, peptide_method = "top", peptide_topx = 2, 
peptide_strictness = "strict",peptide_summary = "mean", transition_topx = 3, 
transition_strictness = "strict",transition_summary = "sum", fasta = NA, 
apex_model = NA, combine_precursors = FALSE, combine_peptide_sequences = FALSE, 
consensus_proteins = TRUE, consensus_peptides = TRUE, 
consensus_transitions = TRUE, ...)
}
\arguments{
  \item{data}{a mandatory data frame containing the columns \code{"run_id"}, \code{"protein_id"}, \code{"protein_intensity"}, and \code{"concentration"} for quantification on the protein level. For quantification on the peptide level, the columns \code{"run_id"}, \code{"protein_id"}, \code{"peptide_id"}, \code{"peptide_sequence"}, \code{"precursor_charge"}, \code{"peptide_intensity"} and \code{"concentration"} are required. For quantification on the transition level, the columns \code{"protein_id"}, \code{"peptide_id"}, \code{"transition_id"}, \code{"peptide_sequence"}, \code{"precursor_charge"}, \code{"transition_intensity"} and \code{"concentration"} are required. The id columns can be defined in any format, while the \code{"_intensity"} and \code{"concentration"} columns need to be numeric and in non-log form. The data may contain calibration data (with numeric \code{"concentration"} and test data (with \code{"concentration"} = "?"))}
  \item{peptide_method}{one of \code{"top"}, \code{"all"}, \code{"iBAQ"}, \code{"APEX"} or \code{"NSAF"} peptide to protein intensity estimation methods.}
  \item{peptide_topx}{(\code{"top"} only:) a positive integer value of the top x peptides to consider for \code{"top"} methods.}
  \item{peptide_strictness}{(\code{"top"} only:) whether \code{peptide_topx} should only consider proteins with the minimal peptide number (\code{"strict"}) or all (\code{"loose"}).}
  \item{peptide_summary}{(\code{"top"} and \code{"all"} only:) how to summarize the peptide intensities: \code{"mean"}, \code{"median"}, \code{"sum"}.}
  \item{transition_topx}{a positive integer value of the top x transitions to consider for transition to peptide intensity estimation methods.}
  \item{transition_strictness}{whether \code{transition_topx} should only consider peptides with the minimal transition number (\code{"strict"}) or all (\code{"loose"}).}
  \item{transition_summary}{how to summarize the transition intensities: \code{"mean"}, \code{"median"}, \code{"sum"}.}
  \item{fasta}{(\code{"iBAQ"}, \code{"APEX"} and \code{"NSAF"}:) the path and filename to an amino acid fasta file containing the proteins of interest.}
  \item{apex_model}{(\code{"APEX"} only:) The \code{"APEX"} model to use (see APEX).}
  \item{combine_precursors}{whether to sum all precursors of the same peptide.}
  \item{combine_peptide_sequences}{whether to sum all variants (modifications) of the same peptide sequence.}
  \item{consensus_proteins}{if multiple runs are provided, select identical proteins among all runs.}
  \item{consensus_peptides}{if multiple runs are provided, select identical peptides among all runs.}
  \item{consensus_transitions}{if multiple runs are provided, select identical transitions among all runs.}

  \item{...}{future extensions.}
}

\details{
The ProteinInference module provides functionality to infer protein quantities from the measured precursor or fragment intensities or peptide spectral counts. If the dataset contains targeted MS2-level data, the paired precursor and fragment ion signals, the transitions, are first summarized to the precursor level. Different methods for aggregation can be specified, e.g. sum, mean or median and a limit for the selection of the most intense transitions can be provided. It is further possible to exclude precursors, which do not have sufficient transitions to fulfill this boundary. To summarize precursor intensities or spectral counts to theoretical protein intensities, the mean, TopN (Silva et al., 2006; Malmstrom et al., 2009; Schmidt et al., 2011; Ludwig et al., 2012), APEX (Lu et al., 2006), iBAQ (Schwanhausser et al., 2011) and NSAF (Zybailov et al., 2006) protein intensity estimators are provided. For APEX, iBAQ and NSAF, the protein database in FASTA format needs to be supplied. In terms of targeted data acquisition, for both APEX and iBAQ methods all peptides of a protein must be targeted. The results are reported in the same unified data structure as from the previous step
}
\value{
  A standard aLFQ import data frame on protein level.
}
\references{
Silva, J. C., Gorenstein, M. V., Li, G.-Z., Vissers, Johannes P. C. & Geromanos, S. J. \emph{Absolute quantification of proteins by LCMSE: a virtue of parallel MS acquisition.} Mol. Cell Proteomics 5, 144-156 (2006).

Malmstrom, J. et al. \emph{Proteome-wide cellular protein concentrations of the human pathogen Leptospira interrogans.} Nature 460, 762-765 (2009).

Schmidt, A. et al. \emph{Absolute quantification of microbial proteomes at different states by directed mass spectrometry.} Molecular Systems Biology 7, 1-16 (2011).

Ludwig, C., Claassen, M., Schmidt, A. & Aebersold, R. \emph{Estimation of Absolute Protein Quantities of Unlabeled Samples by Selected Reaction Monitoring Mass Spectrometry.} Molecular & Cellular Proteomics 11, M111.013987-M111.013987 (2012).

Lu, P., Vogel, C., Wang, R., Yao, X. & Marcotte, E. M. \emph{Absolute protein expression profiling estimates the relative contributions of transcriptional and translational regulation.} Nat Biotech 25, 117-124 (2006).

Schwanhausser, B. et al. \emph{Global quantification of mammalian gene expression control.} Nature 473, 337-342 (2011).

Zybailov, B. et al. \emph{Statistical Analysis of Membrane Proteome Expression Changes in Saccharomyces c erevisiae.} J. Proteome Res. 5, 2339-2347 (2006).

Gerster S., Kwon T., Ludwig C., Matondo M., Vogel C., Marcotte E. M., Aebersold R., Buhlmann P. \emph{Statistical approach to protein quantification.} Molecular & Cellular Proteomics 13, M112.02445 (2014).
}
\author{George Rosenberger \email{gr2578@cumc.columbia.edu}}

\seealso{\code{\link{import}}, \code{\link{AbsoluteQuantification}}, \code{\link{ALF}}, \code{\link{APEX}}, \code{\link{apexFeatures}}, \code{\link{proteotypic}}}

\examples{
data(UPS2MS)

data_ProteinInference <- ProteinInference(UPS2_SRM)
print(data_ProteinInference)
}

\keyword{ ProteinInference }
\keyword{ label-free }
\keyword{ quantification }
\keyword{ peptide inference }
\keyword{ protein inference }
