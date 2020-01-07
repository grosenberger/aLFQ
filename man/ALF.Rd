\name{ALF}
\alias{ALF}
\alias{ALF.default}

\title{Generate ALF report}
\description{Estimation of Absolute Protein Quantities of Unlabeled Samples by Targeted Mass Spectrometry.}
\usage{
\method{ALF}{default}(data, report_filename="ALF_report.pdf", 
prediction_filename="ALF_prediction.csv", peptide_methods = c("top"), 
peptide_topx = c(1,2,3), peptide_strictness = "loose", 
peptide_summary = "mean", transition_topx = c(1,2,3), 
transition_strictness = "loose", transition_summary = "sum", fasta = NA, 
apex_model = NA, combine_precursors = FALSE, combine_peptide_sequences = FALSE, 
consensus_proteins = TRUE, consensus_peptides = TRUE, consensus_transitions = TRUE,
cval_method = "boot", cval_mcx = 1000, ...)
}
\arguments{
  \item{data}{a mandatory data frame containing the columns \code{"run_id"}, \code{"protein_id"}, \code{"peptide_id"}, \code{"peptide_sequence"}, \code{"precursor_charge"}, \code{"peptide_intensity"} and \code{"concentration"} are required. For quantification on the transition level, the columns \code{"protein_id"}, \code{"peptide_id"}, \code{"transition_id"}, \code{"peptide_sequence"}, \code{"precursor_charge"}, \code{"transition_intensity"} and \code{"concentration"} are required. The id columns can be defined in any format, while the \code{"_intensity"} and \code{"concentration"} columns need to be numeric and in non-log form. The data may contain calibration data (with numeric \code{"concentration"} and test data (with \code{"concentration"} = "?"))}
  \item{report_filename}{the path and filename of the PDF report.}
  \item{prediction_filename}{the path and filename of the predictions of the optimal model.}
  \item{peptide_methods}{a vecter containing a combination of \code{"top"}, \code{"all"}, \code{"iBAQ"}, \code{"APEX"} or \code{"NSAF"} peptide to protein intensity estimation methods.}
  \item{peptide_topx}{(\code{"top"} only:) a positive integer value of the top x peptides to consider for \code{"top"} methods.}
  \item{peptide_strictness}{(\code{"top"} only:) whether \code{peptide_topx} should only consider proteins with the minimal peptide number (\code{"strict"}) or all (\code{"loose"}).}
  \item{peptide_summary}{(\code{"top"} and \code{"all"} only:) how to summarize the peptide intensities: \code{"mean"}, \code{"median"}, \code{"sum"}.}
  \item{transition_topx}{a positive integer value of the top x transitions to consider for transition to peptide intensity estimation methods.}
  \item{transition_strictness}{whether \code{transition_topx} should only consider peptides with the minimal transition number (\code{"strict"}) or all (\code{"loose"}).}
  \item{transition_summary}{how to summarize the transition intensities: \code{"mean"}, \code{"median"}, \code{"sum"}.}
  \item{fasta}{(\code{"iBAQ"}, \code{"APEX"} and \code{"NSAF"} only:) the path and filename to an amino acid fasta file containing the proteins of interest.}
  \item{apex_model}{(\code{"APEX"} only:) The \code{"APEX"} model to use (see APEX).}
  \item{combine_precursors}{whether to sum all precursors of the same peptide.}
  \item{combine_peptide_sequences}{whether to sum all variants (modifications) of the same peptide sequence.}
  \item{consensus_proteins}{if multiple runs are provided, select identical proteins among all runs.}
  \item{consensus_peptides}{if multiple runs are provided, select identical peptides among all runs.}
  \item{consensus_transitions}{if multiple runs are provided, select identical transitions among all runs.}
  \item{cval_method}{a method for doing crossvalidation: \code{"boot"} (bootstrapping), \code{"mc"} (monte carlo cross-validation), \code{"loo"} (leaving-one-out).}
  \item{cval_mcx}{a positive integer value of the number of folds for cross-validation.}
  \item{...}{future extensions.}
}

\details{
The ALF module enables model selection for TopN transitions and peptides for protein quantification (Ludwig et al., 2012). The workflow is completely automated and a report and prediction (using the best model) is generated.
}
\value{
  The reports specified in the function call.
}
\references{
  Ludwig, C., Claassen, M., Schmidt, A. \& Aebersold, R. \emph{Estimation of Absolute Protein Quantities of Unlabeled Samples by Selected Reaction Monitoring Mass Spectrometry.} Molecular \& Cellular Proteomics 11, M111.013987-M111.013987 (2012).
}
\author{George Rosenberger \email{gr2578@cumc.columbia.edu}}

\seealso{\code{\link{import}}, \code{\link{ProteinInference}}, \code{\link{AbsoluteQuantification}}, \code{\link{APEX}}, \code{\link{apexFeatures}}, \code{\link{proteotypic}}}

\examples{
\dontrun{data(UPS2MS)}

\dontrun{ALF(UPS2_SRM)}

\dontrun{data(LUDWIGMS)}

\dontrun{ALF(LUDWIG_SRM)}
}

\keyword{ ALF }
\keyword{ SRM }
\keyword{ SWATH }
\keyword{ SIS }
\keyword{ AQUA }
