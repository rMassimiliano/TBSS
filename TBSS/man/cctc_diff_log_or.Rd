% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/case_case_time_control_stat.R
\name{cctc_diff_log_or}
\alias{cctc_diff_log_or}
\title{Discrepancy measure for the case case time controll TBSS}
\usage{
cctc_diff_log_or(
  control_reference,
  control_hazard,
  case_reference,
  case_hazard
)
}
\arguments{
\item{control_reference:}{number of control in the reference window}

\item{control_hazard:}{number of control in the hazard window}

\item{case_reference:}{number of cases in the reference window}

\item{case_hazard:}{number of cases in the hazard window}
}
\description{
Compute the discrepancy measure used for the case-case-time-control TBSS. The discrepancty is the log difference of the odds ratios for the case-case-time-control analysis.
}
