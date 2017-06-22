#' Data: Multiattribute Decisions (Hilbig & Moshagen, 2014)
#'
#' Dataset with multiattribute decisions across 3 item types (Hilbig & Moshagen, 2014).
#' @details
#' Each participant made 32 choices for each of 3 item types with four cues (with validities .9, .8, .7, and .6).
#' @format A data frame 3 variables:
#' \describe{
#'   \item{\code{B1}}{Frequency of choosing Option B for Item Type 1}
#'   \item{\code{B2}}{Frequency of choosing Option B for Item Type 2}
#'   \item{\code{B3}}{Frequency of choosing Option B for Item Type 3}
#' }
#' @template ref_hilbig14
#' @examples
#' data(hilbig2014)
#' head(hilbig2014)
#'
#' \dontrun{
#' n <- rep(32, 3)
#' preds <- list(WADDD = c(-1, 1, -1),
#'               TTB =   c(-1, -1, -1),
#'               EQW =   c(-1, 1, 0),
#'               GUESS = c(0, 0, 0))
#' cnmls <- compute_cnml(preds, n, c = .5)
#' cnml_baseline <- compute_cnml(1:3, n, c = 1)
#' select_nml(hilbig2014[1:5,], n,
#'            c(list(baseline=cnml_baseline), cnmls))
#' }
"hilbig2014"


#' Data: Multiattribute Decisions (Heck, Hilbig & Moshagen, 2017)
#'
#' Dataset with multiattribute decisions across 4 item types (Heck, Hilbig & Moshagen, 2017).
#' @details
#' Each participant made 40 choices for each of 4 item types with four cues (with validities .9, .8, .7, and .6).
#' @format A data frame 4 variables:
#' \describe{
#'   \item{\code{B1}}{Frequency of choosing Option B for Item Type 1}
#'   \item{\code{B2}}{Frequency of choosing Option B for Item Type 2}
#'   \item{\code{B3}}{Frequency of choosing Option B for Item Type 3}
#'   \item{\code{B4}}{Frequency of choosing Option B for Item Type 4}
#' }
#' @template ref_hilbig14
#' @examples
#' data(heck2017)
#' head(heck2017)
#'
#' \dontrun{
#' n <- rep(40, 4)
#' preds <- list(WADDD = c(-1, -1, -1, 1),
#'               TTB =   c(-1, -1, -1, -1),
#'               EQW =   c(-1, 0, -1, 1),
#'               GUESS = c(0, 0, 0, 0))
#' cnmls <- compute_cnml(preds, n, c = .5)
#' cnml_baseline <- compute_cnml(1:4, n, c = 1)
#' select_nml(heck2017[1:5,], n,
#'            c(list(baseline=cnml_baseline), cnmls))
#' }
"heck2017"
