#' Data: Multiattribute Decisions (Hilbig & Moshagen, 2014)
#'
#' Choice frequencies of multiattribute decisions across 3 item types (Hilbig & Moshagen, 2014).
#'
#' @details
#' Each participant made 32 choices for each of 3 item types with four cues (with validities .9, .8, .7, and .6).
#'
#' The pattern of cue values of Option A and and B was as follows:
#' \describe{
#'   \item{Item Type 1: }{A = (1, 1, 1, -1) vs. B = (-1, 1, -1, 1)}
#'   \item{Item Type 2: }{A = (1, -1, -1, -1) vs. B = (-1, 1, 1, -1)}
#'   \item{Item Type 3: }{A = (1, 1, 1, -1) vs. B = (-1, 1, 1, 1)}
#' }
#' @format A data frame 3 variables:
#' \describe{
#'   \item{\code{B1}}{Frequency of choosing Option B for Item Type 1}
#'   \item{\code{B2}}{Frequency of choosing Option B for Item Type 2}
#'   \item{\code{B3}}{Frequency of choosing Option B for Item Type 3}
#' }
#' @template ref_hilbig2014
#' @examples
#' data(hilbig2014)
#' head(hilbig2014)
#'
#' # validities and cue values
#' v <- c(.9, .8, .7, .6)
#' cueA <- matrix(
#'   c(
#'     1, 1, 1, -1,
#'     1, -1, -1, -1,
#'     1, 1, 1, -1
#'   ),
#'   ncol = 4, byrow = TRUE
#' )
#' cueB <- matrix(
#'   c(
#'     -1, 1, -1, 1,
#'     -1, 1, 1, -1,
#'     -1, 1, 1, 1
#'   ),
#'   ncol = 4, byrow = TRUE
#' )
#'
#' # get strategy predictions
#' strategies <- c(
#'   "baseline", "WADDprob", "WADD",
#'   "TTB", "EQW", "GUESS"
#' )
#' preds <- strategy_multiattribute(cueA, cueB, v, strategies)
#' c <- c(1, rep(.5, 5)) # upper bound of probabilities
#'
#' # use Bayes factor for strategy classification
#' n <- rep(32, 3)
#' strategy_postprob(k = hilbig2014[1:5, ], n, preds)
"hilbig2014"


#' Data: Multiattribute Decisions (Heck, Hilbig & Moshagen, 2017)
#'
#' Choice frequencies with multiattribute decisions across 4 item types (Heck, Hilbig & Moshagen, 2017).
#'
#' @details
#' Each participant made 40 choices for each of 4 item types with four cues
#' (with validities .9, .8, .7, and .6).
#' The pattern of cue values of Option A and and B was as follows:
#' \describe{
#'   \item{Item Type 1: }{A = (-1, 1, 1, -1) vs. B = (-1, -1, -1, -1)}
#'   \item{Item Type 2: }{A = (1, -1, -1, 1) vs. B = (-1, 1, -1, 1)}
#'   \item{Item Type 3: }{A = (-1, 1, 1, 1) vs. B = (-1, 1, 1, -1)}
#'   \item{Item Type 4: }{A = (1, -1, -1, -1) vs. B = (-1, 1, 1, -1)}
#' }
#' Raw data are available as \code{\link{heck2017_raw}}
#' @format A data frame 4 variables:
#' \describe{
#'   \item{\code{B1}}{Frequency of choosing Option B for Item Type 1}
#'   \item{\code{B2}}{Frequency of choosing Option B for Item Type 2}
#'   \item{\code{B3}}{Frequency of choosing Option B for Item Type 3}
#'   \item{\code{B4}}{Frequency of choosing Option B for Item Type 4}
#' }
#' @template ref_heck2017
#' @examples
#' data(heck2017)
#' head(heck2017)
#' n <- rep(40, 4)
#'
#' # cue validities and values
#' v <- c(.9, .8, .7, .6)
#' cueA <- matrix(
#'   c(
#'     -1, 1, 1, -1,
#'     1, -1, -1, 1,
#'     -1, 1, 1, 1,
#'     1, -1, -1, -1
#'   ),
#'   ncol = 4, byrow = TRUE
#' )
#' cueB <- matrix(
#'   c(
#'     -1, -1, -1, -1,
#'     -1, 1, -1, 1,
#'     -1, 1, 1, -1,
#'     -1, 1, 1, -1
#'   ),
#'   ncol = 4, byrow = TRUE
#' )
#'
#' # get predictions
#' strategies <- c(
#'   "baseline", "WADDprob", "WADD",
#'   "TTBprob", "TTB", "EQW", "GUESS"
#' )
#' strats <- strategy_multiattribute(cueA, cueB, v, strategies)
#'
#' # strategy classification with Bayes factor
#' strategy_postprob(heck2017[1:4, ], n, strats)
"heck2017"



#' Data: Multiattribute Decisions (Heck, Hilbig & Moshagen, 2017)
#'
#' Raw data with multiattribute decisions (Heck, Hilbig & Moshagen, 2017).
#'
#' @details
#' Each participant made 40 choices for each of 4 item types with four cues
#' (with validities .9, .8, .7, and .6).
#' Individual choice freqeuncies are available as \code{\link{heck2017}}
#'
#' @format A data frame with 21 variables:
#' \describe{
#'   \item{\code{vp}}{ID code of participant}
#'   \item{\code{trial}}{Trial index}
#'   \item{\code{pattern}}{Number of cue pattern}
#'   \item{\code{ttb}}{Prediction of take-the-best (TTB)}
#'   \item{\code{eqw}}{Prediction of equal weights (EQW)}
#'   \item{\code{wadd}}{Prediction of  weighted additive (WADD)}
#'   \item{\code{logoddsdiff}}{Log-odds difference (WADDprob)}
#'   \item{\code{ttbsteps}}{Number of TTB steps (TTBprob)}
#'   \item{\code{itemtype}}{Item type as in paper}
#'   \item{\code{reversedorder}}{Whether item is reversed}
#'   \item{\code{choice}}{Choice}
#'   \item{\code{rt}}{Response time}
#'   \item{\code{choice.rev}}{Choice (reversed)}
#'   \item{\code{a1}}{Value of Cue 1 for Option A}
#'   \item{\code{a2}}{Value of Cue 2 for Option A}
#'   \item{\code{a3}}{Value of Cue 3 for Option A}
#'   \item{\code{a4}}{Value of Cue 4 for Option A}
#'   \item{\code{b1}}{Value of Cue 1 for Option B}
#'   \item{\code{b2}}{Value of Cue 2 for Option B}
#'   \item{\code{b3}}{Value of Cue 3 for Option B}
#'   \item{\code{b4}}{Value of Cue 4 for Option B}
#' }
#' @template ref_heck2017
#' @seealso \code{\link{heck2017}} for the aggregated choice frequencies per item type.
#' @examples
#' data(heck2017_raw)
#' head(heck2017_raw)
#'
#' \donttest{
#' # get cue values, validities, and predictions
#' cueA <- heck2017_raw[, paste0("a", 1:4)]
#' cueB <- heck2017_raw[, paste0("b", 1:4)]
#' v <- c(.9, .8, .7, .6)
#' strat <- strategy_multiattribute(
#'   cueA, cueB, v,
#'   c(
#'     "TTB", "TTBprob", "WADD",
#'     "WADDprob", "EQW", "GUESS"
#'   )
#' )
#'
#' # get unique item types
#' types <- strategy_unique(strat)
#' types$unique
#'
#' # get table of choice frequencies for analysis
#' freq <- with(
#'   heck2017_raw,
#'   table(vp, types$item_type, choice)
#' )
#' freqB <- freq[, 4:1, 1] + # reversed items: Option A
#'   freq[, 5:8, 2] # non-rev. items: Option B
#' head(40 - freqB)
#' data(heck2017)
#' head(heck2017) # same frequencies (different order)
#'
#' # strategy classification
#' pp <- strategy_postprob(
#'   freqB[1:4, ], rep(40, 4),
#'   types$strategies
#' )
#' round(pp, 3)
#' }
"heck2017_raw"


# colnames(heck2017_raw)
# heck2017_raw <- heck2017_raw %>% select(-choiceprob)
# save(heck2017_raw, file="data/heck2017_raw.RData")
