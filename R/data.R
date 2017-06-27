#' Data: Multiattribute Decisions (Hilbig & Moshagen, 2014)
#'
#' Choice frequencies of multiattribute decisions across 3 item types (Hilbig & Moshagen, 2014).
#'
#' @details
#' Each participant made 32 choices for each of 3 item types with four cues (with validities .9, .8, .7, and .6).
#'
#' The pattern of cue values of Option A and and B was as follows:
#' \itemize{
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
#' @template ref_hilbig14
#' @examples
#' data(hilbig2014)
#' head(hilbig2014)
#'
#' \dontrun{
#' # validities and cue values
#' v <- c(.9, .8, .7, .6)
#' cueA <- matrix(c(1,  1,  1, -1,
#'                  1, -1, -1, -1,
#'                  1,  1,  1, -1),
#'                ncol = 4, byrow = TRUE)
#' cueB <- matrix(c(-1, 1, -1,  1,
#'                  -1, 1,  1, -1,
#'                  -1, 1,  1,  1),
#'                ncol = 4, byrow = TRUE)
#'
#' # get strategy predictions
#' strategies <- c("baseline", "WADDprob", "WADD",
#'                 "TTB", "EQW", "GUESS")
#' preds <- predict_multiattribute(cueA, cueB, v, strategies)
#' c <- c(1, rep(.5, 5))  # upper bound of probabilities
#'
#' # get NML complexity and perform model selection
#' n <- rep(32, 3)
#' cnmls <- compute_cnml(preds, n, c, cores = 3)
#' select_nml(hilbig2014[1:5,], n, cnmls)
#' }
"hilbig2014"


#' Data: Multiattribute Decisions (Heck, Hilbig & Moshagen, 2017)
#'
#' Choice frequencies with multiattribute decisions across 4 item types (Heck, Hilbig & Moshagen, 2017).
#'
#' @details
#' Each participant made 40 choices for each of 4 item types with four cues
#' (with validities .9, .8, .7, and .6).
#' The pattern of cue values of Option A and and B was as follows:
#' \itemize{
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
#' @template ref_heck17
#' @examples
#' data(heck2017)
#' head(heck2017)
#' n <- rep(40, 4)
#'
#' # cue validities and values
#' v <- c(.9, .8, .7, .6)
#' cueA <- matrix(c(-1,  1,  1, -1,
#'                   1, -1, -1,  1,
#'                  -1,  1,  1,  1,
#'                   1, -1, -1, -1),
#'                ncol = 4, byrow = TRUE)
#' cueB <- matrix(c(-1, -1, -1, -1,
#'                  -1, 1 , -1, 1 ,
#'                  -1, 1 , 1 , -1,
#'                  -1, 1 , 1 , -1),
#'                ncol = 4, byrow = TRUE)
#'
#' # get predictions
#' strategies <- c("baseline", "WADDprob", "WADD",
#'                 "TTBprob", "TTB", "EQW", "GUESS")
#' preds <- predict_multiattribute(cueA, cueB, v, strategies)
#' c <- c(1, rep(.5, 6))  # upper bound of error probabilities
#'
#' # strategy classification with Bayes factor
#' select_bf(heck2017[1:4,], n, preds, c)
#'
#' \dontrun{
#' # strategy classification by NML (can take hours)
#' cnmls <- compute_cnml(preds, n, c = c, cores = 3)
#' select_nml(heck2017[1:4,], n, cnmls)
#' }
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
#' @template ref_heck17
#' @examples
#' data(heck2017_raw)
#' head(heck2017_raw)
#'
#' \dontrun{
#' # check predictions
#' cueA <- heck2017_raw[,paste0("a",1:4)]
#' cueB <- heck2017_raw[,paste0("b",1:4)]
#' v <- c(.9, .8, .7, .6)
#' ttb <- predict_multiattribute(cueA, cueB, v, "TTB")
#' table(ttb, heck2017_raw$ttb)
#' }
"heck2017_raw"


# colnames(heck2017_raw)
# heck2017_raw <- heck2017_raw %>% select(-choiceprob)
# save(heck2017_raw, file="data/heck2017_raw.RData")
