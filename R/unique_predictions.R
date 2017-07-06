# #' Find Unique Predictions
# #'
# #' Allows to find unique item types, which are defined as patterns of cue values
# #' that lead to identical strategy predictions.
# #' @param predictions a matrix or data frame with predictions (one strategy/model per column).
# #'   Negative value = strength of preference for Option A.
# #'   Positive value = strength of preference for Option B.
# #'   0 = guessing with probability .50. See \code{\link{predict_multiattribute}}
# #' @param reversed whether to exclude reversed patterns
# #'    (i.e., those that are identical after switching Option A and B)
# #' @param return whether to return \code{"subset"} of input matrix of predictions or the \code{"index"} of unique lines
# #' @export
# unique_predictions <- function (predictions, reversed = TRUE,
#                                 return = "subset"){
#
#   if (return == "subset"){
#     res <- unique(predictions)
#     # if (reversed){
#     #   apply(res, 1, function(p) any(duplicated()))
#     # }
#   } else {
#     res <- which(!duplicated(predictions))
#   }
#
#   # if (reversed){
#   #   find
#   #   r <- apply(predictions == - predictions, 1, all)
#   #   r <- duplicated(predictions, - predictions)
#   # } else {
#   #   r <- 1:nrow(predictions)
#   # }
#   res
# }
