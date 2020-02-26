## -*- coding: utf-8 -*-

#' Creates the combined response variable
#'
#' \code{combined_response} calculates combined responses given the information
#' about the stimulus, ratings, and accuracy or decision.
#'
#' This function calculates the variable distribution of which is modelled by
#' Signal Detection Theory. If ratings are not provided this becomes a simple
#' binary classification decision: 1 (2) = participant decided that the stimulus
#' belongs to class 1 (2). If ratings are provided the resulting variable
#' indicates how strongly does the participant feel that the stimulus belongs to
#' the second class. For example, if ratings represent confidence on a 1-4 scale
#' and participant decided that stimulus = 1 with maximum confidence then
#' response = 1, if participant decided that stimulus = 1 with confidence 1 then
#' response = 4, if participant decided that stimulus = 2 with confidence 1 then
#' response = 5, if participant decided that stimulus = 2 with confidence = 4
#' then response = 8.
#'
#' @param stimulus a vector encoding the stimulus class. This can be of any type
#'   as long as it contains only two kinds of values, which will be converted to
#'   1 and 2 if necessary.
#' @param rating a vector of ratings. Ratings have to be encoded as intergers
#'   and all the possible values between min(rating) and max(rating) have to be
#'   present in the data. Either ratings or accuracy has to be provided. When
#'   ratings are not provided the resulting combined response is just the binary
#'   classification response calculated based on the provided accuracy and
#'   stimulus information.
#' @param accuracy an optional numeric vector: 1 = a correct, 0 = an incorrect
#'   binary classification response. Either \code{accuracy} or \code{decision}
#'   has to be provided.
#' @param decision an optional numeric vector: i = participant decided that the
#'   stimulus belongs to class 1, j = participant decided that the stimulus
#'   belongs to class 2, where j > i. Either \code{accuracy} or \code{decision}
#'   has to be provided.
#' @return a vector of combined or binary response values.
#' @examples
#' data(gabor)
#' ## Using this function we can obtain binary classification decisions from stimulus and accuracy data
#' gabor$dec = combined_response(gabor$stim, accuracy = gabor$acc)
#' ## ... or we can obtain combined response values
#' gabor$resp = combined_response(gabor$stim, gabor$rating, gabor$acc)
#' aggregate(resp ~ rating + dec, gabor, unique)
#' ## This is false, because gabor$stim encodes stimulus classess as 0 or 1 and dec encodes the decision as 1 or 2
#' all(gabor$acc == (gabor$stim == gabor$dec))
#' ## but this is true so we're fine
#' all(gabor$acc == ((gabor$stim+1) == gabor$dec))
#' @export
combined.response = function(stimulus, rating = NULL, accuracy = NULL, decision = NULL){
    ## 0/1 coding
    stimulus = as.numeric(as.factor(as.character(stimulus))) - 1
    if(is.null(decision)){
        if(is.null(decision) && is.null(accuracy))
            stop("Neither decision nor accuracy data were provided")
        decision = accuracy * stimulus + (1 - accuracy) * (1 - stimulus)
    }else{
        decision = as.numeric(as.factor(as.character(decision))) - 1
        if(!is.null(accuracy))
            if(!all(accuracy = (decision == stimulus)))
                stop("Accuracy does not reflect the agrreement between decision and stimulus variables. Check if accuracy = 1 (0) when decision = (!=) stimulus.")
    }
    if(is.null(rating)){
        decision + 1
    }else{
        if(!all(diff(sort(unique(rating))) == 1))
            stop(sprintf('Not all ratings between min(rating) = %d and max(rating) = %d are present.', min(rating, na.rm = T), max(rating, na.rm = T)))
        rating = rating - min(rating, na.rm = T) + 1
        max.rating = max(rating, na.rm = T)
        (1 - decision) * (max.rating + 1 - rating) + decision * (max.rating + rating)
    }
}
