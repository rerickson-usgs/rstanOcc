#' Probability of detection one
#'
#' This function calculates the probability
#' of detecting eDNA in at least one sample
#' assuming eDNA is present at a site.
#'
#' 
#'
#' @param J number of samples per sampleing event
#' @param K number of molecular replicates
#' @param theta sample-level capture probability
#' @param p_detection molecular-level detection probability 
#' @examples
#' sample_detection_one()
#' 
#' 
#' @export 
sample_detection_one <-
    function(
             J = 50,
             K = 8,
             theta = 0.06,
             p_detection = 0.3
             ){
        
        j_index = J:0
	prob =
            sum( choose( J, j_index) * (1 - theta) ^ j_index * (theta * (1 - p_detection)^K )  ^ rev(j_index))       
	return(1.0 - prob)
    }
