#' Prepare input data
#'
#' This function prepares the raw input data for model fitting. All default values are taken from Morris et al (2019).
#'
#' Steps include:
#' 1. Setting values below the detection threshold to half the detection threshhold (following standard practice).
#' 2. Filtering out subjects who do not suppress viral load below the detection threshold by a certain time.
#' 3. Filtering out subjects who do not have a decreasing sequence of viral load (within some buffer range).
#' 4. Filtering out subjects who do not have enough data for model fitting.
#' 5. Removing last data point of subjects with last two points very close to the detection threshold. This prevents skewing of the model fit
#' @param data raw data set. Must have the following column names: 'id' stating the unique identifier for each subject, 'vl' stating the viral load measurements for each subject; 'time' stating the time at which each measurement was taken.
#' @param detection_threshold detection threshold of the assay used to measure viral load. Measurements below this value will be assumed to represent undetectable viral load levels. Default value is 20.
#' @param censortime the maximum time point to inculde in the analysis. Subjects who do not suppress viral load below the detection threshold within this time will be discarded from model fitting. Default value is 365.
#' @param decline_buffer the maximum allowable deviation of values away from a strictly decreasing sequence in viral load. This allows for e.g. measurement noise and small fluctuations in viral load. Default value is 500.
#' @param n_min_single the minimum number of data points required to be included in the analysis. Defaults to 3. It is highly advised not to go below this threshold.
#' @param threshold_buffer numerical value indicating the range above the detection threshold which represents potential skewing of model fits. Subjects with their last two data points within this range will have the last point removed. Default value is 10.
#' @param nsuppression numerical value (1 or 2) indicating whether suppression is defined as having one observation below the detection threshold, or two sustained observations. Default value is 1.
#' @import dplyr
#' @export
#' @examples
#'
#' set.seed(1234567)
#'
#' simulated_data <- simulate_data(nsubjects = 20)
#'
#' filter_data(simulated_data)
#'
filter_data <- function(data, detection_threshold = 20, censortime = 365, decline_buffer = 500,
                        n_min_single = 3, threshold_buffer = 10, nsuppression = 1){

    # Check that data frame includes columns for 'id', 'time', 'vl'
    if (!(all(c("vl", "time", "id") %in% names(data)))) {
        stop("'data' must be a data frame with named columns for 'id', 'time', and 'vl'")
    }

    if (!is.numeric(data$time)) {
        stop("Column for the time of observations ('time') must be numeric")
    }

    if (is.factor(data$id)) {
        data$id <- as.character(data$id)
    }

    if (any(is.na(data$id))) {
        warning("Some subjects have missing IDs; removing these from the data")

        data <- data %>% filter(!is.na(id))
    }


    if (!(nsuppression %in% c(1,2))) {
        warning("nsuppression must take the numeric value 1 or 2 to define the criteria for reaching suppression; reverting to default nsuppression = 1")

        nsuppression <- 1
    }

    if (nsuppression == 1) {
        # 1. Change everything <= detection_threshhold to 1/2 * detection_threshhold
        data_filtered <- data %>% mutate(vl = case_when(vl <= detection_threshold ~ detection_threshold/2,
                                                        vl >= detection_threshold ~ vl) ) %>%
            # 2. Look at only those who reach control within user defined censortime
            filter(time <= censortime) %>% group_by(id) %>% filter(any(vl <= detection_threshold)) %>% ungroup() %>%
            # 3a. Isolate data from the highest VL measurement (from points 1 - 3) to the first point below detection
            filter(!is.na(vl)) %>% group_by(id) %>%
            slice(which.max(vl[1:3]):Position(function(x) x <= detection_threshold, vl)) %>% ungroup() %>%
            # 3b. Only keep VL sequences that are decreasing with user defined buffer...
            group_by(id) %>% filter(all(vl <= cummin(vl) + decline_buffer)) %>%
            # 4. ...AND have min # dps above the detection threshold
            filter(length(vl[vl > detection_threshold]) >= n_min_single)
    } else if (nsuppression == 2) {
        data_filtered <- data %>% mutate(vl = case_when(vl <= detection_threshold ~ detection_threshold/2,
                                                        vl >= detection_threshold ~ vl) ) %>%
            filter(!is.na(vl)) %>%
            filter(time <= censortime) %>% group_by(id) %>%
            # NOW: must have 2 consecutive measurements below threshold
            mutate(firstbelow = intersect(which(vl <= detection_threshold), which(vl <= detection_threshold) + 1)[1] - 1 ) %>%
            mutate(firstbelow = time[firstbelow]) %>% filter(time <= firstbelow) %>% ungroup() %>%
            # Continue as above
            group_by(id) %>%
            slice(which.max(vl[1:3]):Position(function(x) x <= detection_threshold, vl)) %>% ungroup() %>%
            group_by(id) %>% filter(all(vl <= cummin(vl) + decline_buffer)) %>%
            filter(length(vl[vl > detection_threshold]) >= n_min_single)
    }

    # 5. Remove last data point of subjects with last two points very close to the threshold to prevent skewing model fit
    if (nrow(data_filtered) > 0) {
        data_filtered <- data_filtered %>% group_by(id) %>%
            mutate(n = n(), index = 1:n(), tag = ifelse(vl[n-1] - vl[n] < threshold_buffer, TRUE, FALSE) ) %>%
        filter(!(tag == TRUE & index == n)) %>% ungroup() %>% select(-index, -n, -tag)
}
    return(data_filtered)
}
