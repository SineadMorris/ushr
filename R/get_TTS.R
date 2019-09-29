#' Prepare input data for non-parametric TTS calculations.
#'
#' This function prepares the raw input data for TTS interpolation. All default values are taken from Morris et al (2019).
#'
#' Steps include:
#' 1. Setting values below the suppression threshold to half the suppression threshhold (following standard practice).
#' 2. Filtering out subjects who do not suppress viral load below the suppression threshold by a certain time.
#' 3. Filtering out subjects who do not have a decreasing sequence of viral load (within some buffer range).
#' @param data raw data set. Must have the following column names: 'id' stating the unique identifier for each subject, 'vl' stating the viral load measurements for each subject; 'time' stating the time at which each measurement was taken.
#' @param suppression_threshold numeric value indicating the suppression threshold: measurements below this value will be assumed to represent viral suppression. Default value is 20.
#' @param censortime the maximum time point to inculde in the analysis. Subjects who do not suppress viral load below the detection threshold within this time will be discarded from model fitting. Default value is 365.
#' @param decline_buffer the maximum allowable deviation of values away from a strictly decreasing sequence in viral load. This allows for e.g. measurement noise and small fluctuations in viral load. Default value is 500.
#' @importFrom rlang .data
#' @export
#' @examples
#'
#' set.seed(1234567)
#'
#' simulated_data <- simulate_data(nsubjects = 20)
#'
#' filter_dataTTS(data = simulated_data, suppression_threshold = 10)
#'
filter_dataTTS <- function(data, suppression_threshold = 20,
                           censortime = 365, decline_buffer = 500){

    # Check that data frame includes columns for 'id', 'time', 'vl'
    if (!(all(c("vl", "time", "id") %in% names(data)))) {
        stop("Data frame must have named columns for 'id', 'time', and 'vl'")
    }

    # 1. Change everything <= suppression_threshhold to 1/2 * dsuppression_threshhold
    data_filtered <- data %>% mutate(vl = case_when(vl <= suppression_threshold ~ suppression_threshold/2,
                                                    vl >= suppression_threshold ~ vl) ) %>%
        # 2. Look at only those who reach control within user defined censortime
        filter(.data$time <= censortime) %>% group_by(id) %>%
        filter(any(.data$vl <= suppression_threshold)) %>% ungroup() %>%
        # 3a. Isolate data from the highest VL measurement (from points 1 - 3) to the first point below detection
        filter(!is.na(.data$vl)) %>% group_by(id) %>%
        slice(which.max(.data$vl[1:3]):Position(function(x) x <= suppression_threshold, .data$vl)) %>%
        ungroup() %>%
        # 3b. Only keep VL sequences that are decreasing with user defined buffer...
        group_by(id) %>% filter(all(.data$vl <= cummin(.data$vl) + decline_buffer))

    return(data_filtered)
}


#' Biphasic root function
#'
#' This function expressed the root of the biphasic suppression equation, V(t) = suppression_threshold.
#'
#' @param timevec vector of the times, t, at which V(t) should be calculated
#' @param params named vector of all parameters needed to compute the biphasic model, V(t)
#' @param suppression_threshold suppression threshold: measurements below this value will be assumed to represent viral suppression. Typically this would be the detection threshold of the assay. Default value is 20.
#' @export
#'
biphasic_root <- function(timevec, params, suppression_threshold){
    value <- params["A"] * exp (- timevec * params["delta"]) + params["B"] * exp( - timevec * params["gamma"]) - suppression_threshold
    as.numeric(value)
}


#' Single phaseroot function
#'
#' This function expressed the root of the single phase suppression equation, V(t) = suppression_threshold.
#'
#' @param timevec vector of the times, t, at which V(t) should be calculated
#' @param params named vector of all parameters needed to compute the biphasic model, V(t)
#' @param suppression_threshold suppression threshold: measurements below this value will be assumed to represent viral suppression. Typically this would be the detection threshold of the assay. Default value is 20.
#' @export
#'
single_root <- function(timevec, params, suppression_threshold){
    if (all(c("B", "gamma") %in% params)) {
        value <- params["B"] * exp( - timevec * params["gamma"]) - suppression_threshold
    } else{
        value <- params["Bhat"] * exp( - timevec * params["gammahat"]) - suppression_threshold
    }

    as.numeric(value)
}


#' Parametric TTS function
#'
#' This function computes the parametric form of the time to suppression
#'
#' @param params named vector of all parameters needed to compute the suppression model, V(t)
#' @param rootfunction specifies which function should be used to calculate the root: biphasic or single phase.
#' @param suppression_threshold suppression threshold: measurements below this value will be assumed to represent viral suppression. Typically this would be the detection threshold of the assay. Default value is 20.
#' @param uppertime the maximum time interval to search for the time to suppression. Default value is 365.
#' @export
#'
get_parametricTTS <- function(params, rootfunction, suppression_threshold, uppertime){
    TTS <- rep(NA, nrow(params))

    for (i in 1:nrow(params)){
        TTS[i] = stats::uniroot(rootfunction, lower = 1, upper = uppertime,
                         params = params[i,], suppression_threshold = suppression_threshold)$root
    }
    return(TTS)
}


#' Non-parametric TTS function
#'
#' This function computes the non-parametric form of the time to suppression
#'
#' @param vl vector of viral load measurements.
#' @param suppression_threshold numeric value for the suppression threshold: measurements below this value will be assumed to represent viral suppression. Typically this would be the detection threshold of the assay. Default value is 20.
#' @param time voector indicating the time when vl measurements were taken.
#' @param npoints numeric value of the number of interpolation points to be considered.
#' @export
#'
get_nonparametricTTS <- function(vl, suppression_threshold, time, npoints){

    TTS <- time[which(vl == suppression_threshold)[1]]
    firstbelow <- which(vl < suppression_threshold)[1]

    if(is.na(TTS) | (!is.na(time[firstbelow]) & (time[firstbelow] < TTS)) ){
        lastabove <- time[firstbelow - 1]

        yax <- c(vl[firstbelow - 1], vl[firstbelow])
        xax <- c(lastabove, time[firstbelow])

        interpolation <- stats::approx(xax, yax, n = npoints)
        TTS <- interpolation$x[interpolation$y <= suppression_threshold][1]
    }
    return(TTS)

}


#' Time to suppression (TTS) function
#'
#' This function calculates the time to suppress HIV below a specified threshold.
#'
#' Options include: parametric (i.e. using fitted model) or non-parametric (i.e. interpolating the processed data).
#' @param model_output output from fitting model. Only required if parametric = TRUE.
#' @param data raw data set. Must have the following column names: 'id' stating the unique identifier for each subject, 'vl' stating the viral load measurements for each subject; 'time' stating the time at which each measurement was taken. Only required if parametric = FALSE.
#' @param suppression_threshold suppression threshold: measurements below this value will be assumed to represent viral suppression. Typically this would be the detection threshold of the assay. Default value is 20.
#' @param uppertime the maximum time interval to search for the time to suppression. Default value is 365.
#' @param decline_buffer the maximum allowable deviation of values away from a strictly decreasing sequence in viral load. This allows for e.g. measurement noise and small fluctuations in viral load. Default value is 500.
#' @param parametric logical TRUE/FALSE indicating whether time to suppression shoudl be calculated usingthe parametric (TRUE) or non-parametric (FALSE) method. If TRUE, a fitted model object is required. If FALSE, the raw data frame is required. Defaults to TRUE.
#' @param ARTstart logical TRUE/FALSE indicating whether the time to suppression should be represented as time since ART initiation. Default = FALSE. If TRUE, ART initiation times must be included as a data column named 'ART'.
#' @param npoints numeric value of the number of interpolation points to be considered. Default is 1000.
#' @importFrom rlang .data
#' @export
#' @examples
#'
#' set.seed(1234567)
#'
#' simulated_data <- simulate_data(nsubjects = 20)
#'
#' get_TTS(data = simulated_data, suppression_threshold = 10, parametric = FALSE)
#'
get_TTS <- function(model_output = NULL, data = NULL,
                    suppression_threshold = 20, uppertime = 365, decline_buffer = 500,
                    parametric = TRUE, ARTstart = FALSE, npoints = 1000){

    # 1. Parametric TTS ----------------------------------------------------------------
    if(parametric == TRUE){

        if(is.null(model_output)){
            stop("Model output not found. You must supply the fitted model to calculate parametric TTS values. Try ?get_model_fits.")
        }

        if (length(model_output$biphasicCI) > 0 & length(model_output$singleCI) > 0) {
            # Biphasic
            biphasic_params <- model_output$biphasicCI %>%
                select(-.data$lowerCI, -.data$upperCI) %>% spread(.data$param, .data$estimate) %>%
                mutate(TTS = get_parametricTTS(params = ., rootfunction = biphasic_root, suppression_threshold, uppertime),
                       model = "biphasic", calculation = "parametric")

            # Single phase
            single_params <- model_output$singleCI %>%
                select(-.data$lowerCI, -.data$upperCI) %>% spread(.data$param, .data$estimate) %>%
                mutate(TTS = get_parametricTTS(params = ., rootfunction = single_root, suppression_threshold, uppertime),
                       model = "single phase", calculation = "parametric")

            # All
            TTS_output <- biphasic_params %>% full_join(single_params) %>%
                select(.data$id, .data$TTS, .data$model, .data$calculation)

        } else if (length(model_output$biphasicCI) > 0 & length(model_output$singleCI) == 0) {
            # Biphasic
            biphasic_params <- model_output$biphasicCI %>%
                select(-.data$lowerCI, -.data$upperCI) %>% spread(.data$param, .data$estimate) %>%
                mutate(TTS = get_parametricTTS(params = ., rootfunction = biphasic_root, suppression_threshold, uppertime),
                       model = "biphasic", calculation = "parametric")

            # All
            TTS_output <- biphasic_params %>% select(.data$id, .data$TTS, .data$model, .data$calculation)

        } else if (length(model_output$biphasicCI) == 0 & length(model_output$singleCI) > 0) {
            # Single phase
            single_params <- model_output$singleCI %>%
                select(-.data$lowerCI, -.data$upperCI) %>% spread(.data$param, .data$estimate) %>%
                mutate(TTS = get_parametricTTS(params = ., rootfunction = single_root, suppression_threshold, uppertime),
                       model = "single phase", calculation = "parametric")

            # All
            TTS_output <- single_params %>% select(.data$id, .data$TTS, .data$model, .data$calculation)
        }

    }

    # 2. Non-parametric TTS ----------------------------------------------------------------
    if(parametric == FALSE){

        if(is.null(data)){
            stop("Data not found. You must supply the data to calculate non-parametric TTS values")
        }

        # Check dataframe includes columns for 'id', 'time', 'vl'
        if(!(all(c("vl", "time", "id") %in% names(data)))){
            stop("Data frame must have named columns for 'id', 'time', and 'vl'")
        }

        # Filter out subjects to focus on those who reach suppression below the specified threshold.
        data_filtered <- filter_dataTTS(data, suppression_threshold, uppertime, decline_buffer)

        TTS_output <- data_filtered %>%
            mutate(TTS = get_nonparametricTTS(.data$vl, suppression_threshold, .data$time, npoints)) %>%
            ungroup() %>% distinct(.data$id, .data$TTS) %>% mutate(calculation = "non-parametric")
    }

    if(ARTstart == TRUE){
        print("Calculating TTS as time since ART initiation...")

        if(is.null(data$ART)){
            print("Data frame is missing ART column. Returning original TTS values.")
        } else {
            ARTdata <- data %>% distint(.data$id, .data$ART)

            TTS_output <- TTS_output %>% left_join(ARTdata) %>% mutate(TTS = .data$TTS - .data$ART)
        }
    }
    return(TTS_output)
}
