#' Get plotting theme
#'
#' This function gets the plotting theme for ggplot.
#'
#' @param textsize numeric value for base text size on ggplot. Default is 9.
#'
get_plottheme <- function(textsize){
    mytheme <- theme_bw() + theme(axis.text = element_text(size = textsize),
                                  axis.title = element_text(size = textsize + 2),
                                  legend.text = element_text(size = textsize),
                                  legend.title = element_text(size = textsize + 2),
                                  strip.text.x = element_text(size = textsize))
    return(mytheme)
}

#' Plot data
#'
#' This function plots raw, filtered, or simulated data.
#'
#' @param data data frame of raw, filtered, or simulated data. Must include columns for 'id', 'vl', and 'time'.
#' @param textsize numeric value for base text size on ggplot. Default is 9.
#' @param pointsize numeric value for base point size on ggplot. Default is 1.
#' @param linesize numeric value for line width on ggplot. Default is 0.5.
#' @param facet_col numeric value for number of columns to use when plotting subject panels. Defaults to NULL (i.e. ggplot default).
#' @param detection_threshold numeric value indicating the detection threshold of the assay used to measure viral load. Default value is 20.
#' @import ggplot2
#' @importFrom rlang .data
#' @export
#' @examples
#'
#' set.seed(1234567)
#'
#' simulated_data <- simulate_data(nsubjects = 20)
#'
#' model_output <- ushr(data = simulated_data, detection_threshold = 10)
#'
#' plot_data(simulated_data, detection_threshold = 10)
#'
plot_data <- function(data, textsize = 9, pointsize = 1, linesize = 0.5,
                      facet_col = NULL, detection_threshold = 20){

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package \"ggplot2\" is required for automated plotting.
             Either install it, or plot manually.")
    }

    mytheme <- get_plottheme(textsize)

    data %>% ggplot(aes(x = .data$time, y = .data$vl)) + geom_point(size = pointsize) +
        geom_hline(aes(yintercept = detection_threshold), size = linesize, linetype = "dashed") +
        facet_wrap(~ id, ncol = facet_col) + mytheme +
        scale_y_log10("HIV viral load") + scale_x_continuous("Time")
}



#' Plot model fits
#'
#' This function plots the output from model fitting.
#'
#' @param model_output output from model fitting using ushr().
#' @param type character string indicating whether the biphasic or single phase fits should be plotted. Must be either "biphasic" or "single". Defaults to "biphasic".
#' @param detection_threshold numeric value indicating the detection threshold of the assay used to measure viral load. Default value is 20.
#' @param textsize numeric value for base text size on ggplot. Default is 9.
#' @param pointsize numeric value for base point size on ggplot. Default is 1.
#' @param linesize numeric value for line width on ggplot. Default is 0.5.
#' @param facet_col numeric value for number of columns to use when plotting subject panels. Defaults to NULL (i.e. ggplot default).
#' @export
#' @examples
#'
#' set.seed(1234567)
#'
#' simulated_data <- simulate_data(nsubjects = 20)
#'
#' model_output <- ushr(data = simulated_data, detection_threshold = 20)
#'
#' plot_model(model_output, type = "biphasic", detection_threshold = 20)
#'
plot_model <- function(model_output, type = "biphasic", detection_threshold = 20,
                          textsize = 9, pointsize = 1, linesize = 0.5,
                          facet_col = NULL){

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package \"ggplot2\" is required for automated plotting.
             Either install it, or plot manually.")
    }
    if (is.null(model_output)) {
        stop("Please include model output from 'ushr()' for plotting. If you want to plot data, use 'plot_data()'.")
    }
    # get desired fits for plotting
    if (type == "biphasic") {
        fits <- model_output$biphasic_fits
    } else {
        fits <- model_output$single_fits
    }

    if (nrow(fits) == 0) {
        stop("There are no fits of the type you have chosen. Try plotting the other type of fit ('biphasic', or 'single').")
    }

    filtered_data <- model_output$data_filtered %>% filter(id %in% unique(fits$id))

    mytheme <- get_plottheme(textsize)

    filtered_data %>% ggplot() + geom_point(aes(x = time, y = vl, group = id), size = pointsize) +
        geom_line(data = fits, aes(x = time, y = fit, group = id), lty = 1, col = "black", size = linesize) +
        facet_wrap(~id, ncol = facet_col) + mytheme +
        xlab("Time") + scale_y_log10("HIV viral load") +
        geom_hline(aes(yintercept = detection_threshold), linetype = "dashed")
}


#' Plot time to suppression distribution
#'
#' This function plots a histogram of the time to suppression estimates.
#'
#' @param TTS_output output from estimating time to suppression (TTS) values using get_TTS()..
#' @param textsize numeric value for base text size on ggplot. Default is 9.
#' @param bins numeric value indicating the number of bins for the histogram. Default is 20.
#' @export
#' @examples
#'
#' set.seed(1234567)
#'
#' simulated_data <- simulate_data(nsubjects = 20)
#'
#' TTSestimates <- get_TTS(data = simulated_data, suppression_threshold = 10, parametric = FALSE)
#'
#' plot_TTS(TTSestimates, bins = 5)
#'
plot_TTS <- function(TTS_output, textsize = 9, bins = 20){
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package \"ggplot2\" is required for automated plotting.
             Either install it, or plot manually.")
    }
    mytheme <- get_plottheme(textsize)

    ggplot(data = TTS_output, aes(x = TTS)) +
        geom_histogram(bins = bins, fill = "grey", colour = "black") +
        mytheme + ylab("Frequency") + xlab("Time to suppression")
}


#' Summarize model output
#'
#' This function summarizes the output of model fitting..
#'
#' @param model_output output from model fitting using ushr().
#' @param data dataframe of original data used for model fitting. Must include named 'id' column as a subject identifier
#' @param stats logical indicator: should the median and sd lifespans also be returned? Default is FALSE.
#' @import dplyr
#' @import tidyr
#' @importFrom rlang .data
#' @export
#' @examples
#'
#' set.seed(1234567)
#'
#' simulated_data <- simulate_data(nsubjects = 20)
#'
#' model_output <- ushr(data = simulated_data, detection_threshold = 10)
#'
#' summarize_model(model_output, data = simulated_data)
#'
summarize_model <- function(model_output, data, stats = FALSE){

    # get biphasic parameter estimates
    biphasicfits <- model_output$biphasicCI %>%
        select(- .data$lowerCI, - .data$upperCI) %>% spread(.data$param, .data$estimate) %>%
        mutate(ShortLifespan = 1/.data$delta, LongLifespan = 1/.data$gamma, Model = "Biphasic")

    singlephasefits <- model_output$singleCI  %>%
        select(- .data$lowerCI, - .data$upperCI) %>% spread(.data$param, .data$estimate) %>%
        mutate(SingleLifespan = 1/.data$gammahat, Model = "Single phase")

    allfits <- biphasicfits %>% full_join(singlephasefits)

    allinfo <- data %>% distinct(id) %>%
        mutate(Included = ifelse(id %in% allfits$id, "Yes", "No")) %>%
        left_join(allfits) %>%
        select(.data$id, .data$Included, #Reason = ExclusionReason,
               .data$Model, .data$ShortLifespan, .data$LongLifespan, .data$SingleLifespan) %>%
        mutate_if(is.numeric, round, 2) %>% replace(., is.na(.), "")


    biphasicstats <- model_output$biphasicCI %>%
        select(- .data$lowerCI, - .data$upperCI) %>% spread(.data$param, .data$estimate) %>%
        mutate(ShortLifespan = 1/.data$delta, LongLifespan = 1/.data$gamma) %>%
        gather(Param, estimate, A:LongLifespan) %>%
        group_by(.data$Param) %>%
        summarize(Median = median(.data$estimate), SD = sd(.data$estimate)) %>%
        mutate(Median = signif(.data$Median, 3), SD = signif(.data$SD, 3), Model = "Biphasic")


   singlestats <- model_output$singleCI %>%
        select(- .data$lowerCI, - .data$upperCI) %>% spread(.data$param, .data$estimate) %>%
        mutate(SingleLifespan = 1/.data$gammahat) %>%
        gather(Param, estimate, Bhat:SingleLifespan) %>%
        group_by(.data$Param) %>%
        summarize(Median = median(.data$estimate), SD = sd(.data$estimate)) %>%
        ungroup() %>%
        mutate(Median = signif(.data$Median, 3), SD = signif(.data$SD, 3), Model = "Single phase")


    if (stats) {
        return(list(summary = allinfo, biphasicstats = biphasicstats, singlestats = singlestats))
    } else {
        return(allinfo)
    }
}
