# Load a series of convenient & helpful packages
require(ggplot2)  # for pretty multivariate plotting
require(xcms)  # the main package for metabolomics data analysis
require(xlsx)  # for reading xlsx files
require(stringr)  # for easy consistent handling of strings and regular expressions in R
require(repr)  # for resizing plots in the Jupyter notebook
require(plyr)  # for use of the . function
require(viridis)  # for nice colormaps
require(MASS)
require(dplyr)  # for easy dataframe manipulation, must be imported after plyr
require(tidyr)  # for easy dataframe manipulation
require(magrittr)  # for the %<>% operator
require(reshape2)  # for more difficult and confusing dataframe manipulation
require(ecipex)  # calculates isotopic fine structure & exact masses for chemical formulae
require(roxygen2)  # automatically generate documentation


#' Find median and maximum intensities of features in a data frame output by 
#' \code{xcms::diffreport} or \code{xcms::peakTable()}.
#'
#' @param xtab a data frame, usually output by \code{xcms::diffreport} or \code{xcms::peakTable()}.
#' @param col_vector (optional) a vector of integers pointing to the columns of xtab with intensity values.
#' @param log.scale a boolean, if TRUE, then returned values will be log10 transformed.
#' 
#' @export
find_row_medians_and_maxes <- function(xtab, col_vector, log.scale=T){
    # TODO: implement default choosing of col_vector    
    data <- xtab[, col_vector]
    if(log.scale == T){
        xtab$median_intensity <- log10(apply(data, 1, FUN=median))
        xtab$max_intensity <- log10(apply(data, 1, FUN=max))
        return(xtab)
        } else{
            xtab$median_intensity <- apply(data, 1, FUN=median)
            xtab$max_intensity <- apply(data, 1, FUN=max)
            return(xtab)
              }
    }


#' Find the mass defect of m/z values in a data frame that usually arises from  
#' \code{xcms::diffreport} or \code{xcms::peakTable()}.
#'
#' @param xtab a data frame, usually output by \code{xcms::diffreport} or \code{xcms::peakTable()}.
#' @col_vector (optional) a vector of integers pointing to the columns of xtab with intensity values.
#' @log.scale a boolean, if TRUE, then returned values will be log10 transformed.
#' 
#' @export
find_mz_defect <- function(xtab){
    if(!any(names(xtab) %in% c('mz', 'mzmed'))){
        cat('Error, mz column not found.')
        return()
        }
    possible_names <- c('mz', 'mzmed')
    col_name <- possible_names[possible_names %in% names(xtab)]
    data <- xtab[, col_name]
    integer_mz <- round(data)
    decimal_mz <- data - integer_mz
    xtab$mz_defect <- decimal_mz
    return(xtab)
    }
    

#' Retrieve a single TIC from a single \code{xcms}-readable LC-MS data file.
#'
#' @param file path & filename for a \code{xcms}-readable LC-MS data file (e.g. *.mzML).
#' @param rtcor whether to use corrected RT data or raw RT data.
#' @return a dataframe containing RT and intensity information.
#' 
#' @export
get_tic_from_file <- function(file, rtcor=NULL) {
     object <- xcmsRaw(file)
     cbind(if (is.null(rtcor)) object@scantime else 
     		  rtcor, rawEIC(object, mzrange=range(object@env$mz))$intensity) 
}

#' Retrieve many TICs and return a long-format data frame
#' 
#'  Based on code from http://metabolomics-forum.com/viewtopic.php?f=26&t=122
#' @param xcmsSet the \code{xcms::xcmsSet} object from which to get TICs.
#' @param files the full-path filenames from which to get TICs (redundant if \code{\xcmsSet} is NULL).
#' @param rt whether to use raw or corrected RT data.
#' @return a dataframe containing RT, intensity, and file information.
#' 
#' @export
get_tics_from_files <- function(xcmsSet=NULL, files=NULL, rt=c("raw", "corrected")) {
  if (is.null(xcmsSet)) {
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                      "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
    if (is.null(files))
        files <- getwd()
    info <- file.info(files)
    listed <- list.files(files[info$isdir], pattern = filepattern,
                          recursive = TRUE, full.names = TRUE)
    files <- c(files[!info$isdir], listed)
  } else {
    files <- filepaths(xcmsSet)
  }

  N <- length(files)
  TIC <- vector("list",N)
  
      
  # initialize dataframe
  TIC_df <- data.frame(rt=c(), intensity=c(), file=c())
  for (i in 1:N) {
      cat(files[i],"\n")
      if (!is.null(xcmsSet) && rt == "corrected")
        rtcor <- xcmsSet@rt$corrected[[i]] else 
          rtcor <- NULL
      TIC[[i]] <- getTIC(files[i], rtcor=rtcor)
      this_file = rep(files[i], length(TIC[[i]]))
      TIC_df <- rbind(TIC_df, data.frame(rt=TIC[[i]][, 1], intensity=TIC[[i]][, 2], file=this_file))
      # 
  }

  return(TIC_df)
}


#' A dispatching function for \code{getEICdf}. 
#' 
#' @export
getEICdf <- function(xobj, ...){
	type <- class(xobj)
	if(!(type %in% c('xcmsSet', 'xcmsRaw'))){
		cat('Unknown object type.')
		return()
		}
	if(type == 'xcmsSet'){
		return(getEICdf.set(xobj, ...))
		}
	
	if(type == 'xcmsRaw'){
		return(getEICdf.raw(xobj, ...))
		}
}

#' A wrapper around the \code{xcms::getEIC()} function to return results in a data frame.
#'
#' @param xset An S4 object as returned by \code{xcms::xcmsSet()}.
#' @param ... other parameters passed through to \code{xcms::getEIC)}.
#' @return A \code{data.frame} with columns 'rt', 'intensity', 'sample', 'mzrange.1',
#'										    'mzrange.2', 'groupname', 'class', 
#' 											'relative_rt' 
#' 
#' @export
getEICdf.set <- function(xset, ...)
                    {
                    xeic <- getEIC(xset, ...)
    
                    df_out <- data.frame(rt=c(), intensity=c(), sample=c(), mzrange=c(), 
                                         groupname=c(), class=c(), relative_rt=c())
                    for(sample in names(xeic@eic))
                        {
                        current_class <- sampclass(xset)[sampnames(xset) == sample]
                        for(idx in 1:length(xeic@groupnames))
                            {
                            groupname = xeic@groupnames[idx]
                            mz_range = c(xeic@mzrange[idx, 'mzmin'], xeic@mzrange[idx, 'mzmax'])
                            df_size = length(xeic@eic[[sample]][[idx]])
                            this_sample <- rep(sample, df_size)
                            this_mz_range <- matrix(rep(mz_range, df_size), ncol=2, byrow=T)
                            this_groupname <- rep(groupname, df_size)
                            this_class <- rep(current_class, df_size)
                            this_rt <- xeic@eic[[sample]][[idx]][, 1]
                            this_relative_rt <- this_rt - median(this_rt)
                            this_df <- data.frame(rt=this_rt, 
                                                  intensity=xeic@eic[[sample]][[idx]][, 2], 
                                                  sample=this_sample,
                                                  mzrange=this_mz_range,
                                                  groupname=this_groupname,
                                                  class=this_class,
                                                  relative_rt=this_relative_rt
                                                  )
                            
                            df_out <- rbind(df_out, this_df)
                            }
                        }
    
                    options(repr.plot.width=6, repr.plot.height=24)

                    p <- ggplot(data=df_out, aes(x=relative_rt, y=intensity, color=class, group=sample)) + 
                        geom_line() + 
                        facet_grid(groupname~., scale='free') + 
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 90))
            
                    return(list(eic_data=df_out, plot=p))
                    }
                    
#' A wrapper around the \code{xcms::getEIC()} function to return results in a data frame.
#'
#' @param xraw An S4 object as returned by \code{xcms::xcmsRaw()}.
#' @param ... other parameters passed through to \code{xcms::getEIC()}.
#' @return A \code{data.frame} with columns 'rt', 'intensity', 'sample', 'mzrange.1',
#'										    'mzrange.2', 'groupname', 'class', 
#' 											'relative_rt' 
#'
#' This version handles supplying \code{mzrange} and \code{rtrange}.
#' 
#' @export
getEICdf.raw <- function(xraw, ...){
                    xeic <- getEIC(xraw, ...)
    
                    df_out <- data.frame(rt=c(), intensity=c(), sample=c(), mzrange=c())
                    for(idx in 1:length(xeic@eic[[1]]))
                        {
                        mz_range = c(xeic@mzrange[idx, 1], xeic@mzrange[idx, 2])

                        df_size = length(xeic@eic[[1]][[idx]])
                        this_sample <- rep(xraw@filepath, df_size)
                        this_mz_range <- matrix(rep(mz_range, df_size), ncol=2, byrow=T)

                        this_rt <- xeic@eic[[1]][[idx]][, 1]
                        intensity <- xeic@eic[[1]][[idx]][, 2]

                        this_df <- data.frame(rt=this_rt, 
                                              intensity=xeic@eic[[1]][[idx]][, 2], 
                                              sample=this_sample,
                                              mzrange=this_mz_range
                                              )

                        df_out <- rbind(df_out, this_df)
                    }
                    df_out$mz <- as.factor((df_out$mzrange.1 + df_out$mzrange.2)/2)
                    return(df_out)
                }

          
          
          
#' Call \code{xcms::diffreport()} and make a plot of the results.
#'
#' @param xset The xcmsSet object to analyze.
#' @param class1 The first class to include in the comparison, 
#'    			 must be a character found in \code{sampclass(xcms)}.
#' @param class2 The second class to include in the comparison.
#' @param eicwidth width (in seconds) of calculated EICs.
#' @param max_intensity_cutoff log10-scaled intensity threshold below which features will not be plotted.
#' @param label_top_n the number n of the top features by _both_ p-value and _fold-change_ to number.
#' @param plot_height the plot height
#' @param plot_width the plot width
#' @param rt_tol, the allowed difference between retention times for the same feature, in seconds.
#' @param mz_tol, the allowed difference between m/z values for the same feature, in Da.
#' @param rt_range, the allowed time window in which features must fall.  Useful for avoiding features in flow-through or column reequilibration stage.
#' @return A named list with item \code{'plot'} and \code{'reporttab'}
#' 
#' @export
do_diffreport_and_plot <- function(xset, 
                                   class1, 
                                   class2, 
                                   eicwidth=60,
                                   max_intensity_cutoff=5,
                                   label_top_n=15,
                                   plot_height=6,
                                   plot_width=6,
                                   rt_tol=12,
                                   mz_tol=0.05,
                                   rt_range=c(0, 3000),
                                   ...){
    
    # make report and add desired max/median/mass_defect annotations
    reporttab <- diffreport(xset, class1=class1, class2=class2, eicwidth=eicwidth, ...)
    data_cols <- names(reporttab) %in% sampnames(xset)
    reporttab <- find_row_medians_and_maxes(reporttab, col_vector = data_cols)
    reporttab <- find_mz_defect(reporttab)
    
    # remove "missing" pvalues/fold-changes that are actually infinite because of 0 intensity in one sample
    reporttab[is.na(reporttab$pvalue), 'pvalue'] <- min(reporttab$pvalue, na.rm=T)
    reporttab[is.na(reporttab$fold), 'fold'] <- max(reporttab$fold, na.rm=T)
    
    
    # prepare for plotting
    options(repr.plot.width=8, repr.plot.height=8)
    
    # make labels for all points on plot 
    reporttab$lbl <- paste(sprintf('%.4f', reporttab$mzmed),
                           sprintf('%s', reporttab$name), 
                           sep='\n'
                          )

    # plot only features above a maximum intensity cutoff, with well-conserved rt and mz, and in good rt range
    keeps <- reporttab$max_intensity >= max_intensity_cutoff &
                reporttab$rtmax - reporttab$rtmin < rt_tol &
                reporttab$mzmax - reporttab$mzmin < mz_tol &
                findInterval(reporttab$rtmed, rt_range) == 1

    # keep labels only for top n features by fold change or by p-value
    fold_cutoff <- log10(reporttab$fold[order(-reporttab$fold)])[keeps][label_top_n]
    pval_cutoff <- -log10(reporttab$pvalue[order(reporttab$pvalue)])[keeps][label_top_n]
    discards <- (-log10(reporttab$pvalue) < pval_cutoff & 
                log10(reporttab$fold) < fold_cutoff ) 
                        
    reporttab[discards, 'lbl'] <- ''
    
    p <- ggplot(data=reporttab[keeps, ], aes(x=-log10(pvalue), y=log10(fold), size=max_intensity, color=rtmed)) + 
        geom_point() + 
        theme_bw() +
        scale_color_gradientn(colours=viridis(20)) + 
        geom_text(data=reporttab[keeps & !discards, ], 
                   aes(label=lbl), size=3.5, color='black', alpha=0.9, check_overlap=T)
                      

    reporttab_filtered <- reporttab[keeps & !discards, ]
    return(list(reporttab=reporttab_filtered, plot=p))
}



#' Find and retain values from a long query that are within a tolerance of a shorter target vector.
#'
#' @param query_mzs A vector of mzs from which values that are close to target_mzs will be identified.
#' @param target_mzs A vector of mzs of interest with which to search query.
#' @param ppm, a mass tolerance in parts per million.
#' @param mz_tol, a mass tolerance in Da.  Superfluous if ppm is supplied.
#' @return A boolean vector of length(query_mzs) indicating matching positions.
#' 
#' @export
filter_by_mz <- function(query_mzs, target_mzs, ppm = 10, mz_tol = NULL){
    # ensure that users only supply one mass tolerance argument (either ppm or mz_tol)
    if(!missing(ppm) & !missing(mz_tol)){
        stop("Only a single mz tolerance parameter {mz_tol or ppm} can be supplied, not both.")
    }
    
    # find min and max mzs
    if(missing(mz_tol)){
        min_mzs <- target_mzs * (1 - ppm / 1e6)
        max_mzs <- target_mzs * (1 + ppm / 1e6)
    } else{
        min_mzs <- target_mzs - mz_tol
        max_mzs <- target_mzs + mz_tol
    }
    
    # allocate boolean matrix
    n_rows <- length(query_mzs)
    n_cols <- length(target_mzs)
    bool_mat <- matrix(FALSE, nrow = n_rows, ncol = n_cols)
    
    # loop over target mzs to create columns of bool_mat
    for(idx in seq_along(target_mzs)){
        min_mz <- min_mzs[idx]
        max_mz <- max_mzs[idx]
        bool_mat[, idx] <- query_mzs >= min_mz & query_mzs <= max_mz
    }
    
    # apply any() to rows of matrix
    return(apply(bool_mat, FUN = any, MARGIN = 1))
}


#' Parse the oddly formatted .csv files exported by Agilent's MassHunter software into an R data frame.
#'
#' @param filename A string with the path of the relevant .csv file.
#' @return A data frame with columns "V2", "V3", "file" and "intensity"
#' 
#' @export
parse_agilent_cgram_csv <- function(filename){
    raw.text <- readLines(filename)
    n.lines <- length(raw.text)
    
    # header lines start with '#' and come in pairs;
    # the first line has metadata and the second has actual headers
    header.lines <- str_detect(raw.text, '[#]') %>% which
    data.lines <- !str_detect(raw.text, '[#]')
    
    # make the dataframe of all data rows, unannotated
    raw.df <- read.csv(text = raw.text[data.lines], 
                       sep = ",", 
                       header = F, 
                       stringsAsFactors = F,
                       colClasses = c('NULL', 'numeric', 'numeric'),
                       quote = '"'
                       )
    
    
    # keep only odd elements of header_lines http://stackoverflow.com/a/13462110/4480692
    header.starts <- header.lines[c(TRUE, FALSE)]
    num.headers <- length(header.starts)
    data.starts <- header.starts + 2
    data.ends <- c(header.starts[2:num.headers] - 1, n.lines)
    
    # now must adjust for line numbering in filtered data
    data.lengths <- data.ends - data.starts + 1
    new.ends <- cumsum(data.lengths)
    new.starts <- c(1, new.ends[1:(num.headers-1)] + 1)
    
    # initialize annotation column
    raw.df$signal <- NA
    raw.df$file <- NA
    
    # add metadata by block
    for(idx in seq_along(new.starts)){
        metadata <- raw.text[header.starts[idx]]
        file <- str_extract(metadata, '[:graph:]+[.]d')
        signal <- str_extract(metadata, '[-]ESI|[/+]ESI|DAD')
	mzs <- str_extract(metadata, 'EIC[:12 ]*[(].*[)]')
        raw.df[new.starts[idx]:new.ends[idx], 'file'] <- file
        raw.df[new.starts[idx]:new.ends[idx], 'signal'] <- signal
	raw.df[new.starts[idx]:new.ends[idx], 'mzs'] <- mzs
        }

    # rename data colums
    names(raw.df)[1:2] <- c('time', 'intensity')	
    return(raw.df)
}

#' Parse the oddly formatted .csv files exported by Agilent's MassHunter software into an R data frame.
#'
#' @param filename A string with the path of the relevant .csv file.
#' @return A data frame with columns "V2", "V3", "file" and "intensity"
#' 
#' @export
parse_agilent_spectrum_csv <- function(filename){
    raw.text <- readLines(filename)
    n.lines <- length(raw.text)
    
    # header lines start with '#' and come in pairs;
    # the first line has metadata and the second has actual headers
    header.lines <- str_detect(raw.text, '[#]') %>% which
    data.lines <- !str_detect(raw.text, '[#]')
    
    # make the dataframe of all data rows, unannotated
    raw.df <- read.csv(text = raw.text[data.lines], 
                       sep = ",", 
                       header = F, 
                       stringsAsFactors = F,
                       colClasses = c('NULL', 'numeric', 'numeric'),
                       quote = '"'
                       )
    
    
    # keep only odd elements of header_lines http://stackoverflow.com/a/13462110/4480692
    header.starts <- header.lines[c(TRUE, FALSE)]
    num.headers <- length(header.starts)
    data.starts <- header.starts + 2
    data.ends <- c(header.starts[2:num.headers] - 1, n.lines)
    
    # now must adjust for line numbering in filtered data
    data.lengths <- data.ends - data.starts + 1
    new.ends <- cumsum(data.lengths)
    new.starts <- c(1, new.ends[1:(num.headers-1)] + 1)
    
    # initialize annotation column
    raw.df$signal <- NA
    raw.df$file <- NA
    
    # add metadata by block
    for(idx in seq_along(new.starts)){
        metadata <- raw.text[header.starts[idx]]
        file <- str_extract(metadata, '[:graph:]+[.]d')
        signal <- str_extract(metadata, '[-]ESI|[/+]ESI|DAD')
	precursor <- str_extract(metadata, '[(][0-9.]*[^(]* -> .*[)]')
	rts <- str_extract(metadata, '[(][0-9.]*[ ]sec[)]')
        raw.df[new.starts[idx]:new.ends[idx], 'file'] <- file
        raw.df[new.starts[idx]:new.ends[idx], 'signal'] <- signal
	raw.df[new.starts[idx]:new.ends[idx], 'rts'] <- rts
	raw.df[new.starts[idx]:new.ends[idx], 'precursor'] <- precursor
        }

    # rename data colums
    names(raw.df)[1:2] <- c('mz', 'intensity')	
    return(raw.df)
}


#' Convert a molecular formula to a protonated m/z value using ecipex.
#'
#' @param formula a string parsable by ecipex as a molecular formula
#' @param num_protons, the number of protons gained (if negative, lost) during ionization.
#' @return mz, an m/z value in Daltons.
#' 
#' @export
find_protonated_mz <- function(formula, num_protons = 1){
    ELECTRON_MASS <- 0.000548579909
    HYDROGEN_MASS <- 1.00782503207
    proton_mass <- HYDROGEN_MASS - ELECTRON_MASS
    mass = ecipex(formula, limit = 1e-2)[[1]]$mass[1] + num_protons*proton_mass 
    charge = abs(num_protons)
    return(mass / charge)
}



#' Make it easier to make consistently colored plots with ggplot.
#' http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#' @param n an integer
#' @return hcl, the colors for ggplot
#' 
#' @export
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}