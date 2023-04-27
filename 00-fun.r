extract.time <- function(mydata, datafield, datafield_name, keep.time = TRUE, keep.date=TRUE){
  if (keep.time == TRUE | keep.date == TRUE){
          raw_datafields <- colnames(mydata)[which(startsWith(colnames(mydata), datafield))]
          nr_cols <- length(raw_datafields)
          my_cols <- c(sprintf(paste0(datafield_name, "_", "datetime", "_%1d"), seq(1, nr_cols)))
          colnames(mydata)[which(startsWith(colnames(mydata), datafield))] <- my_cols
	
          for(x in 1:nr_cols) {
            x.col <- my_cols[x]

            if (keep.date == TRUE & keep.time == TRUE) {
                x.date <- paste0(datafield_name, "_", "date", "_", x)
                x.time <- paste0(datafield_name, "_", "time", "_", x)
                        } else if (keep.date == TRUE & keep.time == FALSE) {
                            x.date <- paste0(datafield_name, "_", "date", "_", x)
                            x.time <- NA
                          } else if (keep.date == FALSE & keep.time == TRUE) {
                              x.date <- NA
                              x.time <- paste0(datafield_name, "_", "time", "_", x)
                          }

            mydata <- separate(data = mydata,
			col = x.col,
			into = c(x.date, x.time),
			sep = "T",
			remove = TRUE)

                        x.month <- paste0(datafield_name, "_", "month")
                        setDT(mydata)[, (paste0(x.month, "_", x)) := month(as.Date(mydata[[x.date]]))]
                }

          set(mydata, , grep(paste0(datafield_name, "_date"), colnames(mydata)), NULL)
          x.month.cols <- c(sprintf(paste0(datafield_name, "_", "month", "_%1d"), seq(1, nr_cols)))
          for(col in x.month.cols) {
            set(mydata, j = col, value = as.numeric(mydata[[col]]))
					}
          mydata[[x.month]] <- rowMeans(mydata[, ..x.month.cols], na.rm=TRUE)

	  			for(y.col in x.month.cols) {print(paste0("Nr NA in ", y.col, ":", nrow(mydata[is.na(mydata[[y.col]])])))}
          mydata[[x.month]] <- round(mydata[[x.month]], digits = 0)
          set(mydata, , (x.month.cols), NULL)

          return(mydata)

  } else {
      stop("keep.time and keep.date can't be both FALSE at the same time.")
    }
}

exclude.id <- function(mydata, datafield, datafield_name, my_exclusions, my_cols, exclude = TRUE) {
  nr_cols <- length(colnames(mydata)[which(startsWith(colnames(mydata), datafield))])
  my_cols <- c(sprintf(paste0(datafield_name, "_%1d"), seq(1, nr_cols)))
  df_cols <- colnames(mydata)[which(startsWith(colnames(mydata), datafield))]
  
  #Assign 1 if individual has at least one type of specified string in my_exclusions
  mydata[, (my_cols) := lapply(.SD, function(x) ifelse(x %in% my_exclusions, 1, 0)), .SDcols = df_cols]
  
  #Sum up the instances
  #Remove cols
  mydata[[datafield_name]] <- rowSums(mydata[, ..my_cols], na.rm=TRUE)
  mydata[[datafield_name]] <- ifelse(mydata[[datafield_name]] > 0, 1, 0)
  mydata[, (my_cols) := NULL]
  
  if (exclude) {
    mydata <- mydata[mydata[[datafield_name]] == 0, ]
    mydata[, c(datafield_name) := NULL]
  }
  
  return(mydata)
}

get.means <- function(mydata, datafield, datafield_name, log = FALSE) {
  nr_cols <- length(colnames(mydata)[which(startsWith(colnames(mydata), datafield))])
  my_cols <- c(sprintf(paste0(datafield_name, "_%1d"), seq(1, nr_cols)))
  colnames(mydata)[which(startsWith(colnames(mydata), datafield))] <- my_cols
  
  #Sum up the instances
  #Remove cols
  mydata[[datafield_name]] <- rowSums(mydata[, ..my_cols], na.rm=TRUE)
  
  #If log = TRUE, retrun logged values of selected datafield
  if(log) mydata[, (datafield_name) := lapply(.SD, function(x) log(x + 1)), .SDcols = datafield_name]
  set(mydata, , my_cols, NULL)
  
  return(mydata)
}

last.nm <- function(mydata, datafield, datafield_name, cols.rm = TRUE) {
  my_cols <- grep(datafield, names(mydata), value= TRUE)
  x <- mydata[, ..my_cols]
  x[, res := NA_character_]
  wh = x[, .I]
  for (v in (length(x)-1):1){
    if (!length(wh)) break
    set(x, j="res", i=wh, v = x[[v]][wh])
    wh = wh[is.na(x$res[wh])]
  }
  mydata[, eval(datafield_name) := x$res]
  if(cols.rm) {
    return(mydata[, c(my_cols) := NULL])
  } else return(mydata)
}

mad.cutoff <- function(data, column, cutoff = 2, remove = TRUE) {
  if (cutoff <= 0) stop("Cutoff must be an unsigned positive numeric.")
  if (!class(data)[1] %in% c("data.table", "data.frame")) stop("data must be a data.frame or simliar object.")
  if (!column %in%  colnames(data)) stop("column must be a valid column of data")
  x <- data[[column]]
  ov <- x[which((abs(x - median(x)) / mad(x)) > cutoff)]
  or <- which(data[[column]] %in% ov)
  if (remove == TRUE) {
    return(data[!or, ])
  } else {
    return(or)
  }
}

sd.cutoff <- function(data, column, sd_const = 2, remove = TRUE) {
  if (sd_const <= 0) stop("Cutoff must be an unsigned positive numeric.")
  if (!class(data)[1] %in% c("data.table", "data.frame")) stop("data must be a data.frame or simliar object.")
  if (!column %in%  colnames(data)) stop("column must be a valid column of data")
  x <- data[[column]]
  ov <- x[which(abs(x - mean(x)) > sd_const*sd(x))]
  or <- which(data[[column]] %in% ov)
  if (remove == TRUE) {
    return(data[!or, ])
  } else {
    return(or)
  }
}

shapiro.dt <- function(x, t = c(NA), s = 5000, iter = 100) {
  w <- c()
  p <- c()
  for(i in 1:iter) {
    tv <- x[sample(length(x), s)]
    st <- shapiro.test(tv)
    w <- c(w, st$statistic)
    p <- c(p, st$p.value)
  }
  w <- mean(w)
  p <- mean(p)
  dt <- data.table(
    Trait = t,
    W = w,
    P = p,
    iterations = iter
  )
  return(dt)
}

