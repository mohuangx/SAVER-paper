library(sm)

boxplot.tj <- function(y, xloc = 1, width.box = 0.25, lwd.box = 2, width.hor = 0.25, 
                       lwd.hor = 2, range.wisk = 1.5, lwd.wisk = 2, pch.box = 16, cex.boxpoint = 2, 
                       plot.outliers = FALSE, pch.out = 1, cex.out = 1, color = "black", 
                       fill = "white") {
  
  # makes boxplot with dot as median and solid whisker Interquartile range =
  # (.75 quantile) - (.25 quantile).  Note: Wiskers are not always symmetrical;
  # top wisker extends up to max(y) constrained by y <= (.75 quantile) +
  # range.wisk*Interquartile range bottom whisker is determined by min(y)
  # constrained by y >= (.25 quantile) - range.wisk*Interquartile range
  
  Q <- quantile(y, c(0.25, 0.5, 0.75), na.rm = TRUE)
  names(Q) <- NULL  # gets rid of percentages
  IQ.range <- Q[3] - Q[1]
  low <- Q[1] - range.wisk * IQ.range
  high <- Q[3] + range.wisk * IQ.range
  index <- which((y >= low) & (y <= high))
  wisk.low <- min(y[index])
  wisk.high <- max(y[index])
  outliers <- y[which((y < low) | (y > high))]
  
  
  
  # plot box:
  xleft <- xloc - width.box/2
  xright <- xloc + width.box/2
  ybottom <- Q[1]
  ytop <- Q[3]
  rect(xleft, ybottom, xright, ytop, lwd = lwd.box, border = color, col = fill)
  
  # plot whiskers:
  segments(xloc, wisk.low, xloc, Q[1], lwd = lwd.wisk, col = color)
  segments(xloc, Q[3], xloc, wisk.high, lwd = lwd.wisk, col = color)
  
  # plot horizontal segments:
  x0 <- xloc - width.hor/2
  x1 <- xloc + width.hor/2
  # plot median:
  segments(x0, Q[2], x1, Q[2], lwd = lwd.hor, col = color)
  # segments(x0, wisk.low, x1, wisk.low, lwd = lwd.hor, col = color)
  # segments(x0, wisk.high, x1, wisk.high, lwd = lwd.hor, col = color)
  
  # plot outliers:
  if (plot.outliers == TRUE) {
    xloc.p <- rep(xloc, length(outliers))
    points(xloc.p, outliers, pch = pch.out, cex = cex.out, col = color)
  }
}

boxplot.ej <- function(y, xloc = 1, width.box = 0.25, lwd.box = 2, width.hor = 0.25, 
                       lwd.hor = 2, range.wisk = 1.5, lwd.wisk = 2, pch.box = 16, cex.boxpoint = 2, 
                       plot.outliers = FALSE, pch.out = 1, cex.out = 1, color = "black", 
                       fill = "white") {
  
  # makes boxplot with dot as median and solid whisker Interquartile range =
  # (.75 quantile) - (.25 quantile).  Note: Wiskers are not always symmetrical;
  # top wisker extends up to max(y) constrained by y <= (.75 quantile) +
  # range.wisk*Interquartile range bottom whisker is determined by min(y)
  # constrained by y >= (.25 quantile) - range.wisk*Interquartile range
  
  Q <- quantile(y, c(0.25, 0.5, 0.75), na.rm = TRUE)
  names(Q) <- NULL  # gets rid of percentages
  IQ.range <- Q[3] - Q[1]
  low <- Q[1] - range.wisk * IQ.range
  high <- Q[3] + range.wisk * IQ.range
  index <- which((y >= low) & (y <= high))
  wisk.low <- min(y[index])
  wisk.high <- max(y[index])
  outliers <- y[which((y < low) | (y > high))]
  
  
  
  # plot box:
  xleft <- xloc - width.box/2
  xright <- xloc + width.box/2
  ybottom <- Q[1]
  ytop <- Q[3]
  rect(xleft, ybottom, xright, ytop, lwd = lwd.box, border = color, col = fill)
  
  # plot whiskers:
  segments(xloc, wisk.low, xloc, Q[1], lwd = lwd.wisk, col = color)
  segments(xloc, Q[3], xloc, wisk.high, lwd = lwd.wisk, col = color)
  
  # plot horizontal segments:
  x0 <- xloc - width.hor/2
  x1 <- xloc + width.hor/2
  # plot median:
  segments(x0, Q[2], x1, Q[2], lwd = lwd.hor, col = color)
  segments(x0, wisk.low, x1, wisk.low, lwd = lwd.hor, col = color)
  segments(x0, wisk.high, x1, wisk.high, lwd = lwd.hor, col = color)
  
  # plot outliers:
  if (plot.outliers == TRUE) {
    xloc.p <- rep(xloc, length(outliers))
    points(xloc.p, outliers, pch = pch.out, cex = cex.out, col = color)
  }
}


vioplot.singmann <- function(x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, fill = "white",
                             horizontal = FALSE, col = NULL, border = "black", lty = 1, lwd = 1, rectCol = "black", 
                             colMed = "white", pchMed = 19, at, add = FALSE, wex = 1, mark.outlier = TRUE, 
                             pch.mean = 4, ids = NULL, drawRect = TRUE, yaxt = "s", width.box = 0.25,
                             width.hor = 0.25, na.rm = FALSE) {
  
  # process multiple datas
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  # pass 1 - calculate base range - estimate density setup parameters for
  # density estimation
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  outliers <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  
  # global args for sm.density function-call
  args <- list(display = "none")
  
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    if (!is.null(ids)) 
      names(data) <- ids
    if (is.null(names(data))) 
      names(data) <- as.character(1:(length(data)))
    
    # calculate plot parameters 1- and 3-quantile, median, IQR, upper- and
    # lower-adjacent
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25, na.rm = na.rm)
    q3[i] <- quantile(data, 0.75, na.rm = na.rm)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    
    # strategy: xmin = min(lower, data.min)) ymax = max(upper, data.max))
    est.xlim <- c(min(lower[i], data.min), max(upper[i], data.max))
    
    # estimate density curve
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), args))
    
    # calculate stretch factor the plots density heights is defined in range 0.0
    # ... 0.5 we scale maximum estimated point to 0.4 per data
    hscale <- 0.4/max(smout$estimate) * wex
    
    # add density curve x,y pair to lists
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
    min.d <- boxplot.stats(data)[["stats"]][1]
    max.d <- boxplot.stats(data)[["stats"]][5]
    # height[[i]] <- height[[i]][(base[[i]] > min.d) & (base[[i]] < max.d)]
    # height[[i]] <- c(height[[i]][1], height[[i]], height[[i]][length(height[[i]])])
    # base[[i]] <- base[[i]][(base[[i]] > min.d) & (base[[i]] < max.d)]
    # base[[i]] <- c(min.d, base[[i]], max.d)
    outliers[[i]] <- list(data[(data < min.d) | (data > max.d)], names(data[(data < 
                                                                               min.d) | (data > max.d)]))
    
    # calculate min,max base ranges
  }
  # pass 2 - plot graphics setup parameters for plot
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5) else range(at) + min(diff(at))/2 * c(-1, 1)
    
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  } else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  
  # setup plot
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    
    box()
    for (i in 1:n) {
      # plot left/right density curve
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), c(base[[i]], 
                                                                  rev(base[[i]])), col = col, border = border, lty = lty, lwd = lwd)
      
      if (drawRect) {
        # plot IQR
        boxplot.tj(datas[[i]], xloc = at[i], width.box = width.box, width.hor = width.hor, fill = fill)
        
        # lines( at[c( i, i)], c(lower[i], upper[i]) ,lwd=lwd, lty=lty)
      }
    }
  } else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    
    box()
    for (i in 1:n) {
      # plot left/right density curve
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], rev(at[i] + 
                                                                         height[[i]])), col = col, border = border, lty = lty, lwd = lwd)
      
      if (drawRect) {
        # plot IQR
        boxplot.ej(datas[[i]], xloc = at[i], width.box = width.box, width.hor = width.hor, fill = fill)
        
        # lines( at[c( i, i)], c(lower[i], upper[i]) ,lwd=lwd, lty=lty)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, q1 = q1, q3 = q3))
}