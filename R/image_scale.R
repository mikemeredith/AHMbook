# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# Helper function to draw scale for image (from SCR book)
# (introduced somewhere in AHM1 Chapter 9)

# cex.legend added 12 July 2019, v.0.1.4.9083

image_scale <-
function (z, col, x, y = NULL, size = NULL, digits = 2, labels = c("breaks",
    "ranges"), cex.legend=1)
{
    n <- length(col)
    usr <- par("usr")
    mx <- mean(usr[1:2])
    my <- mean(usr[3:4])
    dx <- diff(usr[1:2])
    dy <- diff(usr[3:4])
    if (missing(x))
        x <- mx + 1.05 * dx/2
    else if (is.list(x)) {
        if (length(x$x) == 2)
            size <- c(diff(x$x), -diff(x$y)/n)
        y <- x$y[1]
        x <- x$x[1]
    }
    else x <- x[1]
    if (is.null(size))
        if (is.null(y)) {
            size <- 0.618 * dy/n
            y <- my + 0.618 * dy/2
        }
        else size <- (y - my) * 2/n
    if (length(size) == 1)
        size <- rep(size, 2)
    if (is.null(y))
        y <- my + n * size[2]/2
    i <- seq(along = col)
    rect(x, y - i * size[2], x + size[1], y - (i - 1) * size[2],
        col = rev(col), xpd = TRUE)
    rng <- range(z, na.rm = TRUE)
    bks <- seq(from = rng[2], to = rng[1], length = n + 1)
    bks <- formatC(bks, format = "f", digits = digits)
    labels <- match.arg(labels)
    if (labels == "breaks")
        ypts <- y - c(0, i) * size[2]
    else {
        bks <- paste(bks[-1], bks[-(n + 1)], sep = " - ")
        ypts <- y - (i - 0.5) * size[2]
    }
    text(x = x + 1.2 * size[1], y = ypts, labels = bks, adj = ifelse(size[1] >
        0, 0, 1), xpd = TRUE, cex=cex.legend)
}
