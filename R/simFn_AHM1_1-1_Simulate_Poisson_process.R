# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# sim.fn - AHM1 section 1.1 p4

# A function to help to understand the relationship between point patterns,
# abundance data and occurrence data (also called presence/absence or distribution data)
# (introduced in AHM1 Section 1.1)

sim.fn <- function(quad.size = 10, cell.size = 1, intensity = 1, show.plot = TRUE){
#
# Function that simulates animal or plant locations in space according
# to a homogenous Poisson process. This process is characterized by the
# intensity, which is the average number of points per unit area.
# The resulting point pattern is then discretized to obtain abundance data and
# presence/absence (or occurrence) data. The discretization of space is
# achieved by choosing the cell size.
# Note that you must choose cell.size such that the ratio of quad.size
# to cell.size is an integer.
# Argument show.plot should be set to FALSE when running simulations
# to speed things up.

# Compute some preliminaries
exp.M <- intensity * quad.size^2       # Expected population size in quadrat
breaks <- seq(0, quad.size, cell.size) # boundaries of grid cells
n.cell <- (quad.size / cell.size)^2    # Number of cells in the quadrat
mid.pt <- breaks[-length(breaks)] + 0.5 * cell.size # cell mid-points

# Simulate three processes: point process, cell abundance summary and cell occurrence summary
# (1) Generate and plot the mother of everything: point pattern
M <- rpois(1, exp.M)          # Realized population size in quadrat is Poisson
u1 <- runif(M, 0, quad.size)  # x coordinate of each individual
u2 <- runif(M, 0, quad.size)  # y coordinate of each individual

# (2) Generate abundance data
# Summarize point pattern per cell: abundance (N) is number of points per cell
N <- as.matrix(table(cut(u1, breaks=breaks), cut(u2, breaks= breaks)))
lambda <- round(mean(N),2)    # lambda: average realized abundance per cell
var <- var(c(N))              # Spatial variance of N

# (3) Generate occurrence (= presence/absence) data
# Summarize point pattern even more:
# occurrence (z) is indicator for abundance greater than 0
z <- N   ;  z[z>1] <- 1     # Convert abundance info to presence/absence info
psi <- mean(z)              # Realized occupancy in sampled sites


# Visualisation
if(show.plot){
  op <- par(mfrow = c(2, 2), mar = c(5,5,5,2), cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.3)
  on.exit(par(op))
  tryPlot <- try( {
    # (1) Visualize point pattern
    plot(u1, u2, xlab = "x coord", ylab = "y coord", cex = 1, pch = 16, asp = 1,
       main = paste("Point pattern: \nIntensity =", intensity, ", M =", M, "inds."),
       xlim = c(0, quad.size), ylim = c(0, quad.size), frame = FALSE, col = "red") # plot point pattern
    polygon(c(0, quad.size, quad.size, 0), c(0,0, quad.size, quad.size), lwd = 3,
       col = NA, border = "black")   # add border to grid boundary

    # (2) Visualize abundance pattern
    # Plot gridded point pattern with abundance per cell
    plot(u1, u2, xlab = "x coord", ylab = "y coord", cex = 1, pch = 16, asp = 1,
       main = paste("Abundance pattern: \nRealized mean density =", lambda, "\nSpatial variance =", round(var,2)),
       xlim = c(0, quad.size), ylim = c(0, quad.size), frame = FALSE, col = "red") # plot point pattern
    # Overlay grid onto study area
    for(i in 1:length(breaks)){
       for(j in 1:length(breaks)){
       segments(breaks[i], breaks[j], rev(breaks)[i], breaks[j])
       segments(breaks[i], breaks[j], breaks[i], rev(breaks)[j])
       }
    }
    # Print abundance (N) into each cell
    for(i in 1:length(mid.pt)){
      for(j in 1:length(mid.pt)){
       text(mid.pt[i],mid.pt[j],N[i,j],cex =10^(0.8-0.4*log10(n.cell)),col="blue")
      }
    }
    polygon(c(0, quad.size, quad.size, 0), c(0,0, quad.size, quad.size), lwd = 3, col = NA, border = "black")   # add border to grid boundary

    # (3) Visualize occurrence (= presence/absence) pattern
    # Summarize point pattern even more:
    # occurrence (z) is indicator for abundance greater than 0
    plot(u1, u2, xlab = "x coord", ylab = "y coord", cex = 1, pch = 16, asp = 1,
       main = paste("Occurrence pattern: \nRealized occupancy =", round(psi,2)), xlim = c(0, quad.size),
       ylim = c(0, quad.size), frame = FALSE, col = "red") # plot point pattern
    # Overlay grid onto study area
    for(i in 1:length(breaks)){
       for(j in 1:length(breaks)){
          segments(breaks[i], breaks[j], rev(breaks)[i], breaks[j])
          segments(breaks[i], breaks[j], breaks[i], rev(breaks)[j])
       }
    }
    # Shade occupied cells (which have abundance N > 0 or occurrence z = 1)
    for(i in 1:(length(breaks)-1)){
       for(j in 1:(length(breaks)-1)){
          polygon(c(breaks[i], breaks[i+1], breaks[i+1], breaks[i]),
          c(breaks[j], breaks[j], breaks[j+1], breaks[j+1]),
          col = "black", density = z[i,j]*100)
       }
    }
    polygon(c(0, quad.size, quad.size, 0), c(0,0, quad.size, quad.size), lwd = 3, col = NA, border = "black")   # add border to grid boundary

    # (4) Visualize abundance distribution across sites
    # plot(table(N), xlab = "Abundance (N)", ylab = "Number of cells",
    # col = "black", xlim = c(0, max(N)), main = "Frequency of N with mean density (blue)", lwd = 3, frame = FALSE)
    histCount(N, NULL, xlab = "Abundance (N)", ylab = "Number of cells",
      color = "grey", main = "Frequency of N with mean density (blue)")
    abline(v = lambda, lwd = 3, col = "blue", lty=2)
    }, silent = TRUE )
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
}

# Numerical output
return(list(quad.size = quad.size, cell.size = cell.size, intensity = intensity, exp.N = exp.M,
   breaks = breaks, n.cell = n.cell, mid.pt = mid.pt, M = M, u1 = u1, u2 = u2, N = N, z = z, psi = psi))
}

