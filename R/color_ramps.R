
# Convenience wrappers for colorRampPalette which take care of specifying colors.
# Colors are from colorbrewer2.org and all are color-blind friendly

if(FALSE) {
# Some palettes we can use
# Sequential 9-class ramp, like heat colors, use for raster maps
YlOrRd9 <- c('#ffffcc', '#ffeda0', '#fed976', '#feb24c', '#fd8d3c', '#fc4e2a', '#e31a1c', '#bd0026', '#800026')
# Instead of terrain.colors
YlGn9 <- c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529')
Greys9 <- c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000')
Greens9 <- c('#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b')

# Diverging 11-class ramps
RdYlBu <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')
RdBu <- c('#67001f', '#b2182b', '#d6604d', '#f4a582', '#fddbc7', '#d1e5f0', '#92c5de', '#4393c3', '#2166ac', '#053061')
BrBG <- c('#543005', '#8c510a', '#bf812d', '#dfc27d', '#f6e8c3', '#f5f5f5', '#c7eae5', '#80cdc1', '#35978f', '#01665e', '#003c30')

# Special 10-class, but middle 2 same (if small, don't care if + or -), try this for crosscorrPlot
RdYlYlBu <- c('#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf', '#ffffbf', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4')
}

fixRange <- function(rng, max) {
  rng <- unique(round(rng))
  rngmin <- max(1, min(rng, max))
  rngmax <- min(max, max(rng, 1))
  if(rngmin == rngmax)
    warning("You specified a single colour; this will not be pretty.", call.=FALSE)
  return(rngmin:rngmax)
}

rampYOR <- function(n=5, range=1:9, bias=1, ...) {
  n <- max(round(n))
  range <- fixRange(range, max=9)
  cols <- c('#ffffcc', '#ffeda0', '#fed976', '#feb24c', '#fd8d3c', '#fc4e2a', '#e31a1c', '#bd0026', '#800026')
  colorRampPalette(cols[range], bias=bias, ...)(n)
}

rampGreys <- function(n=5, range=1:9, bias=1, ...) {
  n <- max(round(n))
  range <- fixRange(range, max=9)
  cols <- c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000')
  colorRampPalette(cols[range], bias=bias, ...)(n)
}

rampGreens <- function(n=5, range=1:9, bias=1, ...) {
  n <- max(round(n))
  range <- fixRange(range, max=9)
  cols <- c('#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b')
  colorRampPalette(cols[range], bias=bias, ...)(n)
}

rampBYR <- function(n=5, range=1:11, bias=1, ...) {
  n <- max(round(n))
  range <- fixRange(range, max=11)
  cols <- rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090',
      '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))
  colorRampPalette(cols[range], bias=bias, ...)(n)
}

rampGBr <- function(n=5, range=1:11, bias=1, ...) {
  n <- max(round(n))
  range <- fixRange(range, max=11)
  cols <- rev(c('#543005', '#8c510a', '#bf812d', '#dfc27d', '#f6e8c3',
      '#f5f5f5', '#c7eae5', '#80cdc1', '#35978f', '#01665e', '#003c30'))
  colorRampPalette(cols[range], bias=bias, ...)(n)
}

if(FALSE) {
# examples
image(volcano, col=rampYOR(225))
image(volcano, col=rampBYR(225))
image(volcano, col=rampGBr(225))
points(runif(500), runif(500), pch=16, cex=0.7)
image(volcano, col=rampGBr(225, range=2:10))
points(runif(500), runif(500), pch=16, cex=0.7)

# Try nonsense values
image(volcano, col=rampBYR(225, range=2.3)) # Gives warning
image(volcano, col=rampBYR(225, range=23)) # Gives warning
image(volcano, col=rampBYR(225, range=c(-3, 7.1))) # 'range' becomes 1:7
image(volcano, col=rampBYR(225, alpha=0.7)) # ignores 'alpha'

# Try with rasters
library(raster)
r <- raster(system.file("external/test.grd", package="raster"))
plot(r)
plot(r, col=rampGBr(225))
plot(r, col=rampGBr(225, bias=3))
plot(r, col=rampGBr(225, bias=0.5))

}