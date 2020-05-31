
# Convenience wrappers for colorRampPalette which take care of specifying colors.
# Colors are from colorbrewer2.org and all are color-blind friendly

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
