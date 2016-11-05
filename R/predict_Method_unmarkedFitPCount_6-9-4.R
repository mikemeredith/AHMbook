# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kéry & Andy Royle, Academic Press, 2016.

# predict Method for unmarkedFitPCount - section 6.9.4 p265

# Revised predict function for "unmarkedFitPCount" (in Section 6.9.4)
# -------------------------------------------------------------------
# (1) this has not been tested.
# (2) Only gives 95% confidence interval.
# (introduced in Section)
setMethod("predict", "unmarkedFitPCount",
    function(object, type, newdata, backTransform = TRUE, na.rm = TRUE,
        appendData = FALSE, level=0.95, ...)
{
    if(type %in% c("psi", "alpha"))
        stop(type, " is scalar, so use backTransform instead")
    if(missing(newdata) || is.null(newdata))
        newdata <- getData(object)
    formula <- object@formula
    detformula <- as.formula(formula[[2]])
    stateformula <- as.formula(paste("~", formula[3], sep=""))
    if(inherits(newdata, "unmarkedFrame"))
        class(newdata) <- "unmarkedFrame"
    cls <- class(newdata)[1]
    if(!cls %in% c("unmarkedFrame", "data.frame", "RasterStack"))
        stop("newdata should be an unmarkedFrame, data.frame, or RasterStack", call.=FALSE)
    if(identical(cls, "RasterStack"))
        if(!require(raster))
            stop("raster package is required")
    switch(cls,
    unmarkedFrame = {
        designMats <- unmarked:::getDesign(newdata, formula, na.rm = na.rm)
        switch(type,
            state = {
                X <- designMats$X
                offset <- designMats$X.offset
                },
            det = {
                X <- designMats$V
                offset <- designMats$V.offset
                })
        },
    data.frame = {
        switch(type,
            state = {
                mf <- model.frame(stateformula, newdata)
                X <- model.matrix(stateformula, mf)
                offset <- model.offset(mf)
                },
            det = {
                mf <- model.frame(detformula, newdata)
                X <- model.matrix(detformula, mf)
                offset <- model.offset(mf)
                })
            },
    RasterStack = {
        cd.names <- names(newdata)
        npix <- prod(dim(newdata)[1:2])
        isfac <- is.factor(newdata)
        if(any(isfac))
            stop("This method currently does not handle factors")
        z <- as.data.frame(matrix(getValues(newdata), npix))
        names(z) <- cd.names
        switch(type,
               state = {
                   varnames <- all.vars(stateformula)
                   if(!all(varnames %in% cd.names))
                       stop("At least 1 covariate in the formula is not in the raster stack.\n   You probably need to assign them using\n\t 'names(object) <- covariate.names'")
                   mf <- model.frame(stateformula, z, na.action="na.pass")
                   X.terms <- attr(mf, "terms")
                   X <- model.matrix(X.terms, mf)
                   offset <- model.offset(mf)
               },
               det= {
                   varnames <- all.vars(detformula)
                   if(!all(varnames %in% cd.names))
                       stop("At least 1 covariate in the formula is not in the raster stack.\n   You probably need to assign them using\n\t 'names(object) <- covariate.names'")
                   mf <- model.frame(detformula, z, na.action="na.pass")
                   X.terms <- attr(mf, "terms")
                   X <- model.matrix(X.terms, mf)
                   offset <- model.offset(mf)
               })
    })
    out <- data.frame(matrix(NA, nrow(X), 4,
        dimnames=list(NULL, c("Predicted", "SE", "lower", "upper"))))
    mix <- object@mixture
    lam.mle <- coef(object, type="state")
    if(identical(mix, "ZIP") & identical(type, "state")) {
        psi.hat <- plogis(coef(object, type="psi"))
        if(is.null(offset))
            offset <- rep(0, nrow(X))
#warning("Method to compute SE for ZIP model has not been written. Scratch that.
#Method has been written but not tested/evaluated.
#Also, you only get a 95% confidence interval for the ZIP model. ")
    }
    for(i in 1:nrow(X)) {
        if(nrow(X) > 5000) {
            if(i %% 1000 == 0)
                cat("  doing row", i, "of", nrow(X), "\n")
        }
        if(any(is.na(X[i,])))
            next
        if(identical(mix, "ZIP") & identical(type, "state")) {
## for the ZIP model the predicted values on the log scale have us add log(1-psi.hat) to
### the normal linear prediction
            out$Predicted[i] <-                X[i,] %*% lam.mle + offset[i] + log(1 - psi.hat)
## to compute the approximate SE, I compute the variance of the usual linear part -- that is easy
## and to that I add the variance of log(1-psi.hat) obtained by the delta approximation
logit.psi<-coef(object,type="psi")
#  To do that I took derivative of log(1-psi.hat) using application of chain rule.... hopefully correctly.
delta.approx.2ndpart<-   ( ((1/(1-psi.hat))*(exp(logit.psi)/((1+exp(logit.psi))^2)))^2 ) * (SE(object)["psi(psi)"]^2)
## now the SE is the sqrt of the whole thing
out$SE[i]<- sqrt(  t(X[i,])%*%vcov(object)[1:ncol(X),1:ncol(X)]%*%X[i,] + delta.approx.2ndpart   )
## Here I use a 95% confidence interval b/c I'm not sure how to use "confint"!!!
ci <- c(out$Predicted[i]-1.96*out$SE[i],out$Predicted[i] + 1.96*out$SE[i])
out$lower[i]<- ci[1]
out$upper[i]<- ci[2]
            if(backTransform){
                out$Predicted[i] <- exp(out$Predicted[i])
### If back-transform, delta approx says var = (exp(linear.predictor)^2)*Var(linear.predictor)
### also I exponentiate the confidence interval.....
out$SE[i]<- out$Predicted[i]*out$SE[i]
ci<-exp(ci)
# formula from Goodman 1960 JASA.  This is the se based on "lambda*(1-psi)"
## not sure how well it compares to what I did above.
#part2<-  coef(object,type="psi")
#var.psi.part<- (exp(part2)/((1+exp(part2))^2))*(SE(object)["psi(psi)"]^2)
#part1<- X[i,]*exp(X[i,]%*%lam.mle)
#var.lambda.part<- t(part1)%*%vcov(object)[1:ncol(X),1:ncol(X)]%*%(part1)
#out$SE[i]<-out$Predicted[i]*out$Predicted[i]*var.psi.part + (1-psi.hat)*(1-psi.hat)*var.lambda.part - var.psi.part*var.lambda.part
#ci<- c( NA, NA)
}

        } else {
            lc <- linearComb(object, X[i,], type, offset = offset[i])
            if(backTransform)
                lc <- backTransform(lc)
            out$Predicted[i] <- coef(lc)
            out$SE[i] <- SE(lc)
            ci <- confint(lc, level=level)
        }
        out$lower[i] <- ci[1]
        out$upper[i] <- ci[2]
    }
    if(appendData) {
        if(!identical(cls, "RasterStack"))
            out <- data.frame(out, as(newdata, "data.frame"))
        else
            out <- data.frame(out, z)
    }
    if(identical(cls, "RasterStack")) {
        E.mat <- matrix(out[,1], dim(newdata)[1], dim(newdata)[2],
                        byrow=TRUE)
        E.raster <- raster(E.mat)
        extent(E.raster) <- extent(newdata)
        out.rasters <- list(E.raster)
        for(i in 2:ncol(out)) {
            i.mat <- matrix(out[,i], dim(newdata)[1], dim(newdata)[2],
                            byrow=TRUE)
            i.raster <- raster(i.mat)
            extent(i.raster) <- extent(newdata)
            out.rasters[[i]] <- i.raster
        }
        out.stack <- stack(out.rasters)
        names(out.stack) <- colnames(out)
        out <- out.stack
    }
    return(out)
})

