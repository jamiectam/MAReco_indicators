library(devtools)
install_github("slarge/gradientForest2/pkg/gradientForest")

set.seed(627)
#########
require(reshape2)
require(gradientForest)
require(ggplot2)
require(extendedForest)

###########

combined.cumulative.importance.plot <-
  function (obj,weight=c("uniform","species","rsq.total","rsq.mean","site","site.species","site.rsq.total","site.rsq.mean")[3],
            use.diff=FALSE, prednames=names(obj$X)[-1], show.weights=FALSE, show.gf.names=TRUE, sort=TRUE, ...)
  {
    if ((nw <- length(weight)) > 1) show.gf.names <- show.weights <- FALSE
    gearnames <- names(obj$dens)[-1]
    CU <- if(!show.gf.names) list() else 
      lapply(prednames, function(predictor) do.call("rbind",lapply(obj$gf.names[[predictor]], function(gear) {
        cu <- cumimp(obj,predictor,gf.name=gear)
        dens <- density(obj,predictor,gf.name=gear,gridded=T)
        if (use.diff) cu$y <- diff(c(0,cu$y))
        data.frame(predictor=predictor,gf.name=gear,value=cu$x,CU=cu$y,dens=dens$y)
      })))      
    CU <- c(CU,
            lapply(prednames, function(predictor) do.call("rbind",lapply(weight, function(w) {
              cu <- cumimp(obj,predictor,weight=w)
              if (use.diff) cu$y <- diff(c(0,cu$y))
              data.frame(predictor=predictor,gf.name=paste("combined",w,sep="."),value=cu$x,CU=cu$y,dens=-1)
            })))
    )
    CU <- do.call("rbind",CU)
    imp <- importance(obj)
    o <- order(-imp)
    CU$predictor <- ordered(CU$predictor, levels=if(sort) names(sort(-imp)) else prednames)
    sps <- trellis.par.get("superpose.line")
    n <- nlevels(CU$gf.name)
    ng <- n - nw
    CU <- transform(CU, maxdens=tapply(dens,interaction(value,predictor),max,na.rm=T)[interaction(value,predictor)])
    CU <- transform(CU, hue=as.numeric(unclass(gf.name))/n, saturation=pmax(dens/maxdens,0.1))
    sps$lwd <- rep(c(3,3),c(ng,nw))
    sps$col <- hsv(h=1:n/n,s=1,v=1:n < n) #rainbow(n,end=1)
    trellis.par.set("superpose.line",sps)
    pgfun <- function(x,y,subscripts,groups,group.value,...){
      n <- length(x)
      panel.segments(x[-1],y[-1],x[-n],y[-n],lwd=with(subset(CU[subscripts,],gf.name==group.value),ifelse(dens != -1,3,1)),
                     col=with(subset(CU[subscripts,],gf.name==group.value),hsv(h=hue[-1],s=ifelse(is.na(saturation[-1]),0,saturation[-1]),v=dens != -1)))
    }
    print(xyplot(
      CU~value|predictor,
      data=CU,
      groups=gf.name,      
      type='l',
      scales=list(y=list(relation="same"),x=list(relation="free"),cex=0.5),
      par.strip.text=list(cex=0.5),
      ylab=if (use.diff) "Importance" else "Cumulative Importance",
      xlab="Predictor value",
      as.table=T,
      panel="panel.superpose",
      panel.groups=if(show.weights) pgfun else panel.xyplot,
      key=list(
        space="top",
        columns=min(n,3),
        text=list(levels(CU$gf.name)),
        lines=list(type='l',col=trellis.par.get("superpose.line")$col[1:n],
                   lwd=trellis.par.get("superpose.line")$lwd[1:n])
      )
    ))
  }

###########
combined.performance.plot <-
  function(obj,horizontal = FALSE, show.names = FALSE, las = 2, cex.axis = 0.7, ...)
  {
    Ylab = expression(R^2)
    Box.temp <- boxplot(split(obj$rsq,factor(rep(names(obj$nspec),obj$nspec))),plot=FALSE)
    Bxp.temp <- bxp(Box.temp, show.names = show.names, horizontal = horizontal, width=obj$nspec, las = las, cex.axis = cex.axis)
    axis(labels=Box.temp$names,side=1+horizontal,at=Bxp.temp,cex.axis=0.45,padj=0,las=2)
    mtext(Ylab,side=2-horizontal,line=2.5)
    title("Overall performance of random forests by gradientForest objects")
  }

#################


combinedGradientForest <-
  function(..., nbin=101, method=2, standardize=c("before","after")[1])
  {
    std.options <- c("before","after")
    if (is.na(std.option <- pmatch(standardize,std.options)))
      stop(paste('Unmatched standardize value "',standardize,'". Expecting "before" or "after"',sep=""))
    
    fList <- list(...)
    ngear <- length(fList)
    if(!all(sapply(fList,inherits,"gradientForest")))
      stop("Every argument must be a gradientForest")
    #
    #   assign forest names
    if(is.null(gearnames <- names(fList)))
      gearnames <- paste("F",1:ngear,sep="")
    if (any(empty <- gearnames==""))
      gearnames[empty] <- paste("F",1:ngear,sep="")[empty]
    names(fList) <- gearnames
    #
    #   check the predictor names are the same
    npred <- sapply(fList, function(obj) ncol(obj$X))
    nsite <- sapply(fList, function(obj) nrow(obj$X))
    #    if(!all(npred == npred[1]))
    #      stop("Every forest must have the same number of predictors")
    prednames <- lapply(fList, function(obj) sort(names(obj$X)))  # TO DO: allow for sorting
    allpreds <- unique(sort(unlist(prednames)))
    #   find gradientForest objects that support each predictor
    gf.names <- lapply(namenames(allpreds), function(predictor) gearnames[sapply(prednames, is.element, el=predictor)])
    #    if(!all(prednames == prednames[,1]))
    #      stop("Every forest must have the same predictors")
    #    prednames <- prednames[,1]
    #    npred <- npred[1]
    #
    #   combined predictor matrix and common bins for importance curve grid
    create.df.aux <- function(X) as.data.frame(do.call("cbind",lapply(namenames(allpreds), function(pred) {res <- X[[pred]]; if(is.null(res)) rep(0,nrow(X)) else res})))
    create.df <- function(X,transpose=F) {
      if (transpose) X <- as.data.frame(t(X))
      X <- create.df.aux(X)
      if (transpose) X <- as.data.frame(t(X))
      X
    }
    
    X <- do.call("rbind",lapply(gearnames, function(a) {cbind(gf.name=a,create.df(fList[[a]]$X))}))
    #    X <- do.call("rbind",lapply(gearnames, function(a) {cbind(gf.name=a,fList[[a]]$X)}))
    bins <- do.call("cbind",lapply(X[allpreds], function(x) bin(x,nbin=nbin)))
    imp.rsq.list <- lapply(fList, function(x) create.df(x$imp.rsq,transpose=T))
    imp.rsq <- do.call("cbind",imp.rsq.list)
    rsq <- unlist(lapply(fList, function(x) x$result))
    #
    #   combined density calculation
    X_r <- cbind(X[,1], stack(X[allpreds]))
    names(X_r) <- c("gf.name","value","predictor")
    nspec <- sapply(fList,"[[","species.pos.rsq")
    X_r$nspec <- nspec[X_r$gf.name]
    X_r <- na.omit(X_r)
    dens <- with(X_r,tapply(1:nrow(X_r),predictor,function(sub) {
      whiten(density(value[sub],weight=nspec[sub]/sum(nspec[sub])),lambda=0.95)
    }))
    dens <- c(list(Combined=dens),lapply(fList, function(x) x$dens))
    #
    #   Gather the overall cumulative importances from each gradientForest
    #   Combine cumulative importances
    #   Normalize relative to combined importance
    #
    gridded.cumulative.importance <- function(obj, predictor) {
      cu <- cumimp(obj, predictor=predictor, standardize_after=(std.options[std.option]=="after"))
      grid <- bins[,predictor]
      if (length(cu$x)==1)
        y <- approx(cu$x,cu$y,grid,rule=2,method="constant",yleft=0)$y
      else y <- approx(cu$x,cu$y,grid,rule=2,method="linear")$y
      list(x=grid, y=y)
    }
    #
    #   Linear interpolation to grid
    #   Density outside survey is zero
    #
    interpolate <- function(xy, grid) {
      res <- approx(xy$x,xy$y,grid,rule=1,method="linear")$y
      res[is.na(res)] <- 0
      res
    }
    
    CU <- lapply(namenames(allpreds), function(predictor)
      lapply(fList[gf.names[[predictor]]], gridded.cumulative.importance, predictor=predictor))
    rsq.total <- sapply(lapply(fList,"[[","result"),sum)
    imp.rsq.total <- sapply(imp.rsq.list,rowSums,na.rm=TRUE)
    for (predictor in allpreds) {
      g <- gf.names[[predictor]]
      weight <- rbind(
        uniform = rep(1,length(g)), 
        species = nspec[g], 
        rsq.total = imp.rsq.total[predictor,g],
        rsq.mean = imp.rsq.total[predictor,g]/nspec[g],
        site = nsite[g], 
        site.species = nsite[g]*nspec[g], 
        site.rsq.total = nsite[g]*imp.rsq.total[predictor,g],
        site.rsq.mean = nsite[g]*imp.rsq.total[predictor,g]/nspec[g]
      )
      densList <- lapply(dens[c("Combined",g)],"[[",predictor) # list of densities, combined version first
      grid <- bins[,predictor]
      densMat <- sapply(densList, interpolate, grid=grid)
      CU[[predictor]][["density"]] <- list(x=grid,y=densMat)
      if (method==2) {
        CUmat <- combine.cumulative.importance(CU[[predictor]][g], densMat, grid, weight)
      } else if (method==1) {
        imp <- rowMeans(imp.rsq)[predictor]
        CUmat <- combine.cumulative.importance.method1(CU[[predictor]][g], densMat, grid, weight, imp)
      } else stop(paste("Unknown method:",method))
      for (i in rownames(weight))
        CU[[predictor]][[paste("combined",i,sep=".")]] <- list(x=grid,y=CUmat[,i])
    }
    
    out <- list(
      call = match.call(),
      X = X,
      dens = dens,
      imp.rsq = imp.rsq,
      rsq = rsq,
      nspec = nspec,
      CU = CU,
      gf.names = gf.names
    )
    class(out) <- c("combinedGradientForest","list")
    out
  }
###########################
cumimp.combinedGradientForest <-
  function (x, predictor, weight=c("uniform","species","rsq.total","rsq.mean","site","site.species","site.rsq.total","site.rsq.mean")[3], gf.name, ...)
  {
    if (!inherits(x,"combinedGradientForest"))
      stop(paste("'x' must be a combinedGradientForest object"))
    if (length(predictor) != 1)
      stop(paste("'predictor' must be a single string"))
    if (!is.element(predictor,names(x$X)[-1]))
      stop(paste("Predictor",predictor,"does not belong to combinedGradientForest object"))
    if (is.na(option <- pmatch(weight,c("uniform","species","rsq.total","rsq.mean","site","site.species","site.rsq.total","site.rsq.mean"))))
      stop(paste('Unmatched weight "',weight,'". Expecting one of "uniform", "species", "rsq.total", "rsq.mean", "site", "site.species", "site.rsq.total" or "site.rsq.mean"',sep=""))
    
    if (missing(gf.name)) {
      res <- x$CU[[predictor]][[paste("combined",weight,sep=".")]]
    } else {
      res <- x$CU[[predictor]][[gf.name]]
    }
    res
  }
########################
density.combinedGradientForest <-
  function(x,predictor,gridded=F,gf.name,...)
  {
    if (!inherits(x,"combinedGradientForest"))
      stop(paste("'x' must be a combinedGradientForest object"))
    if(!gridded) {
      if(missing(gf.name))
        x$dens$Combined[[predictor]]
      else x$dens[[gf.name]][[predictor]]
    } else {
      if(missing(gf.name))
        with(x$CU[[predictor]]$density, list(x=x, y=y[,"Combined"]))
      else with(x$CU[[predictor]]$density, list(x=x, y=y[,gf.name]))
    }  
  }
#############################
density.gradientForest <-
  function(x,predictor,...)
  {
    if (!inherits(x,"gradientForest"))
      stop(paste("'x' must be a gradientForest object"))
    x$dens[[predictor]]
  }

###########################
`importance.combinedGradientForest` <-
  function (x, type=c("Weighted","Raw","Species")[1], sort=TRUE, ...)
  {
    if (!inherits(x,"combinedGradientForest"))
      stop(paste("'x' must be a combinedGradientForest object"))
    weighted <- rowSums(x$imp.rsq, na.rm=TRUE)/ncol(x$imp.rsq)
    if (sort)
      o <- order(-weighted)
    else o <- 1:length(weighted)
    nam <- rownames(x$imp.rsq)
    res <- switch(pmatch(type,c("Weighted","Raw","Species")),
                  weighted[o],
                  rowSums(sweep(x$imp.rsq,2,x$rsq,"/"), na.rm=TRUE)[o]/ncol(x$imp.rsq),
                  if (sort) sort(x$rsq,decreasing=T) else x$rsq
    )
    if (is.null(res))
      stop(paste('Unmatched type "',type,'". Expecting one of "Weighted", "Raw" or "Species"',sep=""))
    else res
  }
##################
`plot.combinedGradientForest` <-
  function (x, plot.type = c("Overall.Importance","Predictor.Ranges",
                             "Predictor.Density","Cumulative.Importance","Performance")[1],
            par.args=NULL,plot.args=NULL,...)
    
  {   plot.options <- c("Overall.Importance","Predictor.Ranges","Predictor.Density","Cumulative.Importance","Performance")
  if (!inherits(x,"combinedGradientForest"))
    stop(paste("'x' must be a combinedGradientForest object"))
  if (is.na(plot.option <- pmatch(plot.type,plot.options)))
    stop(paste('Unmatched plot.type "',plot.type,'". Expecting one of "Overall.Importance", "Predictor.Ranges", "Predictor.Density", "Cumulative.Importance" or "Performance")',sep=""))
  
  old.par<-par(no.readonly=TRUE)
  on.exit(par(old.par))
  
  amend.args <- function(default.args, new.args) {
    # replace those that match
    for(arg in intersect(names(default.args), names(new.args))) 
      default.args[[arg]] <- new.args[[arg]]
    # append those that don't match
    extra <- new.args[is.na(match(names(new.args),names(default.args)))]
    c(default.args,extra)
  }
  
  
  if(plot.options[plot.option]=="Overall.Importance"){  
    plot.args.def <- amend.args(list(las = 1, cex.axis = 0.7, cex.names = 0.7), plot.args)
    plot.args.def<- amend.args(plot.args.def,list(...))     
    par.args.def <- amend.args(list(mfrow = c(1, 2), mar = c(4, 6, 2, 1)), par.args)
    par(par.args.def)    
    do.call("overall.importance.combinedGradientForest.plot",c(list(obj=quote(x)),plot.args.def))
  }
  
  
  if(plot.options[plot.option]=="Predictor.Ranges"){	
    plot.args.def<- amend.args(plot.args,list(...))     
    if(!is.null(par.args)) par(par.args)
    do.call("predictor.ranges.plot",c(list(obj=quote(x)),plot.args.def))
  }
  
  
  
  if(plot.options[plot.option]=="Predictor.Density"){	
    plot.args.def<- amend.args(plot.args,list(...))     
    if(!is.null(par.args)) par(par.args)  
    do.call("predictor.density.plot",c(list(obj=quote(x)),plot.args.def))
  }
  
  
  
  if(plot.options[plot.option]=="Cumulative.Importance"){	
    plot.args.def <- amend.args(list(weight="rsq.total", use.diff=FALSE, prednames=names(x$X)[-1], 
                                     show.weights=FALSE, show.gf.names=TRUE, sort=TRUE), plot.args)
    plot.args.def<- amend.args(plot.args.def,list(...))     
    if(!is.null(par.args)) par(par.args)
    do.call("combined.cumulative.importance.plot",c(list(obj=quote(x)),plot.args.def)) 
  }
  
  if(plot.options[plot.option]=="Performance"){	
    plot.args.def <- amend.args(list(horizontal = FALSE, show.names = FALSE, las = 2, cex.axis = 0.7), plot.args)
    plot.args.def<- amend.args(plot.args.def,list(...))  
    old.mar<-par()$mar   
    par.args.def <- amend.args(list(mfrow=c(1,1),mar=old.mar+c(0,2.5,0,0)), par.args)
    par(par.args.def)    
    do.call("combined.performance.plot",c(list(obj=quote(x)),plot.args.def))
    par(mar=old.mar)
  }         
  
  invisible()
  }
############################3
`plot.gradientForest` <-
  function(x,  plot.type=c("Overall.Importance","Split.Density","Cumulative.Importance","Performance")[1], 
           par.args=NULL, plot.args=NULL, ...) 
  {
    if (!inherits(x,"gradientForest"))
      stop(paste("'x' must be a gradientForest object"))
    plot.options <- c("Overall.Importance","Split.Density","Cumulative.Importance","Performance")
    if (is.na(plot.option <- pmatch(plot.type,plot.options)))
      stop(paste('Unmatched plot.type "',plot.type,'". Expecting one of "Overall.Importance", "Split.Density", "Cumulative.Importance" or "Performance"',sep=""))
    
    old.par<-par(no.readonly=TRUE)
    on.exit(par(old.par))
    
    amend.args <- function(default.args, new.args) {
      # replace those that match
      for(arg in intersect(names(default.args), names(new.args))) 
        default.args[[arg]] <- new.args[[arg]]
      # append those that don't match
      extra <- new.args[is.na(match(names(new.args),names(default.args)))]
      c(default.args,extra)
    }
    
    if(plot.options[plot.option]=="Overall.Importance"){	
      plot.args.def <- amend.args(list(cex.axis = 0.7, cex.names = 0.7, las=1, horiz = TRUE), plot.args)
      plot.args.def<- amend.args(plot.args.def,list(...))     
      par.args.def <- amend.args(list(mfrow = c(1, 2), mar = c(4, 6, 2, 1)), par.args)
      par(par.args.def)    
      do.call("overall.importance.plot",c(list(obj=quote(x)),plot.args.def))
    }
    
    if(plot.options[plot.option]=="Split.Density"){      
      plot.args.def <- amend.args(list(leg.posn="topright",bin=F, nbin=101, leg.panel=1, barwidth=1, cex.legend=0.8, line.ylab=1.5), plot.args)
      plot.args.def<- amend.args(plot.args.def,list(...))     
      par.args.def <- amend.args(list(mar =c(4.5, 1.5, 0.5, 4.5), omi = c(0.1, 0.25, 0.1, 0.1)), par.args)
      par(par.args.def)    
      do.call("Split.density.plot.method2",c(list(obj=quote(x)),plot.args.def))
    }
    
    if(plot.options[plot.option]=="Cumulative.Importance") {
      plot.args.def <- amend.args(list(leg.posn="topleft",legend=TRUE, common.scale=F, line.ylab=1.0, cex.legend=0.75, show.species=TRUE, show.overall=TRUE, leg.nspecies=10), plot.args)
      plot.args.def<- amend.args(plot.args.def,list(...))     
      par.args.def <- amend.args(list(mar=c(0.0,2.1,1.1,0), omi=c(0.75, 0.75, 0.1, 0.1)), par.args)
      par(par.args.def)    
      do.call("species.cumulative.plot",c(list(obj=quote(x)),plot.args.def))
    }
    
    if(plot.options[plot.option]=="Performance")  {                               
      plot.args.def <- amend.args(list(horizontal = FALSE, show.names = FALSE, las=2, cex.axis = 0.7,cex.labels=0.7,line=2), plot.args)
      plot.args.def<- amend.args(plot.args.def,list(...))     
      par.args.def <- amend.args(list(mfrow=c(1,1)), par.args)
      par(par.args.def)    
      do.call("performance.plot",c(list(obj=quote(x)),plot.args.def))   
    }
    
    invisible()	
    
  }

########################3
`predict.combinedGradientForest` <-
  function (object, newdata, extrap=TRUE, ...)
  {
    if (!inherits(object,"combinedGradientForest"))
      stop(paste("'object' must be a combinedGradientForest object"))
    linfun <- function(xold,yold,xnew)
      yold[1] + (xnew-xold[1])*diff(yold)/diff(xold)
    if (missing(newdata))
      newdata <- object$X[,-1]
    if(!inherits(newdata,"data.frame"))
      stop("newdata must be a data.frame")
    newnames <- names(newdata)
    if(!all(ok <- newnames %in% names(object$X)[-1])) {
      badnames <- paste(newnames[!ok], collapse=", ")
      stop(paste("the following predictors are not in the gradientForest:\n\t",badnames,sep=""))
    }
    for (varX in newnames) {
      ci <- cumimp(object, varX, ...)
      xold <- range(ci$x)
      yold <- range(ci$y)
      xnew <- range(newdata[,varX],na.rm=T)
      if (extrap)
        ynew <- linfun(xold, yold, xnew)
      else 
        ynew <- yold
      if (xnew[1] < xold[1]) {
        ci$x <- c(xnew[1],ci$x)
        ci$y <- c(ynew[1],ci$y)
      }
      if (xnew[2] > xold[2]) {
        ci$x <- c(ci$x,xnew[2])
        ci$y <- c(ci$y,ynew[2])
      }
      f <- approxfun(ci, rule = 2)  
      newdata[,varX] <- f(newdata[,varX])     
    }
    class(newdata) <- c("predict.gradientForest", "data.frame")
    newdata
  }
#####################
`predict.gradientForest` <-
  function (object, newdata, extrap=TRUE, ...)
  {
    if (!inherits(object,"gradientForest"))
      stop(paste("'object' must be a gradientForest object"))
    linfun <- function(xold,yold,xnew)
      yold[1] + (xnew-xold[1])*diff(yold)/diff(xold)
    if (missing(newdata))
      newdata <- object$X
    if(!inherits(newdata,"data.frame"))
      stop("newdata must be a data.frame")
    newnames <- names(newdata)
    if(!all(ok <- newnames %in% names(object$X))) {
      badnames <- paste(newnames[!ok], collapse=", ")
      stop(paste("the following predictors are not in the gradientForest:\n\t",badnames,sep=""))
    }
    for (varX in newnames) {
      ci <- cumimp(object, varX, ...)
      xold <- range(ci$x)
      yold <- range(ci$y)
      xnew <- range(newdata[,varX],na.rm=T)
      if (extrap)
        ynew <- linfun(xold, yold, xnew)
      else 
        ynew <- yold
      if (xnew[1] < xold[1]) {
        ci$x <- c(xnew[1],ci$x)
        ci$y <- c(ynew[1],ci$y)
      }
      if (xnew[2] > xold[2]) {
        ci$x <- c(ci$x,xnew[2])
        ci$y <- c(ci$y,ynew[2])
      }
      f <- approxfun(ci, rule = 2)  
      newdata[,varX] <- f(newdata[,varX])     
    }
    class(newdata) <- c("predict.gradientForest", "data.frame")
    newdata
  }
####################
predictor.density.plot <-
  function (obj, ...)
  {
    bind.varXY <- function(obj) {do.call("rbind", lapply(names(obj), function(predictor) data.frame(predictor=predictor,value=obj[[predictor]]$x,density=obj[[predictor]]$y))) }
    dens <- do.call("rbind", lapply(names(obj$dens), function(gear) cbind(gf.name=gear,bind.varXY(obj$dens[[gear]]))))
    o <- order(-importance(obj))
    dens$predictor <- ordered(dens$predictor, levels=names(sort(-importance(obj))))
    spl <- trellis.par.get("superpose.line")
    n <- length(levels(dens$gf.name))
    spl$lwd[1] <- 2
    trellis.par.set("superpose.line",spl)
    print(xyplot(
      density~value|predictor,
      data=dens,
      groups=gf.name,
      type='l',
      scales=list(relation="free",cex=0.5),
      par.strip.text=list(cex=0.5),
      ylab="Density",
      xlab="Predictor value",
      as.table=T,
      key=list(
        space="top",
        columns=3,
        text=list(levels(dens$gf.name)),
        lines=list(type='l',col=spl$col[1:n],
                   lwd=spl$lwd[1:n])
      )
    ))
  }
#####################
predictor.ranges.plot <-
  function (obj, ...)
  {
    o <- order(-importance(obj))
    X_r <- cbind(obj$X[,1], stack(obj$X[-1]))
    names(X_r) <- c("gf.name","value","predictor")
    X_r$predictor <- ordered(X_r$predictor, levels=names(sort(-importance(obj))))
    print(bwplot(gf.name ~ value|predictor, X_r, scales=list(x=list(relation="free",cex=0.5)),
                 par.strip.text=list(cex=0.6),xlab="Predictor value", as.table=T, ...))
  }
########################
`print.combinedGradientForest` <-
  function(x,...)
  {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = " ")
    cat("\ngradientForest objects:\n")
    cat(paste(levels(x$X[,1]), sep = " "), "\n", sep = " ")
    cat("\nNumber of Species:\n")
    print(x$nspec)
    cat("\nPredictors:\n")
    cat(paste(names(x$X)[-1], sep = " "), "\n", sep = ", ", fill=T)
    invisible(x)
  }

##########################
