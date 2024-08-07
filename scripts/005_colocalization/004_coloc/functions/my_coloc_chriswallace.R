
prior.adjust <- function(summ,newp12,p1=1e-4,p2=1e-4,p12=1e-6) {
  if(is.list(summ) && "summary" %in% names(summ))
    summ <- summ$summary
  if(!identical(names(summ), c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")))
    stop("not a coloc summary vector")
  ## back calculate likelihoods
  f <- function(p12)
    prior.snp2hyp(summ["nsnps"],p12=p12,p1=p1,p2=p2)
  pr1 <- f(newp12)
  pr0 <- matrix(f(p12),nrow=nrow(pr1),ncol=ncol(pr1),byrow=TRUE)
  ## make everything conformable
  ## if(is.matrix(summ) && nrow(summ)==1) summ <- as.vector(summ)
  ## if(nrow(pr1)==1) pr1 <- as.vector(pr1)
  ## if(nrow(pr1)>1) pr1 <- t(pr1)
  newpp <- matrix(summ[-1],nrow=nrow(pr1),ncol=ncol(pr1),byrow=TRUE) * pr1/pr0 # prop to, not equal to
  newpp/rowSums(newpp)
}


prior.snp2hyp <- function(nsnp,p12=1e-6,p1=1e-4,p2=1e-4) {
  if(any(p12<p1*p2) || any(p12 > p1) || any(p12 > p2))
    return(NULL)
  tmp <- cbind(nsnp * p1,
               nsnp * p2,
               nsnp * (nsnp-1) * p1 * p2,
               nsnp * p12)
  tmp <- cbind(1-rowSums(tmp),tmp)
  ## if(nrow(tmp)==1) {
  ##     tmp <- c(tmp)
  ##     names(tmp) <- paste0("H",0:4)
  ## } else 
  colnames(tmp) <- paste0("H",0:4)
  tmp
}

manh.plot <- function(df, wh,
                      position=if("position" %in% names(df)) {
                        df$position
                      } else {
                        1:nrow(df)
                      }) {
  znm <- if(wh==1) { "z.df1" } else {"z.df2" }
  ## print(znm)
  ## print(head(df))
  logp <- - ( pnorm(-abs(df[[znm]]),log.p=TRUE) + log(2) ) / log(10)
  ## mycol <- ifelse(A$snp %in% nCV, "red","black")
  Pal <- colorRampPalette(c('white','blue'))
  
  ##This adds a column of color values
  ## based on the y values
  Col <- Pal(100)[ceiling(100*df$SNP.PP.H4)]
  plot(position,logp,col="gray20",
       bg = Col, # Fill colour
       pch = 21, # Shape: circles that can filed
       frame.plot = FALSE, # Remove the frame 
       xlab=if("position" %in% names(df)) {
         "Chromosome position"
       } else {
         "SNP number"
       },
       ylab="-log10(p)",
       xaxt='n')
  ## main=paste("Trait",wh))
  axis(side=1,labels=FALSE) 
}

##' Shows how prior and posterior per-hypothesis probabilities change as a function of p12
##'
##' Function is called mainly for plotting side effect.  It draws two plots, showing how prior and posterior probabilities of each coloc hypothesis change with changing p12.  A decision rule sets the values of the posterior probabilities considered acceptable, and is used to shade in green the region of the plot for which the p12 prior would give and acceptable result.  The user is encouraged to consider carefully whether some prior values shown within the green shaded region are sensible before accepting the hypothesis.  If no shading is shown, then no priors give rise to an accepted result.
##' @title Prior sensitivity for coloc
##' @param obj output of coloc.detail or coloc.process
##' @param rule a decision rule.  This states what values of posterior probabilities "pass" some threshold.  This is a string which will be parsed and evaluated, better explained by examples.  "H4 > 0.5" says post prob of H4 > 0.5 is a pass.  "H4 > 0.9 & H4/H3 > 3" says post prob of H4 must be > 0.9 AND it must be at least 3 times the post prob of H3."
##' @param dataset1 optional the dataset1 used to run SuSiE. This will be used to make a Manhattan plot if plot.manhattans=TRUE.
##' @param dataset2 optional the dataset2 used to run SuSiE. This will be used to make a Manhattan plot if plot.manhattans=TRUE.
##' @param npoints the number of points over which to evaluate the prior values for p12, equally spaced on a log scale between p1*p2 and min(p1,p2) - these are logical limits on p12, but not scientifically sensible values.
##' @param doplot draw the plot. set to FALSE if you want to just evaluate the prior and posterior matrices and work with them yourself
##' @param plot.manhattans if TRUE, show Manhattans of input data
##' @param preserve.par if TRUE, do not change par() of current graphics device - this is to allow sensitivity plots to be incoporated into a larger set of plots, or to be plot one per page on a pdf, forexample
##' @param row when coloc.signals() has been used and multiple rows are returned in the coloc summary, which row to plot
##' @return list of 3: prior matrix, posterior matrix, and a pass/fail indicator (returned invisibly)
##' @export
##' @author Chris Wallace
my_sensitivity <- function(obj,rule="", trait1_title, trait2_title,
                           dataset1=NULL,dataset2=NULL,
                           npoints=100,doplot=TRUE,plot.manhattans=TRUE,
                           preserve.par=FALSE,
                           row=1) {
  stopifnot("list" %in% class(obj))
  stopifnot("priors" %in% names(obj))
  stopifnot("summary" %in% names(obj))
  if(rule=="")
    stop("please supply a rule to define colocalisation, eg 'H4 > thr' where thr is some probability of H4 that you accept as colocalisation")
  rule.init <- rule
  rule <- gsub("(H.)","PP.\\1.abf",rule,perl=TRUE)
  
  ## massage results object
  results=obj$results
  ## multiple signals?
  multiple=FALSE
  if(is.data.table(obj$summary)) { # we're not in coloc.abf anymore
    if(!(row %in% 1:nrow(obj$summary)))
      stop("row must be between 1 and ",nrow(obj$summary))
    pp <- unlist(c(obj$summary[row,grep("PP|nsnp",names(obj$summary)),with=FALSE]))
    if(paste0("SNP.PP.H4.row",row) %in% names(results)) {
      multiple=TRUE
      results[["SNP.PP.H4"]]  <- results[[paste0("SNP.PP.H4.row",row)]]
    }
    if(paste0("z.df1.row",row) %in% names(results)) { # might be passed here or in separate dataset objects
      results[["z.df1"]]  <- results[[paste0("z.df1.row",row)]]
      results[["z.df2"]]  <- results[[paste0("z.df2.row",row)]]
    } else {
      pp <- unlist(c(obj$summary[row,grep("PP|nsnp",names(obj$summary)),with=FALSE]))
    }
  } else {
    pp <- obj$summary
  }
  ## need to add z score from datasets?
  if(!is.null(dataset1) && !is.null(dataset2)) {
    df1=with(dataset1,data.table(snp=snp,position=position,z.df1=beta/sqrt(varbeta)))
    df2=with(dataset2,data.table(snp=snp,position=position,z.df2=beta/sqrt(varbeta)))
    df=merge(df1,df2,by=c("snp","position"),all=TRUE)
    results=merge(results,df,by="snp")
  }
  
  p12 <- obj$priors["p12"]
  p1 <- obj$priors["p1"]
  p2 <- obj$priors["p2"]
  check <- function(pp) { with(as.list(pp),eval(parse(text=rule))) }
  pass.init <- check(pp)
  message("Results ",if(check(pp)) { "pass" } else { "fail" }, " decision rule ",rule.init)
  
  testp12 <- 10^seq(log10(p1*p2),log10(min(p1,p1)),length.out=npoints)
  testH <- prior.snp2hyp(pp["nsnps"],p12=testp12,p1=p1,p2=p2)
  testpp <- as.data.frame(prior.adjust(summ=pp,newp12=testp12,p1=p1,p2=p2,p12=p12))
  colnames(testpp) <- gsub("(H.)","PP.\\1.abf",colnames(testpp),perl=TRUE)
  pass <- check(testpp)
  w <- which(pass)
  
  if(doplot) {
    H <- as.character(0:4)
    op <- par('mfcol', 'mar', 'mfrow','mar','mgp','las','tck')
    on.exit(par(op))
    if(!preserve.par) {
      if(plot.manhattans)
        layout(mat = matrix(1:4,2,2),
               heights = c(1, 1), # Heights of the two rows
               widths = c(2, 3)) # Widths of the two columns
      else
        par(mfcol=c(1,2))
    }
    par(mar = c(3, 3, 2, 1) # Dist' from plot to side of page
        ,mgp = c(2, 0.4, 0) # Dist' plot to label
        ,las = 1 # Rotate y-axis text
        ,tck = -.01 # Reduce tick length
    )
    if(plot.manhattans) {
      if(!("z.df1" %in% colnames(results)) || !("z.df2" %in% colnames(results)))
        stop("please supply dataset1, dataset2, if you want to view Manhattan plots. Otherwise set plot.manhattans=FALSE.")
      manh.plot(results,1)
      title(main=paste(trait1_title, if(multiple) { paste("row",row) } else { "" }))
      manh.plot(results,2)
      title(main=paste(trait2_title, if(multiple) { paste("row",row) } else { "" }))
    }
    m <- list(testH,as.matrix(testpp))
    ti <-   list("Prior probabilities", "Posterior probabilities")
    for(i in 1:2) {
      ym <- if(i==1) { max(m[[i]][,-1]) } else { max(m[[i]]) }
      matplot(testp12,m[[i]],log="x",xlab="p12",ylab="Prob",
              type="b",
              bg = 1:5, # Fill colour
              pch = 21, # Shape: circles that can filed
              col="gray20",
              frame.plot = FALSE, # Remove the frame
              panel.first = abline(h = seq(0, 1, 0.2), col = "grey80"),
              ylim=c(0,ym))
      title(main=ti[[i]],adj=0)
      title(sub=paste("shaded region:",rule.init),adj=0)
      ## title(main=paste("Acceptance rule (shaded region):",rule.init))
      ## legend("topleft",pch=rep(21,5),pt.bg=1:5,legend=paste0("H",0:4))
      if(i==1)
        legend("left",inset=c(0.1,0),bg="white",pch=rep(21,5),pt.bg=1:5,pt.cex=2,legend=paste0("H",0:4))
      abline(v=p12,lty="dashed",col="gray")
      text(p12,0.5,"results",srt=90,col="gray40")
      if(any(pass))
        rect(xleft=testp12[min(w)],ybottom=0,
             xright=testp12[max(w)],ytop=1,
             col=rgb(0,1,0,alpha=0.1), border="green")
      ## add text showing rule
      ## mtext(paste("shaded region:",rule.init),side=3,adj=1)
    }
  }
  
  invisible(cbind(testpp,p12=testp12,pass=pass))
}

