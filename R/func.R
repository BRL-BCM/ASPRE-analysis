require(CTD)
require(R.utils)

make.met.rowname = function(mets){
  rown=tolower(trimws(mets))
  rown = sapply(rown, function(x) gsub("\\*","",x, perl = FALSE))
  rown = sapply(rown, function(x) gsub("(?<=[\\s])\\s*|^\\s+|\\s+$","",x, perl = TRUE))
  rown = trimws(rown)
  return(rown)
}

getFill = function(data_mx){
  fill = apply(data_mx, 1, function(i) 1-(sum(is.na(i))/length(i)))
  fill = as.vector(fill)
  names(fill)=rownames(data_mx)
  return(fill)
}

imputeMissingValues = function(data, ref) {
  data = data[which(rownames(data) %in% rownames(ref)),]
  ref = ref[which(rownames(ref) %in% rownames(data)),]
  data = data[sort(rownames(data)),]
  ref = ref[sort(rownames(ref)),]
  imputed.data = data
  for (met in 1:nrow(ref)) {
    rowData = ref[met,]
    if (any(is.na(rowData))) {
      rowData = as.numeric(rowData[-which(is.na(rowData))])
    } else {
      rowData = as.numeric(rowData)
    }
    # Impute using uniform random variable, where a = 0.99*observed minimum, and b = observed minimum
    min_row = min(rowData)
    if (min_row<0) {
      min_row = -1*min_row
      imputed.data[met, is.na(data[met,])] = tryCatch(-1*runif(sum(is.na(data[met,])), min = 0.99*min_row, max= min_row), 
                                                      error = function(e) e, warning=function(w) print(sprintf("%s: met%d", w, met)))
    } else {
      imputed.data[met, is.na(data[met,])] = tryCatch(runif(sum(is.na(data[met,])), min = 0.99*min(rowData), max= min(rowData)), 
                                                      error = function(e) e, warning=function(w) print(sprintf("%s: met%d", w, met)))
    }
  }
  return(imputed.data)
}


zscoreData = function(data, ref) {
  print("zscoreData() called.")
  
  # Only metabolites that also occur in the reference population can be z-scored
  data = data[which(rownames(data) %in% rownames(ref)),]
  ref = ref[which(rownames(ref) %in% rownames(data)),]
  data = data[sort(rownames(data)),]
  ref = ref[sort(rownames(ref)),]
  
  # Log transform data
  data = log(data)
  ref = log(data.matrix(ref))
  
  zscore.data = matrix(NA, nrow=nrow(data), ncol=ncol(data))
  rownames(zscore.data) = rownames(data)
  colnames(zscore.data) = colnames(data)
  for (met in 1:nrow(data)) {
    met_data = as.numeric(ref[met,])
    rmSamples = unique(c(which(is.na(met_data)), which(is.infinite(met_data))))
    if (length(rmSamples)>0) {
      x = met_data[-rmSamples]
    } else {
      x = met_data
    }
    if (all(is.na(x))) {
      
    } else {
      #x = x[intersect(which(x>quantile(x, 0.025)), which(x<quantile(x, .975)))]
      d = qqnorm(x, plot.it = FALSE);
      x = as.numeric(d$y)
      z = as.numeric(d$x)
      df = data.frame(x=x,z=z)
      t = lm(x~z, data=df)
      mn.est = as.numeric(t$coefficients[1])
      sd.est = as.numeric(t$coefficients[2])
      rm(d,x,z,df,t)
      zscore.data[met,] = (data[met, ]-mn.est)/sd.est
    }
  }
  
  return(zscore.data)
}

combineDatasets = function(ref, research) {
  ref = as.matrix(ref)
  research = as.matrix(research)
  
  unionMets = unique(c(rownames(ref), rownames(research)))
  data = matrix(0, nrow=length(unionMets), ncol=ncol(ref)+ncol(research))
  rownames(data) = unionMets
  colnames(data) = c(colnames(ref), colnames(research))
  for (r in 1:length(unionMets)) {
    if (unionMets[r] %in% rownames(ref)) {
      data[r,colnames(ref)] = as.numeric(ref[unionMets[r], ])
    } else {
      data[r,colnames(ref)] = rep(NA, ncol(ref))
    }
    
    if (unionMets[r] %in% rownames(research)) {
      data[r,colnames(research)] = as.numeric(research[unionMets[r], ])
    } else {
      data[r,colnames(research)] = rep(NA, ncol(research))
    }
  }
  
  return(data)
}


#### UMAP -GGPLOT2 WRAPPER ####
require(umap)
umap = function (mydata, n_neighbors=15, min_dist=0.1, labels = FALSE, printres = FALSE, seed = FALSE, 
          axistextsize = 18, legendtextsize = 18, dotsize = 5, textlabelsize = 4, 
          legendtitle = "Group", controlscale = FALSE, scale = 1, low = "grey", 
          high = "red", colvec = c("skyblue", "gold", "violet", "darkorchid", 
                                   "slateblue", "forestgreen", "violetred", "orange", "midnightblue", 
                                   "grey31", "black"), printheight = 20, printwidth = 22, 
          text = FALSE) 
{
  config=umap.defaults
  config$n_neighbors=n_neighbors
  config$min_dist=min_dist
  if (controlscale == TRUE && class(labels) %in% c("character", 
                                                   "factor") && scale %in% c(1, 2)) {
    stop("when categorical labels, use scale=3")
  }
  if (controlscale == TRUE && class(labels) %in% c("numeric") && 
      scale %in% c(3)) {
    stop("when continuous labels, use scale=1 or scale=2")
  }
  if (controlscale == FALSE && scale %in% c(2, 3)) {
    warning("if your trying to control the scale, please set controlscale=TRUE")
  }
  if (sum(is.na(labels)) > 0 && class(labels) %in% c("character", 
                                                     "factor")) {
    warning("there is NA values in the labels vector, setting to unknown")
    labels <- as.character(labels)
    labels[is.na(labels)] <- "Unknown"
  }
  if (sum(is.na(text)) > 0 && class(text) %in% c("character", 
                                                 "factor")) {
    warning("there is NA values in the text vector, setting to unknown")
    text <- as.character(text)
    text[is.na(text)] <- "Unknown"
  }
  message("***UMAP wrapper function***")
  message("running...")
  if (seed != FALSE) {
    set.seed(seed)
  }
  if (labels[1] == FALSE && text[1] == FALSE) {
    umap <- umap::umap(t(as.matrix(mydata)),config = config)
    scores <- data.frame(umap$layout)
    p <- ggplot(data = scores, aes(x = X1, y = X2)) + geom_point(colour = "skyblue", 
                                                                 size = dotsize) + theme_bw() + theme(legend.position = "none", 
                                                                                                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                                                                      axis.text.y = element_text(size = axistextsize, colour = "black"), 
                                                                                                      axis.text.x = element_text(size = axistextsize, colour = "black"), 
                                                                                                      axis.title.x = element_text(size = axistextsize), 
                                                                                                      axis.title.y = element_text(size = axistextsize)) + 
      scale_colour_manual(values = colvec)
    if (printres == TRUE) {
      message("printing UMAP to current directory...")
      png("UMAP.png", height = printheight, width = printwidth, 
          units = "cm", res = 900, type = "cairo")
      print(p)
      dev.off()
    }
  }
  else if (labels[1] != FALSE && text[1] == FALSE) {
    umap <- umap::umap(t(as.matrix(mydata)),config = config)
    scores <- data.frame(umap$layout)
    if (controlscale == TRUE) {
      if (scale == 1) {
        p <- ggplot(data = scores, aes(x = X1, y = X2)) + 
          geom_point(aes(colour = labels), size = dotsize) + 
          theme_bw() + theme(panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                            colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                          colour = "black"), axis.title.x = element_text(size = axistextsize), 
                             axis.title.y = element_text(size = axistextsize), 
                             legend.title = element_text(size = legendtextsize), 
                             legend.text = element_text(size = legendtextsize)) + 
          labs(colour = legendtitle) + scale_colour_distiller(palette = "Spectral")
      }
      else if (scale == 2) {
        p <- ggplot(data = scores, aes(x = X1, y = X2)) + 
          geom_point(aes(colour = labels), size = dotsize) + 
          theme_bw() + theme(panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                            colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                          colour = "black"), axis.title.x = element_text(size = axistextsize), 
                             axis.title.y = element_text(size = axistextsize), 
                             legend.title = element_text(size = legendtextsize), 
                             legend.text = element_text(size = legendtextsize)) + 
          labs(colour = legendtitle) + scale_colour_gradient(low = low, 
                                                             high = high)
      }
      else if (scale == 3) {
        p <- ggplot(data = scores, aes(x = X1, y = X2)) + 
          geom_point(aes(colour = labels), size = dotsize) + 
          theme_bw() + theme(panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                            colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                          colour = "black"), axis.title.x = element_text(size = axistextsize), 
                             axis.title.y = element_text(size = axistextsize), 
                             legend.title = element_text(size = legendtextsize), 
                             legend.text = element_text(size = legendtextsize)) + 
          labs(colour = legendtitle) + scale_colour_manual(values = colvec)
      }
    }
    else {
      p <- ggplot(data = scores, aes(x = X1, y = X2)) + 
        geom_point(aes(colour = labels), size = dotsize) + 
        theme_bw() + theme(panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                          colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                        colour = "black"), axis.title.x = element_text(size = axistextsize), 
                           axis.title.y = element_text(size = axistextsize), 
                           legend.title = element_text(size = legendtextsize), 
                           legend.text = element_text(size = legendtextsize)) + 
        labs(colour = legendtitle)
    }
    if (printres == TRUE) {
      message("printing UMAP to current directory...")
      png("UMAPlabeled.png", height = printheight, width = printwidth, 
          units = "cm", res = 900, type = "cairo")
      print(p)
      dev.off()
    }
  }
  else if (labels[1] != FALSE && text[1] != FALSE) {
    umap <- umap::umap(t(as.matrix(mydata)),config = config)
    scores <- data.frame(umap$layout)
    scores$label <- text
    if (controlscale == TRUE) {
      if (scale == 1) {
        p <- ggplot(data = scores, aes(x = X1, y = X2, 
                                       label = label)) + geom_point(aes(colour = labels), 
                                                                    size = dotsize) + theme_bw() + theme(panel.grid.major = element_blank(), 
                                                                                                         panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                                                                                                        colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                                                                                                      colour = "black"), axis.title.x = element_text(size = axistextsize), 
                                                                                                         axis.title.y = element_text(size = axistextsize), 
                                                                                                         legend.title = element_text(size = legendtextsize), 
                                                                                                         legend.text = element_text(size = legendtextsize)) + 
          labs(colour = legendtitle) + scale_colour_distiller(palette = "Spectral") + 
          geom_text(vjust = "inward", hjust = "inward", 
                    size = textlabelsize)
      }
      else if (scale == 2) {
        p <- ggplot(data = scores, aes(x = X1, y = X2, 
                                       label = label)) + geom_point(aes(colour = labels), 
                                                                    size = dotsize) + theme_bw() + theme(panel.grid.major = element_blank(), 
                                                                                                         panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                                                                                                        colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                                                                                                      colour = "black"), axis.title.x = element_text(size = axistextsize), 
                                                                                                         axis.title.y = element_text(size = axistextsize), 
                                                                                                         legend.title = element_text(size = legendtextsize), 
                                                                                                         legend.text = element_text(size = legendtextsize)) + 
          labs(colour = legendtitle) + scale_colour_gradient(low = low, 
                                                             high = high) + geom_text(vjust = "inward", 
                                                                                      hjust = "inward", size = textlabelsize)
      }
      else if (scale == 3) {
        p <- ggplot(data = scores, aes(x = X1, y = X2, 
                                       label = label)) + geom_point(aes(colour = labels), 
                                                                    size = dotsize) + theme_bw() + theme(panel.grid.major = element_blank(), 
                                                                                                         panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                                                                                                        colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                                                                                                      colour = "black"), axis.title.x = element_text(size = axistextsize), 
                                                                                                         axis.title.y = element_text(size = axistextsize), 
                                                                                                         legend.title = element_text(size = legendtextsize), 
                                                                                                         legend.text = element_text(size = legendtextsize)) + 
          labs(colour = legendtitle) + scale_colour_manual(values = colvec) + 
          geom_text(vjust = "inward", hjust = "inward", 
                    size = textlabelsize)
      }
    }
    else {
      p <- ggplot(data = scores, aes(x = X1, y = X2, label = label)) + 
        geom_point(aes(colour = labels), size = dotsize) + 
        theme_bw() + theme(panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                          colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                        colour = "black"), axis.title.x = element_text(size = axistextsize), 
                           axis.title.y = element_text(size = axistextsize), 
                           legend.title = element_text(size = legendtextsize), 
                           legend.text = element_text(size = legendtextsize)) + 
        labs(colour = legendtitle) + geom_text(vjust = "inward", 
                                               hjust = "inward", size = textlabelsize)
    }
    if (printres == TRUE) {
      message("printing UMAP to current directory...")
      png("UMAPlabeled.png", height = printheight, width = printwidth, 
          units = "cm", res = 900, type = "cairo")
      print(p)
      dev.off()
    }
  }
  else if (labels[1] == FALSE && text[1] != FALSE) {
    umap <- umap::umap(t(as.matrix(mydata)),config = config)
    scores <- data.frame(umap$layout)
    scores$label <- text
    p <- ggplot(data = scores, aes(x = X1, y = X2, label = label)) + 
      geom_point(aes(colour = factor(rep(1, ncol(mydata)))), 
                 size = dotsize) + theme_bw() + theme(legend.position = "none", 
                                                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                      axis.text.y = element_text(size = axistextsize, colour = "black"), 
                                                      axis.text.x = element_text(size = axistextsize, colour = "black"), 
                                                      axis.title.x = element_text(size = axistextsize), 
                                                      axis.title.y = element_text(size = axistextsize)) + 
      scale_colour_manual(values = colvec) + geom_text(vjust = "inward", 
                                                       hjust = "inward", size = textlabelsize)
    if (printres == TRUE) {
      message("printing UMAP to current directory...")
      png("UMAP.png", height = printheight, width = printwidth, 
          units = "cm", res = 900, type = "cairo")
      print(p)
      dev.off()
    }
  }
  message("done.")
  return(p)
}


####################################################
## Performs t-tests comparing each row of dataGroup1
## with the corresponding row in dataGroup2.
## Uses the apply function to perform the loop.
####################################################
perform_all_t_tests = function(dataGroup1,dataGroup2,paired=FALSE){
    nGroup1 = ncol(dataGroup1)
    nGroup2 = ncol(dataGroup2)
    dataAll = cbind(dataGroup1,dataGroup2)
    tTestWithErrorHandling = function(x){
        testResult = try(t.test(x[1:nGroup1],x[(nGroup1+1):(nGroup1+nGroup2)],paired = paired),silent=TRUE);
        if(is.character(testResult)){
            return(c(NA,NA,NA))
        }else{
            return(c(testResult$p.value,testResult$estimate))
        }
    }
    results = matrix(unlist(apply(dataAll,1,tTestWithErrorHandling)),ncol=3,byrow=TRUE)
    colnames(results) = c("P.value","Mean.group.1","Mean.group.2")
    return(results)
}

perform_all_t_tests.4mx = function(data_mx,group=list(), p.adjust.method = "holm", paired=FALSE, do.FoldChange = FALSE, do.MeanDiff = TRUE) {
  names(group)=gsub(" ","_",names(group))
  if(length(group)==2){
    p.adjust.method=p.adjust.method
    p.df = perform_all_t_tests(data_mx[,group[[1]]], data_mx[,group[[2]]],paired = paired)
    p.adj=p.adjust(p.df[,"P.value"],method = p.adjust.method)
    p.df=cbind(p.adj,p.df)
    p.df=as.matrix(p.df)
    rownames(p.df) = as.character(rownames(data_mx))
    colnames(p.df) = c(sprintf("P.adj.%s.%s",names(group)[1],names(group)[2]),
                       sprintf("P.%s.%s",names(group)[1],names(group)[2]),
                       sprintf("Mean.%s",names(group)[1]),
                       sprintf("Mean.%s",names(group)[2])
    )
    if(do.FoldChange){
      p.df=cbind(p.df,log2FC=log2(p.df[,4]/p.df[,3]))
      colnames(p.df) = c(colnames(p.df)[1:4],sprintf("log2FC.%s.%s",names(group)[1],names(group)[2]))
    }
    if(do.MeanDiff){
      p.df=cbind(p.df,meanDiff=p.df[,4]-p.df[,3])
      colnames(p.df) = c(colnames(p.df)[1:(ncol(p.df)-1)],sprintf("MeanDiff.%s-%s",names(group)[2],names(group)[1]))
    }
  }else{
    print("Number of groups exceeded 2.")
  }
  return(p.df)
}


perform_all_t_tests.pair_wise = function(data_mx,group=list(),p.adjust.method = "holm",paired=FALSE, my_comparison = NULL,do.FoldChange = FALSE, do.MeanDiff = TRUE,p.adjust.method2 = "BH"){
  names(group)=gsub(" ","_",names(group))
  if (is.null(my_comparison)) {
    comparison_pairs = as.list(data.frame(combn(unique(names(group)), 2),stringsAsFactors = F))
  } else {
    comparison_pairs = lapply(my_comparison,function(x) gsub(" ","_",x))
  }
  p.df = lapply(comparison_pairs,
                function(x) perform_all_t_tests.4mx(data_mx,group=group[x],
                                                    p.adjust.method = p.adjust.method,
                                                    paired=paired,
                                                    do.FoldChange = do.FoldChange,
                                                    do.MeanDiff = do.MeanDiff)
  )
  p.adj = data.frame(lapply(p.df,function(x) matrix(x[,1],dimnames = list(dimnames(x)[[1]],dimnames(x)[[2]][1]))))
  p.result = list(
    comparison_pairs = comparison_pairs,
    p.list = p.df,
    p.value = data.frame(lapply(p.df,function(x) matrix(x[,2],dimnames = list(dimnames(x)[[1]],dimnames(x)[[2]][2])))),
    p.adj= p.adj,
    p.adj2 = t(apply(p.adj,c(1),p.adjust,method=p.adjust.method2))
  )
  if (do.FoldChange) {
    p.result[["logFC"]] = data.frame(lapply(p.df,function(x) matrix(x[,grep("log2FC",colnames(x))],
                                                                    dimnames = list(dimnames(x)[[1]],dimnames(x)[[2]][grep("log2FC",colnames(x))]))))
  }
  if (do.MeanDiff){
    p.result[["meanDiff"]] = data.frame(lapply(p.df,function(x) matrix(x[,grep("MeanDiff.",colnames(x))],
                                                                       dimnames = list(dimnames(x)[[1]],dimnames(x)[[2]][grep("MeanDiff.",colnames(x))]))))
  }
  return(p.result)
}

is_outlier <- function(x) {
  y=rep(FALSE,length(x))
  if(length(which(is.na(x)))>0){
    ind=which(is.na(x))
    a=x[-ind]
    for (i in 1:length(y)){
      if(!i %in% ind){
        y[i]= x[i] < quantile(a, 0.25) - 1.5 * IQR(a) | x[i] > quantile(a, 0.75) + 1.5 * IQR(a)
      }else{
        y[i]= FALSE
      }
    }
  }else{
      y= x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x)
    }
  return(y)
}

is_outlier1 <- function(x) {
  y=rep(FALSE,length(x))
  if(length(which(is.na(x)))>0){
    ind=which(is.na(x))
    a=x[-ind]
    for (i in 1:length(y)){
      if(i != ind){
        y[i]= x[i] < quantile(a, 0.25) - 1 * IQR(a) | x[i] > quantile(a, 0.75) + 1 * IQR(a)
      }else{
        y[i]= FALSE
      }
    }
  }else{
    y= x < quantile(x, 0.25) - 1 * IQR(x) | x > quantile(x, 0.75) + 1 * IQR(x)
  }
  return(y)
}

PEASA.getDataMxByInput=function(input){
  #----get data_mx----
  load(file = "/Volumes/Seagate-2C/BaylorMHG/PE_ASA/1st_052820/dataSet_new_compID.RData")
  metReport=list()
  #----input----
  rmX=input$rmX
  level=input$level
  fillThreshold=1-input$fillThreshold
  DiffThr=input$DiffThr
  ASA2Outlier=input$ASA2Outlier
  rmBatchEffect=input$rmBatchEffect
  # input$level
  if (level == "zscore"){
    data_mx = dataset$data_zscore
  }else if ( level == "imputed scaled"){
    data_mx = dataset$data_norm
  }else if(level == "anchor"){
    data_mx = dataset$data_raw
  }else{
    print("not included.")
  }
  data_mx=apply(data_mx,c(1,2),as.numeric)
  fill = getFill(data_mx)
  data_mx=data_mx[!fill==0,]
  # input$level
  if(rmX == "Yes"){
    temp= compAnno$BIOCHEMICAL[match(rownames(data_mx), compAnno$COMP_ID)]
    ind=grep("x - ",temp)
    temp=temp[!is.na(temp)]
    if(length(ind)>0){
      data_mx=data_mx[-ind,]
    }
  }
  data_mx=apply(data_mx,c(1,2),as.numeric)
  
  #---- generate Histogram of metabolites fill rate ----
  fill = getFill(data_mx)
  name = rep("fill.rate",length(fill))
  temp=data.frame(name=name,fill.rate=fill)
  
  # input$fillThreshold remove mets not passing missing rate threshold
  fillThreshold=1-input$fillThreshold
  if(length(which(fill<fillThreshold)) >0 ){
    data_mx = data_mx[-which(fill<fillThreshold),]
  }
  
  #---- filter fill----
  # input$DiffThr batch missing rate difference
  data_mx.placebo = data_mx[,sampleAttr$sIDs[!sampleAttr$ASA]]
  data_mx.ASA = data_mx[,sampleAttr$sIDs[sampleAttr$ASA]]
  
  fill.placebo = getFill(data_mx.placebo)
  fill.ASA = getFill(data_mx.ASA)
  DiffThr=input$DiffThr
  
  #---- filter batchDiff----
  data_mx=data_mx[!abs(fill.placebo-fill.ASA)>DiffThr,]
  rm(list=setdiff(ls(),c(lsf.str(),"rmX","level","fillThreshold","DiffThr","ASA2Outlier","rmBatchEffect","metReport","data_mx","dataset","compAnno","sampleAttr","input")))
  
  #---- remove provided outliers ----
  data_mx=data_mx[,!colnames(data_mx) %in% input$outliers]
  sampleAttr=sampleAttr[!sampleAttr$sIDs %in% input$outliers,]
  

  
  #---- rmOutliersLowFillRateMets ----
  if(length(grep("ASA_outlier",sampleAttr$BErm))>1){
    temp=getFill(data_mx[,grep("ASA_outlier",sampleAttr$BErm)])
    # tmp=temp[order(temp)][1:30]
    # compAnno$BIOCHEMICAL[match(names(tmp),compAnno$COMP_ID)]
    # data_mx=data_mx[temp>0.4,]
  }
  
  #---- rmBatchEffect ----
  
  if(input$batchEffectCorrectMethod=="RUV"){
    require(MetNorm)
    load("/Volumes/Seagate-2C/BaylorMHG/PE_ASA/1st_052820/batcheffect/QCM_old_useraw.RData")
    # spikes=rownames(data_mx)%in%(QCMs[c(2,3,4,8,9,10)])
    # spikes=rownames(data_mx)%in%(c(QCMs[c(2,3,10)],"36808"))
    # QCMs=unique(c(QCMs[c(2,3,10)]))
    # QCMs=unique(c(QCMs[c(2,3,10)],"36808"))
    # QCMs=unique(c(QCMs[c(1,2,3,9,10)]))
    # QCMs=unique(c(QCMs[c(1,2,3)]))
    data_mx=imputeMissingValues(data_mx,data_mx)
    # QCMs=c("1564","36808","37506")
    # QCMs=c("36808")
    # QCMs=unlist(unique(sapply(QCMs,function(x) get_cor_QCM(x,t(data_mx)))))
    spikes=rownames(data_mx)%in%(QCMs)
    # NA.ind=apply(data_mx,c(1,2),is.na)
    
    data_mx.nocorrect=data_mx
    # smpl.ind=sampleAttr$sIDs[sampleAttr$consss_batch=="non-outlier"]
    
    data_mx.corrected=NormalizeRUVRand(t(data_mx),ctl=spikes,k=2,lambda=NULL,plotk=TRUE)
    data_mx.corrected=t(data_mx.corrected$newY)
    data_mx=data_mx.corrected
    # data_mx[,smpl.ind]=data_mx.nocorrect[,smpl.ind]
    # data_mx[NA.ind]=NA
    
  }else if(input$batchEffectCorrectMethod=="harmony"){
    load("/Volumes/Seagate-2C/BaylorMHG/PE_ASA/1st_052820/outlier.ASA_tr2.RData")
    # data_mx.nocorrect=data_mx
    # smpl.ind=sampleAttr$sIDs[sampleAttr$consss_batch=="non-outlier"]
    # NA.ind=apply(data_mx,c(1,2),is.na)
    V = t(data_mx)
    V=t(imputeMissingValues(data_mx,data_mx))
    # sampleAttr$consss_batch=replace(sampleAttr$consss_batch,sampleAttr$consss_batch=="non-outlier"&sampleAttr$ASA,"ASA-non-outlier")
    meta_data = sampleAttr[match(rownames(V),sampleAttr$sIDs),]
    harmony_embeddings <- harmony::HarmonyMatrix(
      V, meta_data, 'consss_batch', do_pca = FALSE, verbose=TRUE
    )
    data_mx=t(harmony_embeddings)
    # data_mx[,smpl.ind]=data_mx.nocorrect[,smpl.ind]
    
    # data_mx[NA.ind]=NA
   
    
  }else if(input$batchEffectCorrectMethod=="ComBat"){
    require(sva)
    data_mx=imputeMissingValues(data_mx,data_mx)
    # mod = model.matrix(~as.factor(ASA_tri)+as.factor(pe), data=sampleAttr)
    x=ComBat(dat = data_mx, batch = sampleAttr$consss_batch, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
    data_mx=x
  }
  #---- ----
  metReport$data_mx=data_mx
  metReport$sampleAttr=sampleAttr
  

  return(metReport)
}


gg_circle <- function(rx, ry, xc, yc, color="black", fill=NA, ...) {
  x <- xc + rx*cos(seq(0, pi, length.out=100))
  ymax <- yc + ry*sin(seq(0, pi, length.out=100))
  ymin <- yc + ry*sin(seq(0, -pi, length.out=100))
  annotate("ribbon", x=x, ymin=ymin, ymax=ymax, color=color, fill=fill, ...)
}

ropls.plot <- function(d, plottype = "score", xvar, yvar, hotelling = FALSE, ellipse = FALSE, col.var = NULL, col.pca = NULL,useMetName = FALSE){
  N <- nrow(d@scoreMN)
  sm <- d@modelDF
  if(length(d@subsetVi>0)){
    y <- d@suppLs$yMCN[d@subsetVi]
  }else{
    y <- d@suppLs$yMCN
    }
  hotFisN <- (N - 1) * 2 * (N^2 - 1) / (N^2 * (N - 2)) * qf(0.95, 2, N - 2)
  
  #Define scores
  if(length(d@orthoScoreMN) >= length(d@scoreMN)){
    score = data.frame(d@scoreMN, d@orthoScoreMN,label=rownames(d@scoreMN))
  }else{
    score = data.frame(d@scoreMN,label=rownames(d@scoreMN))}
  
  #Define loadings
  if(useMetName){
    rownames(d@loadingMN)=compAnno$BIOCHEMICAL[match(rownames(d@loadingMN),compAnno$COMP_ID)]
  }
  if(length(d@orthoLoadingMN) >= length(d@loadingMN)){
    loading = data.frame(d@loadingMN, d@orthoLoadingMN,label=rownames(d@loadingMN))
    }else{loading = data.frame(d@loadingMN,label=rownames(d@loadingMN))}
  
  #plotting scores
  if(plottype == "score"){
    if(d@typeC == "PCA"){
      p <- ggplot(score, aes_string(x = xvar, y = yvar, col = col.pca)) + geom_point(size = 2)
    }else{
      p <- ggplot(score, aes_string(x = xvar, y = yvar, col = "y",label="label")) + geom_point(size = 2) +
        scale_colour_gradient(low = "gold3", high = "royalblue")
      }
    
    p <- p + xlab(paste(xvar,":", sm[xvar, "R2X"] * 100, "%")) + ylab(paste(yvar,":", sm[yvar, "R2X"] * 100, "%"))
    
    if(hotelling){
      p <- p + gg_circle(
        rx = sqrt(as.numeric(var(score %>% select(xvar))) * hotFisN),
        ry = sqrt(as.numeric(var(score %>% select(yvar))) * hotFisN),
        xc = 0, yc = 0)}
    
    if (ellipse) {
      p <- p + stat_ellipse(
        geom = "polygon", alpha = 0.3, linetype = "blank",
        aes_string(fill = "y"), type = "norm")}
  }
  
  #plotting loadings
  if(plottype == "loading"){
    p <- ggplot(loading, aes_string(x = xvar, y = yvar, col = col.var, label="label")) + geom_point(size = 2)
    p <- p + xlab(paste(xvar,":", sm[xvar, "R2X"] * 100, "%")) + ylab(paste(yvar,":", sm[yvar, "R2X"] * 100, "%"))
    p <- p + geom_hline(yintercept = 0, color = "gray") + geom_vline(xintercept = 0, color = "gray")
  }
  
  return(p)
}

# load("ASAonPREG.loadings.df.RData")
# # df=as.matrix(df)
# var=df[,c("BIOCHEMICAL","VIP.PREG","VIP.ASA")]
# 
# 
# ASA.rank=df$VIP.ASA
# names(ASA.rank)=df$BIOCHEMICAL
# ASA.rank=ASA.rank[ASA.rank>1]
# ASA.rank=ASA.rank[order(ASA.rank,decreasing = T)]
# ASA.rank=names(ASA.rank)
# 
# PREG.rank=df$VIP.PREG
# names(PREG.rank)=df$BIOCHEMICAL
# PREG.rank=PREG.rank[PREG.rank>1]
# PREG.rank=PREG.rank[order(PREG.rank,decreasing = T)]
# PREG.rank=names(PREG.rank)
# 
# fgsea.set=list(ASA.rank=ASA.rank,PREG.rank=PREG.rank)
# fgsea.stats=df$Strength
# names(fgsea.stats)=df$BIOCHEMICAL
# fgsea.stats=fgsea.stats[1:30]
# 
# fgseaRes <- fgsea(pathways = fgsea.set, stats = fgsea.stats, minSize=5, maxSize=500, nperm=10000)
# 
# var.name="BIOCHEMICAL"
# value.name="VIP.PREG"
# 
# ropls.enrich <- function(var, var.name, value.name, filterset = FALSE){
#   var.arrange = var %>% dplyr::distinct(!!sym(var.name), .keep_all = TRUE) %>%
#     arrange(-!!sym(value.name))
#   
#   fgsea.test = var.arrange %>%
#     select(!!sym(var.name), !!sym(value.name)) %>%
#     tibble::deframe()
#   
#   fgsea.set = var.arrange %>%
#     select(-!!sym(value.name)) %>%
#     tidyr::gather(key = "collection", value = "value", (var %>% select(-!!sym(var.name), -!!sym(value.name)) %>% colnames)) %>%
#     tidyr::unite("set", collection, value, sep = "_")
#   
#   if(filterset != FALSE){
#     fgsea.set = fgsea.set %>% filter(grepl(filterset,.$set))
#   }
#   
#   fgsea.set = split(fgsea.set %>% select(!!sym(var.name)) %>% as.matrix %>% as.character(), fgsea.set$set)
#   fgseaRes <- fgsea(pathways = fgsea.set, stats = fgsea.test, minSize=5, maxSize=500, nperm=10000)
#   
#   return(fgseaRes)
# }
# 
rampred <<- function(numdata,brk=c("quantile"),rev = FALSE){
  if(brk == "quantile"){
    brks = quantile(numdata, probs = seq(.05, .95, .05), na.rm = TRUE)
  }else if(brk == "range"){
    brks = quantile(seq(min(numdata, na.rm = TRUE),max(numdata, na.rm = TRUE),diff(range(numdata, na.rm = TRUE))/20), probs = seq(.00, 1, .05), na.rm = TRUE)
  }else{
    print("choose quantile or range")
  }
  # brks = quantile(seq(min(numdata),max(numdata),diff(range(numdata))/20), probs = seq(.00, 1, .05), na.rm = TRUE)
  clrs = round(seq(255, 40, length.out = length(brks) + 1), 0) %>% {paste0("rgb(255,", ., ",", ., ")")}
  if(rev){
    clrs = rev(clrs)
  }
  return(list(brks=brks,clrs=clrs))
}

getzscore=function(data_mx){
  zscore=apply(data_mx, 1,function(x) (x-mean(x[!is_outlier(x)]))/sd(x[!is_outlier(x)]))
}

get_cor_QCM=function(IS,Y){
  IS<-Y[,IS]
  r<-numeric(dim(Y)[2])
  for(j in 1:length(r)){
    r[j]<-cor(IS,Y[,j]) }
  ctl<-logical(length(r))
  # ctl[which(r>round(quantile(r,0.9),2))]<-TRUE
  ctl[which(r>round(quantile(r,0.995),2))]<-TRUE
  
  cor_QCM=colnames(Y)[ctl]
  names(cor_QCM)=compAnno$BIOCHEMICAL[match(cor_QCM,compAnno$COMP_ID)]
  return(cor_QCM)
}
