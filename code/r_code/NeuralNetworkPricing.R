LoadNeuralNetwork <- function(json_file_path){
# Description: Loads a neural network stored in a json file.
  tmp <- fromJSON(json_file_path)
  nLayers <- (length(tmp)-4)/2
  
  # Input and output scalings:
  scaleMeanIn <- matrix(tmp[[length(tmp)-3]],ncol=1)
  scaleStdIn <- sqrt(matrix(tmp[[length(tmp)-2]],ncol=1))
  scaleMeanOut <- matrix(tmp[[length(tmp)-1]],ncol=1)
  scaleStdOut <- sqrt(matrix(tmp[[length(tmp)]],ncol=1))
  
  w <- vector(mode = "list", length = nLayers)
  b <- vector(mode = "list", length = nLayers)
  for (i in 1:nLayers){
    w[[i]] <- t(tmp[[2*i-1]])
    b[[i]] <- as.matrix(tmp[[2*i]])
  }
  
  nn <- setNames(list(w,b,scaleMeanIn,scaleStdIn,scaleMeanOut,scaleStdOut),
                 c("weights","biases","scaleMeanIn","scaleStdIn","scaleMeanOut","scaleStdOut"))
  
  return(nn)
}
elu <- function(x){
# Implements the elu activation function.
  val <- x
  val[x<=0] <- exp(x[x<=0])-1
  return(val)
}
EvaluateNeuralNetwork <- function(nn,x){
# Description: Evaluates a neural network.
  nLayers <- length(nn$biases)
  val <- (x - nn$scaleMeanIn)/nn$scaleStdIn
  for (i in 1:(nLayers-1)){
    val <- elu(nn$weights[[i]]%*%val + nn$biases[[i]])
  }
  val <- (nn$weights[[nLayers]]%*%val + nn$biases[[nLayers]])*nn$scaleStdOut + nn$scaleMeanOut
  return(val)
}
EvaluateModelInGrid <- function(model,x){
# Description: Evaluates a model in the fixed contract grid.
  
  # Check bounds:
  if (any(x < model$lb | x > model$ub)){
    stop("EvaluateModelInGrid: Parameter bounds are violated.")
  }
  
  # Evaluate sub-networks:
  iv = numeric(length(model$k))
  for (i in 1:length(model$nn)){
    iv[model$out_idx[[i]]] = EvaluateNeuralNetwork(model$nn[[i]],x[model$in_idx[[i]]])
  }
  return(iv)
}
AreContractsInDomain <- function(model,kq,Tq){
# Description: Checks which contracts are within the domain supported 
# by the neural networks.

  if (length(kq)!=length(Tq)){
	  stop("AreContractsInDomain: Input vectors were not of the same size.")
  }
  uniqT <- unique(Tq)
  uniqTGrid <- unique(model$T)
  minTGrid <- min(uniqTGrid)
  maxTGrid <- max(uniqTGrid)
  idxValid <- logical(length(kq))
  for (i in 1:length(uniqT)){
    idxT <- Tq == uniqT[i]
    if (uniqT[i] > maxTGrid || uniqT[i] < minTGrid)
    {
      idxValid[idxT] <- FALSE
    }
    else
    {
      if (uniqT[i]==maxTGrid){
        idxAbove <- length(uniqTGrid)
      } else {
        idxAbove <- min(which(uniqTGrid > uniqT[i]))
      }
      idxBelow <- idxAbove - 1
      idxGridBelow <- model$T == uniqTGrid[idxBelow]
      idxGridAbove <- model$T == uniqTGrid[idxAbove]
      idxValid[idxT] <- (kq[idxT] >= max(min(model$k[idxGridBelow]),
                                        min(model$k[idxGridAbove]))) & 
                        (kq[idxT] <= min(max(model$k[idxGridBelow]),
                                          max(model$k[idxGridAbove])))
    }
  }
  return(idxValid)
}
GetImpliedVolatility  <- function(model,par,kq,Tq){
# Description: Evaluates a neural network based model for
# arbitrary contracts within the supported domain.
  
  if (any(!AreContractsInDomain(model,kq,Tq))){
    stop("GetImpliedVolatility: Not all contracts are within the contract grid.")
  } else {
    iv_grid <- EvaluateModelInGrid(model,par)
    iv <- numeric(length(kq))
    uniqT <- unique(Tq)
    uniqTGrid <- unique(model$T)
    maxTGrid = max(uniqTGrid)
    for (i in 1:length(uniqT)){
      idxT <- Tq == uniqT[i]
      if (uniqT[i]==maxTGrid){
        idxAbove <- length(uniqTGrid)
      } else {
        idxAbove <- min(which(uniqTGrid > uniqT[i]))
      }
      idxBelow <- idxAbove - 1
      T_above <- uniqTGrid[idxAbove]
      T_below <- uniqTGrid[idxBelow]
      idxGridBelow <- model$T == uniqTGrid[idxBelow]
      idxGridAbove <- model$T == uniqTGrid[idxAbove]
      
      iv_below_grid <- iv_grid[idxGridBelow]
      iv_above_grid <- iv_grid[idxGridAbove]
      k_below_grid <- model$k[idxGridBelow]
      k_above_grid <- model$k[idxGridAbove]
      
      # Fit splines:
      spline_lower <- splinefun(k_below_grid,iv_below_grid,method="natural")
      spline_upper <- splinefun(k_above_grid,iv_above_grid,method="natural")
      
      # Evaluate spline:
      iv_below <- spline_lower(kq[idxT])
      iv_above <- spline_upper(kq[idxT])
      
      # Interpolate in time dimension:
      frac <- (uniqT[i]-T_below)/(T_above - T_below)
      iv[idxT] <- sqrt( ((1-frac)*T_below*iv_below^2 + frac*T_above*iv_above^2)/uniqT[i] )
    }
  }
  return(iv)
}
LoadModel <- function(contract_grid_folder,weights_folder,model_name){
# Description: Loads a neural network pricing model.

  # Load contract grid:
  logM <- read.delim(paste(contract_grid_folder,"//logMoneyness.txt",sep=""), header = FALSE)
  logM <- logM$V1
  expiries <- read.delim(paste(contract_grid_folder,"//expiries.txt",sep=""), header = FALSE)
  expiries <- expiries$V1
  
  # Load neural networks:
  json_files = c("_weights_1.json",
                 "_weights_2.json",
                 "_weights_3.json",
                 "_weights_4.json",
                 "_weights_5.json",
                 "_weights_6.json")
  nns <- idxIn <- idxOut <- list()
  idxOutStart <- 1
  for (i in 1:length(json_files)){
    fullpath_weights <- paste(weights_folder,"//",model_name,json_files[i],sep="")
    nns[[i]] <- LoadNeuralNetwork(fullpath_weights)
    idxIn[[i]] <- seq(from=1,by=1,to=length(nns[[i]]$scaleMeanIn))
    idxOutEnd <- idxOutStart + length(nns[[i]]$scaleMeanOut) - 1
    idxOut[[i]] <- seq(from=idxOutStart,by=1,to=idxOutEnd)
    idxOutStart <- idxOutEnd + 1
  }
  
  # Set the forward variance curve grid-points (in case relevant):
  Txi = c(seq(from=0.0025,by=0.0025,to=0.0175),
          seq(from=0.02,by=0.02,to=0.14),
          seq(from=0.16,by=0.12,to=1),
          seq(from=1.25,by=0.25,to=2),3)
  
  # Set parameter bounds:
  if (model_name == "rheston"){
	  lb = c(0,0.1,-1,rep(0.05^2,length(Txi)));
	  ub = c(0.5,1.25,0,rep(1,length(Txi)));
  } else if (model_name == "rbergomi"){
	  lb = c(0,0.75,-1,rep(0.05^2,length(Txi)));
	  ub = c(0.5,3.50,0,rep(1,length(Txi)));  
  } else if (model_name == "rbergomi_extended"){
	  lb = c(0.75,-1,-0.5,-0.5,rep(0.05^2,length(Txi)));
	  ub = c(3.50,0,0.5,0.5,rep(1,length(Txi)));  
  } else if (model_name == "heston"){
	  lb = c(0,0.05^2,0,-1,0.05^2);
	  ub = c(25,1,10,0,1);  
	  Txi = c()
  }
  # Define model:
  model <- setNames(list(logM,expiries,nns,idxIn,idxOut,lb,ub,Txi),
                    c("k","T","nn","in_idx","out_idx","lb","ub","Txi"))
  
  return(model)
}