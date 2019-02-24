# TODO: Add comment
# Date of Creation: 06/07/2018
# Last Update: 02/23/2019
# Author: Carlos Yasojima
#########################################################
###################### LIBRARIES ########################
library(gstat)
library(sp) 
library(automap)
library(GA)
library(reshape)
library(NISTunits)
library(SearchTrees)
library(RGeostats)
library(fpc)
library(outliers)
library(scales)

############ DATABASE SELECT ###############
# 1 - Meuse, 2 - Wolfcamp, 3 - Broomsbarn, 4 - Coalash, 5 - WalkerLake
databases = c("meuse","wolfcamp","broomsbarn","coalash","walkerlake")
databaseSelected = databases[1]

######### VARIABLES INIT (FIXED) ###########
model1 = "exponential"
model2 = "Exp"
model3 = "Exponential"
trainDataSize = 1  
sillMultiplier = 5
nlags = 20

######### VARIABLES INIT ###########
outputFileName = "TestLog.txt"
nTestsForEachCluster = 5

############ AG PARAMETERS #################
gaPopulation = 50
gaIter = 3

############ KNN PARAMETERS ################
# Neighbours = 4 means 3 neighbours ########
# Cause the closes is the point itself #####
nNeighbours = 4  

############ K-MEANS PARAMTERS #############
kmeansClusters = 3

############################################################################################################################################

## DATABASES TRANSFORMATION AND LOADING ####
switch(databaseSelected,
       meuse={
         data(meuse)
         data = meuse[,c(1,2,6)]
         var = "zinc"},
       
       wolfcamp={
         train = read.table('WOLFCAMP.DAT',header=TRUE, skip=3)
         colnames(train) = c('x','y','plevel')
         train = rbind(c(42.8,127.6,1464.0), train)
         data = train[,c(1:3)]
         var = "plevel"
       },
       
       broomsbarn={
         train = read.table('BroomsBarn.dat', skip=4, dec = ".")
         colnames(train) = c("x","y","K","log10K","pH","P","log10P")
         data = train[,c(1,2,3)]
         var = "K"
       },
       
       coalash={
         data(coalash)
         train = coalash
         colnames(train) = c('x','y','coalash')
         data = train[,c(1:3)]
         var = "coalash"
       },
       
       walkerlake={
         data(walker)
         train = cbind(walker@coords,walker@data$V)
         train = as.data.frame(train)
         colnames(train) = c("x","y","V")
         data = train[,c(1:3)]
         var = "V"})

####################### FUNCTIONS #######################

######## CONVERT DEGREE RADIANS ##########
deg2rad = function(deg) {
  return((pi * deg) / 180)
}

rad2deg = function(rad) {
  return((180 * rad) / pi)
}

######### K-MEANS #########
kmeansClustering = function(data, nCluster) {
  if(nCluster == 1) {
    clusterBase = kmeans(data, nCluster)
    clusteredData = cbind(data, clusterBase$cluster)
    index = ncol(clusteredData)
    colnames(clusteredData)[index] = "cluster"
  } else {
    clusterBase = kmeans(data, nCluster)
    distances = dist(data)
    clusterAnalysis = cluster.stats(distances, clustering = clusterBase$cluster)
    print(clusterBase$betweenss/clusterBase$totss)
    print(clusterAnalysis$dunn)
    print(clusterAnalysis$avg.silwidth)
    print(clusterAnalysis$clus.avg.silwidths)
    clusteredData = cbind(data, clusterBase$cluster)
    index = ncol(clusteredData)
    colnames(clusteredData)[index] = "cluster"
  }
  return(clusteredData)
}

################### DATA PARTITIONING ###################
partitionData = function(data, percentage) {
  smp_size <- floor(percentage * nrow(data))
  train_ind <- sample(seq_len(nrow(data)), size = smp_size)
  train <- data[train_ind, ]
  test <- data[-train_ind, ]
  if(percentage == 1) {
    train = cbind(train, rep(1))
    index = ncol(train)
    colnames(train)[index] = "train"
    train = train[order(as.numeric(row.names(train))), ]
    partitionedData = train
  } else {
    train = cbind(train, rep(1))
    test = cbind(test, rep(0))
    index = ncol(train)
    colnames(train)[index] = "train"
    colnames(test)[index] = "train"
    train = train[order(as.numeric(row.names(train))), ]
    test = test[order(as.numeric(row.names(test))), ]
    partitionedData = rbind(train, test)
  }
  
  return(partitionedData)
}

################### DATA CLUSTER ASSIGNMENT #####################
assignCluster = function(data, nNeighbours) {
  A <- SpatialPoints(cbind(data$x,data$y))
  B <- SpatialPoints(cbind(data$x,data$y))
  
  nNearest = nNeighbours
  tree <- createTree(coordinates(B))
  inds <- knnLookup(tree, newdat=coordinates(A), k=nNearest)
  nROWS = nrow(data)	
  data = cbind(data,inds)
  index = (ncol(data)) - nNeighbours
  nNearest = nNearest + index
  
  for (i in 1:nROWS) {
    for(j in (index+1):nNearest) {
      data[i,j] = data[data[i,j],index]
    }
  }
  
  if(nNearest > (index+1)) {
    clusterSelected = apply(data[,(index+1):nNearest],1,function(x) names(which.max(table(x))))
    data = cbind(data[,c(1:(index-1))], clusterSelected)
  } else {
    clusterSelected = data[,(index+1)]
    data = cbind(data[,c(1:(index-1))], clusterSelected)
  }
  
  return(data)
}

################### GA OPTIMIZATION (GA1) #####################
gaOptim = function(optimData, nCluster, krigVar, popSize, generations, model) {
  coordinates(optimData) = ~x+y
  MSE = vector(length=nCluster)
  MSE2 = vector(length=nCluster)
  detCoefficients = vector(length=nCluster)
  GASolutions = matrix( ,nrow = nCluster, ncol = 5)
  columnIndex = which(colnames(optimData@data)==krigVar)
  
  maxSill = var(optimData@data[,columnIndex])
  maxSill = maxSill * sillMultiplier
  print(maxSill)
  maxNugget = 0
  distances = dist(optimData@coords)
  maxRange = max(distances)
  
  lags = nlags
  w = maxRange/lags
  
  var = paste(krigVar, '~ 1')
  var = as.formula(var)
  
  for(i in 1:nCluster) {
    assign(paste("cluster", i, sep=""), 
           optimData[optimData@data$clusterSelected  %in% i,])
  }
  
  for(i in 1:nCluster) { 
    Object = get(paste0("cluster", i))
    optimBase = Object
    
    fitnessFunction = function(x1)
    {	
      lzn.vgm = variogram(var, optimBase, cutoff=maxRange, width = w) 
      lzn.model = vgm(model, range = x1[3], psill = x1[2], nugget = x1[1], anis = c(x1[4], x1[5]))
      lzn.fit = fit.variogram(object = lzn.vgm, model = lzn.model, fit.sills = FALSE, fit.ranges = FALSE, fit.kappa = FALSE)
      xValid = krige.cv(var, optimBase, model=lzn.fit) 
      
      return (-sum(xValid$residual^2))
    }
    
    if(nrow(optimBase) > 4) {
      GA = ga(type = "real-valued", fitness = fitnessFunction, 
              min = c(0, 0, 0, 1, 0), max = c(maxNugget, maxSill, maxRange, 180, 1), 
              popSize = popSize, maxiter = generations)
      
      GASolutions[i,] = GA@solution
      lzn.vgm = variogram(var, optimBase) 
      lzn.model = vgm(range = GA@solution[1,3], model, 
                      psill = GA@solution[1,2], nugget = GA@solution[1,1], 
                      anis = c( GA@solution[1,4], GA@solution[1,5]))
      lzn.fit = fit.variogram(object = lzn.vgm, model = lzn.model, fit.sills = FALSE, fit.ranges = FALSE, fit.kappa = FALSE)
      
      xValid = krige.cv(var, optimBase, model=lzn.fit)
      
      MSE[i] = sum(xValid$residual^2)
      MSE2[i] = sum(xValid$residual)
      detCoefficients[i] = cor(xValid$observed, xValid$observed - xValid$residual, method = "pearson")^2
    } else {
      MSE[i] = 0
      MSE2[i] = 0
      detCoefficients[i] = 0
    }
  }
  
  colnames(GASolutions) <- c("Nugget", "Sill", "Range", "Angle", "Factor")
  result = list(MSE, MSE2, detCoefficients, GASolutions)
  names(result) <- c("MSE", "MSE2","DetCoefficient", "VariogParams")
  return(result)
}

################### GA OPTIMIZATION (GA2) #####################
gaOptim2 = function(optimData, nCluster, krigVar, popSize, generations, model) {
  coordinates(optimData) = ~x+y
  MSE = vector(length=nCluster)
  MSE2 = vector(length=nCluster)
  detCoefficients = vector(length=nCluster)
  GASolutions = matrix( ,nrow = nCluster, ncol = 3)
  columnIndex = which(colnames(optimData@data)==krigVar)
  
  maxSill = var(optimData@data[,columnIndex])
  maxSill = maxSill * sillMultiplier
  maxNugget = 0
  distances = dist(optimData@coords)
  maxRange = max(distances)
  lags = nlags
  w = maxRange/lags
  
  var = paste(krigVar, '~ 1')
  var = as.formula(var)
  
  for(i in 1:nCluster) {
    assign(paste("cluster", i, sep=""), 
           optimData[optimData@data$clusterSelected  %in% i,])
  }
  
  for(i in 1:nCluster) { 
    Object = get(paste0("cluster", i))
    optimBase = Object
    
    fitnessFunction = function(x1)
    {	
      lzn.vgm = variogram(var, optimBase, cutoff = maxRange, width = w) 
      lzn.model = vgm(model, range = x1[3], psill = x1[2], nugget = x1[1])
      lzn.fit = fit.variogram(object = lzn.vgm, model = lzn.model, fit.sills = FALSE, fit.ranges = FALSE, fit.kappa = FALSE)
      xValid = krige.cv(var, optimBase, model=lzn.fit) 
      
      return (-sum(xValid$residual^2))
    }
    
    if(nrow(optimBase) > 3) {
      GA = ga(type = "real-valued", fitness = fitnessFunction, 
              min = c(0, 0, 0), max = c(maxNugget, maxSill, maxRange), 
              popSize = popSize, maxiter = generations)
      
      GASolutions[i,] = GA@solution
      lzn.vgm = variogram(var, optimBase) 
      lzn.model = vgm(range = GA@solution[1,3], model, 
                      psill = GA@solution[1,2], nugget = GA@solution[1,1])
      lzn.fit = fit.variogram(object = lzn.vgm, model = lzn.model, fit.sills = FALSE, fit.ranges = FALSE, fit.kappa = FALSE)
      
      xValid = krige.cv(var, optimBase, model=lzn.fit)
      
      MSE[i] = sum(xValid$residual^2)
      MSE2[i] = sum(xValid$residual)
      detCoefficients[i] = cor(xValid$observed, xValid$observed - xValid$residual, method = "pearson")^2
    } else {
      MSE[i] = 0
      MSE2[i] = 0
      detCoefficients[i] = 0
    }
  }
  
  colnames(GASolutions) <- c("Nugget", "Sill", "Range")
  result = list(MSE, MSE2, detCoefficients, GASolutions)
  names(result) <- c("MSE", "MSE2","DetCoefficient", "VariogParams")
  return(result)
}

################### LEVENBERG-MARQUADT + LEAST SQUARES (LM-WLS) #####################
lsOptim1 = function(optimData, nCluster, krigVar, model) {
  coordinates(optimData) = ~x+y
  MSE = vector(length=nCluster)
  MSE2 = vector(length=nCluster)
  detCoefficients = vector(length=nCluster)
  LSSolutions = matrix( ,nrow = nCluster, ncol = 3)
  columnIndex = which(colnames(optimData@data)==krigVar)
  
  var = paste(krigVar, '~ 1')
  var = as.formula(var)
  
  for(i in 1:nCluster) {
    assign(paste("cluster", i, sep=""), 
           optimData[optimData@data$clusterSelected  %in% i,])
  }
  
  for(i in 1:nCluster) { 
    Object = get(paste0("cluster", i))
    optimBase = Object
    
    temp = get(paste0("cluster", i))
    initialSill = var(temp@data[,columnIndex])	
    print(initialSill)
    initialNugget = 0
    print(initialNugget)
    distances = dist(temp@coords)
    initialRange = max(distances)
    print(initialRange)
    
    lags = nlags
    w = initialRange/lags
    
    lzn.vgm = variogram(var, optimBase, cutoff = initialRange, width = w) 
    lzn.fit = fit.variogram(lzn.vgm, vgm(psill = initialSill, model = model2, range = initialRange/2, nugget = initialNugget), fit.sills = c(FALSE,TRUE))
    xValid = krige.cv(var, optimBase, model=lzn.fit) 
    
    LSSolutions[i,] = c(lzn.fit$psill[1], lzn.fit$psill[2], lzn.fit$range[2])
    MSE[i] = sum(xValid$residual^2)
    MSE2[i] = sum(xValid$residual)
    detCoefficients[i] = cor(xValid$observed, xValid$observed - xValid$residual, method = "pearson")^2
    
  }
  
  colnames(LSSolutions) = c("nugget", "sill", "range")
  result = list(MSE, MSE2, detCoefficients, LSSolutions)
  names(result) <- c("MSE", "MSE2", "DetCoefficient", "VariogParams")
  return(result)
}

################### ITERATIVE LEAST SQUARES METHOD OPTIMIZATION (GN-ILS2) #####################
iterativeLsOptim = function(optimData, nCluster, krigVar, model, model2) {
  columnIndex = which(colnames(optimData)==krigVar)
  columnNumber = ncol(optimData)
  MSE = vector(length=nCluster)
  MSE2 = vector(length=nCluster)
  ILSSolutions = matrix(, nrow = nCluster, ncol = 5)
  
  var = paste(krigVar, '~ 1')
  var = as.formula(var)
  
  for(i in 1:nCluster) {
    assign(paste("cluster", i, sep=""), optimData[optimData[,columnNumber]  %in% i,])
    assign(paste("cluster", i, sep=""), db.create(get(paste0("cluster", i))[,c(1,2,columnIndex)], flag.grid=FALSE,ndim=2,autoname=F))
    
    temp = get(paste0("cluster", i))
    
    variogram = vario.calc(temp, nlag=nlags)
    
    fit.model = model.auto(variogram, struct=(c(model2)), 
                           maxiter = 10000, title="Anisotropic Variogram to be Fitted", auth.aniso = TRUE)
    resultTEMP = RGeostats::xvalid(temp, fit.model, neigh = neigh.create(type = 0))
    
    ILSparams = cbind(0, fit.model@basics[[1]]@sill)
    ILSparams = cbind(ILSparams, fit.model@basics[[1]]@range)
    ILSparams = cbind(ILSparams, rad2deg(acos(fit.model@basics[[1]]@aniso.rotmat[1])))
    if(((fit.model@basics[[1]]@aniso.coeffs[1,] / fit.model@basics[[1]]@aniso.coeffs[2,]) > 0) & 
       ((fit.model@basics[[1]]@aniso.coeffs[1,] / fit.model@basics[[1]]@aniso.coeffs[2,]) < 1)) {
      ILSparams = cbind(ILSparams, fit.model@basics[[1]]@aniso.coeffs[1,] / fit.model@basics[[1]]@aniso.coeffs[2,])
    } else {
      ILSparams = cbind(ILSparams, fit.model@basics[[1]]@aniso.coeffs[2,] / fit.model@basics[[1]]@aniso.coeffs[1,])
    }
    
    ILSSolutions[i,] = ILSparams	
    
    MSE[i] = sum(resultTEMP@items[,5]^2)
    MSE2[i] = sum(resultTEMP@items[,5])
  }
  
  colnames(ILSSolutions) = c("nugget", "sill", "range", "angle", "factor")
  result = list(MSE, MSE2, ILSSolutions)
  names(result) <- c("MSE", "MSE2", "VariogParams")
  return(result)
}

################### ITERATIVE LEAST SQUARES METHOD OPTIMIZATION (GN-ILS1) #####################
iterativeLsOptim2 = function(optimData, nCluster, krigVar, model, model2) {
  columnIndex = which(colnames(optimData)==krigVar)
  columnNumber = ncol(optimData)
  MSE = vector(length=nCluster)
  MSE2 = vector(length=nCluster)
  ILSSolutions = matrix(, nrow = nCluster, ncol = 5)
  
  var = paste(krigVar, '~ 1')
  var = as.formula(var)
  
  for(i in 1:nCluster) {
    assign(paste("cluster", i, sep=""), optimData[optimData[,columnNumber]  %in% i,])
    assign(paste("cluster", i, sep=""), db.create(get(paste0("cluster", i))[,c(1,2,columnIndex)], flag.grid=FALSE,ndim=2,autoname=F))
    
    temp = get(paste0("cluster", i))
    
    variogram = vario.calc(temp, nlag=nlags, dirvect = c(0,45,90,135))
    
    fit.model = model.auto(variogram, struct=(c(model2)), 
                           maxiter = 10000, title="Anisotropic Variogram to be Fitted", auth.aniso = TRUE, auth.lock2d = TRUE)
    resultTEMP = RGeostats::xvalid(temp, fit.model, neigh = neigh.create(type = 0))
    
    ILSparams = cbind(0, fit.model@basics[[1]]@sill)
    ILSparams = cbind(ILSparams, fit.model@basics[[1]]@range)
    ILSparams = cbind(ILSparams, rad2deg(acos(fit.model@basics[[1]]@aniso.rotmat[1])))
    if(((fit.model@basics[[1]]@aniso.coeffs[1,] / fit.model@basics[[1]]@aniso.coeffs[2,]) > 0) & 
       ((fit.model@basics[[1]]@aniso.coeffs[1,] / fit.model@basics[[1]]@aniso.coeffs[2,]) < 1)) {
      ILSparams = cbind(ILSparams, fit.model@basics[[1]]@aniso.coeffs[1,] / fit.model@basics[[1]]@aniso.coeffs[2,])
    } else {
      ILSparams = cbind(ILSparams, fit.model@basics[[1]]@aniso.coeffs[2,] / fit.model@basics[[1]]@aniso.coeffs[1,])
    }
    
    ILSSolutions[i,] = ILSparams	
    
    MSE[i] = sum(resultTEMP@items[,5]^2)
    MSE2[i] = sum(resultTEMP@items[,5])
  }
  
  colnames(ILSSolutions) = c("nugget", "sill", "range", "angle", "factor")
  result = list(MSE, MSE2, ILSSolutions)
  names(result) <- c("MSE", "MSE2", "VariogParams")
  return(result)
}

#################### MAIN EXEC ##########################
sink()
sink()
scores = scores(data[,var], type = "z", prob=0.99)
print(scores)
print(data[scores==TRUE,])
data = data[scores==FALSE,]

for(i in 1:ncol(data)) {
  data[,i] = rescale(data[,i], to=c(-1,1))
}

partitionedData = partitionData(data, trainDataSize)
train = partitionedData[partitionedData$train %in% 1,]
test = partitionedData[partitionedData$train %in% 0,]

nTests = nTestsForEachCluster
fileName = outputFileName
nCluster = kmeansClusters

sink(fileName, type="output", append = TRUE)
cat("--------- Parameters ---------- ")
cat("\n")
cat("Database Selected: ")
cat(databaseSelected)
cat("\n")
cat("Number of Iterations for each Cluster: ")
cat(nTestsForEachCluster)
cat("\n")
cat("GA Population Size: ")
cat(gaPopulation)
cat("\n")
cat("GA Number of Iterations: ")
cat(gaIter)
cat("\n")
cat("KNN Number of Neighbours: ")
cat(nNeighbours-1)
cat("\n")
cat("Number of Clusters: ")
cat(kmeansClusters)
cat("\n")
cat("\n")
sink()

GA1ResultsMeans = matrix(nrow = nTests, ncol = 2)
GA2ResultsMeans = matrix(nrow = nTests, ncol = 2)
LS1ResultsMeans = matrix(nrow = nTests, ncol = 2)
ILSResultsMeans = matrix(nrow = nTests, ncol = 2)
ILS2ResultsMeans = matrix(nrow = nTests, ncol = 2)

GA1ResultsSd = matrix(nrow = nTests, ncol = 2)
GA2ResultsSd = matrix(nrow = nTests, ncol = 2)
LS1ResultsSd = matrix(nrow = nTests, ncol = 2)
ILSResultsSd = matrix(nrow = nTests, ncol = 2)
ILS2ResultsSd = matrix(nrow = nTests, ncol = 2)

for(z in 1:nCluster) {
  GA1Results = matrix(nrow = nTests, ncol = 2)
  GA2Results = matrix(nrow = nTests, ncol = 2)
  LS1Results = matrix(nrow = nTests, ncol = 2)
  ILSResults = matrix(nrow = nTests, ncol = 2)
  ILS2Results = matrix(nrow = nTests, ncol = 2)
  
  sink(fileName, type="output", append = TRUE)
  cat("CLUSTER ---------- ")
  cat(z)
  cat("\n")
  cat("                          GA --- GA2 -- LM -- ILS -- ILS2           ")
  cat("\n")
  sink()
  
  for(i in 1:nTests) {
    clusteredData = kmeansClustering(train, z)
    
    assignClusterData = assignCluster(clusteredData, nNeighbours)

    result1 = try(gaOptim(optimData = assignClusterData, nCluster = z, krigVar = var, popSize = gaPopulation, generations =  gaIter, model = model2))
    result2 = try(gaOptim2(optimData = assignClusterData, nCluster = z, krigVar = var, popSize = gaPopulation, generations =  gaIter, model = model2))
    result3 = try(lsOptim1(optimData = assignClusterData, nCluster = z, krigVar = var, model = model1))
    result4 = try(iterativeLsOptim2(optimData = assignClusterData, nCluster = z, krigVar = var, model = model2, model2 = model3))
    result5 = try(iterativeLsOptim(optimData = assignClusterData, nCluster = z, krigVar = var, model = model2, model2 = model3))
    
    sink(fileName, type="output", append = TRUE)
    
    cat("Iteration: ")
    cat(i)
    cat("              ")
    variance = var(assignClusterData[,3])
    NMSEfactor = 1/(variance*nrow(assignClusterData))
    try(cat(round(sum(result1$MSE*NMSEfactor),4)))
    cat(" ")
    GA1Results[i,1] = try(sum(result1$MSE*NMSEfactor))

    try(cat(round(sum(result2$MSE*NMSEfactor),4)))
    cat(" ")
    GA2Results[i,1] = try(sum(result2$MSE*NMSEfactor))

    try(cat(round(sum(result3$MSE*NMSEfactor),4)))
    cat(" ")
    LS1Results[i,1] = try(sum(result3$MSE*NMSEfactor))
    
    try(cat(round(sum(result4$MSE*NMSEfactor),4)))
    cat(" ")
    ILSResults[i,1] = try(sum(result4$MSE*NMSEfactor))
    
    try(cat(round(sum(result5$MSE*NMSEfactor),4)))
    cat(" ")
    ILS2Results[i,1] = try(sum(result5$MSE*NMSEfactor))
    cat("               ")
    
    cat("\n")
    if(i == nTests) {
      cat("Mean: ")
      cat("                    ")
      try(cat(round(mean(GA1Results[,1]),4)))
      GA1ResultsMeans[z,1] = try(round(mean(GA1Results[,1]),4))
      cat(" ")
      try(cat(round(mean(GA2Results[,1]),4)))
      GA2ResultsMeans[z,1] = try(round(mean(GA2Results[,1]),4))
      cat(" ")
      try(cat(round(mean(LS1Results[,1]),4)))
      LS1ResultsMeans[z,1] = try(round(mean(LS1Results[,1]),4))
      cat(" ")
      try(cat(round(mean(ILSResults[,1]),4)))
      ILSResultsMeans[z,1] = try(round(mean(ILSResults[,1]),4))
      cat(" ")
      try(cat(round(mean(ILS2Results[,1]),4)))
      ILS2ResultsMeans[z,1] = try(round(mean(ILS2Results[,1]),4))
      
      cat("\n")
      
      cat("Stand. Deviation: ")
      cat("        ")
      try(cat(round(sd(GA1Results[,1]),4)))
      GA1ResultsSd[z,1] = try(round(sd(GA1Results[,1]),4))
      cat(" ")
      try(cat(round(sd(GA2Results[,1]),4)))
      GA2ResultsSd[z,1] = try(round(sd(GA2Results[,1]),4))
      cat(" ")
      try(cat(round(sd(LS1Results[,1]),4)))
      LS1ResultsSd[z,1] = try(round(sd(LS1Results[,1]),4))
      cat(" ")
      try(cat(round(sd(ILSResults[,1]),4)))
      ILSResultsSd[z,1] = try(round(sd(ILSResults[,1]),4))
      cat(" ")
      try(cat(round(sd(ILS2Results[,1]),4)))
      ILS2ResultsSd[z,1] = try(round(sd(ILS2Results[,1]),4))
      
      cat("\n")
    }
    sink()
  }
}

sink(fileName, type="output", append = TRUE)
sink()
