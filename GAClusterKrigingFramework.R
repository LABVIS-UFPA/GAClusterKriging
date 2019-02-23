# TODO: Add comment
# Date: 06/07/2018
# Author: Carlos Yasojima
#########################################################
###################### LIBRARIES ########################
library(gstat)
library(sp) 
library(automap)
library(GA)
library(reshape)
library(ggplot2)
#library(dbscan)
library(geoR)
library(NISTunits)
library(SearchTrees)
library(RGeostats)
library(fpc)
library(outliers)
library(scales)
library(gridExtra)
library(DEoptim)
library(georob)

####################### FUNCTIONS #######################

######## CONVERT DEGREE RADIANS ##########
deg2rad = function(deg) {
  return((pi * deg) / 180)
}

rad2deg = function(rad) {
  return((180 * rad) / pi)
}

############## NORMALIZE ################

normalizeByFactor = function(data, var, normalizationFactor) {
  data[,var] = data[,var]/normalizationFactor
  return(data)
}

######### K-MEANS #########
kmeansKneeGraph = function(data, varName) {
  ss = matrix(nrow = 15,ncol = 2)
  for (i in 1:15) {
    cluster = kmeans(data, i)
    ss[i,1] = sum(cluster$withinss)
  }
  ss[,2] = seq(1:15)
  dfm <- melt(ss, id.vars = "clusters")
  dfm = dfm[1:15,]
  colnames(dfm) = c("clusters","Year","SumOfQquares")
  a = ggplot(dfm, aes(clusters, SumOfQquares)) + 
    geom_line() + labs(x = "Number of Clusters", y = "Sum of Squares", subtitle = varName)
  
  return(a)
}

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

########## DBSCAN ##########
dbscanKneeGraph = function(data, nPoints, epsilon) {
  dbscan::kNNdistplot(data, k =  nPoints)
  abline(h = epsilon, lty = 2)
}

dbscanClustering = function(data, nPoints, epsilon) {
  clusterBase <- dbscan(data, eps = epsilon, minPts = nPoints)
  print(clusterBase)
  clusteredData = cbind(data, clusterBase$cluster)
  index = ncol(clusteredData)
  colnames(clusteredData)[index] = "cluster"
  clusteredData <- clusteredData[ which(clusteredData$cluster != 0) , ]
  
  return(clusteredData)
}

########## TREE CLUSTER ##########
#https://journal.r-project.org/archive/2015/RJ-2015-032/RJ-2015-032.pdf

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

################### GA OPTIMIZATION #####################
gaOptim = function(optimData, nCluster, krigVar, popSize, generations, model) {
  coordinates(optimData) = ~x+y
  MSE = vector(length=nCluster)
  MSE2 = vector(length=nCluster)
  detCoefficients = vector(length=nCluster)
  GASolutions = matrix( ,nrow = nCluster, ncol = 5)
  columnIndex = which(colnames(optimData@data)==krigVar)
  
  maxSill = var(optimData@data[,columnIndex])
  maxSill = maxSill * 10
  print(maxSill)
  maxNugget = 0
  distances = dist(optimData@coords)
  maxRange = max(distances)
  
  lags = 20
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
      
      #return (cor(xValid$observed, xValid$observed - xValid$residual, method = "pearson")^2)
      #return (-(1-(cor(xValid$observed, xValid$observed - xValid$residual, method = "pearson")^2)))
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

################### GA OPTIMIZATION #####################
gaOptim2 = function(optimData, nCluster, krigVar, popSize, generations, model) {
  coordinates(optimData) = ~x+y
  MSE = vector(length=nCluster)
  MSE2 = vector(length=nCluster)
  detCoefficients = vector(length=nCluster)
  GASolutions = matrix( ,nrow = nCluster, ncol = 3)
  columnIndex = which(colnames(optimData@data)==krigVar)
  
  maxSill = var(optimData@data[,columnIndex])
  maxSill = maxSill * 1
  maxNugget = 0
  distances = dist(optimData@coords)
  maxRange = max(distances)
  lags = 20
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
      
      #return (cor(xValid$observed, xValid$observed - xValid$residual, method = "pearson")^2)
      #return (-(1-(cor(xValid$observed, xValid$observed - xValid$residual, method = "pearson")^2)))
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

################### DE OPTIMIZATION #####################
deOptim = function(optimData, nCluster, krigVar, popSize, generations, model) {
  coordinates(optimData) = ~x+y
  MSE = vector(length=nCluster)
  MSE2 = vector(length=nCluster)
  detCoefficients = vector(length=nCluster)
  DESolutions = matrix( ,nrow = nCluster, ncol = 5)
  columnIndex = which(colnames(optimData@data)==krigVar)
  
  maxSill = var(optimData@data[,columnIndex])
  maxSill = maxSill * 10
  maxNugget = 0
  distances = dist(optimData@coords)
  maxRange = max(distances)
  print(maxSill)
  print(maxNugget)
  print(maxRange)
  
  lags = 20
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
      lzn.model = vgm(model, range = x1[3], psill = x1[2], nugget = x1[1], anis = c(x1[4], x1[5]))
      lzn.fit = fit.variogram(object = lzn.vgm, model = lzn.model, fit.sills = FALSE, fit.ranges = FALSE, fit.kappa = FALSE)
      xValid = krige.cv(var, optimBase, model=lzn.fit) 
      
      return (sum(xValid$residual^2))
    }
    
    if(nrow(optimBase) > 3) {
      DE = DEoptim(fitnessFunction, lower = c(0, 0, 0, 1, 0), upper = c(maxNugget, maxSill, maxRange, 180, 1), 
                   DEoptim.control(NP = popSize, itermax = generations, F = 0.8, CR = 0.5))
      
      DESolutions[i,] = DE$optim$bestmem
      lzn.vgm = variogram(var, optimBase) 
      lzn.model = vgm(range = DE$optim$bestmem[3], model, 
                      psill = DE$optim$bestmem[2], nugget = DE$optim$bestmem[1], 
                      anis = c(DE$optim$bestmem[4], DE$optim$bestmem[5]))
      lzn.fit = fit.variogram(object = lzn.vgm, model = lzn.model, fit.sills = FALSE, fit.ranges = FALSE, fit.kappa = FALSE)
      
      xValid = krige.cv(var, optimBase, model=lzn.fit)
      
      MSE[i] = sum(xValid$residual^2)
      MSE2[i] = sum(xValid$residual)
      detCoefficients[i] = cor(xValid$observed, xValid$observed - xValid$residual, method = "pearson")^2
    } else {
      MSE[i] = 0
      detCoefficients[i] = 0
    }
  }
  
  colnames(DESolutions) <- c("Nugget", "Sill", "Range", "Angle", "Factor")
  result = list(MSE, detCoefficients, DESolutions)
  names(result) <- c("MSE", "DetCoefficient", "VariogParams")
  return(result)
}


################### GEOR LEAST SQUARES METHOD OPTIMIZATION #####################
####SIGMA = SILL
####PHI = RANGE
####TAUSQ = NUGGET
lsOptim1 = function(optimData, nCluster, krigVar, model) {
  columnIndex = which(colnames(optimData)==krigVar)
  columnNumber = ncol(optimData)
  MSE = vector(length=nCluster)
  MSE2 = vector(length=nCluster)
  LSSolutions = matrix(, nrow = nCluster, ncol = 3)
  
  var = paste(krigVar, '~ 1')
  var = as.formula(var)
  
  for(i in 1:nCluster) {
    assign(paste("cluster", i, sep=""), optimData[optimData[,columnNumber]  %in% i,])
    assign(paste("cluster", i, sep=""), as.geodata(get(paste0("cluster", i)), coords.col = c(1,2), data.col = c(columnIndex)))
    
    temp = get(paste0("cluster", i))
    initialSill = var(temp$data)	
    print(initialSill)
    initialNugget = 0
    print(initialNugget)
    distances = dist(temp$coords)
    initialRange = max(distances)
    print(initialRange)
    
    e.variog <- variog(temp, uvec = seq(0,initialRange,l=10),  #Método para initialrange (ainda precisa de ajuste dependendo da base)
                       estimator.type = "modulus")
    plot(e.variog)
    initialValues = matrix(nrow = 2, ncol=2)
    initialValues[1,] = c(0.0001,0.1)
    initialValues[2,] = c(initialSill, initialRange/2)
    
    ls.fit = variofit(e.variog, ini.cov.pars = initialValues, cov.model = model, fix.nugget = TRUE)
    
    resultGEOR = geoR::xvalid(temp, model = ls.fit)
    
    LSparams = cbind(ls.fit$nugget, ls.fit$cov.pars[1])
    LSparams = cbind(LSparams, ls.fit$cov.pars[2])
    LSSolutions[i,] = LSparams		
    
    MSE[i] = sum(resultGEOR$error^2)
    MSE2[i] = sum(resultGEOR$error)
  }
  
  colnames(LSSolutions) = c("nugget", "sill", "range")
  result = list(MSE, MSE2, LSSolutions)
  names(result) <- c("MSE", "MSE2", "VariogParams")
  return(result)
}

################### GSTAT LEAST SQUARES METHOD OPTIMIZATION #####################
#Default initial parameter values are chosen from the sample variogram, where:
#the range parameter is taken as 1/3 of the maximum sample variogram distance,
#the nugget parameter is taken as the mean of the first three sample variogram values, and
#the partial sill is taken as the mean of the last five sample variogram values.

lsOptim2 = function(optimData, nCluster, krigVar, model) {
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
    
    lags = 20
    w = initialRange/lags
    
    lzn.vgm = variogram(var, optimBase, cutoff = initialRange, width = w) 
    lzn.fit = fit.variogram(lzn.vgm, vgm(psill = initialSill, model = model, range = initialRange/2, nugget = initialNugget), fit.sills = c(FALSE,TRUE))
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

################### GEOROB LEAST SQUARES METHOD OPTIMIZATION WITH ANISOTROPY #####################
#Default initial parameter values are chosen from the sample variogram, where:
#the range parameter is taken as 1/3 of the maximum sample variogram distance,
#the nugget parameter is taken as the mean of the first three sample variogram values, and
#the partial sill is taken as the mean of the last five sample variogram values.

lsOptim3 = function(optimData, nCluster, krigVar) {
  columnIndex = which(colnames(optimData)==krigVar)
  columnNumber = ncol(optimData)
  MSE = vector(length=nCluster)
  MSE2 = vector(length=nCluster)
  LSASolutions = matrix(, nrow = nCluster, ncol = 5)
  
  var = paste(krigVar, '~ 1')
  var = as.formula(var)
  
  for(i in 1:nCluster) {
    assign(paste("cluster", i, sep=""), optimData[optimData[,columnNumber]  %in% i,])
    
    temp = get(paste0("cluster", i))
    initialSill = var(temp[,3])	
    print(initialSill)
    initialNugget = 0
    print(initialNugget)
    distances = dist(temp[,c(1,2)])
    initialRange = max(distances)
    print(initialRange)
    
    ls.fit = georob(var, temp[,c(1,2,3)], locations = ~x+y, variogram.model = "RMexp",
                    param = c(variance = initialSill, nugget = 0.1, scale = initialRange/3),
                    fit.param = default.fit.param(scale = TRUE, variance = TRUE),
                    aniso = default.aniso(f1 = 0.50, f2 = 1, omega = 148.),
                    fit.aniso = default.fit.aniso(f1 = TRUE, f2 = TRUE, phi = TRUE, omega = TRUE), 
                    tuning.psi = 10000)
    
    ls.fit.xvalid = cv(object = ls.fit, seed = 1, lgn = false, nset = 1)
    summary(ls.fit.xvalid)
    
    LSAparams = cbind(0, ls.fit$variogram.object[[1]]$param[1])
    LSAparams = cbind(LSAparams, ls.fit$variogram.object[[1]]$param[4])
    LSAparams = cbind(LSAparams, ls.fit$variogram.object[[1]]$aniso[3])
    
    if(ls.fit$variogram.object[[1]]$aniso[1] > ls.fit$variogram.object[[1]]$aniso[2]) {
      LSAparams = cbind(LSAparams, ls.fit$variogram.object[[1]]$aniso[2]/ls.fit$variogram.object[[1]]$aniso[1])
    } else {
      LSAparams = cbind(LSAparams, ls.fit$variogram.object[[1]]$aniso[1]/ls.fit$variogram.object[[1]]$aniso[2])
    }
    LSASolutions[i,] = LSAparams
    
    coordinates(temp)=~x+y
    lzn.vgm = variogram(var, temp) 
    lzn.model = vgm("Exp", range = LSAparams[2], psill = LSAparams[3], nugget = LSAparams[1], anis=c(LSAparams[4], LSAparams[5]))
    lzn.fit = fit.variogram(object = lzn.vgm, model = lzn.model, fit.sills = FALSE, fit.ranges = FALSE, fit.kappa = FALSE)
    xValid = krige.cv(var, temp, model=lzn.fit) 
    
    MSE[i] = sum(xValid$residual^2)
    MSE2[i] = sum(xValid$residual)
  }
  
  colnames(LSASolutions) = c("nugget", "sill", "range","Angle","Factor")
  result = list(MSE, MSE2, LSASolutions)
  names(result) <- c("MSE", "MSE2", "VariogParams")
  return(result)
}


######################## LIKELIHOOD METHOD ##########################
lhOptim = function(optimData, nCluster, krigVar, model) {
  columnIndex = which(colnames(optimData)==krigVar)
  columnNumber = ncol(optimData)
  MSE = vector(length=nCluster)
  MSE2 = vector(length=nCluster)
  detCoefficients = vector(length=nCluster)
  LHSolutions = matrix(, nrow = nCluster, ncol = 5)
  
  var = paste(krigVar, '~ 1')
  var = as.formula(var)
  
  for(i in 1:nCluster) {
    assign(paste("cluster", i, sep=""), optimData[optimData[,columnNumber]  %in% i,])
    assign(paste("cluster", i, sep=""), as.geodata(get(paste0("cluster", i)), 
                                                   coords.col = c(1,2), data.col = c(columnIndex)))
    
    temp = get(paste0("cluster", i))
    initialSill = var(temp$data)	
    print(initialSill)
    initialNugget = 0
    print(initialNugget)
    distances = dist(temp$coords)
    initialRange = max(distances)
    print(initialRange)
    
    initialValues = matrix(nrow = 2, ncol=2)
    initialValues[1,] = c(0.0001,0.1)
    initialValues[2,] = c(initialSill, initialRange/2)
    
    #Verificar novamente ângulo e fator da anisotropia
    lh.fit = likfit(temp, ini.cov.pars = initialValues, fix.psiA = FALSE, fix.psiR = FALSE, fix.nugget = TRUE, nugget = 0)
    if(lh.fit$cov.pars[1] == 0) {
      lh.fit$cov.pars[1] = 0.0001
    }
    if(lh.fit$cov.pars[2] == 0) {
      lh.fit$cov.pars[2] = 0.0001
    }
    resultGEOR = geoR::xvalid(temp, mode = lh.fit)
    
    radius = lh.fit$aniso.pars[1]
    angle = NISTradianTOdeg(radius)
    summary(lh.fit)
    
    #Transformando geoR ratio em gStat ratio
    maxDist = lh.fit$max.dist
    tempRatio = maxDist/lh.fit$aniso.pars[2]
    finalRatio = tempRatio/maxDist
    
    LHparams = cbind(lh.fit$nugget, lh.fit$cov.pars[1])
    LHparams = cbind(LHparams, lh.fit$cov.pars[2])
    LHparams = cbind(LHparams, angle)
    LHparams = cbind(LHparams, finalRatio)
    LHSolutions[i,] = LHparams	
    print(LHSolutions)
    
    MSE[i] = sum(resultGEOR$error^2)
    MSE2[i] = sum(resultGEOR$error)
  }
  
  colnames(LHSolutions) = c("nugget", "sill", "range", "angle", "factor")
  result = list(MSE, MSE2, LHSolutions)
  names(result) <- c("MSE", "MSE2", "VariogParams")
  return(result)
}

################### ITERATIVE LEAST SQUARES METHOD OPTIMIZATION #####################
#http://rgeostats.free.fr/doc/Autofit.html
#http://rgeostats.free.fr/doc/2D.R
#degrees = (radians * 180) / pi
#radians = (degrees * pi) / 180
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
    
    variogram = vario.calc(temp, nlag=20)
    
    ########### Modelo por LS
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

################### ITERATIVE LEAST SQUARES METHOD OPTIMIZATION #####################
#http://rgeostats.free.fr/doc/Autofit.html
#http://rgeostats.free.fr/doc/2D.R
#degrees = (radians * 180) / pi
#radians = (degrees * pi) / 180
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
    
    variogram = vario.calc(temp, nlag=20, dirvect = c(0,45,90,135,180))
    
    ########### Modelo por LS
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

##############################################################################
##############################################################################
#############################   COMMENTS   ###################################
##############################################################################
##############################################################################
##############################################################################

#the initial sill for each basic structure is set by 
#default to a value equal to the total variance 
#divided by the number of structures;

#The initial range for each basic structure is set by  
#default to half of the maximum distance divided 
#by the number of basic structures (the nugget effect component is discarded).  
#The maximum distance is obtained as the longest distance for which the experimental 
#variogram has been calculated; 

#Initial anisotropy is set to 1 (Isotropic Hyphotesis)

#Angle is arbitrarily set to 0

# uses a regression tree with a ???xed amount of leaf nodes to partition the data in the objective space
#IMPLEMENT TREE-BASED CLUSTERING package: TREECLUST (Model Tree Cluster Kriging (MTCK))
#MTCK uses only one of the trained Kriging models per unseen record to predict 

####################### PARAMS ##########################

# ####################### MEUSE DATABASE #######################
# data(meuse)
# data(meuse.grid)
# data = meuse[,c(1,2,6)] # Selection of data variables
# data.grid = meuse.grid
# #nCluster = 1 #Kmeans number of clusters
# normalizationFactor = 1 #Divisão da variável resposta
# nPoints = 6  #DBscan number of points
# epsilon = 250  #DBscan epsilon circuference
# trainDataSize = 1  #Percentage of train Data
# nNeighbours = 4  #Number of Neighbours in cluster assignment
# var = "zinc"
# gaPopulation = 100
# gaIter = 5
# model1 = "exponential"
# model2 = "Exp"
# model3 = "Exponential"
# #g1 = kmeansKneeGraph(data, "DATABASE: Meuse")

####################### BROOMSBARN DATABASE #######################
train = read.table('BroomsBarn.dat', skip=4, dec = ".")
colnames(train) = c("x","y","K","log10K","pH","P","log10P")
data = train[,c(1,2,3)]
x = rep(1:18, each = 30)
y = rep(1:30, 18)
data.grid = cbind(x,y)
data.grid = as.data.frame(data.grid)
nCluster = 7 #Kmeans number of clusters
nPoints = 4  #DBscan number of points
normalizationFactor = 1 #Divisão da variável resposta
epsilon = 250  #DBscan epsilon circuference
trainDataSize = 1  #Percentage of train Data
nNeighbours = 4  #Number of Neighbours in cluster assignment
var = "K"
gaPopulation = 50
gaIter = 3
model1 = "exponential"
model2 = "Exp"
model3 = "Exponential"
#g2 = kmeansKneeGraph(data, "DATABASE: Broomsbarn")

####################### COALASH DATABASE #######################
# data(coalash)
# train = coalash
# colnames(train) = c('x','y','coalash')
# data = train[,c(1:3)]
# x = rep(1:16, each = 23)
# y = rep(1:23, 16)
# data.grid = cbind(x,y)
# data.grid = as.data.frame(data.grid)
# nCluster = 7  #Kmeans number of clusters
# normalizationFactor = 1 #Divisão da variável resposta
# trainDataSize = 1  #Percentage of train Data
# nNeighbours = 4  #Number of Neighbours in cluster assignment
# var = "coalash"
# gaPopulation = 100
# gaIter = 5
# model1 = "exponential"
# model2 = "Exp"
# model3 = "Exponential"
# #g3 = kmeansKneeGraph(data, "DATABASE: Coalash")

####################### WOLFCAMP #######################
# train = read.table('WOLFCAMP.DAT',header=TRUE, skip=3)
# colnames(train) = c('x','y','plevel')
# train = rbind(c(42.8,127.6,1464.0), train)
# data = train[,c(1:3)]
# #nCluster = 1 #Kmeans number of clusters
# normalizationFactor = 1 #Divisão da variável resposta
# trainDataSize = 1  #Percentage of train Data
# nNeighbours = 4  #Number of Neighbours in cluster assignment
# var = "plevel"
# gaPopulation = 100
# gaIter = 5
# model1 = "exponential"
# model2 = "Exp"
# model3 = "Exponential"
# # g4 = kmeansKneeGraph(data, "DATABASE: Wolfcamp")


####################### LAKE WALKER DATA #######################
# data(walker)
# train = cbind(walker@coords,walker@data$V)
# train = as.data.frame(train)
# colnames(train) = c("x","y","V")
# data = train[,c(1:3)]
# #nCluster = 3 #Kmeans number of clusters
# trainDataSize = 1  #Percentage of train Data
# nNeighbours = 4  #Number of Neighbours in cluster assignment
# var = "V"
# gaPopulation = 100
# gaIter = 5
# model1 = "exponential"
# model2 = "Exp"
# model3 = "Exponential"
# # g5 = kmeansKneeGraph(data, "DATABASE: Lake Walker")

# ####################### SINTETIC DATA #######################
# train = read.table('sint-data1.csv',sep=',',dec='.',header=TRUE)
# data = train[,c(1,2,5)]
# data.grid = data[,c(1,2)]
# nCluster = 7  #Kmeans number of clusters
# trainDataSize = 0.10  #Percentage of train Data
# nNeighbours = 3  #Number of Neighbours in cluster assignment
# var = "QuadraticY"
# gaPopulation = 50
# gaIter = 2
# model1 = "exponential"
# model2 = "Exp"
# model3 = "Exponential"
# #g6 = kmeansKneeGraph(data, "DATABASE: Sintetic")

####################### AGROMET DATA #######################
# train = read.table('Agromet.dat',skip=3, sep=",", dec = ".")
# train = cbind(train[,c(2,3)], train[,1])
# colnames(train) = c("x","y","height")
# data = train[,c(1:3)]
# data.grid = data[,c(1,2)]
# nCluster = 7  #Kmeans number of clusters
# trainDataSize = 0.05  #Percentage of train Data
# nNeighbours = 3  #Number of Neighbours in cluster assignment
# var = "height"
# gaPopulation = 50
# gaIter = 1
# model1 = "exponential"
# model2 = "Exp"
# model3 = "Exponential"
# #g7 = kmeansKneeGraph(data, "DATABASE: Agromet")

# grid.arrange(g1,g2,g3,g4,g5,g6,g7, ncol=3)
#################### MAIN EXEC ##########################
sink()
sink()
## DATA PREPROCESSING + DETECTING AND REMOVING OUTLIERS
scores = scores(data[,var], type = "z", prob=0.99)
print(scores)
print(data[scores==TRUE,])
data = data[scores==FALSE,]

for(i in 1:ncol(data)) {
  data[,i] = rescale(data[,i], to=c(-1,1))
}

## Train and test data partitioning ##
partitionedData = partitionData(data, trainDataSize)
train = partitionedData[partitionedData$train %in% 1,]
test = partitionedData[partitionedData$train %in% 0,]

nTests = 5
fileName = "broomsbarn 15.txt"
nCluster = 1

GA1ResultsMeans = matrix(nrow = nTests, ncol = 2)
# GA2ResultsMeans = matrix(nrow = nTests, ncol = 2)
# LS1ResultsMeans = matrix(nrow = nTests, ncol = 2)
# LS2ResultsMeans = matrix(nrow = nTests, ncol = 2)
# LHResultsMeans = matrix(nrow = nTests, ncol = 2)
# ILSResultsMeans = matrix(nrow = nTests, ncol = 2)
# ILS2ResultsMeans = matrix(nrow = nTests, ncol = 2)

GA1ResultsSd = matrix(nrow = nTests, ncol = 2)
# GA2ResultsSd = matrix(nrow = nTests, ncol = 2)
# LS1ResultsSd = matrix(nrow = nTests, ncol = 2)
# LS2ResultsSd = matrix(nrow = nTests, ncol = 2)
# LHResultsSd = matrix(nrow = nTests, ncol = 2)
# ILSResultsSd = matrix(nrow = nTests, ncol = 2)
# ILS2ResultsSd = matrix(nrow = nTests, ncol = 2)

for(z in 1:nCluster) {
  GA1Results = matrix(nrow = nTests, ncol = 2)
  # GA2Results = matrix(nrow = nTests, ncol = 2)
  # LS1Results = matrix(nrow = nTests, ncol = 2)
  # LS2Results = matrix(nrow = nTests, ncol = 2)
  # LHResults = matrix(nrow = nTests, ncol = 2)
  # ILSResults = matrix(nrow = nTests, ncol = 2)
  # ILS2Results = matrix(nrow = nTests, ncol = 2)
  
  sink(fileName, type="output", append = TRUE)
  cat("CLUSTER ---------- ")
  cat(z)
  cat("\n")
  cat("            GA --- GA2 --- LS1 --- LS2 --- LH --- ILS --- ILS2     GA --- GA2 --- LS1 --- LS2 --- LH --- ILS --- ILS2")
  cat("\n")
  sink()
  
  ## Clustering Method ##
  for(i in 1:nTests) {
    #set.seed(i)
    
    clusteredData = kmeansClustering(train, z)
    ## Cluster Assignment to each data ##
    # a = ggplot(clusteredData, aes(x=x, y=y, size=K, color = as.factor(cluster))) +
    #   geom_point(alpha=0.5) +
    #   labs(title = "Broomsbarn Database - 6 Clusters Without KNN")
    # a + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan","blue","red","yellow")) + theme(legend.position="none")
    
    assignClusterData = assignCluster(clusteredData, nNeighbours)
    # b = ggplot(assignClusterData, aes(x=x, y=y, size=K, color = clusterSelected)) +
    #   geom_point(alpha=0.5) +
    #   labs(title = "Broomsbarn Database - 6 Clusters with KNN (3 Neigh.)")
    # b + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","cyan","blue","red","yellow"))  + theme(legend.position="none")
    plot(x = assignClusterData$x, y = assignClusterData$y, col = assignClusterData$cluster) 
    
    ## Variogram Optimization ##
    result1 = try(gaOptim(optimData = assignClusterData, nCluster = z, krigVar = var, popSize = gaPopulation, generations =  gaIter, model = model2))
    # result2 = try(gaOptim2(optimData = assignClusterData, nCluster = z, krigVar = var, popSize = gaPopulation, generations =  gaIter, model = model2))
    # result3 = try(lsOptim1(optimData = assignClusterData, nCluster = z, krigVar = var, model = model1))
    # result4 = try(lsOptim2(optimData = assignClusterData, nCluster = z, krigVar = var, model = model2))
    # result5 = try(lhOptim(optimData = assignClusterData, nCluster = z, krigVar = var, model = model2))
    # result6 = try(iterativeLsOptim2(optimData = assignClusterData, nCluster = z, krigVar = var, model = model2, model2 = model3))
    # result7 = try(iterativeLsOptim(optimData = assignClusterData, nCluster = z, krigVar = var, model = model2, model2 = model3))
    
    sink(fileName, type="output", append = TRUE)
    
    cat("Iteração: ")
    cat(i)
    cat("  ")
    
    variance = var(assignClusterData[,3])
    NMSEfactor = 1/(variance*nrow(assignClusterData))
    try(cat(round(sum(result1$MSE*NMSEfactor),4)))
    cat(" ")
    GA1Results[i,1] = try(sum(result1$MSE*NMSEfactor))

    # try(cat(round(sum(result2$MSE*NMSEfactor),4)))
    # cat(" ")
    # GA2Results[i,1] = try(sum(result2$MSE*NMSEfactor))
    # 
    # try(cat(round(sum(result3$MSE*NMSEfactor),4)))
    # cat(" ")
    # LS1Results[i,1] = try(sum(result3$MSE*NMSEfactor))
    # 
    # try(cat(round(sum(result4$MSE*NMSEfactor),4)))
    # cat(" ")
    # LS2Results[i,1] = try(sum(result4$MSE*NMSEfactor))
    # 
    # try(cat(round(sum(result5$MSE*NMSEfactor),4)))
    # cat(" ")
    # LHResults[i,1] = try(sum(result5$MSE*NMSEfactor))
    # 
    # try(cat(round(sum(result6$MSE*NMSEfactor),4)))
    # cat(" ")
    # ILSResults[i,1] = try(sum(result6$MSE*NMSEfactor))
    # 
    # try(cat(round(sum(result7$MSE*NMSEfactor),4)))
    # cat(" ")
    # ILS2Results[i,1] = try(sum(result7$MSE*NMSEfactor))
    # cat("   ")
    
    mean = mean(assignClusterData[,3])
    PAEEfactor = 1/(mean*nrow(assignClusterData))
    try(cat(round(sum(result1$MSE2*PAEEfactor),4)))
    cat(" ")
    GA1Results[i,2] = try(sum(result1$MSE2*PAEEfactor))

    # try(cat(round(sum(result2$MSE2*PAEEfactor),4)))
    # cat(" ")
    # GA2Results[i,2] = try(sum(result2$MSE2*PAEEfactor))
    # 
    # try(cat(round(sum(result3$MSE2*PAEEfactor),4)))
    # cat(" ")
    # LS1Results[i,2] = try(sum(result3$MSE2*PAEEfactor))
    # 
    # try(cat(round(sum(result4$MSE2*PAEEfactor),4)))
    # cat(" ")
    # LS2Results[i,2] = try(sum(result4$MSE2*PAEEfactor))
    # 
    # try(cat(round(sum(result5$MSE2*PAEEfactor),4)))
    # cat(" ")
    # LHResults[i,2] = try(sum(result5$MSE2*PAEEfactor))
    # 
    # try(cat(round(sum(result6$MSE2*PAEEfactor),4)))
    # cat(" ")
    # ILSResults[i,2] = try(sum(result6$MSE2*PAEEfactor))
    # 
    # try(cat(round(sum(result7$MSE2*PAEEfactor),4)))
    # cat(" ")
    # ILS2Results[i,2] = try(sum(result7$MSE2*PAEEfactor))
    
    cat("\n")
    if(i == nTests) {
      cat("Média: ")
      cat("  ")
      try(cat(round(mean(GA1Results[,1]),4)))
      GA1ResultsMeans[z,1] = try(round(mean(GA1Results[,1]),4))
      cat(" ")
      # try(cat(round(mean(GA2Results[,1]),4)))
      # GA2ResultsMeans[z,1] = try(round(mean(GA2Results[,1]),4))
      # cat(" ")
      # try(cat(round(mean(LS1Results[,1]),4)))
      # LS1ResultsMeans[z,1] = try(round(mean(LS1Results[,1]),4))
      # cat(" ")
      # try(cat(round(mean(LS2Results[,1]),4)))
      # LS2ResultsMeans[z,1] = try(round(mean(LS2Results[,1]),4))
      # cat(" ")
      # try(cat(round(mean(LHResults[,1]),4)))
      # LHResultsMeans[z,1] = try(round(mean(LHResults[,1]),4))
      # cat(" ")
      # try(cat(round(mean(ILSResults[,1]),4)))
      # ILSResultsMeans[z,1] = try(round(mean(ILSResults[,1]),4))
      # cat(" ")
      # try(cat(round(mean(ILS2Results[,1]),4)))
      # ILS2ResultsMeans[z,1] = try(round(mean(ILS2Results[,1]),4))
      # cat(" ")
      cat("   ")
      try(cat(round(mean(GA1Results[,2]),4)))
      GA1ResultsMeans[z,2] = try(round(mean(GA1Results[,2]),4))
      cat(" ")
      # try(cat(round(mean(GA2Results[,2]),4)))
      # GA2ResultsMeans[z,2] = try(round(mean(GA2Results[,2]),4))
      # cat(" ")
      # try(cat(round(mean(LS1Results[,2]),4)))
      # LS1ResultsMeans[z,2] = try(round(mean(LS1Results[,2]),4))
      # cat(" ")
      # try(cat(round(mean(LS2Results[,2]),4)))
      # LS2ResultsMeans[z,2] = try(round(mean(LS1Results[,2]),4))
      # cat(" ")
      # try(cat(round(mean(LHResults[,2]),4)))
      # LHResultsMeans[z,2] = try(round(mean(LHResults[,2]),4))
      # cat(" ")
      # try(cat(round(mean(ILSResults[,2]),4)))
      # ILSResultsMeans[z,2] = try(round(mean(ILSResults[,2]),4))
      # cat(" ")
      # try(cat(round(mean(ILS2Results[,2]),4)))
      # ILS2ResultsMeans[z,2] = try(round(mean(ILS2Results[,2]),4))
      
      cat("\n")
      
      cat("DesvPadr: ")
      cat("  ")
      try(cat(round(sd(GA1Results[,1]),4)))
      GA1ResultsSd[z,1] = try(round(sd(GA1Results[,1]),4))
      cat(" ")
      # try(cat(round(sd(GA2Results[,1]),4)))
      # GA2ResultsSd[z,1] = try(round(sd(GA2Results[,1]),4))
      # cat(" ")
      # try(cat(round(sd(LS1Results[,1]),4)))
      # LS1ResultsSd[z,1] = try(round(sd(LS1Results[,1]),4))
      # cat(" ")
      # try(cat(round(sd(LS2Results[,1]),4)))
      # LS2ResultsSd[z,1] = try(round(sd(LS2Results[,1]),4))
      # cat(" ")
      # try(cat(round(sd(LHResults[,1]),4)))
      # LHResultsSd[z,1] = try(round(sd(LHResults[,1]),4))
      # cat(" ")
      # try(cat(round(sd(ILSResults[,1]),4)))
      # ILSResultsSd[z,1] = try(round(sd(ILSResults[,1]),4))
      # cat(" ")
      # try(cat(round(sd(ILS2Results[,1]),4)))
      # ILS2ResultsSd[z,1] = try(round(sd(ILS2Results[,1]),4))
      # cat(" ")
      cat("   ")
      try(cat(round(sd(GA1Results[,2]),4)))
      GA1ResultsSd[z,2] = try(round(sd(GA1Results[,2]),4))
      cat(" ")
      # try(cat(round(sd(GA2Results[,2]),4)))
      # GA2ResultsSd[z,2] = try(round(sd(GA2Results[,2]),4))
      # cat(" ")
      # try(cat(round(sd(LS1Results[,2]),4)))
      # LS1ResultsSd[z,2] = try(round(sd(LS1Results[,2]),4))
      # cat(" ")
      # try(cat(round(sd(LS2Results[,2]),4)))
      # LS2ResultsSd[z,2] = try(round(sd(LS2Results[,2]),4))
      # cat(" ")
      # try(cat(round(sd(LHResults[,2]),4)))
      # LHResultsSd[z,2] = try(round(sd(LHResults[,2]),4))
      # cat(" ")
      # try(cat(round(sd(ILSResults[,2]),4)))
      # ILSResultsSd[z,2] = try(round(sd(ILSResults[,2]),4))
      # cat(" ")
      # try(cat(round(sd(ILS2Results[,2]),4)))
      # ILS2ResultsSd[z,2] = try(round(sd(ILS2Results[,2]),4))
      # 
      cat("\n")
    }
    sink()
  }
}
