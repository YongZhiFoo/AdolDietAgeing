

# A function created to find the outer perimeter over which the surface should be fitted for proportional data
findConvex.prop<-function(x,y,rgnames,res=101){
  hull<-cbind(x,y)[chull(cbind(x,y)),]
  x.new<-seq(0,1,len=res)
  y.new<-seq(0,1,len=res)
  ingrid<-as.data.frame(expand.grid(x.new,y.new))                                                              
  Fgrid<-ingrid
  Fgrid[(point.in.polygon(ingrid[,1], ingrid[,2], hull[,1],hull[,2])==0),]<-NA
  names(Fgrid)<-rgnames
  return(Fgrid)
}



#---------------------------------
#---- special functions
#---------------------------------
#--- function for making subdiagonals (from Bill Venables)
  subdiag <- function (v, k) {
    n <- length(v) + abs(k)
    x <- matrix(0, n, n)
    if (k == 0)
        diag(x) <- v
    else if (k < 0)
      { ## sub-diagonal
        j <- 1:(n+k)
        i <- (1 - k):n
        x[cbind(i, j)] <- v
} 
x } 



# FUnction to work out whether points in nutrient spaace fall within the hull given by the raw data

inhull <- function(testpts, calpts, hull=convhulln(calpts), tol=mean(mean(abs(calpts)))*sqrt(.Machine$double.eps)) { 
  
  require(sp)
  require(geometry)
  
  # https://tolstoy.newcastle.edu.au/R/e8/help/09/12/8784.html
  calpts <- as.matrix(calpts) 
  testpts <- as.matrix(testpts) 
  p <- dim(calpts)[2] 
  cx <- dim(testpts)[1] # rows in testpts
  nt <- dim(hull)[1] # number of simplexes in hull 
  nrmls <- matrix(NA, nt, p)
  
  degenflag <- matrix(TRUE, nt, 1) 
  for (i in 1:nt){ 
    nullsp<-t(Null(t(calpts[hull[i,-1],] - matrix(calpts[hull[i,1],],p-1,p, byrow=TRUE))))
    if (dim(nullsp)[1] == 1){
      nrmls[i,]<-nullsp
      degenflag[i]<-FALSE
    }
  }
  
  if(length(degenflag[degenflag]) > 0) warning(length(degenflag[degenflag])," degenerate faces in convex hull")
  nrmls <- nrmls[!degenflag,] 
  nt <- dim(nrmls)[1] 
  
  center = apply(calpts, 2, mean) 
  a<-calpts[hull[!degenflag,1],] 
  nrmls<-nrmls/matrix(apply(nrmls, 1, function(x) sqrt(sum(x^2))), nt, p)
  
  dp <- sign(apply((matrix(center, nt, p, byrow=TRUE)-a) * nrmls, 1, sum))
  nrmls <- nrmls*matrix(dp, nt, p)
  
  aN <- diag(a %*% t(nrmls)) 
  val <- apply(testpts %*% t(nrmls) - matrix(aN, cx, nt, byrow=TRUE), 1,min) 
  
  val[abs(val) < tol] <- 0 
  as.integer(sign(val)) 
}



ggSurface<-function(GAM, data, XYZ, labels, predict_val, surf_min=NA, surf_max=NA, x.limits=NA, y.limits=NA, z.val=NA, exclude=NULL, subtitle="", lab_size=3, nlevels=5, contour_at=NA, skip=0){
  
  require(ggplot2)
  require(sp)
  require(geometry)
  require(mgcv)
  require(metR)
  
  # This specifies the color scheme for surface
  rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
  map<-rgb.palette(256)
  
  # List for the order of plots
  nutrient.order<-XYZ[c(1,2,3)]
  
  # List for the labels
  labels.order<-labels[c(1,2,3)]
  
  # Values to predict over, if they are unspecified
  if(is.na(x.limits)[1] == T){
    x.limits<-c(floor(min(data[,nutrient.order[1]])), ceiling(max(data[,nutrient.order[1]])))
  }
  if(is.na(y.limits)[1] == T){
    y.limits<-c(floor(min(data[,nutrient.order[2]])), ceiling(max(data[,nutrient.order[2]])))
  }		
  
  # If we do not specify values to slice at, use the 25, 50, and 75 %ile
  if(is.na(z.val) == T){
    z.val<-round(median(data[,nutrient.order[3]]))
  }
  
  # Fitted list to hold some results for later
  x.new<-seq(min(x.limits, na.rm=T), max(x.limits, na.rm=T), len=501)
  y.new<-seq(min(y.limits, na.rm=T), max(y.limits, na.rm=T), len=501)
  z.new<-z.val
  predictors<-as.data.frame(expand.grid(x.new, y.new, z.new))
  names(predictors)<-nutrient.order
  in.poly<-as.numeric(inhull(predictors[,c(1:3)], data[,names(predictors)]) != -1)
  
  # Add the predictors for the additional 'confounders'
  predictors<-cbind(predictors, predict_val)
  predictors<-predictors[-which(in.poly == 0),]
  
  # Unless fitted values are directly provided do the predictions
  predictions<-predict(GAM, newdata=predictors, type="response", exclude=exclude, newdata.guaranteed=T)
  
  # Get the proedictions for the kth trait					
  predictions_k<-predictions
  
  # Find the min and max values across all predictions
  mn<-surf_min
  mx<-surf_max
  if(is.na(mn)==T){
    mn<-min(predictions_k, na.rm=T)
  }
  if(is.na(mx)==T){
    mx<-max(predictions_k, na.rm=T)
  }
  locs<-(range(predictions_k, na.rm=TRUE) - mn) / (mx-mn) * 256	
  
  plot_data<-predictors
  plot_data$fit<-predictions_k
  plot_data$x<-plot_data[,nutrient.order[1]]
  plot_data$y<-plot_data[,nutrient.order[2]]
  
  # Set the contour
  if(is.na(contour_at)[1] == T){
    contour_use<-signif((max(predictions_k, na.rm=T)-min(predictions_k, na.rm=T))/nlevels, 1)
  }else{
    contour_use<-contour_at	
  }
  
  # Make the plot
  plot<-ggplot(plot_data, aes(x=x, y=y)) +
    geom_raster(aes(fill=fit), show.legend=F, interpolate=F, na.rm=T) +
    scale_fill_gradientn(colors=map[locs[1]:locs[2]]) +
    geom_contour(data=plot_data, aes(x=x, y=y, z=fit), na.rm=T, color="black", binwidth=contour_use) +	
    geom_label_contour(data=plot_data, aes(x=x, y=y, z=fit), size=lab_size, binwidth=contour_use, skip=skip) +
    theme_bw() +
    labs(x = labels.order[1], y = labels.order[2], subtitle=subtitle) +
    theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
    theme(title=element_text(size=15)) + 
    xlim(x.limits) +
    ylim(y.limits)
  
  return(plot)
  
}

# Function to make the 9-way surface with gg aesthetics based on inputted model and data
gg9waySurface<-function(GAM, data, slices, XYZ, labels, predict_val, ranges, contour_at, pdf_name, exclude=NULL, XYZ_limits){
  
  require(gridExtra)	
  
  # Plot the surfaces
  plots_list<-list()
  
  # Make the surfaces for fat on the third
  for(i in 1:3){
    plots_list[[i]]<-ggSurface(GAM = GAM, data=data, XYZ=XYZ, labels=labels, predict_val= predict_val, exclude = exclude, z.val = slices[1,i], subtitle = paste0(labels[3], " = ", slices[1,i]), surf_min=ranges[1], surf_max=ranges[2], contour_at=contour_at, x.limits=XYZ_limits[[1]], y.limits=XYZ_limits[[2]])
  }
  # Now for carb
  for(i in 1:3){
    plots_list[[i+3]]<-ggSurface(GAM = GAM, data=data, XYZ=XYZ[c(1,3,2)], labels=labels[c(1,3,2)], predict_val= predict_val, exclude = exclude, z.val = slices[2,i], subtitle = paste0(labels[2], " = ", slices[2,i]), surf_min=ranges[1], surf_max=ranges[2], contour_at=contour_at, x.limits=XYZ_limits[[1]], y.limits=XYZ_limits[[3]])
  }
  # Now for protein
  for(i in 1:3){
    plots_list[[i+6]]<-ggSurface(GAM = GAM, data=data, XYZ=XYZ[c(2,3,1)], labels=labels[c(2,3,1)], predict_val= predict_val, exclude = exclude, z.val = slices[3,i], subtitle = paste0(labels[1], " = ", slices[3,i]), surf_min=ranges[1], surf_max=ranges[2], contour_at=contour_at, x.limits=XYZ_limits[[2]], y.limits=XYZ_limits[[3]])
  }
  
  
  pdf(pdf_name, height=12, width=12)
  
  grid.arrange(plots_list[[1]], 
               plots_list[[2]], 
               plots_list[[3]], 
               plots_list[[4]], 
               plots_list[[5]], 
               plots_list[[6]], 
               plots_list[[7]], 
               plots_list[[8]], 
               plots_list[[9]], layout_matrix=rbind(c(1,2,3), c(4,5,6), c(7,8,9)))
  
  dev.off()
  
  return(plots_list)
  
}

########## Yong's code below

# Function to make one surface that captures the full range of the 2 nutrients. Because the individual nutrient intakes are positively correlated with totak caloric intake, we lose the outer ends of the surface when we slice the plot based on the third nutrient in ggSurface. So here we "flatten" the third nutrient instead. Hence ggSurfaceflat
ggSurfaceflat<-function(GAM, data, XYZ, labels, predict_val, surf_min=NA, surf_max=NA, x.limits=NA, y.limits=NA, z.val=NA, exclude=NULL, subtitle="", lab_size=3, nlevels=5, contour_at=NA, skip=0){
  
  require(ggplot2)
  require(sp)
  require(geometry)
  require(mgcv)
  require(metR)
  
  # This specifies the color scheme for surface
  rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
  map<-rgb.palette(256)
  
  # List for the order of plots
  nutrient.order<-XYZ[c(1,2,3)]
  
  # List for the labels
  labels.order<-labels[c(1,2,3)]
  
  # Values to predict over, if they are unspecified
  if(is.na(x.limits)[1] == T){
    x.limits<-c(floor(min(data[,nutrient.order[1]])), ceiling(max(data[,nutrient.order[1]])))
  }
  if(is.na(y.limits)[1] == T){
    y.limits<-c(floor(min(data[,nutrient.order[2]])), ceiling(max(data[,nutrient.order[2]])))
  }		
  
  # If we do not specify values to slice at, use the 25, 50, and 75 %ile
  if(is.na(z.val) == T){
    z.val<-round(median(data[,nutrient.order[3]]))
  }
  
  # Fitted list to hold some results for later
  x.new<-seq(min(x.limits, na.rm=T), max(x.limits, na.rm=T), len=501)
  y.new<-seq(min(y.limits, na.rm=T), max(y.limits, na.rm=T), len=501)
  z.new<-z.val
  predictors<-as.data.frame(expand.grid(x.new, y.new, z.new))
  names(predictors)<-nutrient.order
  in.poly<-as.numeric(inhull(predictors[,c(1:2)], data[,names(predictors)[1:2]]) != -1)
  
  # Add the predictors for the additional 'confounders'
  predictors<-cbind(predictors, predict_val)
  predictors<-predictors[-which(in.poly == 0),]
  
  # Unless fitted values are directly provided do the predictions
  predictions<-predict(GAM, newdata=predictors, type="response", newdata.guaranteed=T)
  
  # Get the proedictions for the kth trait					
  predictions_k<-predictions
  
  # Find the min and max values across all predictions
  mn<-surf_min
  mx<-surf_max
  if(is.na(mn)==T){
    mn<-min(predictions_k, na.rm=T)
  }
  if(is.na(mx)==T){
    mx<-max(predictions_k, na.rm=T)
  }
  locs<-(range(predictions_k, na.rm=TRUE) - mn) / (mx-mn) * 256	
  
  plot_data<-predictors
  plot_data$fit<-predictions_k
  plot_data$x<-plot_data[,nutrient.order[1]]
  plot_data$y<-plot_data[,nutrient.order[2]]
  
  # Set the contour
  if(is.na(contour_at)[1] == T){
    contour_use<-signif((max(predictions_k, na.rm=T)-min(predictions_k, na.rm=T))/nlevels, 1)
  }else{
    contour_use<-contour_at	
  }
  
  # Make the plot
  plot<-ggplot(plot_data, aes(x=x, y=y)) +
    geom_raster(aes(fill=fit), show.legend=F, interpolate=F, na.rm=T) +
    scale_fill_gradientn(colors=map[locs[1]:locs[2]]) +
    geom_contour(data=plot_data, aes(x=x, y=y, z=fit), na.rm=T, color="black", binwidth=contour_use) +	
    geom_label_contour(data=plot_data, aes(x=x, y=y, z=fit), size=lab_size, binwidth=contour_use, skip=skip) +
    theme_bw() +
    labs(x = labels.order[1], y = labels.order[2]
         #, subtitle=paste0(labels[3], " = ", z.val)
         ) +
    theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
    theme(title=element_text(size=15)) + 
    xlim(x.limits) +
    ylim(y.limits)
  
  return(plot)
  
}

# Function to make an animation of the 3D surface going from the low endc of one nutrient to the high
ggSurfacetomo<-function(GAM, data, XYZ, labels, predict_val, ranges, slicesnum, x.limits=NA, y.limits=NA, exclude=NULL, subtitle="", lab_size=3, nlevels=5, contour_at=NA, skip=0){

# Define the start and end numbers
start_num <- ceiling(min(data[XYZ[3]]))
end_num <- floor(max(data[XYZ[3]]))

# Create the sequence of numbers based on the number of slices specified
my_sequence <- seq(from = start_num, to = end_num, length.out = slicesnum)

# Plot the surfaces
plots_list<-list()

# Make the surfaces for fat on the third
for(i in 1:slicesnum){
  plot_result<-tryCatch({ggSurface(GAM = GAM, 
                                   data=data, 
                                   XYZ=XYZ, 
                                   labels=labels, 
                                   predict_val= predict_val, 
                                   exclude = NA, 
                                   z.val = my_sequence[i], 
                                   surf_min=ranges[1], 
                                   surf_max=ranges[2], 
                                   contour_at=NA, 
                                   x.limits=NA, 
                                   y.limits=NA)}, 
                        error = function(e) {
                          message(paste("Error creating plot", i, ":", e$message))
                        })
  if (!is.null(plot_result)) {
    plots_list[[i]] <- plot_result
  }
  
}

return(plots_list)

plotnum <- length(plots_list)

for(i in 1:plotnum){
  
  print(plots_list[[i]])
  
}

  
  
}

#Getting the peak intake

macropeak<-function(GAM, data, XYZ, labels, predict_val,  exclude=NULL){
  
  require(sp)
  require(geometry)
  require(mgcv)
  require(metR)
  
  # List for the order of plots
  nutrient.order<-XYZ[c(1,2,3)]
  
  # List for the labels
  labels.order<-labels[c(1,2,3)]
  
  # Values to predict over, if they are unspecified
  x.limits<-c(floor(min(data[,nutrient.order[1]])), ceiling(max(data[,nutrient.order[1]])))
  y.limits<-c(floor(min(data[,nutrient.order[2]])), ceiling(max(data[,nutrient.order[2]])))
  z.limits<-c(floor(min(data[,nutrient.order[1]])), ceiling(max(data[,nutrient.order[3]])))
  
  # Fitted list to hold some results for later
  x.new<-seq(min(x.limits, na.rm=T), max(x.limits, na.rm=T), len=200)
  y.new<-seq(min(y.limits, na.rm=T), max(y.limits, na.rm=T), len=200)
  z.new<-seq(min(z.limits, na.rm=T), max(z.limits, na.rm=T), len=200)
  predictors<-as.data.frame(expand.grid(x.new, y.new, z.new))
  names(predictors)<-nutrient.order
  in.poly<-as.numeric(inhull(predictors[,c(1:3)], data[,names(predictors)]) != -1)
  
  # Add the predictors for the additional 'confounders'
  predictors<-cbind(predictors, predict_val)
  predictors<-predictors[-which(in.poly == 0),]
  
  # Unless fitted values are directly provided do the predictions
  predictions<-predict(GAM, newdata=predictors, type="response", exclude=exclude, newdata.guaranteed=T)
  
  # Find peak predicted outcome
  peak<-which(predictions == max(predictions, na.rm=T))
  
  # Find which values of macronutrient gave the peak
  macro_peak<-predictors[peak,]
  
  print(max(predictions, na.rm=T))
  print(macro_peak)
}

macrotrough<-function(GAM, data, XYZ, labels, predict_val,  exclude=NULL){
  
  require(sp)
  require(geometry)
  require(mgcv)
  require(metR)
  
  # List for the order of plots
  nutrient.order<-XYZ[c(1,2,3)]
  
  # List for the labels
  labels.order<-labels[c(1,2,3)]
  
  # Values to predict over, if they are unspecified
  x.limits<-c(floor(min(data[,nutrient.order[1]])), ceiling(max(data[,nutrient.order[1]])))
  y.limits<-c(floor(min(data[,nutrient.order[2]])), ceiling(max(data[,nutrient.order[2]])))
  z.limits<-c(floor(min(data[,nutrient.order[1]])), ceiling(max(data[,nutrient.order[3]])))
  
  # Fitted list to hold some results for later
  x.new<-seq(min(x.limits, na.rm=T), max(x.limits, na.rm=T), len=200)
  y.new<-seq(min(y.limits, na.rm=T), max(y.limits, na.rm=T), len=200)
  z.new<-seq(min(z.limits, na.rm=T), max(z.limits, na.rm=T), len=200)
  predictors<-as.data.frame(expand.grid(x.new, y.new, z.new))
  names(predictors)<-nutrient.order
  in.poly<-as.numeric(inhull(predictors[,c(1:3)], data[,names(predictors)]) != -1)
  
  # Add the predictors for the additional 'confounders'
  predictors<-cbind(predictors, predict_val)
  predictors<-predictors[-which(in.poly == 0),]
  
  # Unless fitted values are directly provided do the predictions
  predictions<-predict(GAM, newdata=predictors, type="response", exclude=exclude, newdata.guaranteed=T)
  
  # Find peak predicted outcome
  trough<-which(predictions == min(predictions, na.rm=T))
  
  # Find which values of macronutrient gave the peak
  macro_trough<-predictors[trough,]
  
  print(min(predictions, na.rm=T))
  print(macro_trough)
}