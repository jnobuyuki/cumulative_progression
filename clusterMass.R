##############
#### Libraries
##############
library(ggplot2)
library(reshape)
library(plyr)
library(grid)

#########################################################################################
################################## progressionAnalysis: #################################
#########################################################################################
# Wrapper function which computes a cluster-mass permutation test for as set of contrasts
#########################################################################################
#### Arguments:
###############
# prog.data: the current data set to analyze
# contrasts: a of named contrasts in which the elements are taken from the 'label' column in prog.data
# group: a specified grouping facator (subject or item), defaults to subjects
# makePlots: should the script generate plots?
# cutoff: the lower-bound on the cluster mass alpha-level (defaults to 0.9), lower values promot smaller, more-sustained effects
# minCluster: what's the smallest cluster that should be considered by the algorithm (defaults to 50ms)?
# B: the number of bootstrap samples
#########################################################################################
#### Output:
############
# Prints a table with observed clusters for each contrast, the test statistic, and associated p-value
# Returns a list containing:
#### summary: a data frame with the observed clusters by contrast, associated test statistics, and p-values
#### bootstrap.data: a data frame with the boot-strap samples for the defined contrasts 
#### plots: a list containing the plots for each contrast
#########################################################################################
progressionAnalysis = function(raw.data, contrasts, condCol='cond', group='subj', makePlots=T, cutoff=.9, minCluster=50, B=1000){
  # Get time/condition means by group
  group.means = getGroupMeans(raw.data, contrasts, group, condCol)
  # Calculate the t-values for the relevant contrasts
  group.Ts = computeT(group.means, contrasts)
  df = length(unique(group.means$group))-1
  
  # Bootstrap data lists
  bootstrap.data = data.frame(row.names = c(1:B))
  plots.summary = c()
  
  # Table variables
  contrast.out = c()
  cluster.start = c()
  cluster.end = c()
  test.statistics = c()
  p.value = c()
  
  # Reshape the data to be in the appropriate form
  reshaped.data = melt(as.data.frame(group.means), id.vars = c('time','group'))
  
  # Collapse across group means to get means by-time and comparison for plotting
  plot.data = ddply(reshaped.data, .(time, variable), summarize, charProg = mean(value, na.rm = 1), se = sd(value, na.rm = 1))
  #plot.data$se = plot.data$se/sqrt(length(unique(reshaped.data$group)))
  plot.data$se = plot.data$se/sqrt(df+1)
  
  # Loop over each of the contrasts, calculating a permutation test for each
  for (comparison in names(contrasts)) {
    # Find time-points at which the alpha-level is satisfied
    critical.data = subset(group.Ts, abs(group.Ts[[comparison]]) >= qt(cutoff, df), 
                           select = c('time', comparison))
    
    # If any time points satisfy the alpha-cutoff
    if(nrow(critical.data) != 0){
      # Find clusters and get the derived test statistic (cluster mass) for each
      clusterMass = calculateMass(critical.data, comparison)

      # Generate a plot for the current comparison
      plot = makeProgressionPlot(plot.data, conditions=comparison, start=clusterMass$start, end=clusterMass$end)
      plots.summary[[comparison]] = plot
      
      # Run the permutation test for the current comparison
      sample = calculatePermutationTest(clusterMass, reshaped.data, contrasts[[comparison]], B)
      bootstrap.data[[comparison]] = sample
    }
    
    # If no time-points satisfied hte alpha-criterion, simply create the plot for the comparison
    else{
      plot = makeProgressionPlot(plot.data, conditions=comparison)
      plots.summary[[comparison]] = plot
    }
    
    # Print the current comparison's plot
    if(makePlots){print(plot)}
    
    # compare the derived test-statistic for each cluster in the contrast
    # against the boot-strap sampled distribution of derived test-statistics
    # for the largest cluster in the contrast (i.e. get the p-values for the clusters)
    for(c in unique(clusterMass$data$cluster)){
      cluster.data = subset(clusterMass$data, cluster == c)
      
      if(max(cluster.data$time)-min(cluster.data$time) >= minCluster){
        ts = cluster.data$test.statistic[1]
        tail = ifelse(ts > 0, sum(sample >= ts), sum(sample <= ts))
        p = tail/length(sample)
        
        contrast.out = c(contrast.out, comparison)
        cluster.start = c(cluster.start, min(cluster.data$time))
        cluster.end = c(cluster.end, max(cluster.data$time))
        test.statistics = c(test.statistics, ts)
        p.value = c(p.value, p)
      }
    }
  }
  
  # Print a summary of all the comparisons
  squares = c(1:10)^2
  numPlots = length(plots.summary)
  numCols = 1
  while(squares[numCols] < numPlots){numCols=numCols+1}
  if(makePlots & numPlots > 1){multiplot(plotlist = plots.summary, cols = numCols)}
  
  # Return a list-of-lists (see description of progressionAnalysis function)
  out.table = data.frame('Contrast'=contrast.out, 
                         'Start'=cluster.start, 
                         'End'=cluster.end, 
                         'Test Statistic'=test.statistics, 
                         'p-value'=p.value)
  print(out.table)
  return(list(summary = out.table, bootstrap.data = bootstrap.data, plots=plots.summary))
}

######################################################################
########################### getGroupMeans ############################
######################################################################
# Function for generating group means for each contrast and time point
######################################################################
#### Arguments:
###############
# prog.data: the current data to analyze
# cur.contrasts: a set of contrasts
# group: the grouping factor (defaults to subjects)
######################################################################
#### Output:
############
######################################################################
# A data frame of means by grouping fact, label-value, and time-point
######################################################################
getGroupMeans = function(raw.data, contrasts, group = 'subj', condCol = 'cond') {
  # Change the name of the grouping-factor to "group" for future reference
  names(raw.data)[which(names(raw.data)==group)] = 'group'
  names(raw.data)[which(names(raw.data)==condCol)] = 'label'
  
  # Get mean charProg by group and condition at each time point
  mean.data = ddply(raw.data, .(group, label, t), summarize, charProg = mean(value, na.rm = 1))
  colnames(mean.data) = c('group', 'label', 'time', 'charProg')
  #mean.data$label = as.character(mean.data$label)
  mean.data = cast(mean.data, time + group ~ label, value = 'charProg')
  
  # Loop over contrasts to generate by-group difference scores
  for (comparison in names(contrasts)) {
    cur.comparison = as.character(contrasts[[comparison]])
    mean.data[[comparison]] = mean.data[[cur.comparison[1]]]-mean.data[[cur.comparison[2]]]
  }
  
  return(mean.data)
}

#################################################################################
################################### computeT ####################################
#################################################################################
# For each contrast, calculates the t-value of the comparison for each time-point
#################################################################################
#### Arguments:
###############
# prog.data: a data set of group-means
# contrasts: a set of contrasts over which to compute t-values
#################################################################################
#### Output:
############
# A table of t-values for each contrast at each time-point
#################################################################################
computeT = function(group.data, contrasts) {
  # Create a data frame with the time-window info from the cProg data being tested
  test.data = data.frame('time' = unique(group.data$time))
  
  # Loop over contrasts, testing each comparison individually
  for(comparison in names(contrasts)) {
    cur.comparison = as.character(contrasts[[comparison]])
    test.data[[comparison]] = 'NA'
    # Loop over times within comparison, finding t-statistic
    for(s in test.data$time) {
      cur.data = subset(group.data, time == s)
      cur.t = t.test(cur.data[[cur.comparison[1]]],
                     cur.data[[cur.comparison[2]]],
                     paired=TRUE)
      test.data[[comparison]][test.data$time == s] = cur.t$statistic
    }
    test.data[[comparison]] = as.numeric(as.character(test.data[[comparison]]))
  }
  return(test.data)
}

#####################################################################
########################### calculateMass ###########################
#####################################################################
# Function which iterates over the critical t-values for a comparison
# sums the t-values in a cluster to create the derived test statistic
#####################################################################
#### Arguments:
###############
# critical.data: the critical t-values for a contrast
# current.comparison: the current contrast of interest
#####################################################################
#### Output:
############
# The set of clusters for a contrast, with derived test statistics
#####################################################################
calculateMass = function(critical.data, current.comparison) {
  critical.data = identifyClusters(critical.data)
  for (c in unique(critical.data$cluster)) {
    test.statistic = sum(subset(critical.data, cluster == c)[[current.comparison]])
	  critical.data$test.statistic[critical.data$cluster == c] = test.statistic
	  clustersize = nrow(subset(critical.data, cluster == c))
	  clustersize = test.statistic
	  if (c == 0) {
	    max = clustersize
	    max.cluster = c
  	} 
	  else {
	    if (abs(clustersize) > abs(max)) {
	      max = clustersize
	      max.cluster = c
    	}
  	}
  }
  
  cur.start <- min(subset(critical.data, cluster == max.cluster)$time)
  cur.end <- max(subset(critical.data, cluster == max.cluster)$time)
  return(list(data = critical.data, max.cluster = max.cluster,start = cur.start,end = cur.end))
}

########################################################################
########################### identifyClusters ###########################
########################################################################
# A function which iterates over t-values >= the set alpha-cutoff
# Contiguous t-values are grouped into numbered clusters (starting at 0)
########################################################################
#### Arguments:
###############
# critical.data: the list of t-values which satisfy the alpha-cutoff
########################################################################
#### Output:
############
# the critical.data table with a new column identifying clusters
########################################################################
identifyClusters = function(critical.data) {
  for(row in c(1:nrow(critical.data))){
    if(row == 1){
      cluster = 0
      critical.data$cluster[row] = cluster
    } 
    else {
      if(critical.data$time[row] == critical.data$time[row-1]+10){
        critical.data$cluster[row] = cluster
      } 
      else {
        cluster = cluster + 1
        critical.data$cluster[row] = cluster
      }
    }
  }
  return(critical.data)
}

###########################################################################################
################################# calculatePermutationTest ################################
###########################################################################################
# The core of the cluster-mass permutation algorithm
# For the largest observed cluster in a contrast
#### Permutes the condition labels on the data
#### Generates a new derived test statistic for the contrast
#### Applies B-many times to create a distribution of derived test statisc for the contrast
###########################################################################################
# Arguments:
############
# critical.data: the identity, and time-span of the largest cluster for the contrast
# bootstrap.data: the group-level means for bootstrap sampling
# cur.contrasts: the current contrast being sampled
# B: the number of bootstrap samples to generate
###########################################################################################
##### Output:
#############
# a list of the sampled derived test-statistics for the current contrast
###########################################################################################
calculatePermutationTest = function(critical.data,bootstrap.data,cur.contrasts,B) {
	start = critical.data$start
	end = critical.data$end
	max.cluster = critical.data$max
	
	bootstrap.data = subset(bootstrap.data, time >= start & time <= end)
	bootstrap.data = cast(bootstrap.data,group+variable~time) 
	sampled.statistics = rep('NA',B)
	for (b in 1:B){
		resample.data = bootstrap.data
		for (cur.group in unique(resample.data$group)){                ### loop over subjects
    		cur.row.nums = which(resample.data$group == cur.group & resample.data$variable %in% cur.contrasts)    ### pull out current subject's data
    		if(rbinom(1,1,.5)==0){                                   ### flip a fair coin, and if it's tails .... 
      			resample.data[cur.row.nums,c(3:ncol(resample.data))] = resample.data[c(cur.row.nums[2],cur.row.nums[1]),c(3:ncol(resample.data))]    ### swap the it.int and i.int data for that subject
    		}  ## otherwise leave the data as-is
  		}

		sample.test = data.frame('time' = colnames(resample.data)[c(-1,-2)])   
  		for(t in 3:ncol(resample.data)){                                               ### t is in columns now, so loop over cols
    		data = resample.data[,c(1,2,t)]                                               ### get current data
   			condition.one.vals = subset(data,variable==cur.contrasts[1])[,3]                        ### extract it.int.vals... 
   			condition.two.vals = subset(data,variable==cur.contrasts[2])[,3]                        ### extract it.int.vals... 
    		sample.test$t.value[sample.test$time == colnames(resample.data)[t]] = t.test(
    			condition.one.vals, condition.two.vals, 
    			paired = T)$statistic        							### compute t statistic
  		}
 	
  		sampled.statistics[b] = sum(sample.test$t.value)
  	}
  	return(as.numeric(as.character(sampled.statistics)))
}

#######################################################
################# makeProgressionPlot #################
#######################################################
# A function for generating a plot for a given contrast
#######################################################
#### Arguments:
###############
# data: the data for plotting
# conditions: the contrast currently being plotted
# edge: how much of the time-window should be represented?
# current.col: what color should the graph be?
# start/end: the boundaries on the largest cluster
#######################################################
#### Output:
############
# A ggplot object
#######################################################
makeProgressionPlot = function(data,conditions = test.contrasts,edge = 4000,current.col='#85274e',start=0,end=0) {
  data = subset(data, variable %in% conditions)
  ymax = max(data$charProg+data$se)
  ymin = min(data$charProg-data$se)
  
  ggplot(subset(data, time<=edge),aes(x=time,y=charProg, colour=variable))+
    geom_line(size = 1.5, show_guide = F) +
    ylim(ymin,ymax) + 
    geom_ribbon(aes(ymax = charProg+se, ymin = charProg-se, fill = variable), alpha = .5, show_guide = F) +
    scale_colour_manual(values = current.col) + 
    scale_fill_manual(values = current.col) + 
    labs(colour = '', fill = '', x = 'Time(ms)', y = 'Characters', title = conditions) + 
    geom_abline(intercept = 0, slope = 0, colour = 'black') +
    theme(text = element_text(size = 10, colour = 'black', face = 'bold'),
          axis.text = element_text(size = 10, colour = 'black', face = 'italic'),
          axis.title = element_text(size = 10, colour = 'black', face = 'bold'),
          axis.line = element_line(colour ='black'),
          legend.title.align = .5, 
          legend.text = element_text(size = 10, colour = 'black', face = 'italic'),
          strip.text = element_text(size = 10, colour = 'black', face = 'bold'),
          strip.background = element_rect(fill = NA, colour = NA),
          legend.position = 'bottom',
          panel.background = element_rect(fill = NA),
          panel.grid = element_blank(),
          panel.margin.x = unit(0, 'cm'),
          panel.border = element_rect(fill = NA, colour = NA),
          legend.box = 'horizontal') +
    annotate("rect", xmin = start, xmax = end, ymin = ymin, ymax = ymax,alpha = .2)
}

#################################################################################
################################### multiplot ###################################
#################################################################################
# A function which prints multiple ggplot objects on the same page
# Source: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#################################################################################
multiplot = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots = c(list(...), plotlist)
  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout = matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } 
  else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
