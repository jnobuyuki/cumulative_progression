#################################################################
##### Sample cluster mass analysis of cumulative progression data
##### Author: Shayne Sloggett 
##### March 2016
#################################################################
##### Analyzes cumulative progression data using the cluster mass 
##### permutation test (Maris & Oostenveld 2007)
##### for within subjects/items designs
#################################################################
#### Sample Item
################
# 1 The actress said that the ballerina horribly misepresented herself at the party.
# 2 The actor said that the ballerina horribly misepresented herself at the party.
# 3 The actress said that the lawyer horribly misepresented herself at the party.
# 4 The actor said that the lawyer horribly misepresented herself at the party.
# 5 The actress heard that the ballerina horribly misepresented herself at the party.
# 6 The actor heard that the ballerina horribly misepresented herself at the party.
# 7 The actress heard that the lawyer horribly misepresented herself at the party.
# 8 The actor heard that the lawyer horribly misepresented herself at the party.
#########
rm(list=ls())

##################################
#### Set directory and load script
##################################
setwd('~/documents/cumulative_progression')
source('clusterMass.r')

##################
#### Factor Coding
##################
# Define factors for the conditions (here, a 2x2x2 design)
condition.key = data.frame(cond = c(1:8),
                           matrix = rep(c('+match', '-match'), times = 4),
                           gram = rep(c('+gram', '-gram'), each = 2, times = 2),
                           verb = rep(c('speech', 'perception'), each = 4))

#############################################
#### Read in the data
#### Add factor labels as individuals columns
#############################################
prog.data = read.csv('cProg.sampleData.txt', sep = '\t')
prog.data = merge(prog.data,condition.key,by='cond')

###############################################################################
############################### Contrast Coding ###############################
###############################################################################
#### Contrasts are lists of pair-wise comparisons
#### The names of the contrasts can be arbitrarily specified
#### Pairwise comparisons always subtract the second item from the first item
#### Pairwise comparisons must refer to:
######## (1) values within a single column in the data, OR:
######## (2) contrasts defined EARLIER in the contrast list 
###############################################################################
##### Main Effects:
###################
matrix.effect = list('Matrix Match' = c('+match', '-match'))
gram.effect = list('Grammaticality' = c('+gram', '-gram'))
verb.effect = list('Verb Type' = c('speech', 'perception'))
###############################################################################
#### Interactions:
##################
# Contrasts testing for a three-way interactin of matrix, verb, and grammaticality
# Resolving the effect of matrix match on +/-gram sentences, and speech/perception verbs
matrix.interaction = list('Speech, +Gram' = c(1,2),
                      'Speech, -Gram' = c(3,4),
                      'Perception, +Gram' = c(5,6),
                      'Perception, -Gram' = c(7,8),
                      'Speech' = c('Speech, +Gram', 'Speech, -Gram'),
                      'Perception' = c('Perception, +Gram', 'Perception, -Gram'),
                      '+Gram' = c('Speech, +Gram', 'Perception, +Gram'),
                      '-Gram' = c('Speech, -Gram', 'Perception, -Gram'),
                      'Verb' = c('Speech', 'Perception'),
                      'Gram' = c('+Gram', '-Gram'))

# Contrasts testing for a three-way interactin of matrix, verb, and grammaticality
# Resolving the effect of grammaticality match on +/-match sentences, and speech/perception verbs
gram.interaction = list('Speech, +Match' = c(1,3),
                        'Speech, -Match' = c(2,4),
                        'Speech' = c('Speech, +Match', 'Speech, -Match'),
                        'Perception, +Match' = c(5,7),
                        'Perception, -Match' = c(6,8),
                        'Perception' = c('Perception, +Match', 'Perception, -Match'),
                        'Verb' = c('Speech', 'Perception'))

################################################################################
################################ Analysis Stream ###############################
################################################################################
#### analyses may be performed by calling inividual functions from clusterMass.R
#### or by calling the progressionAnalysis function (a wrapper function)
################################################################################

############################################################################################
#################################### progressionAnalysis ###################################
############################################################################################
# progressionAnalysis(raw.data, contrasts, condCol, group, makePlots, cutoff, minCluster, B)
############################################################################################
#### Arguments:
###############
# raw.data: the raw cumulative progression data for analysis
# contrasts: the list of contrasts to be tested (contrasts must be named)
# condCol: the column in the data containing the parwise comparison IDs (defaults to "cond")
# group: the column containing the grouping factor (defaults to "subj")
# makePlots: a binary value indicating whether plots should be generated for each contrast
# cutoff: the alpha level set to define clusters (defaults to 0.9)
# minCluster: the minimum cluster size to consider (defaults to 50ms)
# B: the number of bootstrap samples to generate for each contrast (defaults to 1000)
############################################################################################
#### Output:
############
# a list containing:
#### Summary: a summary of all of the observed clusters, by contrast, and associated p-values
#### bootstrap.data: a data frame of the bootsrap samples for each contrast
#### plots: a list of the plots created for each contrast
############################################################################################

########################
#### By-Subject Analyses
########################
##################
#### Main Effects:
##################
subj.gram = progressionAnalysis(prog.data, gram.effect, cutoff=.8, condCol='gram')
subj.matrix = progressionAnalysis(prog.data, matrix.effect, cutoff=.8, condCol='matrix')
subj.verb = progressionAnalysis(prog.data, verb.effect, cutoff=.8, condCol='verb')

##################
#### Interactions:
##################
subj.matrix.interaction = progressionAnalysis(prog.data, matrix.interaction, cutoff=.8)
subj.gram.interaction = progressionAnalysis(prog.data, gram.interaction, cutoff=.8)

#####################
#### By-Item Analyses
#####################
##################
#### Main Effects:
##################
item.matrix = progressionAnalysis(prog.data, matrix.effect, group='item', cutoff=.8, condCol='matrix')
item.gram = progressionAnalysis(prog.data, gram.effect, group='item', cutoff=.8, condCol='gram')
item.verb = progressionAnalysis(prog.data, verb.effect, group='item', cutoff=.8, condCol='verb')

##################
#### Interactions:
##################
subj.matrix.interaction = progressionAnalysis(prog.data, matrix.interaction, group='item', cutoff=.8)
subj.gram.interaction = progressionAnalysis(prog.data, gram.interaction, group='item', cutoff=.8)
