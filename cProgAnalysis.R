######################
######################
###### Permutation analysis of cumulative progression data
###### Author: Shayne Sloggett 
###### Edited by: Brian Dillon
###### October 2015
######################
###### Implements cluster mass test of Maris & Oostenveld (2007)
###### for within subjects or within items designs
######################
######################

#########
### Condition key
### a The actress said that the ballerina horribly misepresented herself at the party.
### b The actor said that the ballerina horribly misepresented herself at the party.
### c The actress said that the lawyer horribly misepresented herself at the party.
### d The actor said that the lawyer horribly misepresented herself at the party.
#########
#########
### e The actress heard that the ballerina horribly misepresented herself at the party.
### f The actor heard that the ballerina horribly misepresented herself at the party.
### g The actress heard that the lawyer horribly misepresented herself at the party.
### f The actor heard that the lawyer horribly misepresented herself at the party.
#########
rm(list=ls())

####### Load libraries and run code
setwd('~/desktop/Project Gemstone/Garnet/eyetracking/eyedatator')

# library(cumulativeProgressionsFxns)  one day :)
source('~/desktop/project gemstone/clusterMass.r')

condition.key = data.frame(cond = c(1:8),
                           matrix = rep(c('+match', '-match'), times = 4),
                           gram = rep(c('+gram', '-gram'), each = 2, times = 2),
                           verb = rep(c('speech', 'perception'), each = 4))

# Contrasts for main effects
matrix.effect = list('Matrix Match' = c('+match', '-match'))
gram.effect = list('Grammaticality' = c('+gram', '-gram'))
verb.effect = list('Verb Type' = c('speech', 'perception'))

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

# Read in the data, merge experiment-specific condition labels in
prog.data = read.csv('garnet-cprog.txt', sep = '\t')
prog.data = merge(prog.data,condition.key,by='cond')

#############################
#### By-Subject Analyses ####
#############################
#################
#### Main Effects
#################
subj.matrix = progressionAnalysis(prog.data, matrix.effect, group='subj', cutoff=.8, condCol='matrix')
subj.gram = progressionAnalysis(prog.data, gram.effect, group='subj', cutoff=.8, condCol='gram')
subj.verb = progressionAnalysis(prog.data, verb.effect, group='subj', cutoff=.8, condCol='verb')
#################
#### Interactions
#################

subj.matrix.interaction.200 = progressionAnalysis(prog.data, matrix.interaction, cutoff=.8, condCol='cond', B=200)
subj.matrix.interaction.500 = progressionAnalysis(prog.data, matrix.interaction, cutoff=.8, condCol='cond', B=500)
subj.matrix.interaction.1000 = progressionAnalysis(prog.data, matrix.interaction, cutoff=.8, condCol='cond', B=1000)
subj.matrix.interaction.2500 = progressionAnalysis(prog.data, matrix.interaction, cutoff=.8, condCol='cond', B=2500)
subj.matrix.interaction.5000 = progressionAnalysis(prog.data, matrix.interaction, cutoff=.8, condCol='cond', B=5000)

subj.gram.interaction = progressionAnalysis(prog.data, gram.interaction, cutoff=.8)
