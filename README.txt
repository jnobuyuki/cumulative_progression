# cumulative_progression

Files:
clusterMass.R: 
  an implementation of the cluster mass permutation test (Maris & Oostenvelt 2007) for cumulative progression data

cProgAnalysis.R: 
  a sample script implementing clusterMass.R for an eye-tracking while reading study of reflexives

garnet-cProg.txt: 
  a sample data set to test the algorithm

Runtime:
Currently, computing 10 contrasts with 2500 bootstrap samples per contrast takes ~30 minutes, adding bootstrap samples seems to be more or less linear (5000 samples takes ~1 hour). Future versions will aim to parallelize bootstrap sampling.
