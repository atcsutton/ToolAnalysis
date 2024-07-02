# TimeGridVertex

A simple vertex finder that minimizes the RMS time to determine where a low energy particle interacted. 
Inspired by T2Ks Neut-Fit, the goal is to minimize $t_{RMS}$. Test vetices are spaced coarsely throughout the detector. The test vertex spacing is gradually decreased until we converge on a location with satisfacotry RMS. The spacing between verticies will decrease by half after each iteration.

ClusterMapName - The name of the map from cluster time to vector of hits
UseMCHits - Boolean to use Hits or MCHits
InitialSpacing - Distance between initial vertices in meters. The number of vetices will be determined on the fly.
MinRMS - Determines when the calculation is terminated

## Data


**TimeGridVertex** `std::map<double, Position>`
* Takes a vector of Hits or MCHits and determines the location within the tank of the generating particle. Creates a map between the local cluster time and the vertex Position


## Configuration


```
param1 value1
param2 value2
```
