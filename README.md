# Spatial Heterogenous Compositional Regression via Geographically Weighted Penalized Approach
This repository contains code for the manuscript *Linking COPD Prevalence with Income Distribution: A Spatial Heterogenous Compositional Regression via Geographically Weighted Penalized Approach.*

## Data
- **`USgdist.rds`** contains the spatial distance
- **`USgdist.RData`** contains the exponentially decayed spatial weight matrix with adjustment of states not geographically connected to USA
- **`state_fips_name.rds`** contains a map of state name to FIPS code
- **`cleans22_19.RData`** contains example real data

## HelpFunction
These are functions that are used in both simulation and real data analysis

## Process
- **`USgdist.R`** is used to create **`USgdist.RData`**, decay rate can be change

## ReadDataAnalysis
- **`RealDataAnalysis_example.R`** contains example code for real data analysis
- **`estimate_full_realus.R`** is a help function

## Simulation
- **`Simulation_example.R`** is an example code for simulation
- **`SimulationData.R`** is function used for simulate data
