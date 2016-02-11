# ABCD - SQUID - Model

## Script to model the reflection of a SQUID at the end of a Transmission line. (under construction)

### Circuit diagram

![alt tag](https://cloud.githubusercontent.com/assets/4573907/12975622/7850b250-d0bd-11e5-8500-1065c85d3a93.png)
#### settings are in plotter_experimental.py 

filename = 'target.hdf5' --- define datafile to be loaded for comparison

matplotlib.use('Qt4Agg') --- macosx, Qt4Agg, WX


#### Run with

ipython -i plotter_experimental.py



## Structure

plotter_experimental.py ---	Main Program which uses the other files

ABCD.py 	---		ABCD Matrixes

parsers.py 	---		Data load and save functions

## Currently working on

interface.py
- [ ] Cleanup poss. use a .uic file (low priority atm)
