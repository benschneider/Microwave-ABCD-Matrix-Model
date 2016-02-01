# ABCD - SQUID - Model

## Script to model the reflection of a SQUID at the end of a Transmission line. (under construction)

### Circuit diagram

![alt tag](https://cloud.githubusercontent.com/assets/4573907/8725273/a0e0dba8-2bd3-11e5-86a0-8202fd2715d2.png)

#### settings are in plotter_experimental.py 

filename = 'target.hdf5' --- define datafile to be loaded for comparison

matplotlib.use('Qt4Agg') --- macosx, Qt4Agg, WX


#### Run with

ipython -i plotter_experimental.py



## Structure

plotter_experimental.py ---	Main Program which uses the other files

ABCD.py 	---		ABCD Matrixes

parsers.py 	---		Data load and save functions

interface.py	---		Loading buttons (contains an update function and a few more)

fitdata.py	---		Fitting algorithm


## Currently working on

fitdata.py
- [ ] Optimize fitting on complex data

interface.py
- [ ] Cleanup poss. use a .uic file (low priority atm)
