# Trajectory Optimization 
## Versions / branches
 * main:   1D motion in 2D phase space

 * 3D:      3D motion in 6D phase space

 the two branches can not be easily merged unfortunately.

## How to use:
0) change to your desired branch
1) edit ctoConfig.txt to change parameters
2) python3 searchOpt.py  -1

command line argument for searchOpt.py selects which trajectory to plot on the grid.  -1 means plot all trajectories.

3) in 3D mode, the plotting is a separate program:
```
> python3 animate3D.py <<FILENAME.csv>>
```

