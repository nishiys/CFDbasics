# Abount CFDbasics

# Code List

## 1D Wave Equation
### Scheme
- Lax scheme

### How to use
The following command gives you the solution.

```shell
# move excecutable to the /run directory

cd run
./waveSolver1d
python wavePlot.py
```

## 1D Shock Tube
### Schemes for convection term
Flux Difference Splitting (FDS)
- Roe

Flux Vector Splitting (FVS)
- Steger Warming
- van Leer
