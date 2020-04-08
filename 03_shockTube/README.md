# Shock Tube Solver

## Convective flux schemes

### Flux Differencing Splitting

- Roe

### Flux Vectoer Splitting

- van Leer
- AUSM

## MUSCL scheme

This code has an option to use MUSCL scheme.

## Accuracy in time

2nd order accuracy in time is default.
The code of 1st order accuracy in time is commented out.

## Results

Results seem to be valid comparing with the following github results.

- [Computational\-Fluid\-Dynamics/31\.Sod Shock Tube at master · xuaoxiqi/Computational\-Fluid\-Dynamics](https://github.com/xuaoxiqi/Computational-Fluid-Dynamics/tree/master/31.Sod%20Shock%20Tube)

van Leer w/ MUSCL has some oscillations in the region between expansion wave (Region 2) and Region 3.

## Bugs

- Exact solution might be wrong.
