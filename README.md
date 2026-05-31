# Habitat 1.0 Heat Distribution

A data-driven heat distribution simulation developed using real temperature measurements collected during the **EMMPOL25 analog astronaut mission** at **Habitat 1.0** operated by the Analog Astronaut Training Centre (AATC).

The project combines:

* Python-based sensor data analysis and parameter estimation,
* C++ implementation of a 2D heat diffusion solver,
* Crank–Nicolson finite-difference simulation,
* heat-map visualization and model validation against real measurements.

## Project Overview

Temperature data recorded by Netatmo sensors were used to estimate:

* heating and cooling time constants,
* effective thermal diffusivity,
* wall heat transfer coefficients,
* heater source strength.

The estimated parameters were then used in a numerical solution of the heat equation:

[
\frac{\partial T}{\partial t}
=============================

D\nabla^2T + S(x,y)
]

to simulate heat propagation inside two rooms of Habitat 1.0:

* **Bedroom** (5.05 m × 5.05 m)
* **Geolab** (2.25 m × 6.35 m)

## Repository Structure

```text
Sensors/                     raw sensor data

heat_distribiution.ipynb     data processing and parameter estimation
heat_maps.ipynb              heat-map visualization

bedroom.cpp                  bedroom simulation
geolab.cpp                   geolab simulation

bedroom/                     simulated temperature fields
geolab/                      simulated temperature fields

paper.tex                    scientific report
references.bib               bibliography
```

## Validation Results

| Room    | RMSE [°C] | MAE [°C] |
| ------- | --------- | -------- |
| Bedroom | 1.158     | 1.112    |
| Geolab  | 8.177     | 7.783    |

The bedroom simulation reproduces the measured thermal behaviour with reasonable accuracy, while larger discrepancies observed in the geolab highlight the limitations of the diffusion-only approach and the importance of convective heat transport.

## Main Technologies

* Python

  * NumPy
  * Pandas
  * SciPy
  * Matplotlib

* C++

  * Finite Difference Method (FDM)
  * Crank–Nicolson scheme
  * Gauss–Seidel iteration

## Author

**Natalia Kowalczyk**
AGH University of Krakow
Faculty of Physics and Applied Computer Science

## License

This repository is intended for educational and research purposes.
