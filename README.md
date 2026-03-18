# Adaptive MPC for nonlinear sytems

## Description

This repository contains the MATLAB code that accompanies the paper 
>  Johannes Köhler "Certainty-equivalent adaptive MPC for uncertain nonlinear systems", 2026.


## Prerequisites

- [MATLAB](https://de.mathworks.com/products/matlab.html)  
- [YALMIP](https://yalmip.github.io/)  
- [Mosek](https://www.mosek.com/)  
- [Casdadi](https://web.casadi.org/)  

## Usage
Run either `main.m` file to setup the model, design the MPC, run closed-loop simulations, generate plots and compare quantative metrics. 
Some settings can be toggled in the beginning of code, e.g., to make a short simulation and not simulate some methods that take longer.

![Adaptive MPC navigating a drone through obstacles while learning the model parameters online](Drone.png)

