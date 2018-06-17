clear all; clc; close all;
% Set grid size
xMax = 2;
yMax = 2;
h = 0.1;

% Set model constants (dt, h, kappa, mu, rho_0, mass, etc.)
N = 100;
dt = 0.005;
kappa = 40;
rho_0 = 1000;
mu = 0.3;
beta = 0.8;

% Running 3 Trials: dam break, object, and different rest densities.
DamBreak(xMax, yMax, h, N, kappa, mu, rho_0, beta, dt);
figure
Obstacle(xMax, yMax, h, N, kappa, mu, rho_0, beta, dt);
figure
MixedDensities(xMax, yMax, h, N, kappa, mu, rho_0, beta, dt);