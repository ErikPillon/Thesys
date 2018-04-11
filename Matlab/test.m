% test fo half filled triangular lattice
clc; clear all; close all;

delta = [1,1]/2;
delta2 = [-1,-1]/2;
site = [1,0];
ratio = sqrt(3);
theta = pi/2;
alpha = 25
% alpha = 0.5;
disp('displacement in')
disp(site)
disp('of delta')
disp(delta2)
[gap1, gap2] = f_two_excitations_energy_gap(delta, delta2, site, ratio, theta, alpha)
