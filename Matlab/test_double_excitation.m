% test for double excitation
clc; close all; clear all; 

theta = pi/2;
site = [1,0];
ratio = 1;
alpha = 1;
delta = [1/2,-1/2];
delta2 = [-1/2,1/2];


% the two energy gap can be evaluated with the function for only one displacement.
% We need only to take into account the correction to this configuration
gap1 = f_one_excitation_energy_gap(delta, ratio, theta, alpha)
gap2 = f_one_excitation_energy_gap(delta2, ratio, theta, alpha)

%% evaluation of the energy gap of the two coupled displacements
% the idea is to subtract two "fake" interactions  
% + add the two real "interactions"

a = [1, 0];
rot = [cos(theta), -sin(theta); sin(theta), cos(theta)];
b = (ratio*rot*a')';

dist = site(1)*a + site(2)*b;

% the displacements must be converted in the lattice coordinate
delta = delta(1)*a+delta(2)*b;
delta2 = delta2(1)*a+delta2(2)*b;

% evaluate directly the energy difference
diff1_f = 1/norm(dist-delta)^alpha;
diff1_r = 1/norm(dist-delta+delta2)^alpha;
diff2_f = 1/norm(-delta2-dist)^alpha;
diff2_r = 1/norm(-delta2-dist+delta)^alpha;

energy_gap = gap1+gap2-diff1_f+diff1_r-diff2_f+diff2_r