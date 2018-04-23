% reference for the work: papers 12 , 10 avril '18

clc; clear all; close all;

% first of all we consider the squared lattice without the possible
% movements along the diagonal

t = 1; %hopping parameter

% PUT ALGORITHM HERE!

%% ... and then we allow movements along the diagonal
% script for one and two excitations on a trinangular lattice, that is 
% theta = pi/2 and ratio = 2*sin(pi/3) 

% initialize the variables for the algorithm
delta = [1/2, 1/2]; % displacement
site = [1,0]; % site of the second displacement
delta2 = [-1/2, -1/2]; % second displacement 
ratio = 2*sin(pi/3); 
theta = pi/2;

[energy1,gap] = f_two_excitations_energy_gap(delta,delta2,site,ratio,theta)
energy2 = energy1+gap;

Energies = [energy1*ones(4,1);energy2*ones(2,1)];

k = linspace(-10*pi, 10*pi, 1000);
v = 1; % hopping parameter along the diagonal
delta_x = 0.5;
delta_y = sqrt(3)/2;

Bands = [-t*exp(1i*k*delta_x/2)-t*exp(-1i*k*delta_x/2);...
    -t*exp(1i*k*delta_y/2)-t*exp(-1i*k*delta_y/2);...
    -t*exp(-1i*k*delta_x/2)-t*exp(1i*k*delta_x/2);...
    -t*exp(-1i*k*delta_y/2)-t*exp(1i*k*delta_y/2);...
    -t*(exp(-1i*k*delta_x/2)+exp(-1i*k*delta_y/2)+exp(1i*k*delta_x/2)+exp(1i*k*delta_y/2));...
    -t*(exp(1i*k*delta_x/2)+exp(1i*k*delta_y/2)+exp(-1i*k*delta_x/2)+exp(-1i*k*delta_y/2))];

Bands = Bands+Energies;


figure
title('Energy Bands')
hold on
for count=1:6
    plot(k,Bands(count,:))
end
ylabel('Energy levels $E_i$','interpreter','latex')
xlabel('$k$', 'interpreter', 'latex')
legend('Energy 1','Energy 2','Energy 3','Energy 4','Energy 5','Energy 6');
