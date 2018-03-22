% script for one and two excitations on a trinangular lattice, that is 
% theta = pi/2 and ratio = 1/(2*sin(pi/3)) 

% the goal is to see how the energy changes wrt one and two displacements
% wrt the ground state energy on a triangular lattice
clc; clear all; close all;

%% initialize the variables for the algorithm
delta = [1/2, 1/2]; % displacement
site = [1,0]; % site of the second displacement
delta2 = [-1/2, -1/2]; % second displacement 
ratio = 1/(2*sin(pi/3)); 
theta = pi/2;

gap1 = 0;
gap2 = 0;
i = 1;
alpha_range = [0.5:0.5:5];
for alpha = alpha_range
    [mov1, mov2] = f_two_excitations_energy_gap(delta, delta2, site, ratio, theta, alpha); 
    gap1(i) = mov1;
    gap2(i) = mov2;
    i = i+1;
    alpha
end

figure
plot(alpha_range, gap1, 'b', alpha_range, gap2, 'r')
grid on
title(['Behaviour of energy displacements wrt \alpha.\newline a/b = ' num2str(ratio) '; \theta = ' num2str(theta)])
xlabel('\alpha')
ylabel('Energy gap')
legend('First displacement energy','Second displacement energy')
% filename = 'Energy comparison: 1 vs 2 displacements';
% saveas(gcf, filename)