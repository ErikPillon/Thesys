% this script is aimed for seeing how the gap energy changes wrt the
% variation in alpha
clc; clear all; close all;

delta = [1/2, -1/2];
site = [1,0];
delta2 = [-1/2, 1/2];
ratio = 1;
theta = pi/2;
gap1 = 0;
gap2 = 0;

i = 1;
alpha_range = (0.5:0.1:1);
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
filename = 'Energy comparison: 1 vs 2 displacements';
saveas(gcf, filename)

