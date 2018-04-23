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
alpha_range = [(0.5:0.1:1),(1.5:0.5:5),(6:1:15)];
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
titlename = ['Behaviour of energy displacements wrt $\alpha$. a/b = ' num2str(ratio) '; $\theta$ = ' num2str(theta)];
title(titlename, 'interpreter','latex')
xlabel('\alpha')
ylabel('Energy gap')

figure
plot(alpha_range(1:5), gap1(1:5), 'b', alpha_range(1:5), gap2(1:5), 'r')
grid on
titlename = ['Behaviour of energy displacements wrt $\alpha$.  a/b = ' num2str(ratio) '; $\theta$ = ' num2str(theta)];
title(titlename, 'interpreter','latex')
xlabel('\alpha')
ylabel('Energy gap')

cd Im
filename = 'Energy_comparison:1vs2displacements.eps';
saveas(gcf, filename,'epsc')
filename2 = 'Energy_comparison:1vs2displacements_zoomed.eps';
saveas(gcf, filename2,'epsc')
cd ..
