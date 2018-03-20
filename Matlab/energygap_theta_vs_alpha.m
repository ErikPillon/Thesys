% the goal is compare the energy changes wrt the variation in alpha and in
% theta

clc; clear all; close all;

delta = [1/2, -1/2];
site = [1,0];
delta2 = [-1/2, 1/2];
ratio = 1;
% theta = pi/2;
gap1 = 0;
gap2 = 0;

theta_range = linspace(pi/4,3*pi/4,10);
alpha_range = (0.5:0.25:3);


j = 1;
for theta = theta_range
    i = 1;
    for alpha = alpha_range
        [mov1, mov2] = f_two_excitations_energy_gap(delta, delta2, site, ratio, theta, alpha); 
        gap1(j,i) = mov1;
        gap2(j,i) = mov2;
        i = i+1
    end
    j = j+1;
end

%% plot the result
figure
[x,y] = meshgrid(alpha_range, theta_range);
pcolor(x, y, gap1-gap2)
title(['Behaviour of energy displacements wrt \alphaa and \theta.\newline a/b = ' num2str(ratio)])
xlabel('\alpha')
ylabel('\theta')
legend('energy difference')
filename = 'theta_vs_alpha';
saveas(gcf, filename)

