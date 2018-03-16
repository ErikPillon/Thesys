% this script is aimed for seeing how the gap energy changes wrt the
% variation in alpha
clc; clear all; close all;

delta = [1/2, 1/2];
ratio = 1;
theta = pi/2;

i = 1;
alpha_range = (0.5:0.1:6);
for alpha = alpha_range
    energy_gap(i) = f_one_excitation_energy_gap(delta, ratio, theta, alpha);
    i = i+1
end
figure
plot(alpha_range, energy_gap)
grid on
title(['Energy gap due to one displacement \delta.\newline a/b = ' num2str(ratio) '; \theta = ' num2str(theta)])

