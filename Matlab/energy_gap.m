clc; clear all; close all;

% algorithm for the implementation of the ewald summation
% following notes on pages 3 and following
% 
% suggested reference "Ewald-Rastelli"

%% initialize all the variables
alpha = 1;
vol = 2;
c = 10; % convergence coefficient
target = 0.62; % value to reach

delta = sqrt(2)*[1/2, 1/2]; % displacement of the electron

G_cut = 10;
R_cut = G_cut;
%R_cut = 20;

% the final summation is given by 1/gam[sum1-sum2+sum3-sum4];
sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0;

%% evaluation of all the terms
% sum1
for Gx = -G_cut:G_cut
    for Gy = -G_cut:G_cut
        G = pi*sqrt(2)*[Gx, Gy];
        sum1 = sum1+((cos(G*delta')-1)*f_phi(alpha/2,G*G'/4/c));
    end
end
sum1 = sum1*pi*c^(alpha/2-1)/vol;

% sum2
sum2 = f_selfinteraction(norm(delta), c, alpha)-f_selfinteraction(0, c, alpha);

% sum3
for Rx = -R_cut:R_cut
    for Ry = -R_cut:R_cut
        R = sqrt(2)*[Rx, Ry];
        sum3 = sum3+(f_phi(alpha/2-1,norm(R+delta)^2*c)-f_phi(alpha/2-1,norm(R)^2*c));
    end
end
sum3 = sum3*c^(alpha/2);

% sum4
sum4 = f_phi(alpha/2-1,norm(delta)^2*c)-f_phi(alpha/2-1,0);
sum4 = c^(alpha/2)*sum4;

gap = (sum1-sum2+sum3-sum4)/gamma(alpha/2);
disp('the final energy gap is:')
disp(gap)
disp('and it is supposed to be as close as possible to:')
disp(target)
