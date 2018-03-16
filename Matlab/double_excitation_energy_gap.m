clc; clear all; close all;

% algorithm for the implementation of the ewald summation
% for the double excitation displacement following notes on pages 4 
% and following
% suggested reference "Ewald-Rastelli"

%% initialize all the variables
alpha = 1;
vol = 2;
c = 5; % convergence coefficient
target = 0.65; % value to reach

delta = sqrt(2)*[1/2, 1/2]; % displacement of the electron

G_cut = 20;
R_cut = G_cut;
%R_cut = 20;

% the final summation is given by 1/gam[sum1-sum2+sum3-sum4];
sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0;

%% evaluation of all the terms for the first excited term
% sum1
for Gx = -G_cut:G_cut
    for Gy = -G_cut:G_cut
        G = pi*sqrt(2)*[Gx, Gy];
        sum1 = sum1+((cos(G*(-delta)')-1)*f_phi(-alpha/2,G*G'/4/c));
    end
end
sum1 = sum1*pi*c^(alpha/2-1)/vol;

% sum2
sum2 = f_selfinteraction(norm(delta), c, alpha)-f_selfinteraction(0, c, alpha);

% sum3
for Rx = -R_cut:R_cut
    for Ry = -R_cut:R_cut
        R = sqrt(2)*[Rx, Ry];
        sum3 = sum3+(f_phi(alpha/2-1,norm(R-delta)^2*c)-f_phi(alpha/2-1,norm(R)^2*c));
    end
end
sum3 = sum3*c^(alpha/2);

% sum4
sum4 = f_phi(alpha/2-1,norm(delta)^2*c)-f_phi(alpha/2-1,0);
sum4 = c^(alpha/2)*sum4;

%% evaluation of the energy necessary for the double excitation
% the idea is to double the gap energy + subtract two "fake" interactions  
% + add the two real "interactions"

%intialize the variables
diff1_f = 0; diff1_r = 0; diff2_f = 0; diff2_r = 0;

%evaluate directly the energy difference
diff1_f = 1/norm([sqrt(2)/2,sqrt(2)/2])^alpha;
diff1_r = 1/norm([sqrt(2),0])^alpha;

diff2_f = 1/norm([sqrt(2)/2,sqrt(2)/2])^alpha;
diff2_r = 1/norm([sqrt(2),0])^alpha;

gap = (sum1-sum2+sum3-sum4)/gamma(alpha/2);
gap2 = 2*gap - diff1_f + diff1_r - diff2_f + diff2_r;

disp('the final energy gap is:')
disp(gap2)
disp('and it is supposed to be as close as possible to:')
disp(target)
