clc; clear all; close all;

% TRIANGULAR LATTICE - SOLVED
%
% algorithm for the implementation of the ewald summation
% following notes on pages 3 and following
% 
% suggested reference "Ewald-Rastelli"
%
% this algorithm aims to plot the energy needed for one displacement letting
% \alpha varying in a (equilateral) triangular lattice half filled

%% initialize all the variables
theta = pi/2;
a = [2*sin(pi/3), 0];
b = [0,1];

vol = norm(a)*norm(b)*sin(theta);
c = 5; % convergence coefficient
% in principle the summation should not depend on the above coefficient

% relative site in which we deplace the second particle by delta2
site = [1,0];  

delta = [1/2, 1/2].*[norm(a), norm(b)]; % displacement of the electron
delta2 = [-1/2, -1/2].*[norm(a), norm(b)];

G_cut = 10;
R_cut = 10;

% evaluation of the reciprocal space vectors
K = [a; b]\(2*pi*eye(2));
kx = K(:,1);
ky = K(:,2);

% the final summation is given by 1/gam[sum1-sum2+sum3-sum4];
sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0;

alpha_range = [(0.1:0.05:1),(1.5:0.25:3),(4:1:10),(15:5:50)];
%alpha_range = [(0.1:0.05:1)];

%% evaluation of all the terms for each alpha
index = 1;
% ewald routine
for alpha = alpha_range
    sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0;
    sum12 = 0; sum22 = 0; sum32 = 0; sum42 = 0;
    % sum1
    for Gx = -G_cut:G_cut+1
        for Gy = -G_cut:G_cut+1
            G = [Gx, Gy].*[norm(kx),norm(ky)];
            sum1 = sum1+((cos(G*(delta)')-1)*f_phi(-alpha/2,G*G'/4/c));
            sum12 = sum12+((cos(G*(delta2)')-1)*f_phi(-alpha/2,G*G'/4/c));
        end
    end
    sum1 = sum1*pi*c^(alpha/2-1)/vol;
    sum12 = sum12*pi*c^(alpha/2-1)/vol;
    % sum2
    sum2 = f_selfinteraction(norm(delta), c, alpha)-f_selfinteraction(0, c, alpha);
    sum22 = f_selfinteraction(norm(delta2), c, alpha)-f_selfinteraction(0, c, alpha);
    % sum3
    for Rx = -R_cut:R_cut+1
        for Ry = -R_cut:R_cut+1
            R = [Rx, Ry].*[norm(a), norm(b)];
            sum3 = sum3+(f_phi(alpha/2-1,norm(R+delta)^2*c)-f_phi(alpha/2-1,norm(R)^2*c));
            sum32 = sum32+(f_phi(alpha/2-1,norm(R+delta2)^2*c)-f_phi(alpha/2-1,norm(R)^2*c));
        end
    end
    sum3 = sum3*c^(alpha/2);
    sum32 = sum32*c^(alpha/2);
    % sum4
    sum4 = f_phi(alpha/2-1,norm(delta)^2*c)-f_phi(alpha/2-1,0);
    sum4 = c^(alpha/2)*sum4;
    sum42 = f_phi(alpha/2-1,norm(delta2)^2*c)-f_phi(alpha/2-1,0);
    sum42 = c^(alpha/2)*sum42;
    % final result stored in the vector gap
    gap1(index) = (sum1-sum2+sum3-sum4)/gamma(alpha/2);
    gap2(index) = (sum12-sum22+sum32-sum42)/gamma(alpha/2);
    %% modification of energy 2
    dist = site.*[norm(a), norm(b)];
    % evaluate directly the energy difference
    diff1_f = 1/norm(dist-delta)^alpha;
    diff1_r = 1/norm(dist-delta+delta2)^alpha;
    diff2_f = 1/norm(-delta2-dist)^alpha;
    diff2_r = 1/norm(-delta2-dist+delta)^alpha;
    % update the real energy
    energy2(index) = gap1(index)+gap2(index)-diff1_f+diff1_r-diff2_f+diff2_r;
    index = index+1; %update counter
    alpha
end
% as alpha grows to infinity, the energy gap grows to 1 as expected since 
% now we've got three electrons as closest neighbourhoods but before we 
% had already two, then we have 3-2 = 1

%% plot the results and save them
figure
plot(alpha_range, gap1, alpha_range, energy2)
grid on
title('Energy gap from ground state')
xlabel('\alpha')
ylabel('Energy gap $(e^2/vol)$', 'interpreter','latex')
legend('One displacement','Two displacements')
% saveas(gcf,'Energy_gap_tringular_lattice.eps','epsc')
% plot the second figure zoomed 
figure(2)
plot(alpha_range(1:15), gap1(1:15), alpha_range(1:15), energy2(1:15))
grid on
title('Energy gap from ground state')
xlabel('\alpha')
ylabel('Energy gap $(e^2/vol)$', 'interpreter','latex')
legend('One displacement','Two displacements')
% saveas(gcf,'Energy_gap_tringular_lattice_zoomed.eps','epsc')