clc; clear all; close all;
% %% TRIANGULAR LATTICE
% delta = [1/2,1/2];
% ratio = sqrt(3);
% theta = pi/2;
% 
% E = 0;
% for iter = 0:10
%     x=iter*0.1;
%     Energy = f_lr_plus_nn(delta, ratio, theta, x) 
%     E(iter+1) = Energy;
% end
% figure
% plot(0:10,E)
% title('Energy gap from LR potential to NN (triangular lattice)')  
% xlabel('$x$','interpreter','latex')
% ylabel('Energy $e^2/vol$','interpreter','latex')
% 
% cd Im/
%     saveas(gcf,'Lr_vs_nn_triangular.eps','epsc')
% cd ..

%% SQUARED LATTICE

%% initialize all the variables
vol = 2;
c = 5; % convergence coefficient
% in principle the summation should not depend on the above coefficient

% relative site in which we deplace the second particle by delta2
site = [1,0];  

delta = sqrt(2)*[1/2, 1/2]; % displacement of the electron
delta2 = sqrt(2)*[-1/2, -1/2];

G_cut = 10;
R_cut = 10;

% the final summation is given by 1/gam[sum1-sum2+sum3-sum4];
sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0;

%% evaluation of all the terms for each alpha
index = 1;
% ewald routine

alpha = 1
sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0;
sum12 = 0; sum22 = 0; sum32 = 0; sum42 = 0;
% sum1
for Gx = -G_cut:G_cut+1
    for Gy = -G_cut:G_cut+1
        G = pi*sqrt(2)*[Gx, Gy];
        sum1 = sum1+((cos(G*(-delta)')-1)*f_phi(-alpha/2,G*G'/4/c));
        sum12 = sum12+((cos(G*(-delta2)')-1)*f_phi(-alpha/2,G*G'/4/c));
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
        R = sqrt(2)*[Rx, Ry];
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
%% modification of energy 17
dist = sqrt(2)*site;
% evaluate directly the energy difference
diff1_f = 1/norm(dist-delta)^alpha;
diff1_r = 1/norm(dist-delta+delta2)^alpha;
diff2_f = 1/norm(-delta2-dist)^alpha;
diff2_r = 1/norm(-delta2-dist+delta)^alpha;
% update the real energy
energy2(index) = gap1(index)+gap2(index)-diff1_f+diff1_r-diff2_f+diff2_r;

%% second time for \alpha=\infty
%% initialize all the variables
% relative site in which we deplace the second particle by delta2
site = [1,0];  

delta = sqrt(2)*[1/2, 1/2]; % displacement of the electron
delta2 = sqrt(2)*[-1/2, -1/2];


% the final summation is given by 1/gam[sum1-sum2+sum3-sum4];
sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0;

%% evaluation of all the terms for each alpha
index = 2;
% ewald routine

alpha = 30
sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0;
sum12 = 0; sum22 = 0; sum32 = 0; sum42 = 0;
% sum1
for Gx = -G_cut:G_cut+1
    for Gy = -G_cut:G_cut+1
        G = pi*sqrt(2)*[Gx, Gy];
        sum1 = sum1+((cos(G*(-delta)')-1)*f_phi(-alpha/2,G*G'/4/c));
        sum12 = sum12+((cos(G*(-delta2)')-1)*f_phi(-alpha/2,G*G'/4/c));
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
        R = sqrt(2)*[Rx, Ry];
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
%% modification of energy 17
dist = sqrt(2)*site;
% evaluate directly the energy difference
diff1_f = 1/norm(dist-delta)^alpha;
diff1_r = 1/norm(dist-delta+delta2)^alpha;
diff2_f = 1/norm(-delta2-dist)^alpha;
diff2_r = 1/norm(-delta2-dist+delta)^alpha;
% update the real energy
energy2(index) = gap1(index)+gap2(index)-diff1_f+diff1_r-diff2_f+diff2_r;

x = linspace(0,1,10);

lr_vs_nn_1 = x*gap1(1)+(1-x)*gap1(2);
lr_vs_nn_2 = x*energy2(1)+(1-x)*energy2(2);
figure
plot(x,lr_vs_nn_1,x,lr_vs_nn_2) 
xlabel('$x$','interpreter','latex')
ylabel('Energy $e^2/vol$','interpreter','latex')
title('Energy gap from LR potential to NN (squared lattice)')  

cd Im/
    saveas(gcf,'Lr_vs_nn_squared.eps','epsc')
cd ..