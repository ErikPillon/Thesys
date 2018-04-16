clc; clear all; close all;

%% TRIANGULAR LATTICE
%
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

v = 0;
V = 1;
dx = 0.5;
dy = 0.5;
E1 = energy1;
E2 = energy2;
t = 1;

Energy = @(v,t,kx,ky) [E1, -2*v*cos(dx/2-dy/2),0,0,-t*exp(1i*kx/2),-t*exp(-1i*kx/2);...
     -2*v*cos(dx/2-dy/2), E1,0,0,-t*exp(1i*ky/2),-t*exp(-1i*ky/2);... 
     0,0,E1,-2*v*cos(dx/2-dy/2),-t*exp(-1i*kx/2),-t*exp(1i*kx/2);...
     0,0,-2*v*cos(dx/2-dy/2),E1,-t*exp(-1i*ky/2),-t*exp(1i*ky/2);...
     -t*exp(-1i*kx/2),-t*exp(-1i*ky/2),-t*exp(1i*kx/2),-t*exp(1i*ky/2),E2,0;...
     -t*exp(1i*kx/2),-t*exp(1i*ky/2),-t*exp(-1i*kx/2),-t*exp(-1i*ky/2),0,E2];

x = 1;
y = 1;
Band = zeros(20,20,6);

for kx = linspace(-pi,pi,20)
    y = 1;
    for ky = linspace(-pi,pi,20)
        [~,D] = eig(Energy(v,t,kx,ky));
        Band(x,y,:) = diag(D);
        y = y+1;
    end
    x = x+1;
end
a=V/t;
figure
titlemessage = ['Electronic Energy bands \break $V/t=$ ' num2str(a) '; $v=$' num2str(v)];
title(titlemessage,'interpreter','latex')
zlabel('Energy $(e^2/vol)$','interpreter','latex')
for i=1:6
    hold on
    mesh(Band(:,:,i))
end
view(7,4)

% figurename=['Bands_Vt_' num2str(a) '_v_'  num2str(v) '.eps'];
% saveas(gcf,figurename,'epsc')

%% try now to express the dispersion relation
Disp_r = zeros(6,1);
count = 1;

k_gm = linspace(0,pi,20);
for k = k_gm
    [~,D] = eig(Energy(v,t,k,k));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

k_mx = linspace(pi,0,20);
for k = k_mx
    [~,D] = eig(Energy(v,t,pi,k));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

k_xg = linspace(pi,0,20);
for k = k_xg
    [~,D] = eig(Energy(v,t,k,0));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

figure
titlemessage = ['Electronic dispersion relation \break $V/t=$ ' num2str(a) '; $v=$' num2str(v)];
title(titlemessage,'interpreter','latex')
ylabel('Energy $(e^2/vol)$','interpreter','latex')
for i=1:6
    hold on
    plot(linspace(0,3,count-1),Disp_r(i,:))
end

cd Im
figurename=['DispR_Vt_' num2str(a) '_v_'  num2str(v) '.eps'];
saveas(gcf,figurename,'epsc')
cd ..