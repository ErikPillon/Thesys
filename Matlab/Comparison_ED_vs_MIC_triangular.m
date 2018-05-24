%clc; clear all; close all;
% put the analytical model, put the exact diagonalization algorithm with
% the values evaluated with minimum image convention and see the
% differences

% run the same code of the script of exact diagonalization

v = 0.05;
t = 0.05;
v = 0;

ed

[E1,E2] = f_two_excitations_energy_gap([1,1]/2,-[1,1]/2,[1,0],1,pi/2);

Energy = @(v,t,kx,ky) [E1, -2*v*cos(kx/2-ky/2),0,0,-t*exp(1i*kx/2),-t*exp(-1i*kx/2);...
     -2*v*cos(kx/2-ky/2), E1,0,0,-t*exp(1i*ky/2),-t*exp(-1i*ky/2);... 
     0,0,E1,-2*v*cos(kx/2-ky/2),-t*exp(-1i*kx/2),-t*exp(1i*kx/2);...
     0,0,-2*v*cos(kx/2-ky/2),E1,-t*exp(-1i*ky/2),-t*exp(1i*ky/2);...
     -t*exp(-1i*kx/2),-t*exp(-1i*ky/2),-t*exp(1i*kx/2),-t*exp(1i*ky/2),E2,0;...
     -t*exp(1i*kx/2),-t*exp(1i*ky/2),-t*exp(-1i*kx/2),-t*exp(-1i*ky/2),0,E2];

% collect data for several points in the reciprocal lattice
x = 1;
y = 1;

Band = zeros(6,1);

count = 0;
for kx = linspace(-pi,pi,100)
    for ky = linspace(-pi,pi,100)
        count = count+1;
        [~,D] = eig(Energy(v,t,kx,ky));
        Band(:,count) = diag(D);
    end
end
disp('Finished collection of values from analytical model')
eigen = reshape(Band,6*count,1);

% plotting data
figure
histogram(eig(Hamiltonian),100,'Normalization','probability')
hold on
histogram(eigen,100,'Normalization','probability')
% title({'Comparison ED values with analytical model in triangular lattice';['$V/t$ = ',num2str(1/t),' $v/t$ = ',num2str(v/t)]},'interpreter','latex')
title('Comparison ED values with analytical model')
xlabel('Energy spectrum')
legend('ED','Analytical model')

% cd Im/
%     saveas(gcf,'ED_vs_AnalyticalModel.eps','epsc')
% cd ..

% cd Im/
%     saveas(gcf,'ED_vs_AnalyticalModel_triangular.eps','epsc')
% cd ..

% cd Im/
%     saveas(gcf,'ED_vs_AnalyticalModel_fugacity.eps','epsc')
% cd ..
