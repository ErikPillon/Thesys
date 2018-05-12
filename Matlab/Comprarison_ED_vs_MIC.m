% put the analytical model, put the exact diagonalization algorithm with
% the values evaluated with minimum image convention and see the
% differences

[E1,E2] = f_two_excitations_energy_gap([1,1]/2,-[1,1]/2,[1,0],1,pi/2);

Energy = @(v,t,kx,ky) [E1, 0,0,0,-t*exp(1i*kx/2),-t*exp(-1i*kx/2);...
     0, E1,0,0,-t*exp(1i*ky/2),-t*exp(-1i*ky/2);... 
     0,0,E1,0,-t*exp(-1i*kx/2),-t*exp(1i*kx/2);...
     0,0,0,E1,-t*exp(-1i*ky/2),-t*exp(1i*ky/2);...
     -t*exp(-1i*kx/2),-t*exp(-1i*ky/2),-t*exp(1i*kx/2),-t*exp(1i*ky/2),E2,0;...
     -t*exp(1i*kx/2),-t*exp(1i*ky/2),-t*exp(-1i*kx/2),-t*exp(-1i*ky/2),0,E2];
 
% collect data for several points in the reciprocal lattice
x = 1;
y = 1;
Band = zeros(6,1);

v = 0;
t = 0.05;
count = 0;
for kx = linspace(-pi,pi,10)
    for ky = linspace(-pi,pi,10)
        count = count+1;
        [~,D] = eig(Energy(v,t,kx,ky));
        Band(:,count) = diag(D);
    end
end
disp('Finished collection of values from analytical model')
eigen = reshape(Band,6*count,1);

% run the same code of the script of exact diagonalization
exact_diagonalization_effective_Hamiltonian

% plotting data
figure
histogram(eig(eff_H),50,'Normalization','probability')
hold on
histogram(eigen,50,'Normalization','probability')
title('Comparison ED values with analytical model')
xlabel('Energy spectrum')
legend('ED','Analytical model')

cd Im/
    saveas(gcf,'ED_vs_AnalyticalModel.eps','epsc')
cd ..
