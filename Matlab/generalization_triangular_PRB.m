clc; clear all; close all;
% this algorithm aims to plot the energy needed for one displacement letting
% \alpha varying

% relative site in which we deplace the second particle by delta2
site = [1,0];  
delta = [1/2, 1/2]; % displacement of the electron
delta2 = [-1/2, -1/2];

count = 1;
X = [0,0];
alpha_range = [(0.4:0.05:3),(4:0.5:7),(8:20)];
for alpha=alpha_range
    [x,y]=f_two_excitations_energy_gap(delta,-delta,[0,1],sqrt(3),pi/2,alpha);
    X(count,:)=[x,y];
    count = count+1;
end

%% plot the results and save them
figure
plot(alpha_range, X(:,1), alpha_range, X(:,2))
grid on
title('Energy gap from ground state, triangular lattice','interpreter','latex')
xlabel('\alpha')
ylabel('Energy gap $(e^2/vol)$', 'interpreter','latex')
legend('One displacement','Two displacements')
cd Im/
    saveas(gcf,'Energy_gap_triangularLattice.eps','epsc')
cd ..

% plot the second figure zoomed 
figure(2)
plot(alpha_range(20:35), X(20:35,1), alpha_range(20:35), X(20:35,2))
grid on
title('Energy gap from ground state, triangular lattice','interpreter','latex')
xlabel('\alpha')
ylabel('Energy gap $(e^2/vol)$', 'interpreter','latex')
legend('One displacement','Two displacements')
cd Im/
    saveas(gcf,'Energy_gap_triangularLattice_zoomed.eps','epsc')
cd ..