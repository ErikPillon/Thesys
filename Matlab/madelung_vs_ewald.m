clc; clear all; close all;
% this script aims at showing how the Madelung energy changes in a finite
% size cluster wrt the cluster size.
% For n=L going to infinity we should recollect the classical values

data = [0,0,0];

count = 0;
n_range = [4:2:30];
for n = n_range
    count = count+1;
    [e1,e2] = f_madelung(n);
    data(count,:) = [1/n, e1, e2];
end
[E1,E2] = f_two_excitations_energy_gap([1,1]/2,-[1,1]/2,[0,1],1,pi/2); 

figure
plot(data(:,1),data(:,2),'*',data(:,1),data(:,3),'*')
xlabel('1/L','interpreter','latex')
ylabel('Energy displacement')
title({'Madelung energies evaluated with' '\it minimum image convention'})
legend('One displacement energy', 'Two displacements energy')
axis([0.03,0.28,0.6,0.65])

% save the image
cd Im/
saveas(gcf,'Minimum_Image_Convention.eps','epsc')
cd ..

figure
loglog(data(:,1),abs(data(:,2)-E1),'*',data(:,1),abs(data(:,3)-E2),'*',data(:,1),data(:,1).^(3)/10);
title('Asymptotical behaviour of the error')
xlabel('cluster size')
ylabel('Error')
legend({'One displacement energy', 'Two displacements energy','$1/n^{3}$'},'interpreter','latex')
axis([0.03,0.3,0.3*10^(-5),0.6*10^(-1)])

cd Im/
saveas(gcf,'Behaviour_error_madelung.eps','epsc')
cd ..
