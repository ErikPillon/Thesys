clc; clear all; close all;
% this script aims at showing how the Madelung energy changes in a finite
% size cluster wrt the cluster size.
% For n=L going to infinity we should recollect the classical values

data = [0,0,0];

count = 0;
n_range = [4:2:20];
for n = n_range
    count = count+1;
    [energy1,energy2] = f_madelung(n);
    data(count,:) = [1/n, energy1, energy2];
end

figure
plot(data(:,1),data(:,2),'*',data(:,1),data(:,3),'*');
xlabel('1/L','interpreter','latex')
ylabel('Energy displacement')
title({'Madelung energies evaluated with' '\it minimum image convention'})
legend('One displacement energy', 'Two displacements energy')
axis([0.03,0.28,0.6,0.65])

% save the image
cd Im/
saveas(gcf,'Minimum_Image_Convention.eps','epsc')
cd ..