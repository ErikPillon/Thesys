clc; clear all, close all;
% this script is intended for evaluating all the energy costs for the same
% displacement changing only the ratio between the two Bravais vectors

delta = [1,1]/2;
site = [1,0];
theta = pi/2;

ratio_range = linspace(1/sqrt(3),sqrt(3),20);
count = 1;
E1 = [0,0];
for ratio = ratio_range
    [x,y] = f_two_excitations_energy_gap(delta,-delta,site,ratio,pi/2);
    E1(count,:) = [x,y];
    count = count+1;
end

cd Im/
    figure
    plot(ratio_range,E1(:,1),ratio_range,E1(:,2))
    title('Energy gaps vs ratio between Bravais vectors')  
    xlabel('ratio $a/b$','interpreter','latex')
    ylabel('Energy displacements')
    saveas(gca,'EnergyVsRatio.eps','epsc')
cd ..


% What follows is absolutely useless since it shows the same behaviour of
% the first figure
%
% %% now we change the configuration, i.e., we choose a different site
% site = [0,1];
% count = 1;
% E2 = [0,0];
% for ratio = ratio_range
%     [x,y] = f_two_excitations_energy_gap(delta,-delta,site,ratio,pi/2);
%     E2(count,:) = [x,y];
%     count = count+1;
% end
% 
% figure
% plot(ratio_range,E2(:,1),ratio_range,E2(:,2))
% title('Energy gaps vs ratio between Bravais vectors')  
% xlabel('ratio $a/b$','interpreter','latex')
% ylabel('Energy displacements')