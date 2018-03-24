% the goal of this script is to show how the energy behaves wrt different
% movements in order to find which are the ones that requires the least
% energy for happening
% 
% ONLY TWO POSSIBLES EXCITATIONS ARE ALLOWED!

clc; clear all; close all;
% in this specific case we will consider a square lattice

delta = [1/2, 1/2]; % first displacement
ratio = 1;
theta = pi/2;

possible_movements = [1,1;1,-1;-1,1;-1,-1]/2;
radius = 2;
count = 1;
possible_sites = [0, 0];
for x= -radius:radius
    for y = -radius:radius
        if ~(x == 0 & y == 0)
            possible_sites(count, :) = [x,y];
            count=count+1;
        end
    end
end

i = 1;
second_excitation = zeros(1,5);
for site = possible_sites'
    site = site';
    for delta2 = possible_movements'
        delta2 = delta2';
        if ~(site+delta2 == delta)
            [energy1, energy2] = f_two_excitations_energy_gap(delta, delta2, site); 
            second_excitation(i,:) = [site, delta2, energy2];
            i = i+1
        end
    end
end
disp('SUCCESFULLY COMPLETED!')
[energy1, energy2] = f_two_excitations_energy_gap(delta, delta2, site); 
%% sort the matrix
% the arrangement is made in such a way that the energy is displayed in the
% last column the energy, in the first two columns the (i,j) site in which
% move the particle of a displacement delta2 (3rd and 4th column)
[~,I] = sort(second_excitation(:,end));
energy_sorted = second_excitation(I,:);
energy_sorted(:,end) = energy_sorted(:,end)-2*energy1

data = fopen('doubleMovementsData.txt','w');
fprintf(data,'%4d %4d %4d %4d %12.8f\n',energy_sorted);
