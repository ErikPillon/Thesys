% the goal of this script is to show how the energy behaves wrt different
% movements in order to find which are the ones that requires the least
% energy for happening
% 
% Scheme: all the two possible excitations are saved in the file
% doubleMovementsData.txt, then the idea is to take thaose configurations
% with the less energy and perturb them
%
% UP TO 3 POSSIBLE EXCITATIONS ARE ALLOWED!

clc; clear all; close all;
% in this specific case we will consider a square lattice
possible_excitations;
A = importdata('doubleMovements_alpha1_theta1.txt')

delta = [1/2, 1/2]; % first displacement
ratio = 1;
theta = pi/2;
alpha = 1;

possible_movements = [1,1;1,-1;-1,1;-1,-1]/2;


radius = 2;
count = 1;

% create the matrix for the third possible displacement
possible_sites2 = [0, 0];
for i= -radius:radius
    for j = -radius:radius
        if ~(i == 0 & j == 0)
            possible_sites2(count, :) = [i,j];
            count=count+1;
        end
    end
end

%% main loop
conf = 1;
state = 1;
while state<= 50
    site1 = A(state,1:2);
    delta1 = A(state,3:4);
    for site2 = possible_sites2'
        site2 = site2';
        if ~(site2 == site1)
            for delta2 = possible_movements'
                    delta2 = delta2';
                    % the following if statements are for checking that there
                    % are not superpositions with the already created displacements
                    if ~(site1+delta1 == delta)
                        if ~(site2+delta2 == delta)
                            if ~(site2+delta2 == site1+delta1)
                                [e1, ~, e3] = f_3_excitations(delta, delta1, site1, delta2, site2,ratio,theta,alpha);
                                third_excitation(conf,:) = [site1, delta1, site2, delta2, e3];
                                conf = conf+1
                            end
                        end
                    end
            end
        end
    end
    state = state+1
end

%% sort the matrix
% the arrangement is made in such a way that the energy is displayed in the
% last column the energy, in the first two columns the (i,j) site in which
% move the particle of a displacement delta2 (3rd and 4th column)
[~,I] = sort(third_excitation(:,end));
Energy_sorted = third_excitation(I,:);
Energy_sorted(:,end) = Energy_sorted(:,end)-e1

% save('tripleMovements.txt', 'Energy_sorted', '-ascii', '-double', '-tabs')

