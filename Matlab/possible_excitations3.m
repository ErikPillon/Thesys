% the goal of this script is to show how the energy behaves wrt different
% movements in order to find which are the ones that requires the least
% energy for happening
% 
% UP TO 3 POSSIBLE EXCITATIONS ARE ALLOWED!

clc; clear all; close all;
% in this specific case we will consider a square lattice

delta = [1/2, 1/2]; % first displacement
ratio = 1;
theta = pi/2;
possible_movements = [1,1;1,-1;-1,1;-1,-1]/2;

radius = 2;
count = 1;
possible_sites1 = [0, 0];
for i= -radius:radius
    for j = -radius:radius
        if ~(i == 0 & j == 0)
            possible_sites1(count, :) = [i,j];
            count=count+1;
        end
    end
end

count = 1;
possible_sites2 = [0, 0];
for i= -radius:radius
    for j = -2*radius:2*radius
        if ~(i == 0 & j == 0)
            possible_sites2(count, :) = [i,j];
            count=count+1;
        end
    end
end

%% main loop
i = 1;
third_excitation = zeros(1,9);
for site2 = possible_sites2'
    site2 = site2';
    for site1 = possible_sites1'
        site1 = site1';
        for delta1 = possible_movements'
            delta1 = delta1';
            for delta2 = possible_movements'
                delta2 = delta2';
                % the followinf if statements are for checking that there
                % are not superpositions with the already created displacements
                if ~(site1+delta1 == delta)
                    if ~(site2+delta2 == delta)
                        if ~(site2+delta2 == delta1)
                            [e1, ~, e3] = f_3_excitations(delta, delta1, site1, delta2, site2);
                            third_excitation(i,:) = [site1, delta1, site2, delta2, e3];
                            i = i+1
                        end
                    end
                end
            end
        end
    end
end

%% sort the matrix
% the arrangement is made in such a way that the energy is displayed in the
% last column the energy, in the first two columns the (i,j) site in which
% move the particle of a displacement delta2 (3rd and 4th column)
[~,I] = sort(second_excitation(:,end));
Energy_sorted = second_excitation(I,:);
Energy_sorted(:,end) = Energy_sorted(:,end)-e1


