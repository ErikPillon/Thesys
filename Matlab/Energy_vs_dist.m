% energy vs distance
clc; clear all; close all;

ratio = 1;
theta = pi/2;
delta = [1/2,1/2];
alpha = 3;

energy1 = f_one_excitation_energy_gap([1/2,1/2], 1, pi/2);


if theta == pi/2
   a = [2*cos(atan(ratio)),0];
   b = [0,2*sin(atan(ratio))];
   vol = norm(cross([a,0],[b,0]));  
else
    a = [1, 0];
    rot = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    b = (ratio*rot*a')';
    vol = ratio*sin(theta);
end

possible_movements = [1,1;1,-1;-1,1;-1,-1]/2;
radius = 6; % TO BE CHANGED

count = 1;
possible_sites = [0, 0];
for x= -radius:radius
    for y = -floor(sqrt(radius^2-x^2)):floor(sqrt(radius^2-x^2))
        if ~(x == 0 & y == 0)
            possible_sites(count, :) = [x,y];
            count=count+1;
        end
    end
end

% the displacements must be converted in the lattice coordinate
delta = delta(1)*a+delta(2)*b;
num_plot = 1;
for alpha = [0.66,1,2,5]
    X = [0,0];
    i = 1;
    for site = possible_sites'
        site = site';
        dist = site(1)*a + site(2)*b;
        for delta2 = possible_movements'
            delta2 = delta2';
            delta2 = delta2(1)*a+delta2(2)*b;
            if ~(site+delta2 == delta)
                % evaluate directly the energy difference
                diff1_f = 1/norm(dist-delta)^alpha;
                diff1_r = 1/norm(dist-delta+delta2)^alpha;
                diff2_f = 1/norm(delta2+dist)^alpha;
                diff2_r = 1/norm(dist)^alpha;
                energy2 = diff1_r+diff2_r-diff1_f-diff2_f;
                X(i,:) = [norm(dist),energy2];
                i = i+1;
            end
        end
    end
    % displace for every case the energy
    subplot(2,2,num_plot)
    plot(X(:,1),X(:,2),'*')
    %xlabel('distance')
    %ylabel('correction term')
    axis([0.5,7,-2,0.5])
    titlem = ['\alpha = ' num2str(alpha)];
    title(titlem)
    figurename = ['Energy_vs_alpha.eps'];
    num_plot = num_plot+1;
end
cd Im/
    saveas(gcf,figurename,'epsc')
cd ..
