clc; clear all; close all;

%% SQUARE LATTICE CUPRATE
%
% script for one and two excitations on a square lattice, that is 
% theta = pi/2 and ratio = 1 with norm(displacement)=1, i.e, |a|=sqrt(2) 


%% initialize all the variables
alpha = 1;
vol = 2;%the more rapid will be the convergence of the  summations.
c = 5; % convergence coefficient

delta = sqrt(2)*[1/2, 1/2]; % displacement of the electron
site = [0,1];
delta2 = sqrt(2)*[-1/2, -1/2];

G_cut = 10;
R_cut = G_cut;

% the final summation is given by 1/gam[sum1-sum2+sum3-sum4];
sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0;

%% evaluation of all the terms
% sum1
for Gx = -G_cut:G_cut
    for Gy = -G_cut:G_cut
        G = pi*sqrt(2)*[Gx, Gy];
        sum1 = sum1+((cos(G*delta')-1)*f_phi(-alpha/2,G*G'/4/c));
    end
end
sum1 = sum1*pi*c^(alpha/2-1)/vol;

% sum2
sum2 = f_selfinteraction(norm(delta),c,alpha)-f_selfinteraction(0,c,alpha);

% sum3
for Rx = -R_cut:R_cut
    for Ry = -R_cut:R_cut
        R = sqrt(2)*[Rx, Ry];
        sum3=sum3+(f_phi(alpha/2-1,norm(R+delta)^2*c)-f_phi(alpha/2-1,norm(R)^2*c));
    end
end
sum3 = sum3*c^(alpha/2);

% sum4
sum4 = f_phi(alpha/2-1,norm(delta)^2*c)-f_phi(alpha/2-1,0);
sum4 = c^(alpha/2)*sum4;

energy1 = (sum1-sum2+sum3-sum4)/gamma(alpha/2);

dist = sqrt(2)*site;

% evaluate directly the energy difference
diff1_f = 1/norm(dist-delta)^alpha;
diff1_r = 1/norm(dist-delta+delta2)^alpha;
diff2_f = 1/norm(-delta2-dist)^alpha;
diff2_r = 1/norm(dist)^alpha;

gap = -diff1_f+diff1_r-diff2_f+diff2_r;

energy2 = 2*energy1 +gap;

%% implementation of the band energy and dispersion relations

V = 1;
E1 = energy1;
E2 = energy2;
t = 1;
v = 0*t;
% t = 0.64/2/sqrt(2);

%this matrix needs to be modified from the previous one by adding 8 terms 
Energy = @(v,t,kx,ky) [E1, -2*v*cos(kx/2-ky/2),0,-2*v*cos(kx/2+ky/2),-t*exp(1i*kx/2),-t*exp(-1i*kx/2);...
     -2*v*cos(kx/2-ky/2), E1,-2*v*cos(kx/2+ky/2),0,-t*exp(1i*ky/2),-t*exp(-1i*ky/2);... 
     0,-2*v*cos(kx/2+ky/2),E1,-2*v*cos(kx/2-ky/2),-t*exp(-1i*kx/2),-t*exp(1i*kx/2);...
     -2*v*cos(kx/2+ky/2),0,-2*v*cos(kx/2-ky/2),E1,-t*exp(-1i*ky/2),-t*exp(1i*ky/2);...
     -t*exp(-1i*kx/2),-t*exp(-1i*ky/2),-t*exp(1i*kx/2),-t*exp(1i*ky/2),E2,0;...
     -t*exp(1i*kx/2),-t*exp(1i*ky/2),-t*exp(-1i*kx/2),-t*exp(-1i*ky/2),0,E2];
 
x = 1;
y = 1;
Band = zeros(20,20,6);
 
for kx = linspace(-pi,pi,50)
    y = 1;
    for ky = linspace(-pi,pi,50)
        [~,D] = eig(Energy(v,t,kx,ky));
        Band(x,y,:) = diag(D);
        y = y+1;
    end
    x = x+1;
end
a=V/t;
figure
titlemessage = ['Electronic Energy bands \break $V/t=$ ' num2str(a) '; $v=$' num2str(v)];
title(titlemessage,'interpreter','latex')
zlabel('Energy $(e^2/vol)$','interpreter','latex')
for i=1:6
    hold on
    mesh(Band(:,:,i))
end
view(7,4)

cd Im
figurename=['Bands_cuprate_Vt_' num2str(a) '_v_'  num2str(v) '.eps'];
saveas(gcf,figurename,'epsc')
cd ..

%% evaluation of the dispersion relation
Disp_r = zeros(6,1);
count = 1;

k_gm = linspace(0,pi,50);
for k = k_gm
    [~,D] = eig(Energy(v,t,k,k));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

k_mx = linspace(pi,0,50);
for k = k_mx
    [~,D] = eig(Energy(v,t,pi,k));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

k_xg = linspace(pi,0,50);
for k = k_xg
    [~,D] = eig(Energy(v,t,k,0));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

figure
titlemessage = ['Electronic dispersion relation \break $V/t=$ ' num2str(a) '; $v=$' num2str(v)];
title(titlemessage,'interpreter','latex')
ylabel('Energy $(e^2/vol)$','interpreter','latex')
for i=1:6
    hold on
    plot(linspace(0,3,count-1),Disp_r(i,:))
end
xlabel('\Gamma                            X                             M                            \Gamma')
set(gca,'xtick',[])

cd Im
figurename=['DispR_cuprate_Vt_' num2str(a) '_v_'  num2str(v) '.eps'];
saveas(gcf,figurename,'epsc')
cd ..

%% Energy level for the gamma point
count = 1;
V = 1;

Gamma = [0];
interval = linspace(-2,2,50);
t_int = [linspace(2,1,50),linspace(1,0.05,50)];
figure
hold on
plot(interval,zeros(length(interval),1));
for i = interval
    for t = t_int
        v = i*t;
        [~,D] = eig(Energy(v,t,0,0));
        d = sort(diag(D));
        if abs(d(1))<0.01
            plot(v/t,1/t,'r*')
        end
    end
end
title('Zero energy level of the $\Gamma$ point on the cuprates','interpreter','latex')
xlabel('v/t','interpreter','latex')
ylabel('V/t','interpreter','latex')

cd Im
figurename=['zero_energy_level_cuprates.eps'];
saveas(gcf,figurename,'epsc')
cd ..
