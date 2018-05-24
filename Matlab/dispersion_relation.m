clc; clear all; close all;
%% TRIANGULAR LATTICE
%
% script for one and two excitations on a trinangular lattice, that is 
% theta = pi/2 and ratio = 2*sin(pi/3)=sqrt(3) 

% initialize the variables for the algorithm
delta = [1/2, 1/2]; % displacement
site = [1,0]; % site of the second displacement
delta2 = [-1/2, -1/2]; % second displacement 
ratio = sqrt(3); 
theta = pi/2;

[E1,E2] = f_two_excitations_energy_gap(delta,delta2,site,ratio,theta)

v = 0.5;
V = 1;
t = 0.5;
sr = 1;

Energy = @(v,t,kx,ky) [E1, -v*exp(1i*(kx/2-ky/2)),0,0,-t*exp(1i*kx/2),-t*exp(-1i*kx/2);...
     -v*exp(-1i*(kx/2-ky/2)), E1,0,0,-t*exp(1i*ky/2),-t*exp(-1i*sr*ky/2);... 
     0,0,E1,-v*exp(1i*(kx/2-ky/2)),-t*exp(-1i*kx/2),-t*exp(1i*kx/2);...
     0,0,-v*exp(-1i*(kx/2-ky/2)),E1,-t*exp(-1i*ky/2),-t*exp(1i*ky/2);...
     -t*exp(-1i*kx/2),-t*exp(-1i*ky/2),-t*exp(1i*kx/2),-t*exp(1i*ky/2),E2,0;...
     -t*exp(1i*kx/2),-t*exp(1i*ky/2),-t*exp(-1i*kx/2),-t*exp(-1i*ky/2),0,E2];
% Energy = @(v,t,kx,ky) [E1, -2*v*cos(kx/2-ky/2),0,0,-t*exp(1i*kx/2),-t*exp(-1i*kx/2);...
%      -2*v*cos(kx/2-ky/2), E1,0,0,-t*exp(1i*ky/2),-t*exp(-1i*sr*ky/2);... 
%      0,0,E1,-2*v*cos(kx/2-ky/2),-t*exp(-1i*kx/2),-t*exp(1i*kx/2);...
%      0,0,-2*v*cos(kx/2-ky/2),E1,-t*exp(-1i*ky/2),-t*exp(1i*ky/2);...
%      -t*exp(-1i*kx/2),-t*exp(-1i*ky/2),-t*exp(1i*kx/2),-t*exp(1i*ky/2),E2,0;...
%      -t*exp(1i*kx/2),-t*exp(1i*ky/2),-t*exp(-1i*kx/2),-t*exp(-1i*ky/2),0,E2];

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
figurename=['Bands_Vt_' num2str(a) '_v_'  num2str(v) '.eps'];
saveas(gcf,figurename,'epsc')
cd ..

%% evaluation of the dispersion relation
Disp_r = zeros(6,1);
count = 1;

k_gm = linspace(0,pi,100);
for k = k_gm
    [~,D] = eig(Energy(v,t,k,k));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

k_mx = linspace(pi,0,100);
for k = k_mx
    [~,D] = eig(Energy(v,t,pi,k));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

k_xg = linspace(pi,0,100);
for k = k_xg
    [~,D] = eig(Energy(v,t,k,0));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

figure
subplot(2,2,1)
titlemessage = ['Electronic dispersion relation \break $V/t=$ ' num2str(a) '; $v=$' num2str(v)];
title(titlemessage,'interpreter','latex')
ylabel('Energy $(e^2/vol)$','interpreter','latex')
for i=1:6
    hold on
    plot(linspace(0,3,count-1),Disp_r(i,:))
end
xlabel('\Gamma                            X                             M                            \Gamma')
set(gca,'xtick',[])

%% evaluation of the dispersion relation following the second possible path
Disp_r = zeros(6,1);
count = 1;

k_gm = linspace(0,-pi,100);
for k = k_gm
    [~,D] = eig(Energy(v,t,k,-k));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

k_mx = linspace(pi,0,100);
for k = k_mx
    [~,D] = eig(Energy(v,t,-pi,k));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

k_xg = linspace(-pi,0,100);
for k = k_xg
    [~,D] = eig(Energy(v,t,k,0));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

subplot(2,2,2)
titlemessage = ['Second possible path \break Electronic dispersion relation \break $V/t=$ ' num2str(a) '; $v=$' num2str(v)];
title(titlemessage,'interpreter','latex')
ylabel('Energy $(e^2/vol)$','interpreter','latex')
for i=1:6
    hold on
    plot(linspace(0,3,count-1),Disp_r(i,:))
end
xlabel('\Gamma                            X                             M                            \Gamma')
set(gca,'xtick',[])

%% evaluation of the dispersion relation
Disp_r = zeros(6,1);
count = 1;

k_gm = linspace(0,pi,100);
for k = k_gm
    [~,D] = eig(Energy(v,t,k,k));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

k_mx = linspace(pi,0,100);
for k = k_mx
    [~,D] = eig(Energy(v,t,pi,k));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

k_xg = linspace(pi,0,100);
for k = k_xg
    [~,D] = eig(Energy(v,t,k,0));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

subplot(2,2,3)
titlemessage = ['Electronic dispersion relation \break $V/t=$ ' num2str(a) '; $v=$' num2str(v)];
title(titlemessage,'interpreter','latex')
ylabel('Energy $(e^2/vol)$','interpreter','latex')
for i=1:6
    hold on
    plot(linspace(0,3,count-1),Disp_r(i,:))
end
xlabel('\Gamma                            X                             M                            \Gamma')
set(gca,'xtick',[])

%% evaluation of the dispersion relation following the second possible path
Disp_r = zeros(6,1);
count = 1;

k_gm = linspace(0,pi,100);
for k = k_gm
    [~,D] = eig(Energy(v,t,k,-k));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

k_mx = linspace(pi,0,100);
for k = k_mx
    [~,D] = eig(Energy(v,t,pi,-k));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

k_xg = linspace(pi,0,100);
for k = k_xg
    [~,D] = eig(Energy(v,t,k,0));
    Disp_r(:,count) = diag(D);
    count = count+1;
end

subplot(2,2,4)
titlemessage = ['Second possible path \break Electronic dispersion relation \break $V/t=$ ' num2str(a) '; $v=$' num2str(v)];
title(titlemessage,'interpreter','latex')
ylabel('Energy $(e^2/vol)$','interpreter','latex')
for i=1:6
    hold on
    plot(linspace(0,3,count-1),Disp_r(i,:))
end
xlabel('\Gamma                            X                             M                            \Gamma')
set(gca,'xtick',[])

% cd Im
% figurename=['DispR_TriangL_Vt_' num2str(a) '_v_'  num2str(v) '.eps'];
% saveas(gcf,figurename,'epsc')
% cd ..

% target = E1-2*t1*(cos(2*kx)+cos(2*ky))-2*t2*(cos(kx+ky)+cos(kx-ky));

% %% Energy level for the gamma point
% count = 1;
% V = 1;
% Gamma = [0];
% interval = linspace(-2,2,50);
% t_int = [linspace(2,1,50),linspace(1,0.05,50)];
% figure
% hold on
% plot(interval,zeros(length(interval),1));
% for i = interval
%     for t = t_int
%         v = i*t;
%         [~,D] = eig(Energy(v,t,0,0));
%         d = sort(diag(D));
%         if abs(d(1))<0.005
%             plot(v/t,1/t,'r*')
%         end
%     end
% end
% title('Zero energy level of the $\Gamma$ point on the triangular lattice','interpreter','latex')
% xlabel('v/t','interpreter','latex')
% ylabel('V/t','interpreter','latex')
% 
% cd Im
% figurename=['zero_energy_level_TriangL.eps'];
% saveas(gcf,figurename,'epsc')
% cd ..
% 
% %% Minimum & Maximum Energy level for the gamma point
% % the graph we want is E/t vs V/t 
% count = 1;
% V = 1;
% t_int = linspace(2,0.05,20);
% data = [0,0,0];
% for t = t_int
%     v = 0;
%     [~,D] = eig(Energy(v,t,0,0));
%     d = sort(diag(D));
%     data(count,:) = [V/t,d(1)/t,d(end)/t]
%     count = count+1;
% end
% figure
% plot(data(:,1),data(:,2),data(:,1),data(:,3))
% xlabel('V/t','interpreter','latex')
% ylabel('E/t','interpreter','latex')
% title('something')
% grid on
% 
% cd Im
% figurename=['zero_energy_level_TriangL.eps'];
% saveas(gcf,figurename,'epsc')
% cd ..
% 

