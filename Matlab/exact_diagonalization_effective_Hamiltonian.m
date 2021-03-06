% this script is intended to be the implementation of the exact
% diagonalization for the effective Hamiltonian for our system

% clc; clear all; close all;

% dimension of the lattice to study. Needs to be and even number, in order
% to preserve the periodic boundary conditions
n = 10;

% For the case of n=4 the lattice at the beginning will byou can obtain the precise symbolic form of the 3-by-3 Hilbert matrix:
% 
%     0   1   0   1
%     1   0   1   0
%     0   1   0   1
%     1   0   1   0   
% 
%% Matrix of the ground state
% we create it in such a way that there is a particle every two sites (half
% filled)

A = zeros(n);
for i = 1:n
    for j = 1:n
        if mod(i+j,2)==0
            A(i,j) = 1;
        end
    end
end

ground_state = reshape(A,[1,n^2]);
% ground_state =
% 
%   1   0   1   0   0   1   0   1   1   0   1   0   0   1   0   1

States = zeros(1,n^2+2);

%% First movements state
% the matrix will be made by vectors of the following form
%
% [         (n*n)         ,   x   ,   y   ]
%             ^               ^       ^
%             |               |       |
%       state described     label   position
%       by 0's and 1's
% 
% label    : gives the type of the state we are dealing
% position : gives the position of the particle created wrt the ground
%            state

possible_movements = [0,-1; 1,0; 0,1; -1,0];

number_of_states = 0;
for i=1:n^2
    if ground_state(i) == 1
        label = 1;
        for d = possible_movements'
            number_of_states = number_of_states + 1;
            s = f_swap(ground_state,i,d);
            States(number_of_states,:) = [s, label, i];
            label = label+1; % every time I create a different conf
        end
    end    
end
disp('Matrix of the First displacements finished')
     
%% second movements states
%
% the algorithm goes as follows: 
% take a first excited state, create the second excited state and check
% if it yet exists; if not, append it at the end and modify the effective
% Hamiltonian
% 

% useful Shorthand notations
nSt = number_of_states;

% Hubbard model parameter
t_int = linspace(0.2,0.05,20);
% t_int=1;

% ------------- Begin of the for loop for all the t's


States2 = zeros(1,n^2+2);
eff_H = zeros(nSt);
for fds = 1:nSt
    % fds = first displaced state
    index = States(fds,end);
    x = 1+mod(index-1,n);
    y = ceil(index/n);
    s = States(fds,1:n^2);
    switch States(fds,end-1)
        case 1
            d = [0,1];
            % d = displacement
            rc = [-1,-1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 5;
            % nc = new configuration
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,fds);
            %% second possible displacement
            rc = [1,-1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 6;
            % nc = new configuration
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,fds);
            %% third possible dislpacement
            dv = [1,1];
            % dv = diagonal displacement
            rc = [-1,-1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 2;
            % nc = new configuration
            [eff_H] = f_disp_diag(s,site,dv,States,v,eff_H,nSt,fds,nc);
            %% fourth possible dislpacement
            rc = [0,-1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 2;
            % nc = new configuration
            [eff_H] = f_disp_diag(s,site,dv,States,v,eff_H,nSt,fds,nc);
        case 2
            d = [-1,0];
            % d = displacement
            rc = [1,1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 5;
            % nc = new configuration
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,fds);
            %% second possible displacement
            rc = [1,-1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 6;
            % nc = new configuration
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,fds);
            %% third possible dislpacement
            dv = -[1,1];
            % dv = diagonal displacement
            rc = [1,1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 1;
            % nc = new configuration
            [eff_H] = f_disp_diag(s,site,dv,States,v,eff_H,nSt,fds,nc);
            %% fourth possible dislpacement
            rc = [1,0];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 1;
            % nc = new configuration
            [eff_H] = f_disp_diag(s,site,dv,States,v,eff_H,nSt,fds,nc);
        case 3
            d = [0,-1];
            % d = displacement
            rc = [1,1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 5;
            % nc = new configuration
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,fds);
            %% second possible displacement
            rc = [-1,1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 6;
            % nc = new configuration
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,fds);
            %% third possible dislpacement
            dv = -[1,1];
            % dv = diagonal displacement
            rc = [1,1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 4;
            % nc = new configuration
            [eff_H] = f_disp_diag(s,site,dv,States,v,eff_H,nSt,fds,nc);
            %% fourth possible dislpacement
            rc = [0,1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 4;
            % nc = new configuration
            [eff_H] = f_disp_diag(s,site,dv,States,v,eff_H,nSt,fds,nc);
        case 4
            d = [1,0];
            % d = displacement
            rc = -[1,1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 5;
            % nc = new configuration
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,fds);
            %% second possible displacement
            rc = [-1,1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 6;
            % nc = new configuration
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,fds);
            %% third possible dislpacement
            dv = [1,1];
            % dv = diagonal displacement
            rc = [-1,-1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 3;
            % nc = new configuration
            [eff_H] = f_disp_diag(s,site,dv,States,v,eff_H,nSt,fds,nc);
            %% fourth possible dislpacement
            rc = [-1,0];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 3;
            % nc = new configuration
            [eff_H] = f_disp_diag(s,site,dv,States,v,eff_H,nSt,fds,nc);
    end
end
% ignoring the fake line created accidentally
eff_H = [eff_H(:,1:nSt), eff_H(:,nSt+2:end)];
% making the effective Hamiltonian symmetric
for i = 1:length(eff_H(:,1))
    for j = nSt+1:length(eff_H(1,:))
        eff_H(j,i) = eff_H(i,j);
    end
end

% potential part of the Hamiltonian
[E1,E2] = f_madelung(n); 
potential = diag([E1*ones(1,nSt),E2*ones(1,length(States2(:,1))-1)]);
eff_H = eff_H + potential;

[~,D] = eig(eff_H);
d = sort(diag(D));
% ------------ end of the for loop 

disp('Effective Hamiltonian finished')

%% Plot the data

% figure
% plot(data(:,1),data(:,2),data(:,1),data(:,3))
% xlabel('V/t','interpreter','latex')
% ylabel('E/t','interpreter','latex')
% title('something')
% grid on

figure
histogram(d,50)
title('exact diagonalization eigenvalues effective Hamiltonian')
xlabel('energy spectrum')
ylabel('frequency')
% dim = [0.2 0.5 0.3 0.3];
% str = {'Square lattice','minimum image convention used'};
% annotation('textbox',dim,'String',str,'FitBoxToText','on');

% cd Im
% figurename=['valori_provvisori_ED.eps'];
% saveas(gcf,figurename,'epsc')
% cd ..
