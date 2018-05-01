% this script is intended to be the implementation of the exact
% diagonalization for the effective Hamiltonian for our system
clc; clear all; close all;

% dimension of the lattice to study. Needs to be and even number, in order
% to preserve the periodic boundary conditions
n = 16;

% Hubbard model parameter
t = 1;

% For the case of n=4 the lattice at the beginning will be 
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
% [         (n*n)         ,   1   ,   1   ]
%             ^               ^       ^
%             |               |       |
%       state described     label   position
%       by 0's and 1's
% 
% label    : gives the type of the state we are dealing
% position : gives the position of the particle created wrt the ground
%            state

possible_movements = [0,-1; 1,0; 0,1; -1,0];

number_of_states = 1;
tic
for i=1:n^2
    if ground_state(i) == 1
        label = 1;
        for d = possible_movements'
            s = f_swap(ground_state,i,d);
            States(number_of_states,:) = [s, label, i];
            label = label+1;
            number_of_states = number_of_states + 1;
        end
    end    
end
disp('time for creating the matrix of all first possible displacements')
toc
     
%% second movements states
%
% the algorithm goes as follows: 
% take a first excited state, create the second excited state and check
% if it yet exists; if not, append it at the end and modify the effective
% Hamiltonian
% 

count2 = 1;
States2 = zeros(1,n^2+2);

% useful Shorthand notations
nSt = number_of_states-1;

tic
eff_H = zeros(2);
for fds = 1:nSt
    % fds = first displaced state
    index = States(fds,end);
    x = ceil(index/n);
    y = mod(index,n);
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
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,2*fds-1);
            %% second possible displacement
            rc = [1,-1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 6;
            % nc = new configuration
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,2*fds);
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
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,2*fds-1);
            %% second possible displacement
            rc = [1,-1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 6;
            % nc = new configuration
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,2*fds);
        case 3
            d = [0,1];
            % d = displacement
            rc = [1,1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 5;
            % nc = new configuration
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,2*fds-1);
            %% second possible displacement
            rc = [-1,1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 6;
            % nc = new configuration
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,2*fds);
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
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,2*fds-1);
            %% second possible displacement
            rc = [-1,1];
            % rc = relative coordinates
            xn = 1+mod(x+rc(2)-1,n);
            yn = 1+mod(y+rc(1)-1,n);
            site = n*(yn-1)+xn;
            nc = 6;
            % nc = new configuration
            [States2,eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,2*fds);
    end
end
time2 = toc

States2 = States2(2:end,:);