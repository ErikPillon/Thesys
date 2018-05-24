% we try once more to develop the algorithm for the Exact Diagonalization;
% we will make use of the function *switch*, *Madelung* and *Slater*
%
% the only parameters to fix are 
% n : cluster size
% t : hopping parameter
% v : fugacity parameter
%

%clc; clear all; close all;

n = 30; 

% t = 0.05;
% v = 0.05;

% *n* dimension of the lattice to study. In order to preserve the Periodic 
% Boundary Conditions, n must be even
if mod(n,2) == 1
    error(message('the cluster size *n* must be even'));
end

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

States = zeros(1,n^2+3);

%% First movements state
% the matrix will be made by vectors of the following form
%
% [         (n*n)        ,  [x,y]    ,  l  ]
%             ^               ^         ^
%             |               |         |
%       state described    position   label
%       by 0's and 1's     modified
% 
% label    : gives the configuration of the state we are dealing (1 to 6)
% position : gives the position of the particle created wrt the ground
%            state

s = ground_state;

% possible movements
moves = [-1,0; 0,1; 1,0; 0,-1];
nSt = 0;

for i=1:n^2
    if s(i) == 1
        x = 1+mod(i-1,n);
        y = ceil(i/n);
        conf = 0; % type of configuration
        for delta = moves'
            nSt = nSt +1;
            conf = conf + 1;
            States(nSt,:) = f_switch([x,y],delta,s,conf);
        end
    end
end

disp('First movement States completed')

% initialize to zero the Hamiltonian
Hamiltonian = zeros(2);

disp('Started the second movement States Matrix')
Seconds = zeros(1,n^2+3);
d = 0;

for i = 1:nSt
    pos = [States(i,end-2),States(i,end-1)];
    switch States(i,end)
        case 1
            d = [1,0];
            site1 = 1+mod(pos+[0,-1]-[1,1],n);
            site2 = 1+mod(pos+[0,1]-[1,1],n);
            dd = [1,1];
            site3 = 1+mod(pos+[0,0]-[1,1],n);
            site4 = 1+mod(pos+[0,-1]-[1,1],n);
        case 2
            d = [0,-1];
            site1 = 1+mod(pos+[-1,0]-[1,1],n);
            site2 = 1+mod(pos+[1,0]-[1,1],n);
            dd = -[1,1];
            site3 = 1+mod(pos+[0,0]-[1,1],n);
            site4 = 1+mod(pos+[1,0]-[1,1],n);
        case 3
            d = [-1,0];
            site1 = 1+mod(pos+[0,1]-[1,1],n);
            site2 = 1+mod(pos+[0,-1]-[1,1],n);
            dd = -[1,1];
            site3 = 1+mod(pos+[0,0]-[1,1],n);
            site4 = 1+mod(pos+[0,1]-[1,1],n);
        case 4
            d = [0,1];
            site1 = 1+mod(pos+[1,0]-[1,1],n);
            site2 = 1+mod(pos+[-1,0]-[1,1],n);
            dd = [1,1];
            site3 = 1+mod(pos+[0,0]-[1,1],n);
            site4 = 1+mod(pos+[-1,0]-[1,1],n);
    end
    % first
    r = f_switch(site1,d,States(i,1:end-3),9);
    [wtf,Seconds] = check_if_exist(r,Seconds);
    sign = f_Slater(States(i,1:end-3),site1,r,d);
    Hamiltonian(i,nSt+wtf) = -sign*t;
    % second
    r = f_switch(site2,d,States(i,1:end-3),9);
    [wtf,Seconds] = check_if_exist(r,Seconds);
    sign = f_Slater(States(i,1:end-3),site1,r,d);
    Hamiltonian(i,nSt+wtf) = -sign*t;
    % third
    r = f_switch(site3,dd,States(i,1:end-3),9);
    [wtf] = check_if_exist(r,States);
    sign = f_Slater(States(i,1:end-3),site3,r,dd);
    Hamiltonian(i,wtf) = -sign*v;
    % fourth
    r = f_switch(site4,dd,States(i,1:end-3),9);
    [wtf] = check_if_exist(r,States);
    sign = f_Slater(States(i,1:end-3),site4,r,dd);
    Hamiltonian(i,wtf) = -sign*v;
end
disp('Finished the second movements matrix states')

Hamiltonian = Hamiltonian(:,[1:nSt,nSt+2:end]);
disp('Finished fixing the Hamiltonian')
% making the effective Hamiltonian symmetric
for i = 1:length(Hamiltonian(:,1))
    for j = nSt+1:length(Hamiltonian(1,:))
        Hamiltonian(j,i) = Hamiltonian(i,j);
    end
end

% potential part of the Hamiltonian
[E1,E2] = f_madelung(n); 
potential = diag([E1*ones(1,nSt),E2*ones(1,length(Seconds(:,1))-1)]);
Hamiltonian = Hamiltonian + potential;

[~,D] = eig(Hamiltonian);
histogram(diag(D),50)