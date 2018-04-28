% this script is intended to be the implementation of the exact
% diagonalization for the effective Hamiltonian for our system
clc; clear all; close all;

% dimension of the lattice to study. Needs to be and even number, in order
% to preserve the periodic boundary conditions
n = 4;

% The lattice at the beginning will be 
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

possible_movements = [1,0; -1,0; 0,1; 0,-1];

number_of_states = 1;
tic
for i=1:n^2
    if ground_state(i) == 1
        label = 1;
        for displacement = possible_movements'
            s = f_swap(ground_state,i,displacement);
            States(number_of_states,:) = [s, label, i];
            label = label+1;
            number_of_states = number_of_states + 1;
        end
    end    
end
disp('time for creating the matrix of all first possible displacements')
toc
     
%% second movements states







