function [States2, eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,fds)
% This function is intended to modify the effective Hamiltonian
% 
% A significative part of the code is dedicated to the implementation of 
% the Slater determinant
%
% Idea:
% [ 1 0 1 0 0 1 0 1 1 0 1 0 0 1 0 1]
% is described by 
% c1 c3 c6 c8 c9 c11 c14 c16
% and then if we want annihilate 6 and create in 10 we have to skip 
% initially 3, ie, (-1)^3 and then 4 sites, ie, (-1)^4 
% for obtaining
% [ 1 0 1 0 0 0 0 1 1 1 1 0 0 1 0 1]
% 
%
% first part of the Slater determinant
num = sum(s(1:site-1));
dims = sqrt(length(s));
x = 1+mod(site-1,dims);
y = ceil(site/dims);
xn = 1+mod(x+d(2)-1,dims);
yn = 1+mod(y+d(1)-1,dims);
new_site = (yn-1)*dims+xn;

% create the new state
r = f_swap(s,site,d);   

num = num+sum(r(1:new_site-1)); % second part os Slater determinant

%% add t ot -t
pos = 0;
flag = 0;
while flag == 0 && pos<length(States2(:,1))
    pos = pos+1;
    if r == States2(pos,1:end-2)
        eff_H(fds,nSt+pos) = t*(-1)^num;
        flag = 1;
    end
end

if flag == 0
    States2(end+1,:) = [r,new_site,nc];
    eff_H(fds,nSt+pos+1) = t*(-1)^num;
end