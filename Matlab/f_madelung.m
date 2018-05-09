function [energy1, energy2] = f_madelung(n,alpha)
% INPUT 
% 
% OUTPUT
% 
% NB
% 
% 
% Check inputs
if nargin < 2
  alpha = 1;
end

A = zeros(n);
for i = 1:n
    for j = 1:n
        if mod(i+j,2)==0
            A(i,j) = 1;
        end
    end
end

pot = 0;
range = -n+1:n-2;
% ref = [1,1];
for i=range
    for j=range
        if A(1+mod(i+1-1,n),1+mod(j+1-1,n)) == 1
            if ~(i==0 & j==0)
                pot = pot+1/norm([i,j])^(alpha);
            end
        end
    end
end

A(1,1) = 0;
A(2,1) = 1;
pot1 = 0;
ref = [2,1];

for i=range
    for j=range
        if A(1+mod(i+2-1,n),1+mod(j+1-1,n)) == 1
            if ~(i==0 & j==0)
                pot1 = pot1+1/norm([i,j])^(alpha);
            end
        end
    end
end
energy1 = pot1-pot;

diff1_f = 1/norm([0,1])^alpha;
diff1_r = 1/norm([1,1])^alpha;
diff2_f = 1/norm([1,0])^alpha;
diff2_r = 1/norm([1,1])^alpha;

energy2 = 2*energy1-diff1_f+diff1_r-diff2_f+diff2_r;