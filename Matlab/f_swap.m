function final_state = f_swap(state, index, displacement)
% This function take a generic state of zeros and ones and create a new
% state with 1 in position (index+displacement) and 0 in position (index)
%
% INPUT
% displacement : 2x1 vector
%
% OUTPUT
%
d = displacement;

% foolproof check
if state(index) == 0
    error('f_swap:the desired state is empty');
end

n = length(state);

final_state = state;
final_state(index) = 0;
% in the following 1+mod(-1,...) is for managing properly the zeros in the
% arrays

dims = sqrt(n);

x = 1+mod(index-1,dims);
y = ceil(index/dims);
M = reshape(final_state,[dims,dims]);
if M(1+mod(x+d(2)-1,dims),1+mod(y+d(1)-1,dims)) == 1
    error('f_swap:found 1 but 0 expected');
else 
    M(1+mod(x+d(2)-1,dims),1+mod(y+d(1)-1,dims)) = 1;
end

final_state = reshape(M,[n,1])';
