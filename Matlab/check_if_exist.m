function [position, Matrix, found] = check_if_exist(new_state,Matrix)

M = Matrix;
s = new_state;
found = 0;
pos = 0;
for i = 1:length(M(:,1))
    if s(1:end-3) == M(i,1:end-3)
        found = found+1;
        pos(found) = i;
    end
end
position = pos;

if pos == 0
    position = i+1;
    Matrix(position,:) = s;
end    
   
   
    