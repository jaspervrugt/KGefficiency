function [A] = replace(A,B,L,n);	
% Replace parents by offspring,
for w=1:length(L),
   A(L(w,1),1:n+1,:) = B(w,1:n+1); 
end;
Temp = A(:,:,:);
A = -sortrows(-Temp,[n+1]);