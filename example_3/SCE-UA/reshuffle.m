function [D,x] = reshuffle(A,p,n,m);  
% Reshuffles points back into original population
counter = 1;
for i=1:p,
   for q=1:m
      x(counter,1:n) = A(q,1:n,i);
      D(counter,1:2) = [A(q,n+1,i) counter];
      counter = counter + 1;
   end;
end;
D = -sortrows(-D,[1]);
