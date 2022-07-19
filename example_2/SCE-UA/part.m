function [A] = part(D,p,x,n);
% Partition D into P complexes
counter = 1;			
for i=1:p:length(D),  											
   for j=1:p,														
      A(counter,1:n+1,j) = [x(D(i,2),1:n) D(i,1)];      		
      i = i+1; 												
   end;															
   counter = counter + 1;									
end;																
