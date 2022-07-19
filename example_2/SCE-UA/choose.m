function [B,L] = choose(rho,A,n,q);				
% Selects q elements from A using trapezoidal probability distribution

rho = cumsum(rho); Selected = []; counter = 1;
while length(Selected)<q,
   U = rand;  				% Draw random number U between 0 and 1 using a uniform distribution
   R = find(U<rho); R = R(1,1);  % Combine labelled U with trapezoidal probability 	 
     Test = mod(find(Selected==R),0);
        while length(Test) >0,
           U = rand;  				% Draw random number U between 0 and 1 using a uniform distribution
           R = find(U<rho); R = R(1,1);   % Combine labelled U with trapezoidal probability 	 
           Test = mod(find(Selected==R),0);
        end;
      Selected = [Selected;R];
      B(counter,1:n+1) = A(R,:,:);
      L(counter,1:2) = [R A(R,n+1,:)];
      counter = counter+1;
end;

   
   
   

