function [T] = checkpars(NewPars,minn,maxn);
% Function checks whether NwePars are within bounds or not
% T = 0 => parameter set is feasible,
% T = 1 => parameter set is infeasible
T = 1; [m,n] = size(NewPars);
for i = 1:n,
   if NewPars(1,i) < minn(1,i),
      T = -1;
   end;
   if NewPars(1,i) > maxn(1,i),
      T = -1;
   end;
end;
    