function [z] = mutation(A,n);
% Mutation within predefined hypercube set forth by the boundaries of A
minn = min(A(:,1:n)); maxn = max(A(:,1:n));

for qq = 1:n,
    z(1,qq) = minn(qq) + rand*(maxn(qq)-minn(qq));
end