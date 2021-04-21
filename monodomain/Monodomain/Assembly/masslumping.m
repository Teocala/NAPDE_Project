function [ML] = masslumping(M)
%MASSLUMPING
% This fiùunction performs the mass lumping of M

n=size(M,1);
ML=zeros(n,n);
for ii=1:n
    ML(ii,ii)=sum(M(ii,:));
end

