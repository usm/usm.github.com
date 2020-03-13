function C=USM_kcumdist(X,N,T)

%USM_CDIST cumulative distrbution for a USM dimesntion
%Syntax: C=USM_KCUMDIST(X,N,T)
%Description:
% calculates the cumulative distribution of the USM kernel for a collection of
% USM coordinates X, an order N and a scaling T.
% This function uses USM_KHEIGHT to determine the density distribution
%
%Jonas Almeida, almeidaj@musc.edu, Nov 2004

[H,U]=USM_kheight(X,N,T);
n=length(H)+1;
C=[[0 0];[U+U(1),H]];
S=C(2,1);%integration step
for i=2:n
    C(i,2)=C(i-1,2)+C(i,2)*S;
end
