function [H,U]=USM_kdensity(X,N,T,U)

%USM_KHEIGHT density of USM kernel
%Syntax: [H,U]=USM_kdensity(X,N,T,U)
%Description:
% calculates the height of the USM kernel at position U for a coolection of
% USM coordinates X, an order N and a scaling T. If U is not provided it
% will be calculated at 1/2^N intervals, for subsequent canculation by
% interpolations.
%
%Jonas Almeida, almeidaj@musc.edu, Nov 2004

if nargin<4 %calculation at 1/2^N intervals
    U=[1/2^(N+1):1/2^N:1];
end
[n,m]=size(X);if m==1;X=X';end
[n,m]=size(U);if n==1;U=U';end

n=length(X);
LB=zeros(n,N+1); %Lower boundaries for each point at each segment length
UB=zeros(n,N+1); %Upper boundaries for each point at each segment length
for i=0:N
    K=floor(X.*2^i);
    LB(:,i+1)=K./(2^i);
    UB(:,i+1)=(K+1)./(2^i);
end

%H.LB=LB;
%H.UB=UB;

Hi=((2*T).^([0:N]))./sum(T.^([0:N]));

m=length(U);
Hz=zeros(n,N+1);
H=U;
for i=1:m
    Hz=repmat(Hi,n,1).*(LB<U(i)).*(UB>U(i));
    H(i)=sum(Hz(:));
end
H=H./n;