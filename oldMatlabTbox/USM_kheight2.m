function [H,U]=USM_kheight2(X,N,T,U)

%USM_KHEIGHT2 height of 2D USM kernel
%Syntax: [H,U]=USM_kheight(X,N,T,U)
%Description:
% calculates the height of the USM kernel at position U for a collection of
% 2D USM coordinates X, an order N and a scaling T. If U is not provided it
% will be calculated at 1/2^N intervals, for subsequent canculation by
% interpolations. D is the # dimensions of the USM map
%
%Jonas Almeida, almeidaj@musc.edu, Nov 2004

D=2; %<--- 2D USM expected, generalize this later

if nargin<4 % if no positions are provided calculate the all surface at 1/2^N intervals
    u=[1/2^(N+1):1/2^N:1]; %calculation at 1/2^N intervals
    nu=length(u);
    U=[u(ceil([1/nu:1/nu:nu]))',u(repmat([1:nu],1,nu))'];
end

[n,m]=size(X);if m==1;X=X';[n,m]=size(X);end

%n=length(X);
LB=zeros(n,N+1); %Lower boundaries for each point at each segment length
UB=zeros(n,N+1); %Upper boundaries for each point at each segment length
for j=1:D
    for i=0:N
        K=floor(X(:,j).*2^i);
        LB(:,i+1,j)=K./(2^i);
        UB(:,i+1,j)=(K+1)./(2^i);
    end
end

%Distribution of heights
Hi=(((2^D)*T).^([0:N]))./sum(T.^([0:N]));

%Hz=zeros(n,N+1);

%[n,m]=size(U);
%if n==1;U=U';[n,m]=size(U);end
m=length(U(:,1));

H=zeros(m,1); 

for i=1:m %calculate one height for each set of coordinates
    for j=1:D %for each dimension
        Hz(:,:,j)=repmat(Hi,n,1).*(LB(:,:,j)<U(i,j)).*(UB(:,:,j)>U(i,j));
    end
    Hz=min(Hz,[],3);
    H(i)=sum(Hz(:));
end
H=H./n; % <--- this is where the distribution is normalized for the area to be 1