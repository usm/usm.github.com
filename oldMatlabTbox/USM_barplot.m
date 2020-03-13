function [x,y]=USM_barplot(x,y,opt)

%USM_BARPLOT plots densities and cumulative distributions as "bar lines"
%Syntax [x,y]=USM_barplot(x,y)

if nargin<3;opt='plot';end
if nargin<2;
    y=x;x=[1:length(y(1,:))]./length(y(1,:))-(0.5/length(y(1,:)));
end
[n,m]=size(y);
if m>n;y=y';x=x';[n,m]=size(y);end
[nx,mx]=size(x);
if mx==1;x=repmat(x,1,m);end

yy=[];xx=[];
for i=1:m
    % y
    z=[y(:,i),y(:,i)]';
    z=z(:);
    z=[0;z];
    yy=[yy,z];
    % x
    S=x(1,i);
    z=x(:,i)-S;z=[z,z]';
    z=z(:);
    z=[z;z(end)+2*S];
    xx=[xx,z];
end
x=xx;y=yy;

eval([opt,'(x,y)']);
