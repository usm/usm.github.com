function x=USM_grid(x)

%USM_GRID plots grid on USMap
%Syntax: x=USM_grid(x)
%Description: x is a vector of resolutions, for example 2.^(-[1:4])
%
%Jonas Almeida, almeidaj@musc.edu, Feb 2005

for i=length(x):-1:1
    xx=[0:x(i):1];
    for j=1:length(xx)
        %horizontal lines
        s=(j-1)*x(i);%hold on
        c=1-ones(1,3).*x(i)-0.1;
        plot([0,1],[s,s],'Color',c);
        plot([s,s],[0,1],'Color',c);
    end
end

