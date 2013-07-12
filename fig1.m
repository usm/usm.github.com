function fig1(n)
% figure 1 - Sierpinski's triangle + CGR plane, filled randomly
% n is the number of points, for example fig1(1000)

if nargin==0;n=1000;end

% Spiersinski
subplot(1,2,1)
h = sqrt(1-(1/2)^2);
IM({[0,0],[1/2,h],[1,0]},n,[0,0]);
text(0,0,'A ','HorizontalAlignment','right','VerticalAlignment','top');
text(1/2,h,'B ','HorizontalAlignment','right','VerticalAlignment','bottom');
text(1,0,'  C','HorizontalAlignment','left','VerticalAlignment','top');

% CGR
subplot(1,2,2)
IM({[0,0],[0,1],[1,1],[1,0]},n);
text(0,0,'A ','HorizontalAlignment','right','VerticalAlignment','top');
text(0,1,'C ','HorizontalAlignment','right','VerticalAlignment','bottom');
text(1,1,'  G','HorizontalAlignment','left','VerticalAlignment','bottom');
text(1,0,'   T','HorizontalAlignment','left','VerticalAlignment','top');

% iterated map plot function 
function n = IM(edg,n,y)
% Iterated map of n points in a the space defined by an edges cell array
% edg is s cell array of vector of [x,y] coordinates

if nargin<3;y=[1/2,1/2];end
axis off
axis square
hold on
m = length(edg);
disp(['filling ',num2str(n),' points in a polygon with ',num2str(m),' edges'])
% plot edges
for i=1:m
    plot(edg{i}(1),edg{i}(2),'ko','MarkerFaceColor','k')
end
% iterate
for i=1:n
    x=edg{ceil(rand()*m)}; % pick one edge randomly
    y=y+0.5*(x-y); % <-- the Iterated Map !
    plot(y(1),y(2),'k.','MarkerSize',3);
end





