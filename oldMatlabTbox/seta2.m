function s=seta2(s)

%SETA2 draws 2D arrows
%
%Syntax: s=seta2(s)
% s.xy0: coordinates of beggining of arrow
% s.xyf: coordinates of end of arrow
% s.k: length of ortogonal projection of arrow wing into arrow stem
% s.L: distance between wing and stem

% defaults
% s.LineWidth
% s.r % exclusion radius

warning off
if ~isfield(s,'LineWidth');s.LineWidth=1;end
if ~isfield(s,'r');s.r=0;end
if ~isfield(s,'L');s.L=sqrt(s.LineWidth)*s.k/4;end
if ~isfield(s,'cor_da_linha');s.cor_da_linha='k';end
if ~isfield(s,'LineStyle');s.LineStyle='-';end

dx=s.xyf(1)-s.xy0(1);
dy=s.xyf(2)-s.xy0(2);
h=sqrt(dx^2+dy^2);

if s.r>0   % exclusion radius
    s2=s;
    dd=[dx,dy].*((s.r)./h);
    s2.xy0=s.xy0+dd;
    s2.xyf=s.xyf-dd;
    s2.r=0;
    s2=seta2(s2);
else
    aayx=atan(dy/dx);%if dx==0;aayx=sign(dy)*pi/2;end
    aayx=aayx*(dx~=0)+pi*(dy>0)*(dx==0);
    aayx=aayx*(dy~=0)+((pi/2)*(dx>0)-(pi/2)*(dx<0))*(dy==0);
    kxy=s.xy0+[dx,dy].*((h-s.k)./h);    
    Lxy_esq(1:2)=kxy(1:2)+[cos(aayx-pi/2),sin(aayx-pi/2)].*s.L;    
    Lxy_dto(1:2)=kxy(1:2)+[cos(aayx+pi/2),sin(aayx+pi/2)].*s.L;
    if s.k<h
        plot([s.xy0(1);kxy(1)],[s.xy0(2);kxy(2)],'Color',s.cor_da_linha,'LineWidth',s.LineWidth,'LineStyle',s.LineStyle);
        hold on
    end
    patch([Lxy_esq(1);Lxy_dto(1);s.xyf(1)],[Lxy_esq(2);Lxy_dto(2);s.xyf(2)],s.cor_da_linha,'LineStyle','none');
    %plot(kxy(1),kxy(2),'o','MarkerFaceColor','k')
    %plot(s.xyf(1),s.xyf(2),'o')
    %plot(Lxy_esq(1),Lxy_esq(2),'o')
    %plot(Lxy_dto(1),Lxy_dto(2),'o')
end
warning on
