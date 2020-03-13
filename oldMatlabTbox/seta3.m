function s=seta3(s)

%SETA3 draws 3D arrows
%
%Syntax: s=seta3(s)
% s.xyz0: coordinates of beggining of arrow
% s.xyzf: coordinates of end of arrow
% s.k: length of ortogonal projection of arrow wing into arrow stem
% s.L: distance between wing and stem

% defaults
% s.LineWidth
% s.r % exclusion radius

warning off
if ~isfield(s,'LineWidth');s.LineWidth=1;end
if ~isfield(s,'r');s.r=0;end
if ~isfield(s,'L');s.L=sqrt(s.LineWidth)*s.k/4;end

dx=s.xyzf(1)-s.xyz0(1);
dy=s.xyzf(2)-s.xyz0(2);
dz=s.xyzf(3)-s.xyz0(3);
h=sqrt(dx^2+dy^2+dz^2);


if s.r>0   % exclusion radius
    s2=s;
    dd=[dx,dy,dz].*((s.r)./h);
    s2.xyz0=s.xyz0+dd;
    s2.xyzf=s.xyzf-dd;
    s2.r=0;
    s2=seta3(s2);
else
    aayx=atan(dy/dx);if dx==0;aayx=sign(dy)*pi/2;end
    aayx=aayx*(dx~=0)+pi*(dy>0)*(dx==0);
    aayx=aayx*(dy~=0)+((pi/2)*(dx>0)-(pi/2)*(dx<0))*(dy==0);
    kxyz=s.xyz0+[dx,dy,dz].*((h-s.k)./h);    
    Lxyz_esq(3)=kxyz(3);Lxyz_dto(3)=kxyz(3);
    Lxyz_esq(1:2)=kxyz(1:2)+[cos(aayx-pi/2),sin(aayx-pi/2)].*s.L;
    Lxyz_dto(1:2)=kxyz(1:2)+[cos(aayx+pi/2),sin(aayx+pi/2)].*s.L;
    if h>s.k
        plot3([s.xyz0(1);kxyz(1)],[s.xyz0(2);kxyz(2)],[s.xyz0(3);kxyz(3)],'LineWidth',s.LineWidth);
        hold on
    end
    %plot3(kxyz(1),kxyz(2),kxyz(3),'o','MarkerFaceColor','k')
    %plot3(s.xyzf(1),s.xyzf(2),s.xyzf(3),'o')
    %plot3(Lxyz_esq(1),Lxyz_esq(2),Lxyz_esq(3),'o')
    %plot3(Lxyz_dto(1),Lxyz_dto(2),Lxyz_dto(3),'o')
    patch([Lxyz_esq(1);Lxyz_dto(1);s.xyzf(1)],[Lxyz_esq(2);Lxyz_dto(2);s.xyzf(2)],[Lxyz_esq(3);Lxyz_dto(3);s.xyzf(3)],'b','LineStyle','none');    
end
warning on
