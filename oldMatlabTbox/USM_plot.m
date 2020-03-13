function [y,more]=USM_plot(y,opt,mis)

%USM_PLOT collection of plot functions for usm structures
%Syntax: y=USM_plot(y,opt)
%Index:
%  0 - special variation on plot to make it look like a bar line
%  1 - density of of each dimension, for dna sequences
%  2 - density for each USM dimension
%

if ischar(opt);opt_.case=opt;opt=opt_;end  % in case the option case is being submitted as a character input
if isnumeric(opt);opt_.case=opt;opt=opt_;end  % in case the option case is being submitted as a numeric input
%if ischar(y) %then usm was not calculated yet
%    y=USM_main(y,'usm');
%end

switch opt.case
    case '0'
        
    case '1'
        if ~isfield(y.usm,'K') % then calculate Kernel first with default parameters
            y=USM_main(y,'K');
        end
        [n,m]=size(y.usm(1).K);
        if n~=4;error('this option is only for 2D USM such as compact representation of DNA');end
        L=length(y.usm);
        j=0;
        for i=1:L% for each sequence
            %forward
            U(:,:,1)=repmat(y.usm(i).K(1,:)',1,m);
            U(:,:,2)=repmat(y.usm(i).K(2,:),m,1);
            Uf=min(U,[],3);
            subplot(3,L,i);
            %j=j+1;subplot(3,L,j);
            imagesc(Uf);colorbar
            %bar3(Uf);
            title([y.seq(i).Header,', forward,\newline A=',num2str(sum(Uf(:)))])
            %backward
            U(:,:,1)=repmat(y.usm(i).K(3,:)',1,m);
            U(:,:,2)=repmat(y.usm(i).K(4,:),m,1);
            Ub=min(U,[],3);
            subplot(3,L,L+i);
            %j=j+1;subplot(3,L,j);
            imagesc(Ub);colorbar
            title([y.seq(i).Header,', backward,\newline A=',num2str(sum(Uf(:)))])
            %bidirectional
            U=Uf+Ub;
            subplot(3,L,2*L+i);
            %j=j+1;subplot(3,L,j);
            imagesc(U);colorbar
            title([y.seq(i).Header,', bidirectional,\newline A=',num2str(sum(U(:)))])
            y.usm(i).U=U;
        end
    case '2'
        L=length(y.usm);
        N=5;T=2;
        U=[1/2^(N+1):1/2^N:1]';
        [n,m]=size(y.usm(1).coord);
        LEG='{';for j=1:n/2;LEG=[LEG,'''',y.units(j),''','];end;LEG(end)='}';LEG=eval(LEG);
        for i=1:L
            figure
            [n,m]=size(y.usm(i).coord);  
            for j=1:n/2
                Hf(j,:)=USM_kheight(y.usm(i).coord(j,:)',N,T);
                Cj=USM_kcumdist(y.usm(i).coord(j,:)',N,T);
                C(:,j)=Cj(:,2);
            end
            subplot(2,1,1);
            Hf=Hf';
            %plot(U,Hf);
            USM_barplot(U,Hf);
            legend(LEG)
            hold on
            Hf=min(Hf,[],2);
            USM_barplot(U,Hf,'area');
            A(i)=sum(Hf).*(2*U(1));
            title(['Area:',num2str(A(i))])
            subplot(2,1,2);
            USM_barplot(U,C(2:end,:));
            legend(LEG)
            axis([0 1 0 1])
            clear Hf C 
            %hold on
            %area([0;U+U(1)],min(C,[],2));
        end
    case '3' %plots forward and backward USM 2D maps - this assumes a USM dimesnion 2
        
        % FORWARD MAP
        figure
        hold on        
        USM_grid(2.^(-[1:6]));
        plot(y.usm.coord(1,:),y.usm.coord(2,:),':','Color',0.5*ones(1,3))
        plot(y.usm.coord(1,1),y.usm.coord(2,1),'o','Color',[0 1 0],'MarkerSize',22)
        plot(y.usm.coord(1,end),y.usm.coord(2,end),'o','Color',[1 0 0],'MarkerSize',22)
        N=length(y.usm(1).coord(1,:));
        for i=N:-1:1
            %c=[(i/N)*3/4 (1-i/N)*3/4 0];
            c=[(i/N) (1-i/N) 0.3];
            text(y.usm(1).coord(1,i)+0.01,y.usm(1).coord(2,i),[y.seq.Sequence(1:i-1),'[',y.seq.Sequence(i),']'],'FontSize',8,'Color',c);
            plot(y.usm.coord(1,i),y.usm.coord(2,i),'.','Color',c,'MarkerSize',18)
        end
        axis([0 1 0 1])
        c=[0 0 0];
        H=text(0,0,[y.units(1),' '],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','top');
        
        H=text(1,0,[' ',y.units(3)],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','left');
        set(H,'VerticalAlignment','top');
        
        H=text(0,1,[y.units(2),' '],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','bottom');
        
        if length(y.units)>3;
        H=text(1,1,[' ',y.units(4)],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','left');
        set(H,'VerticalAlignment','bottom');
        end        
        hold off
        title('FORWARD MAP')
        
        
        % BACKWARD MAP
        figure
        hold on        
        USM_grid(2.^(-[1:6]));
        plot(y.usm.coord(3,:),y.usm.coord(4,:),':','Color',0.5*ones(1,3))
        plot(y.usm.coord(3,1),y.usm.coord(4,1),'o','Color',[0 1 0],'MarkerSize',22)
        plot(y.usm.coord(3,end),y.usm.coord(4,end),'o','Color',[1 0 0],'MarkerSize',22)
        N=length(y.usm(1).coord(3,:));
        for i=N:-1:1
            %c=[(i/N)*3/4 (1-i/N)*3/4 0];
            c=[(i/N) (1-i/N) 0.3];
            text(y.usm(1).coord(3,i)+0.01,y.usm(1).coord(4,i),['[',y.seq.Sequence(i),']',y.seq.Sequence(i+1:end)],'FontSize',8,'Color',c);
            plot(y.usm.coord(3,i),y.usm.coord(4,i),'.','Color',c,'MarkerSize',18)
        end
        axis([0 1 0 1])
        c=[0 0 0];
        H=text(0,0,[y.units(1),' '],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','top');
        
        H=text(1,0,[' ',y.units(3)],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','left');
        set(H,'VerticalAlignment','top');
        
        H=text(0,1,[y.units(2),' '],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','bottom');
        
        if length(y.units)>3;
        H=text(1,1,[' ',y.units(4)],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','left');
        set(H,'VerticalAlignment','bottom');
        end        
        hold off
        title('BACKWARD MAP')
        
    case '3b' %as 3 but with arrows and balck and white only
        if ischar(y) %then usm was not calculated yet
            y=USM_main(y,'usm');
        end
        % FORWARD MAP
        figure
        hold on        
        USM_grid(2.^(-[1:6]));
        %plot(y.usm.coord(1,:),y.usm.coord(2,:),':','Color',0.5*ones(1,3))
        plot(y.usm.coord(1,1),y.usm.coord(2,1),'o','Color',[0.5 0.5 0.5],'MarkerSize',22)
        %plot(y.usm.coord(1,end),y.usm.coord(2,end),'o','Color',[1 0 0],'MarkerSize',22)
        N=length(y.usm(1).coord(1,:));
        for i=1:N
            %c=[(i/N) (1-i/N) 0.3];
            c=[0 0 0];
            if i>1;
                s.xy0=y.usm.coord(1:2,i-1)'; % coordinates of beggining of arrow
                s.xyf=y.usm.coord(1:2,i)'; % coordinates of end of arrow
                s.k=0.02;% length of ortogonal projection of arrow wing into arrow stem
                s.L=0.01;% distance between wing and stem
                s.r=0.02;
                s.cor_da_linha=[0.5 0.5 0.5];
                s.LineStyle=':';
                seta2(s);
            end
            text(y.usm(1).coord(1,i)+0.01,y.usm(1).coord(2,i),['...',y.seq.Sequence(1:i-1),'[',y.seq.Sequence(i),']'],'FontSize',12,'FontName','Verdana','FontAngle','italic','Color',c);
            % defaults
            % s.LineWidth
            plot(y.usm.coord(1,i),y.usm.coord(2,i),'.','Color',c,'MarkerSize',18)
        end
        axis([0 1 0 1])
        c=[0 0 0];
        H=text(0,0,[y.units(1),' '],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','top');
        
        H=text(1,0,[' ',y.units(3)],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','left');
        set(H,'VerticalAlignment','top');
        
        H=text(0,1,[y.units(2),' '],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','bottom');
        
        if length(y.units)>3;
        H=text(1,1,[' ',y.units(4)],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','left');
        set(H,'VerticalAlignment','bottom');
        end        
        hold off
        title('FORWARD MAP')
        
        
        % BACKWARD MAP
        figure
        hold on        
        USM_grid(2.^(-[1:6]));
        %plot(y.usm.coord(3,:),y.usm.coord(4,:),':','Color',0.5*ones(1,3))
        plot(y.usm.coord(3,end),y.usm.coord(4,end),'o','Color',[0.5 0.5 0.5],'MarkerSize',22)
        %plot(y.usm.coord(3,end),y.usm.coord(4,end),'o','Color',[1 0 0],'MarkerSize',22)
        N=length(y.usm(1).coord(3,:));
        for i=N:-1:1
            %c=[(i/N) (1-i/N) 0.3];
            c=[0 0 0];
            if i<N;
                s.xy0=y.usm.coord(3:4,i+1)'; % coordinates of beggining of arrow
                s.xyf=y.usm.coord(3:4,i)'; % coordinates of end of arrow
                s.k=0.02;% length of ortogonal projection of arrow wing into arrow stem
                s.L=0.01;% distance between wing and stem
                s.r=0.02;
                s.cor_da_linha=[0.5 0.5 0.5];
                s.LineStyle=':';
                seta2(s);
            end
            text(y.usm(1).coord(3,i)+0.01,y.usm(1).coord(4,i),['[',y.seq.Sequence(i),']',y.seq.Sequence(i+1:end),'...'],'FontSize',12,'FontName','Verdana','FontAngle','italic','Color',c);
            plot(y.usm.coord(3,i),y.usm.coord(4,i),'.','Color',c,'MarkerSize',18)
        end
        axis([0 1 0 1])
        c=[0 0 0];
        H=text(0,0,[y.units(1),' '],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','top');
        
        H=text(1,0,[' ',y.units(3)],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','left');
        set(H,'VerticalAlignment','top');
        
        H=text(0,1,[y.units(2),' '],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','bottom');
        
        if length(y.units)>3;
        H=text(1,1,[' ',y.units(4)],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','left');
        set(H,'VerticalAlignment','bottom');
        end        
        hold off
        title('BACKWARD MAP')
        
    case '4' %plot 2D kernel
        %FORWARD
        plot3(y.usm(1).K.coord(:,1),y.usm(1).K.coord(:,2),y.usm(1).K.forward);
        %mesh(y.usm(1).K.coord(:,1),y.usm(1).K.coord(:,2),y.usm(1).K.forward);
        c=[0 0 0];
        H=text(0,0,0,[y.units(1),' '],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','top');
        
        H=text(1,0,0,[' ',y.units(3)],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','left');
        set(H,'VerticalAlignment','top');
        
        H=text(0,1,0,[y.units(2),' '],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','top');
        
        if length(y.units)>3;
        H=text(1,1,0,[' ',y.units(4)],'FontSize',18,'Color',c);
        set(H,'HorizontalAlignment','left');
        set(H,'VerticalAlignment','top');
        end 
case '5' %plot 2D kernel as 3D bar
        %FORWARD
        %plot3(y.usm(1).K.coord(:,1),y.usm(1).K.coord(:,2),y.usm(1).K.forward+y.usm(1).K.coord(:,2),y.usm(1).K.backward);
        s=1*(1/2^(y.usm.K_parms.N+1)); %step half side
        zz=y.usm(1).K.forward+y.usm(1).K.backward;
        xx=y.usm(1).K.coord(:,1);
        yy=y.usm(1).K.coord(:,2);
        xsquare=[-1 -1 1 1].*s;ysquare=[-1 1 1 -1].*s;zdata=[1 1 1 1];
        figure;
        hold on
        mm=max(zz)/63;c=colormap(1-gray);
        cc=interp1([1:64]',c,1+zz./mm);%disp(max(1+zz./mm))
        for i=1:length(zz)
            %TOP
            patch('XData',[-1 -1 1 1].*s+xx(i),'YData',[1 -1 -1 1].*s+yy(i),'ZData',[1 1 1 1].*zz(i),'FaceColor',cc(i,:));
            %SIDE 1
            patch('XData',[-1 -1 -1 -1].*s+xx(i),'YData',[-1 -1 1 1].*s+yy(i),'ZData',[0 1 1 0].*zz(i),'FaceColor',cc(i,:));
            %SIDE 2
            patch('XData',[1 1 1 1].*s+xx(i),'YData',[-1 -1 1 1].*s+yy(i),'ZData',[0 1 1 0].*zz(i),'FaceColor',cc(i,:));
            %SIDE 3
            patch('XData',[-1 -1 1 1].*s+xx(i),'YData',[-1 -1 -1 -1].*s+yy(i),'ZData',[0 1 1 0].*zz(i),'FaceColor',cc(i,:));
            %SIDE 3
            patch('XData',[-1 -1 1 1].*s+xx(i),'YData',[1 1 1 1].*s+yy(i),'ZData',[0 1 1 0].*zz(i),'FaceColor',cc(i,:));
            %BOTTOM
            patch('XData',[-1 -1 1 1].*s+xx(i),'YData',[-1 1 1 -1].*s+yy(i),'ZData',[0 0 0 0].*zz(i),'FaceColor',cc(i,:));
            
        end
        c=[0.5 0.5 0.5];
        H=text(0,0,0,[y.units(1),' '],'FontSize',12,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','bottom');
        
        H=text(1,0,0,[' ',y.units(3)],'FontSize',12,'Color',c);
        set(H,'HorizontalAlignment','center');
        set(H,'VerticalAlignment','top');
        
        H=text(0,1,0,[y.units(2),' '],'FontSize',12,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','bottom');
        
        if length(y.units)>3;
        H=text(1,1,0,[' ',y.units(4)],'FontSize',12,'Color',c);
        set(H,'HorizontalAlignment','left');
        set(H,'VerticalAlignment','bottom');
        end         
        
        hold off
        
        view(30,70);axis off
case '5b' %plot 2D kernel as 3D bar
        %FORWARD
        %plot3(y.usm(1).K.coord(:,1),y.usm(1).K.coord(:,2),y.usm(1).K.forward+y.usm(1).K.coord(:,2),y.usm(1).K.backward);
        s=1*(1/2^(y.usmK.parms.N+1));%s=1*(1/2^(y.usm.K_parms.N+1)); %step half side
        zz=y.usmK(1).forward+y.usmK(1).backward;
        xx=y.usmK(1).coord(:,1);
        yy=y.usmK(1).coord(:,2);
        xsquare=[-1 -1 1 1].*s;ysquare=[-1 1 1 -1].*s;zdata=[1 1 1 1];
        figure;
        hold on
        mm=max(zz)/63;c=colormap(1-gray);
        cc=interp1([1:64]',c,1+zz./mm);%disp(max(1+zz./mm))
        for i=1:length(zz)
            %TOP
            patch('XData',[-1 -1 1 1].*s+xx(i),'YData',[1 -1 -1 1].*s+yy(i),'ZData',[1 1 1 1].*zz(i),'FaceColor',cc(i,:));
            %SIDE 1
            patch('XData',[-1 -1 -1 -1].*s+xx(i),'YData',[-1 -1 1 1].*s+yy(i),'ZData',[0 1 1 0].*zz(i),'FaceColor',cc(i,:));
            %SIDE 2
            patch('XData',[1 1 1 1].*s+xx(i),'YData',[-1 -1 1 1].*s+yy(i),'ZData',[0 1 1 0].*zz(i),'FaceColor',cc(i,:));
            %SIDE 3
            patch('XData',[-1 -1 1 1].*s+xx(i),'YData',[-1 -1 -1 -1].*s+yy(i),'ZData',[0 1 1 0].*zz(i),'FaceColor',cc(i,:));
            %SIDE 3
            patch('XData',[-1 -1 1 1].*s+xx(i),'YData',[1 1 1 1].*s+yy(i),'ZData',[0 1 1 0].*zz(i),'FaceColor',cc(i,:));
            %BOTTOM
            patch('XData',[-1 -1 1 1].*s+xx(i),'YData',[-1 1 1 -1].*s+yy(i),'ZData',[0 0 0 0].*zz(i),'FaceColor',cc(i,:));
            
        end
        c=[0.5 0.5 0.5];
        H=text(0,0,0,[y.units(1),' '],'FontSize',12,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','bottom');
        
        H=text(1,0,0,[' ',y.units(3)],'FontSize',12,'Color',c);
        set(H,'HorizontalAlignment','center');
        set(H,'VerticalAlignment','top');
        
        H=text(0,1,0,[y.units(2),' '],'FontSize',12,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','bottom');
        
        if length(y.units)>3;
        H=text(1,1,0,[' ',y.units(4)],'FontSize',12,'Color',c);
        set(H,'HorizontalAlignment','left');
        set(H,'VerticalAlignment','bottom');
        end         
        
        hold off
        
        view(30,70);axis off

case '6' %plot 2D kernel as 3D bar
        %FORWARD
        %plot3(y.usm(1).K.coord(:,1),y.usm(1).K.coord(:,2),y.usm(1).K.forward+y.usm(1).K.coord(:,2),y.usm(1).K.backward);
        s=1*(1/2^(y.usm.K_parms.N+1)); %step half side
        zz=y.usm(1).K.forward;
        xx=y.usm(1).K.coord(:,1);
        yy=y.usm(1).K.coord(:,2);
        xsquare=[-1 -1 1 1].*s;ysquare=[-1 1 1 -1].*s;zdata=[1 1 1 1];
        figure;
        hold on
        mm=max(zz)/63;c=colormap(1-gray);
        cc=interp1([1:64]',c,1+zz./mm);%disp(max(1+zz./mm))
        for i=1:length(zz)
            %TOP
            patch('XData',[-1 -1 1 1].*s+xx(i),'YData',[1 -1 -1 1].*s+yy(i),'ZData',[1 1 1 1].*zz(i),'FaceColor',cc(i,:));
            %SIDE 1
            patch('XData',[-1 -1 -1 -1].*s+xx(i),'YData',[-1 -1 1 1].*s+yy(i),'ZData',[0 1 1 0].*zz(i),'FaceColor',cc(i,:));
            %SIDE 2
            patch('XData',[1 1 1 1].*s+xx(i),'YData',[-1 -1 1 1].*s+yy(i),'ZData',[0 1 1 0].*zz(i),'FaceColor',cc(i,:));
            %SIDE 3
            patch('XData',[-1 -1 1 1].*s+xx(i),'YData',[-1 -1 -1 -1].*s+yy(i),'ZData',[0 1 1 0].*zz(i),'FaceColor',cc(i,:));
            %SIDE 3
            patch('XData',[-1 -1 1 1].*s+xx(i),'YData',[1 1 1 1].*s+yy(i),'ZData',[0 1 1 0].*zz(i),'FaceColor',cc(i,:));
            %BOTTOM
            patch('XData',[-1 -1 1 1].*s+xx(i),'YData',[-1 1 1 -1].*s+yy(i),'ZData',[0 0 0 0].*zz(i),'FaceColor',cc(i,:));
            
        end
        c=[0.5 0.5 0.5];
        H=text(0,0,0,[y.units(1),' '],'FontSize',12,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','bottom');
        
        H=text(1,0,0,[' ',y.units(3)],'FontSize',12,'Color',c);
        set(H,'HorizontalAlignment','center');
        set(H,'VerticalAlignment','top');
        
        H=text(0,1,0,[y.units(2),' '],'FontSize',12,'Color',c);
        set(H,'HorizontalAlignment','right');
        set(H,'VerticalAlignment','bottom');
        
        if length(y.units)>3;
        H=text(1,1,0,[' ',y.units(4)],'FontSize',12,'Color',c);
        set(H,'HorizontalAlignment','left');
        set(H,'VerticalAlignment','bottom');
        end         
        
        hold off
        
        view(30,70);%axis off
    case '7' % plot forward and backward kernel intensity for a set of usm coordinates submited as mis
        %USM dimensions
        D=y.usmK.parms.D
        if nargin<3 %then use kernel sorce sequences as the coordinates
            Fxy=[];%colect forward USM coordinates here
            Bxy=[];%colect backward USM coordinates here
            for i=1:length(y.usm) %for all sequences
                Fxy=[Fxy,y.forward(i).coord];
                Bxy=[Bxy,y.backward(i).coord];
            end
        else
            %find out if the coordinates are only forward or both
            if length(mis(1,:))==2*D %both being submitted
                Fxy=mis(:,1:D)';
                Bxy=mis(:,D+1:end)';
            else %only forward is being submitted
                Fxy=mis(:,1:D)';
            end            
        end
        %calculate forward coordinates
        %this is the iportant part - each coordinate has to be replaced by
        %the right index in the Kernel
        
end
            