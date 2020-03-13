function [y,opt]=USM_study(opt,x)

%USM_STUDY collection of studies on USM properties
%Syntax: [y,opt]=USM_study(x,opt)
%Description:
% this collection of studies is indexed to a switch stricture on opt
% the Index:
% 'random dna' - generates radom dna sequences


switch opt
    case 'random dna' %generates radom dna sequences
        units='AGCT';
        fid=fopen('study.fas','w');
        n=x(1);m=x(2); %x is a vector with two elements, the number of sequences and their length
        for i=1:n
            fprintf(fid,'%s\n',['>sequence #',num2str(i)]);
            fprintf(fid,'%s\n',units(ceil(rand(1,m)*4)));
        end
        fclose(fid);
        %also generates usm variable
        y=USM_main('study.fas','K');
    case 'random abc' %generates radom dna sequences
        units='abc';LU=length(units);
        fid=fopen('study_abc.fas','w');
        n=x(1);m=x(2); %x is a vector with two elements, the number of sequences and their length
        for i=1:n
            fprintf(fid,'%s\n',['>sequence #',num2str(i)]);
            fprintf(fid,'%s\n',units(ceil(rand(1,m)*LU)));
        end
        fclose(fid);
        %also generates usm variable
        y=USM_main('study_abc.fas','usm');
    case 'unidens' %computing densities for forward and backward coordinates separately
        %e.g. USM_study('random dna',[1,100])
        if nargin<2;x=[1,20];end
        y=USM_study('random dna',x);% compute usm structure
        % ... incomplete
        
    case 'fig1' %figure 1 for Bioinformatics Kernel paper
        y=USM_main('study.fas','K'); %read sequences from study.fas
        z=USM_plot(y,'3b'); %plos forward and backward map positions
    case 'fig2' %figure 2 for Bioinformatics Kernel paper
        y=USM_main('ABABAAA.fas','usm'); %read sequences from ABABAAA.fas 
        % and now calculates Kernel density
        K=USM_kheight(y.usm.coord(1,:),3,1);
        n=length(K);s=1/n;
        KK=[K';K'];XX=[0:s:1-s;s:s:1];
        area(XX(:),KK(:))
        G=get(gcf);G1=get(G.Children(1));
        set(G1.Children,'FaceColor',[0 0 0]);
        hold on
        ax=axis;
        for i=1:length(y.seq.Sequence)
            plot(y.usm.coord(1,i),ax(4),'vk','MarkerSize',5,'MarkerFaceColor','k')
            text(y.usm.coord(1,i),ax(4)+1,num2str(i),'HorizontalAlignment','center')
        end
        hold off
    case 'fig3'
        opt.parms.N=5;opt.case='K';opt.parms.T=0.5;y=USM_plot(USM_main('study.fas',opt),'5b');
end
            
        