function [y,opt]=USM_main(x,opt)

%USM_MAIN Main function of the USM toolbox (Universal Sequence Maps)
%Syntax: [y,opt]=USM_main(x,opt)
%Description:
% the nature of x and the options set what y is. This is the main function
% of the toolbox, which gets to be used to do almost everything the toolbox supports.
%
% The data structure of y is
% .USMseq  -> symbolic sequence
% .usmf    -> forward USM coordinates
% .usmb    -> backward USM coordinates
%
% the structure of opt (the options variable) is
% .case   -> option for the switch structure
% .mis    -> miscelaneous options / data
%
%Jonas ALmeida, almeidaj@musc.edu, Nov 2004


if nargin<1;   % if opt is not provided
    opt.case='read sequence';    % then create it, it is assumed that x is the sequence, as a character or a character cell array
end

if ischar(opt);opt_.case=opt;opt=opt_;end  % in case the option case is being submitted as a character input

switch opt.case
    case 'read sequence'
        if ischar(x)
            error('not coded')
        elseif iscell(x)
            error('not coded')
        else
            error('not coded')
        end
    case 'read fasta'
        y.filename=[x,', ',datestr(now)];
        y.seq=fastaread(x);
        y.units=unique(strcat(y.seq.Sequence)); % list unique units
    case 'usm' %calculate forward and backward USM coordinates
        if ischar(x) %then a fasta file name is being submitted
            y=USM_main(x,'read fasta');
        else
            y=x;
        end
        % For each sequence
        n=length(y.units);
        for i=1:length(y.seq)
            % First define matrix of hits
            m=length(y.seq(i).Sequence);
            MH=zeros(n,m);
            % Look for each possible symbol at each position in the sequence
            for j=1:n
                MH(j,:)=(y.seq(i).Sequence==y.units(j));
            end
            % Compacts hit matrix
            MH=USM_compact(MH);   %<----- CONTROL COMPACTATION HERE !
            % Calculate forward coordinates;
            %y.forward(i).MH=MH;
            forward(i).coord=USM_CGR(MH,'loop');
            % Calculate forward coordinates;
            MH=MH(:,end:-1:1); % inversed hit matrix
            %y.backward(i).MH=MH;
            backward(i).coord=USM_CGR(MH,'loop'); % do forward CGR with reversed hit matrix
            backward(i).coord=backward(i).coord(:,end:-1:1); % invert the corrdinates too, and this is backward CGR
            y.usm(i).coord=[forward(i).coord;backward(i).coord];
        end
    case 'usm05' %calculate forward and backward USM coordinates
        if ischar(x) %then a fasta file name is being submitted
            y=USM_main(x,'read fasta');
        else
            y=x;
        end
        % For each sequence
        n=length(y.units);
        for i=1:length(y.seq)
            % First define matrix of hits
            m=length(y.seq(i).Sequence);
            MH=zeros(n,m);
            % Look for each possible symbol at each position in the sequence
            for j=1:n
                MH(j,:)=(y.seq(i).Sequence==y.units(j));
            end
            % Compacts hit matrix
            MH=USM_compact(MH);   %<----- CONTROL COMPACTATION HERE !
            % Calculate forward coordinates;
            %y.forward(i).MH=MH;
            y.forward(i).coord=USM_CGR(MH,0.5);
            % Calculate forward coordinates;
            MH=MH(:,end:-1:1); % inversed hit matrix
            %y.backward(i).MH=MH;
            y.backward(i).coord=USM_CGR(MH,0.5); % do forward CGR with reversed hit matrix
            y.backward(i).coord=y.backward(i).coord(:,end:-1:1); % invert the corrdinates too, and this is backward CGR
            y.usm(i).coord=[y.forward(i).coord;y.backward(i).coord];
        end
    case 'K' %Kernel Density distribution
        if ischar(x) %then do usm first
            y=USM_main(x,'usm');
        else
            y=x;
        end
        if ~isfield(opt,'parms')
            opt.parms.N=5;
            opt.parms.T=2;
        end
        Fxy=[];%colect forward USM coordinates here
        Bxy=[];%colect backward USM coordinates here
        [n,m]=size(y.usm(1).coord);
        for i=1:length(y.usm) %for all sequences
            Fxy=[Fxy,y.usm(i).coord(1:n/2,:)];
            Bxy=[Bxy,y.usm(i).coord(n/2+1:end,:)];
        end
        
        y.usmK.parms=opt.parms;
        y.usmK.parms.D=length(y.usm(1).coord(:,1))/2;
        [y.usmK.forward,y.usmK.coord]=USM_kheight(Fxy',opt.parms.N,opt.parms.T);
        y.usmK.backward=USM_kheight(Bxy',opt.parms.N,opt.parms.T);    
end

