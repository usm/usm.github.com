function M=USM_CGR(M,opt)

%USM_CGR Chaos Game Representation (CGR) iteration for a Hit matrix
% Syntax: y=USM_CGR(M,opt)
% Description: using a binary hit matrix, M, this function iterates the
% forward CGR coordinates. By doing M=M(:,end:-1:1) and reversing the y
% too: y=y(:,end:-1:1), 
% the options, opt, allows the use choice on seeding:
%    opt='random'
%    opt='loop'     (default choice, uses the last coordinates as seed)
%    opt=0.5        (or anyother value, it assigns that value to all
%                    coordinates.A seed vector instead of a seed scalar can
%                    also be used)
%
%Jonas Almeida, almeidaj@musc.edu, Nov 2004


if nargin<2;opt='loop';end

M=double(M);
[n,m]=size(M);%y=zeros(n,m);

if isnumeric(opt)
    CGRo=opt;
    if length(CGRo)==1;CGRo=ones(n,1).*CGRo;end
    
elseif ischar(opt)
    switch opt
        case 'random'
            CGRo=rand(n,1);
        case 'loop'
            CGRo=USM_CGR(M(:,end-min([32,m-1]):end),0.5); % do CGR of the last 32 symbols or as many as available if less
            CGRo=CGRo(:,end);
        otherwise
            error('CGR seeding option not defined')
    end
else
    error(['opt variable needs to be a number or string, not a ',class(opt)])
end

% CGR iteration
M(:,1)=CGRo+(M(:,1)-CGRo).*0.5;
for i=2:m
    M(:,i)=M(:,i-1)+(M(:,i)-M(:,i-1)).*0.5;
end
