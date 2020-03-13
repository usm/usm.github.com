function M=USM_compact(M)

%USM_COMPACT compacts hi matrix from sparce to compact USM
%Syntax: M=USM_compact(M)
%Description:
%  compacts a binary hit matrix, M. For example, compacting a
%  unidirectiopnal USM of a nucleotide sequence will produce the original
%  CGR hit matrix.
%
%See also: USM_CGR, which is where M is used to calculate the coordinates
%
%Jonas Almeida, almeidaj@musc.edu, Nov 2004

[n,m]=size(M);
K=repmat([1:n]',1,m);
M=M.*K;
M=sum(M);
B=dec2bin(0:n-1)';
M=B(:,M)=='1';