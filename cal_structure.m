function [nodf,qb,Nm] = cal_structure(matrix_PA)

% Matlab library "BiMat" is used for analyzing the network structures, 
% which is publicly available via the link: https://bimat.github.io/. 

fp = Bipartite(matrix_PA);

% Nestedness 
fp.nestedness = NestednessNODF(fp.matrix);
fp.nestedness.Detect();
nodf = fp.nestedness.N;

% Modularity
fp.community = LeadingEigenvector(fp.matrix);  % Newman
fp.community.DoKernighanLinTunning = true;  
fp.community.Detect();
qb = fp.community.Qb;
Nm = fp.community.N;
