% PRank: A Fast Analytical Rank Estimation
%
% Return the rank estimation of a given key k*={k1,...,kd} and
% d probability lists {P1,...,Pd}.
%
% The program reads d inputs files: "i.csv", for each 0 <= i <= d-1.
% Each file contains P_i, i.e., list of N probabilities in non-increasing order, each in one line.
% The correct subkey index k_i appears in "i.csv" as the last line after the N probabilities.
%
%
% Copyright (C) Avishai Wool and Liron David, 2019. All rights reserved.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% If you find this software useful, we ask that you cite our paper:
%
%    Liron David and Avishai Wool.
%    PRank: Fast Analytical Rank Estimation via Pareto Distributions.
%    COSADE 2019
%    Technical report https://eprint.iacr.org/2018/550.
%
% Send comments / requests / patches to:
%   Liron David <lirondavid@gmail.com>



d=16;
N=256;

% reads the 16 probability lists + corresponding index
probMatrix=zeros(d,N);
n=zeros(1,d);
index=zeros(1,d);

for i=1:d
    prob = csvread(strcat(num2str(i-1),'.csv'));
    n(i)=N;
    probMatrix(i,1:n(i))=(prob(1:n(i),1))';
    index(i)=prob(N+1,1)+1;
end


%merges the d=16 into d=8
mergeProb=zeros(1,1);
mergeIndex=zeros(1,1);
for i=1:d/2
    [y,k]=merge(probMatrix(2*i-1,1:n(2*i-1)),n(2*i-1),probMatrix(2*i,1:n(2*i)),n(2*i),index(2*i-1),index(2*i));
    mergeProb(i,1:n(1)^2)=y;
    mergeIndex(i)=k;
end
d=d/2;
N=n(1)^2;
n=N*(ones(1,d));
probMatrix=mergeProb;
index=mergeIndex;

probIndex=zeros(1,d);

%calculates the Pi(k_i) and p*
for i=1:d
   probIndex(i)=probMatrix(i,index(i));
end 
pFinal=prod(probIndex);

%calculates the Pareto-like of each {P1,..,P8}
stime=tic;
alpha=zeros(1,d);
a=zeros(1,d);

for i=1:d
    [a(i),alpha(i)] = upper_estimation_fast(probMatrix(i,1:n(i)),n(i),index(i),probIndex(i));
end
preParetoTime=toc(stime);

%calculates the upper for d=8
[pareto_upper,pareto_time] = pareto_estimation(probMatrix,index,d,n,a,alpha);
paretoUpper=log2(pareto_upper)
paretoTime=pareto_time+preParetoTime



