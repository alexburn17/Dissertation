function [I_mat]=Immigration(immigration_rate_vec, n)
%
% DESCRIPTION: Immigration rate Param Sweep returns a matrix that includes 
% the means and standard deveviations for the proportions of infected and 
% Susceptable for the model runs for with a sample size of n
% 
% INPUTS:
% immigration_rate_vec is a vector or immigration param values
% n is the number of model runs
%
% OUTPUTS:
% I_mat is a matrix of means and standard deviations for n runs 
%
% standard paramaters that do not change in this model:
random = 1;
initmap = 0;
dim=100;
maxt= 150;
p = .15;
a=7;
g=30;
method = 'synchronous';


for i=1:n
% take first value from passed param list and use in model    
immigration_rate = immigration_rate_vec(1);
[SIR]=greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);
% store the mean proportions for infected and susceptable in vecotors of length n
Infected = SIR(1,:)./(SIR(1,:) + SIR(2,:));
Susceptable = SIR(2,:)./(SIR(1,:) + SIR(2,:));
% calc means for I and S
INF_I1(i) = mean(Infected);
SUS_I1(i) = mean(Susceptable);

immigration_rate = immigration_rate_vec(2);
[SIR]=greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);
% store the mean proportions for infected and susceptable in vecotors of length n
Infected = SIR(1,:)./(SIR(1,:) + SIR(2,:));
Susceptable = SIR(2,:)./(SIR(1,:) + SIR(2,:));
% calc means for I and S
INF_I2(i) = mean(Infected);
SUS_I2(i) = mean(Susceptable);
% calc means for I and S
immigration_rate = immigration_rate_vec(3);
[SIR]=greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);
% store the mean proportions for infected and susceptable in vecotors of length n
Infected = SIR(1,:)./(SIR(1,:) + SIR(2,:));
Susceptable = SIR(2,:)./(SIR(1,:) + SIR(2,:));
% calc means for I and S
INF_I3(i) = mean(Infected);
SUS_I3(i) = mean(Susceptable);

immigration_rate = immigration_rate_vec(4);
[SIR]=greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);
% store the mean proportions for infected and susceptable in vecotors of length n
Infected = SIR(1,:)./(SIR(1,:) + SIR(2,:));
Susceptable = SIR(2,:)./(SIR(1,:) + SIR(2,:));
% calc means for I and S
INF_I4(i) = mean(Infected);
SUS_I4(i) = mean(Susceptable);


end

% create a matrix of all vectors
I_mat = [INF_I1; INF_I2; INF_I3; INF_I4; SUS_I1; SUS_I2; SUS_I3; SUS_I4];

% find the mean, standard deviation and create a catagory variable for 
% infected (1) or suscetalbe (0) and merge into one output matrix. 
mn = mean(I_mat, 2);
st = std(I_mat, 0, 2);
cat = [1,1,1,1,0,0,0,0];
cat = cat';

% final output matrix
I_mat = [mn st cat];


end

