function [D_mat]=Duration(duration_vec, n)
%
% DESCRIPTION: Duration of Infection Param Sweep returns a matrix that 
% includes the means and standard deveviations for the proportions of 
% infected and Susceptable for the model runs for with a sample size of n
% 
% INPUTS:
% duration_vec is a vector of duration (a) param values
% n is the number of model runs
%
% OUTPUTS:
% D_mat is a matrix of means and standard deviations for n runs 
%
% standard paramaters that do not change in this model:
random = 1;
initmap = 0;
dim=100;
maxt= 150;
p = .15;
immigration_rate=0;
g = 30;
method = 'synchronous';


for i=1:n
% take first value from passed param list and use in model    
a = duration_vec(1);
[SIR]=greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);
% store the mean proportions for infected and susceptable in vecotors of length n
Infected = SIR(1,:)./(SIR(1,:) + SIR(2,:));
Susceptable = SIR(2,:)./(SIR(1,:) + SIR(2,:));
% calc means for I and S
INF_D1(i) = mean(Infected);
SUS_D1(i) = mean(Susceptable);

a = duration_vec(2);
[SIR]=greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);
% store the mean proportions for infected and susceptable in vecotors of length n
Infected = SIR(1,:)./(SIR(1,:) + SIR(2,:));
Susceptable = SIR(2,:)./(SIR(1,:) + SIR(2,:));
% calc means for I and S
INF_D2(i) = mean(Infected);
SUS_D2(i) = mean(Susceptable);
% calc means for I and S

a = duration_vec(3);
[SIR]=greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);
% store the mean proportions for infected and susceptable in vecotors of length n
Infected = SIR(1,:)./(SIR(1,:) + SIR(2,:));
Susceptable = SIR(2,:)./(SIR(1,:) + SIR(2,:));
% calc means for I and S
INF_D3(i) = mean(Infected);
SUS_D3(i) = mean(Susceptable);

a = duration_vec(4);
[SIR]=greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);
% store the mean proportions for infected and susceptable in vecotors of length n
Infected = SIR(1,:)./(SIR(1,:) + SIR(2,:));
Susceptable = SIR(2,:)./(SIR(1,:) + SIR(2,:));
% calc means for I and S
INF_D4(i) = mean(Infected);
SUS_D4(i) = mean(Susceptable);


end

% create a matrix of all vectors
D_mat = [INF_D1; INF_D2; INF_D3; INF_D4; SUS_D1; SUS_D2; SUS_D3; SUS_D4];

% find the mean, standard deviation and create a catagory variable for 
% infected (1) or suscetalbe (0) and merge into one output matrix. 
mn = mean(D_mat, 2);
st = std(D_mat, 0, 2);
cat = [1,1,1,1,0,0,0,0];
cat = cat';

% final output matrix
D_mat = [mn st cat];


end
