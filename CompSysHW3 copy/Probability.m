function [P_mat]=Probability(probability_vec, n)
%
% DESCRIPTION: Probabilty of Infection Param Sweep returns a matrix that 
% includes the means and standard deveviations for the proportions of 
% infected and Susceptable for the model runs for with a sample size of n
% 
% INPUTS:
% probability_vec is a vecotor of probabilities of infection (p) param values
% n is the number of model runs
%
% OUTPUTS:
% P_mat is a matrix of means and standard deviations for n runs 
%
% standard paramaters that do not change in this model:
random = 1;
initmap = 0;
dim=100;
maxt= 150;
immigration_rate=0;
method = 'synchronous';
a = 7;
g = 30;


for i=1:n
% take first value from passed param list and use in model    
p = probability_vec(1);
[SIR]=greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);
% store the mean proportions for infected and susceptable in vecotors of length n
Infected = SIR(1,:)./(SIR(1,:) + SIR(2,:));
Susceptable = SIR(2,:)./(SIR(1,:) + SIR(2,:));
% calc means for I and S
INF_P1(i) = mean(Infected);
SUS_P1(i) = mean(Susceptable);

p = probability_vec(2);
[SIR]=greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);
% store the mean proportions for infected and susceptable in vecotors of length n
Infected = SIR(1,:)./(SIR(1,:) + SIR(2,:));
Susceptable = SIR(2,:)./(SIR(1,:) + SIR(2,:));
% calc means for I and S
INF_P2(i) = mean(Infected);
SUS_P2(i) = mean(Susceptable);
% calc means for I and S

p = probability_vec(3);
[SIR]=greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);
% store the mean proportions for infected and susceptable in vecotors of length n
Infected = SIR(1,:)./(SIR(1,:) + SIR(2,:));
Susceptable = SIR(2,:)./(SIR(1,:) + SIR(2,:));
% calc means for I and S
INF_P3(i) = mean(Infected);
SUS_P3(i) = mean(Susceptable);

p = probability_vec(4);
[SIR]=greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);
% store the mean proportions for infected and susceptable in vecotors of length n
Infected = SIR(1,:)./(SIR(1,:) + SIR(2,:));
Susceptable = SIR(2,:)./(SIR(1,:) + SIR(2,:));
% calc means for I and S
INF_P4(i) = mean(Infected);
SUS_P4(i) = mean(Susceptable);


end

% create a matrix of all vectors
P_mat = [INF_P1; INF_P2; INF_P3; INF_P4; SUS_P1; SUS_P2; SUS_P3; SUS_P4];

% find the mean, standard deviation and create a catagory variable for 
% infected (1) or suscetalbe (0) and merge into one output matrix. 
mn = mean(P_mat, 2);
st = std(P_mat, 0, 2);
cat = [1,1,1,1,0,0,0,0];
cat = cat';

% final output matrix
P_mat = [mn st cat];


end

