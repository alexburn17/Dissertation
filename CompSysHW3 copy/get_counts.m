function [S, I, R] = get_counts(map, dim, a)
% Function to return counts of susceptible (S), infected (I) and Recovered individuals(R)
% by Maike Holthuijzen
% and Alex Burnham
% updated March 2018
%
% INPUTS:
% map : 2D matrix of species numbers
% dim : dimension of map
% a : value of a in GH CA
%
% OUTPUTS:
% 3 vectors with counts of cells in S, I, and R states 
%
SIRmap=map;
for i = 1:dim
    for j = 1:dim
        if SIRmap(i,j) == 0
            SIRmap(i, j) = SIRmap(i,j);
        elseif SIRmap(i,j) > a
            SIRmap(i,j) = 2;
        else
            SIRmap(i,j) = 1;
        end
    end
end

S=sum(sum(SIRmap < 1)); %sum values for S, I, and R and store in vectors
I=sum(sum(SIRmap == 1));
R=sum(sum(SIRmap > 1));
end

