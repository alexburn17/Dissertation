function [neighbors,indeces,rows,cols]=getneighbors(map,r,c,randomize)
% Original function by Maggie Eppstein
% modified by Maike Holthuijzen and Alex Burnham
%
%Updated March 2018
% 
% DESCRIPTION: returns nearest neighors of locations <r,c> according to
% the following "stencil" (where <r,c> is the "o" and "x" represents a
% neighbor) -- uses toroidal wraparound on map
%         x 
%       x o x
%         x 
% 
% INPUTS:
% map == 2D matrix of species numbers with one extra row and column
% r,c == scalar or vector of n row and col indeces of cells to find neighbors for
% randomize == optional boolean input flag (true ==> randomize order of neighbors)
%
% OUTPUTS:
% neighbors == n by 4: each row corresponds to all 4 neighbors of one cell
% indeces == 4 by 4: the 1D indeces into the map of the neighbors
% rows, cols == each are n by 4: the 2D subscripts corresponding to the indeces

dim=length(map)-1; %ignore extra row and col in map

rb=addwrap(r,-1,dim); %row behind
cb=addwrap(c,-1,dim); %column behind
rf=addwrap(r,1,dim); %row in front
cf=addwrap(c,1,dim); %column in front

% list of all 4 nearest neighbors, relative to position <r,c>
% for von neuman neighborhood
rows=[rb  r rf r];
cols=[c cf c cb];

if nargin==4 && randomize %if desired, randomize the order of the neighbors
    if randomize
        for guy=1:size(rows,1) % if r and c were passed in as vectors
            order=randperm(8);
            rows(guy,:)=rows(guy,order);
            cols(guy,:)=cols(guy,order);
        end
    end
end

% convert 2D subscripts in 1D indeces into map
indeces=sub2ind(size(map),rows,cols);
neighbors=map(indeces); %these are the actual neighbors