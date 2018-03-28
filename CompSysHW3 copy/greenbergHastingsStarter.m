 function [SIR]=greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random)
% Implement SIR model in a Cellular Automata using Toroidal Topology and
% Von Neumann neighborhoods with options for matrix initialization,
% immigration, and updating
% based on: J. M. Greenberg and S. P. Hastings, "Spatial Patterns for Discrete Models of 
% Diffusion in Excitable Media", SIAM J. Appl. Math., 34:515-523, 1978
% 
% Modified by Maike Holthuijzen and Alex Burnham
% Updated March 2018
%
% Details: Von Nuemann neighborhood for  a cell 'x' is defined as the 4
% neighbors 'o' to the right, left, top and bottom of the cell:
%   o
% o x o
%   o
% Note: This function calls addwrap.m, getneighbors.m, and get_counts.m
% 
% FUNCTION INPUTS:
% dim: (positive integer) dimensions of square map
% 
% a: (positive integer) duration of infection in an individual
% 
% g: (positive integer) duration of immunity in an individual
% 
% maxt: (positive integer) maximum number of timesteps to run (user can abort early)
% 
% method: string specifying synchronous 'synchronous' or asynchronous
% 'asynchronous' updating. In synchronous updating, cells are updated
% simultaneously after all initial cell conditions at time t-1 are taken
% into account. In the case of asynchronous updating, cells are updated
% randomly.
% 
% p: disease transmission probability (0-1 value). If p = 1, then susceptible cells
% will be updated if they are surrounded by an infected neighbor within the
% Von Neumann neighborhood
%
% immigration_rate: probability of immigration. In this function
% immigration is implemented by randomly changing a percentage of cells to
% infected status (e.g. a number between 1 and a). If immigration rate is
% 0, no immigration will take place. Immigration takes place only after all
% cells have been updated for both syncrhonous and asynchronous
% implementations
%
% random: implements stochastic initialization of the SIR matrix. 
% If random is 1, then the initial matrix will be randomly
% populated (note: if random = 1, then initmap MUST be set to 0.)
%
% initmap: a square matrix of initial integer starting conditions, used for
% deterministic simulation. If initmap is passed as a matrix, 'random' must
% be set to 0.
%
% Example function calls:
% deterministic implementation with synchronous updating, no immigration:
% random = 0;
% initmap = [1 0; 0 3]; 
% initmap = padarray(initmap, [7 7], 0, 'both');
% dim=length(initmap);
% maxt= 18;
% p = 1;
% a = 1;
% g = 2;
% method = 'synchronous';
% immigration_rate = 0;
% greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);
% 
% stochastic implementation with probability of infection = 0.2: 
% random = 1;
% initmap = 0;
% dim=16;
% maxt= 18;
% p = .2;
% a=3;
% g=6;
% method = 'synchronous';
% immigration_rate = 0;
% greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random); 
% 
% OUTPUT: 
% SIR: A 3 x maxt column matrix of [nS,nI,nR] (counts of S, I , and R at each timestep)
%
% HIGH-LEVEL ABSTRACT STATES
% S = Susceptible to infection
% I = Infected/Infectious 
% R = Recovered (and immune from re-infection)
% 
% TRANSITION RULES:
% A cell remains I for exactly “a” timesteps then becomes R
% A cell remains R for exactly “g” timesteps then becomes S
% With probability p, an I cell can infect a neighbor that is S (if p==1,deterministic)
% thus, S turns into I with an overall prob of 1-(1-p)^n, if it has n neighbors that are I
%
% LOW-LEVEL STATES FOR IMPLEMENTATION
%
%	A low-level value of 0 is interpreted as S
%	Values from 1 to a are interpreted as I
%	Values from a+1 to a+g are interpreted as R
%

%In deterministic implementation, set map equal to the initial map and
%create an identical matrix to keep track of cell updates
if ~random
    map = initmap;
    map_after = repmat(map, 1);

else
    %In stochastic implementation, the initial matrix is populated with the
    %values between 0 and (a+g), inclusive
    map = randi([0 a+g],dim,dim);
    map_after = repmat(map, 1);
end

%plot initial configuration of SIR matrix
[fighandle,plothandle] = plotMapInNewFigure(map,a,g); % PLOT INITIAL MAP (see function below)


%preallocate vectors for NR, nS, and NI 
nR=zeros(maxt, 1); %number of recovered
nS = zeros(maxt,1); %number of susceptible
nI=zeros(maxt, 1); %number of infected


for t=1:maxt 
    
    %get counts of Ss, Is, and Rs and store in preallocated vectors
    %with function 'get_counts'
    [S, I, R] = get_counts(map, dim, a);

    nS(t)=S; %sum values for S, I, and R and store in vectors
    nI(t)=I;
    nR(t)=R;
    
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%SYNCHRONOUS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(method, 'synchronous')
    
    for row = 1:dim
        for col=1:dim
            if map(row, col) == 0
                %obtain neighbors for susceptible cells according to a Von
                %Neumann neighborhood
                [neighbors,~,~,~]=getneighbors(map,row,col);
                num_neighbors = sum((neighbors >=1) & (neighbors <(a+1)));
                %If p is 1 and there are neighbors who are infected, cell
                %automatically becomes infected
                    if (p == 1) && (num_neighbors >=1)
                        disp('p =1!')
                        map_after(row, col) = 1;
                    
                    %if p is not 1, we calculate the probability of         
                    else 
                        prob_trans = 1-(1-p)^num_neighbors;
                        is_infected = rand(1);
                            %update a susceptible cell to infected with
                            %'prob_trans'
                            if is_infected < prob_trans 
                            map_after(row,col) = 1; 
                            end
                    end
                    
            end
        end
    end

    %perform updates to cells that are NOT 0
    %we add one to each cell that is not 0 and then mod by (a+g+1)
     for i=1:dim
        for j=1:dim
            %if cell is not 0, add 1, then mod by (a+g+1)
            if map(i,j) ~= 0
                map_after(i,j) = mod(map(i,j)+1, (a+g+1));
            end
        end
     end
     
     %take care of case of immigration
         if immigration_rate ~= 0
            total_cells = dim*dim;
            num_infected = randsample(total_cells, floor((immigration_rate*(dim*dim))));           
            [x, y] = meshgrid(1:dim, 1:dim);
            index_grid =  [x(:) y(:)];
            random_indices = index_grid(num_infected,:);
            for i=1:length(random_indices)
                %change some cells in timemap to 'infected' (a number between 1
                %and a
                
                r=random_indices(i, 1);
                c=random_indices(i,2);
              
                randomnumber = randi([1 a],1,1);
            
                map_after(r,c) = randomnumber;
                
            end
        end
    
    
    %UPDATE PLOT WITH NEW MAP DATA
    set(plothandle,'cdata',map_after); 
    figure(gcf),drawnow; 
    
   % Bail out if the user closed the figure 
    if ~ishandle(fighandle)
        plotMapInNewFigure(map_after,a,g); %plot final map and exit
        title('Final configuration')
        break
    end
    
    %reset the map from current step to be the starting map for next
    %timestep
     map = map_after;
     
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%ASYNCHRONOUS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(method, 'asynchronous')
    
    [x, y] = meshgrid(1:dim, 1:dim);
    index_grid =  [x(:) y(:)];
    random_indices = randperm(dim*dim);
    
    %Must consider (dim*dim) number of cells
        for i=1:(dim*dim) 
            my_index = random_indices(i);
            cell_index= index_grid(my_index,:);
            row = cell_index(1);
            col = cell_index(2);
                if map(row, col) == 0
                    [neighbors,~,~,~]=getneighbors(map,row,col);
                    num_neighbors = sum((neighbors >=1) & (neighbors <(a+1)));
                    %If p is 1 and there are neighbors who are infected, cell
                    %automatically becomes infected. We update the actual map
                    %in the case of asynchronous updating
                        if (p == 1) && (num_neighbors >=1)
                            disp('p =1!')
                            map(row, col) = 1;
                    
                    %if p is not 1, we calculate the probability of         
                        else 
                            %disp('p is not 1!');
                            prob_trans = 1-(1-p)^num_neighbors;
                            is_infected = rand(1);
                                %update a susceptible cell to infected with
                                %'prob_trans'
                                if is_infected < prob_trans 
                                map(row,col) = 1; 
                                disp('update a 0 value')
                   
                                end
                        end
                    
                 else 
                %perform updates to cells that are NOT 0
                %we add one to each cell that is not 0 and then mod by (a+g+1)
                %if cell is not 0, add 1, then mod by (a+g+1)
            
                map(row,col) = mod(map(row,col)+1, (a+g+1));
                %disp('map before immigration')

                end
        end
        
        
        if immigration_rate ~= 0
            total_cells = dim*dim;
            num_infected = randsample(total_cells, floor((immigration_rate*(dim*dim))));
            
            [x, y] = meshgrid(1:dim, 1:dim);
            index_grid =  [x(:) y(:)];
            random_indices = index_grid(num_infected,:);
            for i=1:length(random_indices)
                %change some cells in timemap to 'infected' (a number between 1
                %and a
                disp('made a change')
                r=random_indices(i, 1);
                c=random_indices(i,2);

                randomnumber = randi([1 a],1,1);
             
                map(r,c) = randomnumber;
              
            end
        end
 
        
        set(plothandle,'cdata',map); 
        figure(gcf),drawnow; 
        %disp('updated map')
    
       % Bail out if the user closed the figure (I suggest you also add a
        % condition to bail out if no change in the map)
        if ~ishandle(fighandle)
            plotMapInNewFigure(map,a,g); %plot final map and exit
            title('Final configuration')
            break
        end
    end
    

  SIR=vertcat(nI',nS',nR');
   
end



%%% HELPER FUNCTION TO PLOT MAP USING IMAGESC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fighandle,plothandle] = plotMapInNewFigure(map, a, g)
fighandle=figure;
set(fighandle,'position',[42   256   560   420]); % specify location of figure
plothandle=imagesc(map); % doesn't truncate a row and col like pcolor does
colormap(jet);
set(gca,'clim',[0 a+g]); % make sure the color limits down change dynamically
ch=colorbar;
set(ch,'Ytick',[0 a a+g],'Yticklabel',{'S','I up to here','R up to here'})
title(['Figure by MFH and AB'])
axis('square') %make sure aspect ratio is equal

