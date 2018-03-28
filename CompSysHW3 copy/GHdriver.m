%Driver function to reproduce figures in report
%By Maike Holthuijzen and Alex Burnham
%Updated March 5, 2018

%Fig.1. Reproduce top row of fig. 6.7
%No immigration, synchronous updating, deterministic 

random = 0;
initmap = [1 0; 0 3]; 
initmap = padarray(initmap, [7 7], 0, 'both');
dim=length(initmap);
maxt= 18;
p = 1;
a = 1;
g = 2;
method = 'synchronous';
immigration_rate = 0;
greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);

%Fig 2. Reproduce second row of fig. 6.7
%No immigration, synchronous updating, deterministic 
random = 0;
initmap = [1 2; 0 3]; 
initmap = padarray(initmap, [7 7], 0, 'both');
dim=length(initmap);
maxt= 18;
p = 1;
a=1;
g=2;
method = 'synchronous';
immigration_rate = 0;
greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);

%Fig3. %Top row of fig. 6.7
%No immigration, ASYNCHRONOUS updating, deterministic , no immmigration

random = 0;
initmap = [1 0; 0 3]; 
initmap = padarray(initmap, [7 7], 0, 'both');
dim=length(initmap);
maxt= 18;
p = 1;
a = 1;
g = 2;
method = 'asynchronous';
immigration_rate = 0;
greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);


%Fig 4. %using initial matrix for fig. 6.7:
%immigration rate =0.5, Synchronous updating, deterministic 
random = 0;
initmap = [1 0; 0 3]; 
initmap = padarray(initmap, [7 7], 0, 'both');
dim=length(initmap);
maxt= 18;
p = 1;
a = 1;
g = 2;
method = 'synchronous';
immigration_rate = 0.5;
greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);


%fig 5. Stochastic model
%Stochastic, synchronous updating, immigration=0.05, p=0.2
random = 1;
initmap = 0;
dim=16;
maxt= 18;
p = .2;
a=3;
g=6;
method = 'synchronous';
immigration_rate = 0.05;
greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);


%fig 6. Stochastic model
%Stochastic, synchronous updating, no immigration, p=0.2
random = 1;
initmap = 0;
dim = 16;
maxt = 18;
p = .2;
a = 3;
g = 6;
method = 'synchronous';
immigration_rate = 0;
greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);

%Fig 7. Stochastic model 
% Stochastic, synchronous updating, immigration rate = 0.05, p = 0.2
random = 1;
initmap = 0;
dim=90;
maxt= 300;
p = .2;
a=2;
g=4;
method = 'synchronous';
immigration_rate = 0.0001;
greenbergHastingsStarter(dim,a,g, maxt, p, method, immigration_rate, initmap, random);

