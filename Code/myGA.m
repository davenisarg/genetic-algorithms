function [ fxbest, xbest ] = myGA( func, xrange , tolerance, popSize, generationsNum, crossProb, mutateProb )
%genetic algorithm
% 
% output:
% fxbest: the minimum of the func
% xbest: the point that the func gets the minimum
% 
% input:
% func: the name of the func, e.g. 'afunc'
% xmin: the low limitation of the variables, e.g. [-10,-pi]
% xmax: the high limitation of the variables, e.g. [3,exp(1)]
% tolerance: the tolerance during calculation, e.g. 1e-4
% popSize: the size of population, must be a positive even integer, e.g. 100
% generationsNum: the number of generations, the terminal criteria, e.g. 300
% crossProb: the probability of crossing over, e.g. 0.6
% mutateProb: the probability of mutation, e.g. 0.001
% 
% Let there be afunc.m file which contains the function. I did upload that
% with the submission.
%  function [ y ] = afunc( x)
%    %De Jong's
%    y = sum(x.^2);
%  end
% so run the following in command window:
%  
%  [fxbest, xbest] = myGA('afunc',100)

xmin = -xrange/2
xmax = xrange/2
if ~exist('tolerance','var')
    tolerance = 1e-4;
end
if ~exist('popSize','var')
    popSize = 100;
end
if ~exist('generationsNum','var')
    generationsNum = 300;
end
if ~exist('crossProb','var')
    crossProb = 0.6;
end
if ~exist('mutateProb','var')
    mutateProb = 0.001;
end

[initPop, chromosomeEachSizes, chromosomeSize] = init(popSize, tolerance, xmin, xmax);
pop = initPop;
i=1;
tic;
while i <= generationsNum
    [ fvals, varsTransformed ] = decodeandcal (func ,pop ,chromosomeEachSizes ,xmin ,xmax);
    [ popNewSub, elite, fxbest ] = selectchromo( pop, fvals );
    [ popNewSub ] = crossover(popNewSub,crossProb);
    [ popNewSub ] = mutate( popNewSub, mutateProb);
    chromosomeRandom = randi([0,1], [1,chromosomeSize]);
    pop = [popNewSub;chromosomeRandom;elite];
    i = i + 1;
end
toc;
[ fxbest, xbest ] = decodeandcal (func ,elite ,chromosomeEachSizes ,xmin ,xmax);
toc;
end

function [ initPop, chromosomeEachSizes, chromosomeSize ] = init( popSize, tolerance, xmin, xmax)
%INIT Generate the initial populations
% output:
% initPop: the initial population, each row is a chromosome
% chromosomeEachSizes: the vector of bit lengths of each variable
% chromosomeSize: the summarized bits length
% input:
% popSize: how many populations there are
% tolerance: the tolerance of the result, e.g. 1e-4
chromosomeEachSizes = ceil( log2( (xmax-xmin)./tolerance ) );
chromosomeSize = sum(chromosomeEachSizes);
initPop = randi([0,1], [popSize,chromosomeSize]);
end

function [ fvals, varsTransformed ] = decodeandcal (func ,pop ,chromosomeEachSizes ,xmin ,xmax)
%DECODEANDCAL Check the func and variables values
% output:
% fvals: the vector of value of func
% varsTransformed: the variables transformed from chromosomes
% input:
% func: the name of the func to be optimizated
% pop: the previous population generated
% chromosomeEachSizes: the vector of bit length of each variable
% xmin: the low limitation of the x, e.g. [-10, -100]
% xmax: the high limitation of the x, e.g. [10, 50]
varNum = length(chromosomeEachSizes); % the number of variables
popSize = size(pop ,1);
transformed = (xmax-xmin) ./ (2.^chromosomeEachSizes-1); % the scale of values
chromosomeEachSizes = [0 cumsum(chromosomeEachSizes)];
for i = 1:varNum
    
    popVar{i} = pop(:, chromosomeEachSizes(i) + 1:chromosomeEachSizes(i + 1) );
    % the subpop with subchromosomes_i
    var{i} = sum(ones(popSize ,1)*2.^(size(popVar{i},2)-1:-1 :0).*popVar{i} ,2) .* transformed(i) + xmin(i) ;
end
varsTransformed = [var{1,:}];
for i = 1 :popSize
    fvals(i) = eval([func ,'(varsTransformed(i , :) )']);
end
end

function [ popNewSub, elite, fxbest ] = selectchromo( popPrev, fvals )
%SELECTCHROMO Select, choose sub new population and elite
% output:
% popNewSub: the new sub population (size = size-2)
% elite: the elite chromosome of previous population
% input:
% popPrev: the previous population
% fvals: the func values of previous population
fitness = (max(fvals)-fvals)';
popSize = size(popPrev,1);
[fitnessMin, indexMin] = min(fitness);
[fitnessMax, indexMax] = max(fitness);
elite = popPrev(indexMax,:); % the min-fval aka the best individial
fxbest = fvals(indexMax); % the func value on elite
listTemp = [1:popSize];
listTemp(indexMin) = 0;
listTemp(indexMax) = 0;
listTemp = nonzeros(listTemp);
popTemp = popPrev(listTemp,:);
fitnessTemp = fitness(listTemp,:);
popSizeToBeRenew = popSize - 2;
probAdded = cumsum(fitnessTemp / sum(fitnessTemp));
chromosomeSelected = sum(probAdded*ones(1,popSizeToBeRenew)<ones(popSizeToBeRenew,1)*rand(1,popSizeToBeRenew))+1;
popNewSub = popTemp(chromosomeSelected,:);
end

function [ popNew ] = crossover( popPrev, crossProb )
%CROSSOVER Generate the new population by crossing over
% output:
% popNew: the new generated population
% input:
% popPrev: the previous population
% crossProb: the probability of crossing over
[new, sortIndex] = sort(rand(size(popPrev ,1) ,1));
popSorted = popPrev(sortIndex, :);
pairsNum = size(popSorted, 1)/2;
chromosomeEachSizes = size(popSorted, 2);
parisToCross = rand(pairsNum, 1) < crossProb;
pointsToCross = parisToCross.*randi([1,chromosomeEachSizes],[pairsNum, 1]);

for i=1:pairsNum
    popNew([2*i-1,2*i],:)=[popSorted([2*i-1,2*i],1:pointsToCross(i)), popSorted([2*i,2*i-1],pointsToCross(i)+1:chromosomeEachSizes)];
end
end

function [ popNew ] = mutate( popPrev, mutateProb )
%MUTATE Generate the new population by mutation
% output:
% popNew: the new generated population
% input:
% popPrev: the previous population
% mutateProb: the probability of mutation
popNew = popPrev;
pointsToMutate = find(rand(size(popPrev))< mutateProb);
popNew(pointsToMutate) = 1 - popPrev(pointsToMutate);
end