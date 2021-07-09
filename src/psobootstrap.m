% run PSO a bunch of times and record the best solution

results = zeros(1000, 7);

powers = [-2, -2];
beta = [1,1,1];
sens = false;
pts = 3;
algo = {@PSO, 100000, 1000, 0.4};

step = 0; % current number of iterations

for i = 1:1000
    
   [obj_i, xi_i] = DoptFracPoly(beta, powers, pts, sens, algo);
   results(i, :) = [obj_i, xi_i];
   
   step = step + 1;
    
end

results = sortrows(results, 1);

% answer I found:
% 2.5970    0.3799   14.5349    0.2927    0.3390    0.3452    0.3158
% still doesn't beat the 2.5772 design I get using genetic algorithm with
% high variances of crossover and mutation