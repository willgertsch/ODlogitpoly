% find locally optimal designs for logistic model

% design options
beta = [2, 0, -4]; % nominal values, enter 3 values to add squared term
lower = -Inf; % design interval
upper = Inf;
numpts = 3; % number of design ponts
powers = [-0.5, 1]; % powers to use for fractional polynomials, empty => standard

% algorithm options
method = @GA; % algorithm to use
iterations = 100; % number of iterations
swarm = 1000; % size of the swarm
proC = 1; % prob. of crossover
disC = 20; % distr. index of cross-over: variance around parents
proM = 1; % prob. of mutation
disM = 20; % distr. index of mutation: variance of mutation

% sensitivity functions
% matrix inverse can cause problems in computing sensitivity function
% can disable if set to false
sens = false;

% work in progress
% matrix singularity gives issues with A and E opt.
opt = 'D';

% function evaluations
maxFE = iterations * swarm; 

% guess the interval the 2 design points are in
% this doesn't matter so much since the repair function will move the
% points automatically
init_low = -5; 
init_up = 5;
init = @(n, d) [unifrnd(init_low, init_up, n/2, numpts),...
    unifrnd(0, 1, n/2, numpts)]; % init swarm in interval

% create upper and lower bounds for design points and weights
% length of these vectors determine the number of design points
% need an equal number of weights
upperbounds = [upper*ones([1,numpts]), Inf*ones([1,numpts])];
lowerbounds = [lower*ones([1,numpts]), zeros(1,numpts)];

if isequal(method, @GA)
    algo = {method, proC, disC, proM, disM};
else
    algo = method;
end

% find the optimal design
[Dec,Obj,Con] = platemo('objFcn', @(x,d) logistic(x, beta, opt, powers),...
    'lower', lowerbounds, 'upper', upperbounds, ...
    'algorithm', algo, ...
    'decFcn', @(x, d) sum_constraints(x, upper, lower),...
    'maxFE', maxFE, ...
    'N', swarm, ...
    'initFcn', init);

% find the optimal design
out = sortrows([Dec, Obj], numpts*2 + 1); % sort by objective function
xi = out(1, :); % first row has the point with the lowest objective value
disp("Design found:")
disp(xi(1:numpts))
disp(xi(numpts+1:end-1))
fprintf("Objective value: %f\n", xi(end));

if sens
    % check if design points are optimal
    fprintf("Sensitivity function values:\n");
    for i = 1:numpts
        ch_i = ch_logistic(xi(i), xi(1:end-1), beta, opt, powers);
        fprintf("%f ", ch_i);
    end
    fprintf("\n");
    
    % plot
    % set bounds for plot
    % will follow design interval if bounded
    % otherwise will be derived from solutions
    if upper < Inf
        up = upper + 1;
    else
        up = max(xi(1:numpts)) + 1;
    end
    
    if lower > -Inf
        low = lower - 1;
    else
        low = min(xi(1:numpts)) - 1;
    end
    
    plot_logistic(xi(1:end-1), beta, opt, low, up, powers);
end

