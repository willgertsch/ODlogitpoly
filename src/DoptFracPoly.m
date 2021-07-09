% function for finding D-optimal designs for 2nd degree fractional
% polynomials for logistic regression using GA
% code streamlines the original project.m file
% beta: length 3 vector of nominal coefficients
% powers: powers for x in linear predictor with p1 < p2
% pts: number of design points
% sens: calculate and plot sensitivity function for final design
% algo: {@METHOD, maxFE, swarm size, args}
function [obj, design] = DoptFracPoly(beta, powers, pts, sens, algo)

    % sanity check inputs
    if length(beta) ~= 3
       disp("beta needs to be length 3.");
       return
    elseif length(powers) ~= 2
        disp("powers needs to be length 2");
        return
    elseif powers(1) > powers(2)
        disp("Order powers such that p1 <= p2");
        return
    elseif pts < 1
        disp("Invalid number of design points");
        return
    end
    
    % set design interval to [0, inf]
    lower = 0;
    upper = Inf;
    
%     % algorithm options
%     method = @GA; % algorithm to use
%     iterations = 100; % number of iterations
%     swarm = 1000; % size of the swarm
%     proC = 1; % prob. of crossover
%     disC = 20; % distr. index of cross-over: variance around parents
%     proM = 1; % prob. of mutation
%     disM = 20; % distr. index of mutation: variance of mutation
%     algo = {method, proC, disC, proM, disM};
    
    % convert iterations to function evals
%     maxFE = iterations * swarm;

    method = [algo(1), algo(4:end)];
    maxFE = algo{2};
    swarm = algo{3};
    
    % initialize swarm in [0,5]
    init_low = 0; 
    init_up = 5;
    init = @(n, d) [unifrnd(init_low, init_up, n, pts),...
        unifrnd(0, 1, n, pts)];
    
    % create upper and lower bounds for design points and weights
    % length of these vectors determine the number of design points
    % need an equal number of weights
    upperbounds = [upper*ones([1,pts]), Inf*ones([1,pts])];
    lowerbounds = [lower*ones([1,pts]), zeros(1,pts)];
    
    
    % call PlatEMO to find the optimal design
    [Dec,Obj,~] = platemo('objFcn', @(x,d) logistic(x, beta, 'D', powers),...
        'lower', lowerbounds, 'upper', upperbounds, ...
        'algorithm', method, ...
        'decFcn', @(x, d) sum_constraints(x, upper, lower),...
        'maxFE', maxFE, ...
        'N', swarm, ...
        'initFcn', init);  
    
    % find the optimal design
    out = sortrows([Dec, Obj], pts*2 + 1); % sort by objective function
    xi = out(1, :); % first row has the point with the lowest objective value
    disp("Design found:")
    disp(xi(1:pts))
    disp(xi(pts+1:end-1))
    fprintf("Objective value: %f\n", xi(end));
    obj = xi(end);
    design = xi(1:end-1);
    
    if sens
        % check if design points are optimal
        fprintf("Sensitivity function values:\n");
        for i = 1:pts
            ch_i = ch_logistic(xi(i), xi(1:end-1), beta, 'D', powers);
            fprintf("%f ", ch_i);
        end
        fprintf("\n");

        % plot
        % special plot bounds for fractional polynomials
        if ~isempty(powers)
            low = 0;
            up = max(xi(1:pts)) * 1.5;
        end

        plot_logistic(xi(1:end-1), beta, 'D', low, up, powers);
    end
    
    return
end