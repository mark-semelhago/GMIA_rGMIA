function [beta, theta] = maximum_likelihood_estimation(simulation_output,simulation_output_variances,bound_constraints,design_points_index,sobol_bounds,pardiso_flag)
%maximum_likelihood_estimation estimates the sample mean (beta) and
%conditional precisions (theta) associated with the simulation)


    %Generate the Sobol point set used for searching the box defined by the
    %bounds
    d = size(bound_constraints,2);
    p = sobolset(d);
    thetas = -1*(sobol_bounds(1,2)*p(2:2001,:)+sobol_bounds(1,1));


    %Test the precision matrices parameterized by generated thetas for positive
    %definiteness
    positive_definite = ones(size(thetas,1),1);
    for i=1:size(thetas,1)
        Q = prec_mat_constr([1, -thetas(i,:)]', bound_constraints(:,2)-bound_constraints(:,1)+1);
        [~, positive_definite(i)] = chol(Q,'lower');
    end
    
    %Evaluate each set of beta and theta and choose the one that has the
    %lowest negative log likelihood
    positive_definite = (positive_definite==0);
    thetas = thetas(positive_definite,:);
    theta_0s = rand(size(thetas,1),1);
    betas = zeros(size(thetas,1),1);
    neglogL = zeros(size(thetas,1),1);
    for i = 1:size(thetas,1)
        [neglogL(i), betas(i)] = negloglikelihood([theta_0s(i) thetas(i,:)]', bound_constraints, design_points_index,simulation_output,simulation_output_variances, pardiso_flag);
    end
    [~,opt] = min(neglogL);
    theta = [theta_0s(opt) -thetas(opt,:)]';
    beta = betas(opt);


end



