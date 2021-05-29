function [neglogL, beta] = negloglikelihood(theta, bound_constraints, design_points_index, simulation_output, simulation_output_variances, pardiso_flag)
%negloglikelihood evaluates the negative log likelihood given the mean and
%conditional precision parameters assuming a GMRF structure

Xi1 = setdiff(linspace(1,prod(bound_constraints(:,2)-bound_constraints(:,1)+1),prod(bound_constraints(:,2)-bound_constraints(:,1)+1))',design_points_index);
n2 = length(design_points_index);
Q = prec_mat_constr(theta, bound_constraints(:,2)-bound_constraints(:,1)+1);
Q_11 = Q(Xi1,Xi1);
Q_12 = Q(Xi1,design_points_index);
Q_22 = Q(design_points_index,design_points_index);

if pardiso_flag
    global pardiso1 id1
    id1 = 1;
    pardiso1 = pardiso(id1,-2);
    pardiso1.factorize(id1,tril(Q_11));
    Sigma_22 = pardiso1.solve(id1,tril(Q_11),full(Q_12));
else
    Sigma_22 = Q_11\Q_12;
end
Sigma_22 = full(inv(Q_22 - Q_12'*Sigma_22));
T = (Sigma_22 + diag(simulation_output_variances))*theta(1); %This step is necessary to avoid numerical error.

beta = (ones(n2,1)'*((Sigma_22 + diag(simulation_output_variances))\ones(n2,1)))\(ones(n2,1)'*((Sigma_22+diag(simulation_output_variances))\simulation_output));
neglogL =-0.5*n2*log(theta(1)) +0.5*(log(det(T))) + 0.5*(simulation_output-beta*ones(n2,1))'*((Sigma_22+diag(simulation_output_variances))\(simulation_output-beta*ones(n2,1)));




end


