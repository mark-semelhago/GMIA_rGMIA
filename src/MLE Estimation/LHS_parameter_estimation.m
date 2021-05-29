function eDO = LHS_parameter_estimation(simulator, bound_constraints, num_processors, rng_seed, k, rep, sobol_bounds, experimental_design_object_filename, pardiso_flag)
%LHS_parameter_estimation estimates the MLE parameters and stores the
%simulation information associated with an experimental design.
%Function inputs are the following:
%simulator: function handle associated with the simulator
%bound_constraints: matrix containing the box constraints of the feasible
%region (e.g., [4 6;3 6] implies a [4, 6] x [3, 6] feasible region
%intersected with the integer lattice
%num_processors: number of processors available for use
%rng_seed: random number generator seed
%k: number of design points to use
%rep: 2-element vector, where the first element specifies how many
%simulation replications to perform on a first visit, and the second how
%many simulation replications to perform on subsequent visits
%sobol_bounds: bounds to be used for the sobol random number generator to
%generate points (for a 2D feasible region [0.2, 0.11] seems to work well
%empirically)
%experimental_design_object_filename: string containing the filename of the
%workspace where the ExperimentalDesignObject will be stored
%pardiso_flag: binary (0/1) flag to indicate whether PARDISO should be used
%(1) or not (0)


%Initialize rng and number of processors to be used
maxNumCompThreads(num_processors);
rng(rng_seed);



%Initialize feasible solutions
d = size(bound_constraints,1);

X = [repmat(linspace(bound_constraints(1,1),bound_constraints(1,2),bound_constraints(1,2)-bound_constraints(1,1)+1)',prod(bound_constraints(2:end,2)-bound_constraints(2:end,1)+1),1), zeros(prod(bound_constraints(:,2)-bound_constraints(:,1)+1),d-1)];
for i = 2:d
X(:,i) = repmat(repelem(linspace(bound_constraints(i,1),bound_constraints(i,2),bound_constraints(i,2)-bound_constraints(i,1)+1)',prod(bound_constraints(1:i-1,2)-bound_constraints(1:i-1,1)+1,1)),prod(bound_constraints(i+1:end,2)-bound_constraints(i+1:end,1)+1,1));
end

n = size(X,1);


%Initialize sample statistics
all_cumulative_reps = zeros(n,1);
all_simulation_output = zeros(n,1);
all_simulation_output_variances = zeros(n,1);





%Select design points following LHS methodology
design_points = ceil((repmat((bound_constraints(:,2)-bound_constraints(:,1))'+1,k,1)).*lhsdesign(k,d,'criterion','maximin')) + repmat(bound_constraints(:,1)'-1,k,1);
design_points = sortrows(design_points);
design_points_index = cumprod([1;bound_constraints(:,2)-bound_constraints(:,1) + 1]');
design_points_index = sum(((design_points-repmat(bound_constraints(:,1)',k,1))+1-[zeros(size(design_points,1),1) ones(size(design_points,1),size(design_points,2)-1)]).*repmat(design_points_index(1:end-1),size(design_points,1),1),2);

%Collect simulation output and compute sample statistics
simulation_output = simulator(design_points);
simulation_output_variances = var(simulation_output,0,2);
response_sum(design_points_index) = sum(simulation_output,2);
response_sum_of_square_differences(design_points_index) = rep(1)*var(simulation_output,1,2);






[beta, theta] = maximum_likelihood_estimation(all_simulation_output(design_points_index),response_sum_of_square_differences/rep(1),bound_constraints,design_points_index,sobol_bounds, pardiso_flag);


eDO = ExperimentalDesign(bound_constraints,sparse(response_sum),sparse(response_sum_of_square_differences),sparse(all_cumulative_reps),design_points_index,beta,theta);



save(strcat(experimental_design_object_filename,'.mat'),'eDO');


end

