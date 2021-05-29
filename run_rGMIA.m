function run_info = run_rGMIA(it, N, p, q_par, max_iterations)
global id1 id2 id3 pardiso1 pardiso2 pardiso3;

rng(it);

id1 = 1;
id2 = 2;
id3 = 3;

pardiso1 = pardiso(id1,-2);
pardiso2 = pardiso(id2,-2);
pardiso3 = pardiso(id3,-2);

user_input = struct;
simulator = 'restaurant';
simulator = str2func(simulator);
c = [1 3]';
T = 1;
burnin = 0;
lambda = 2000*ones(3,1);
mu = 10*ones(3,1);
rev = 0.01*linspace(1,3,3)';
l = 0.005*[1 3]';
simulator = @(dp, r) simulator(dp, c, T, burnin, lambda, mu, rev, l, r);
user_input.simulator = simulator;


user_input.bound_constraints = [1 1000;1 1000];
user_input.replications = [10 2]';
user_input.max_num_iterations = 120/(q_par+1);
user_input.time_budget = Inf;
user_input.max_cei_threshold = 0;
user_input.pardiso_flag = 1;
user_input.search_set_size = 200;
user_input.num_screening_solutions = N;
user_input.num_simulation_workers = q_par+1;
user_input.verbose = 1;
user_input.num_rapid_search_iterations = p-1;


%Define run_info
run_info = struct;
run_info.max_cei_solutions = [];
run_info.current_optimal_solutions = [];
run_info.max_ceis = [];
run_info.iteration_timings = [];
%run_info.bound_constraints = user_input.bound_constraints;
%run_info.design_points = [];
%run_info.replications = user_input.replications;
%run_info.num_rapid_search_iterations = user_input.num_rapid_search_iterations;
%run_info.num_parallel_workers = user_input.num_parallel_workers;
%run_info.max_num_iterations = user_input.max_num_iterations;
%run_info.time_budget = user_input.time_budget;
%run_info.max_cei_threshold = user_input.max_cei_threshold;
%run_info.mean_MLE = user_input.mean_MLE;
%run_info.precision_MLE = user_input.precision_MLE;
%run_info.iteration_timing = [];
%run_info.current_optimal_solutions = [];
%run_info.max_cei_solutions = [];
%run_info.parallel_solutions = [];
%run_info.response_sample_means = [];
%run_info.response_sample_variances = [];
%run_info.max_cei = [];
%run_info.search_sets = [];
%run_info.random_seed = [];
%run_info.pardiso_flag = [];

%Estimate GMRF MLEs if not already estimated
%if isempty(mean_MLE) && isempty(precision_MLE)
%    [GlobalMetamodelInformation.mean_MLE, GlobalMetamodelInformation.precision_MLE, GlobalMetamodelInformation.design_points] = MLE_estimation(user_input);
%    run_info.mean_MLE = GlobalMetamodelInformation.mean_MLE;
%    run_info.precision_MLE = GlobalMetamodelInformation.precision_MLE;
%    run_info.design_points = GlobalMetamodelInformation.design_points;
%elseif isempty(mean_MLE) && ~isempty(precision_MLE)
%    error('Missing mean_MLE information.');
%elseif ~isempty(mean_MLE) && isempty(precision_MLE)
%    error('Missing precision_MLE information.');
%end




%Create PARDISO objects
%if user_input.pardiso_flag == 1
%    global id1 pardiso1;
%    id1 = 1;
%    pardiso1 = pardiso(id1,-2);
%    if ~isempty(user_input.search_set_size)
%        global id2 pardiso2;
%        id2 = 2;
%        pardiso2 = pardiso(id2,-2);
%    end
%end


%Construct precision matrix
%precision_matrix = construct_precision_matrix(precision_MLE, bound_constraints);


load(strcat('restaurant1000_it',num2str(1),'.mat'));


run_info = rGMIA(user_input, eDO, run_info);

save(strcat('restaurant1000_it',num2str(it),'N',num2str(N),'p',num2str(p),'q',num2str(q_par),'run_info.mat'));

end
