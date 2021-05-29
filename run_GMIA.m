function run_info = run_GMIA(it)

global pardiso1 id1
id1 = 1;
pardiso1 = pardiso(id1,-2);

user_input = struct;
simulator = 'gmrf';
simulator = str2func(simulator);
grid_dimensions = [401 401]';
load('griewank_401_trend.mat');
Z = reshape(Z,401^2,1);
var_noise = 0.0001*ones(401^2,1);
simulator = @(dp, r) simulator(dp, grid_dimensions,Z,var_noise, r);


user_input.simulator = simulator;
user_input.replications = [10 2]';
user_input.max_num_iterations = 500;
user_input.time_budget = Inf;
user_input.max_cei_threshold = 0;
user_input.pardiso_flag = 1;
user_input.search_set_size = [];
user_input.num_screening_solutions = [];
user_input.num_simulation_workers = 1;
user_input.verbose = 1;
user_input.num_rapid_search_iterations = [];
user_input.bound_constraints = [1 401;1 401];


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

load(strcat('griewank_401_trendit',num2str(it),'.mat'));


run_info = GMIA(user_input, eDO, run_info);

save(strcat('run_info_griewank_401_trendit',num2str(it),'.mat'));

end
