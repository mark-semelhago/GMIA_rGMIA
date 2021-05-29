%Sample run file

%Define simulator
simulator = 'restaurant';
simulator = str2func(simulator);
c = [1 3]';
T = 1;
burnin = 0;
lambda = 100*ones(3,1);
mu = 10*ones(3,1);
rev = 0.01*linspace(1,3,3)';
l = 0.005*[1 3]';
simulator = @(dp, r) simulator(dp, c, T, burnin, lambda, mu, rev, l, r);
user_input.simulator = simulator;

%Construct ExperimentalDesignObject
LHS_parameter_estimation(simulator, [1, 20;1, 20], 1, 1, 20,[10, 2], [0.2,0.11], 'eDO', 0)

%Define user_input struct fields
user_input.simulator = simulator;
user_input.replications = [10 2]';
user_input.max_num_iterations = 500;
user_input.time_budget = Inf;
user_input.max_cei_threshold = 0;
user_input.pardiso_flag = 0;
user_input.search_set_size = [];
user_input.num_screening_solutions = [];
user_input.num_simulation_workers = 1;
user_input.verbose = 1;
user_input.num_rapid_search_iterations = [];
user_input.bound_constraints = [1 10;1 100];

%Define run_info struct fields
run_info = struct;
run_info.max_cei_solutions = [];
run_info.current_optimal_solutions = [];
run_info.max_ceis = [];
run_info.iteration_timings = [];


%Run GMIA
run_info = GMIA(user_input, eDO, run_info);

%Save results
save('restaurant_simulation_GMIA_results.mat');
