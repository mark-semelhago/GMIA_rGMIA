function run_info = rGMIA(user_input, eD, run_info)
%rGMIA is the main function that runs the rGMIA algorithm (see GMIA for
%running GMIA). It takes as input a struct, user_input, to specify run
%parameters; eD, an ExperimentalDesign object (see class definition) and a
%struct, run_info, to store runtime statistics.
%user_input has the following fields:
%simulator: function handle associated with simulation that only takes in
%the solution to simulate and number of simulation replications
%replications: 2x1 vector where the first element is number of reps
%performed upon first visiting solution, and second is number of reps on
%subsequent visits
%max_num_iterations: maximum number of iterations of GMIA to be run
%time_budget: maximum amount of time (in sec) for GMIA to be run (evaluated
%at the beginning of each iteration
%max_cei_threshold: condition to terminate when max CEI value falls below
%threshold
%pardiso_flag: binary (1/0) flag whether to use PARDISO library or not
%search_set_size: number of solutions to be included in search set
%num_simulation_workers: number of processors available for parallel use
%verbose: binary (1/0) flag whether to output iteration number, estimated
%objective function value, max CEI value at iteration and total time
%elapsed at each iteration
%num_rapid_search_iteration: number of rapid search iterations to perform
%before performing global search iteration
%dense dx2 matrix containing inclusive lower and upper-bounds of box
%constraints for each dimension
%run_info has the following fields:
%max_cei_solutions: linear index of max CEI solution at each iteration,
%assuming column-major indexing
%current_optimal_solutions: linear index of solution with lowest sample
%mean at each iteration, assuming column-major indexing
%max_ceis: max CEI value at each iteration
%iteration_timings: time (in seconds) to perform each iteration of GMIA


%Set number of parallel processors to use
maxNumCompThreads(user_input.num_simulation_workers);

%Initialize global metamodel object
gM = GlobalMetamodel(eD, user_input);

iter = 1;

tic;


%Find current optimal solution
gM.compute_current_optimal_solution();


%Compute conditional statistics (means, variances, covariances) for all
%solutions
gM.compute_conditional_statistics_non_partitioned();


%Compute CEI values for all solutions
gM.compute_complete_expected_improvement();

%Compute max CEI solution if single processor is used
if user_input.num_simulation_workers == 1
    gM.compute_max_cei_solution();
end

%Partition feasible region into search set and fixed set
gM.construct_partition();

%Construct rapid metamodel object based on constructed partition
rM = RapidMetamodel(gM);

%Compute conditional statistics for solutions in the rapid metamodel object
%(only solutions in the search set)
rM.compute_conditional_statistics();

%Find max q-CEI solutions if more than 1 processor available
if user_input.num_simulation_workers > 1
    	myCluster = parcluster('local');
	myCluster.NumWorkers = user_input.num_simulation_workers;
	parpool('local',user_input.num_simulation_workers);
	rM.compute_max_q_cei_solutions();    
end


%Record run information for first iteration in run_info
if user_input.num_simulation_workers == 1
    run_info.max_cei_solutions = [run_info.max_cei_solutions; gM.max_cei_solutions];
    
else
    run_info.max_cei_solutions = [run_info.max_cei_solutions; rM.max_cei_solutions'];
end
run_info.max_ceis = [run_info.max_ceis; max(gM.complete_expected_improvement)];
run_info.current_optimal_solutions = [run_info.current_optimal_solutions; gM.current_optimal_solution];
run_info.iteration_timings = [run_info.iteration_timings; toc];

if (user_input.verbose)
    fprintf('Iteration: %d; Est Obj Fcn Val: %f; Max CEI: %f; Time Elapsed: %f\n',iter, full(gM.response_sample_means(gM.current_optimal_solution)),...
        max((gM.complete_expected_improvement(rM.max_cei_solutions))), sum(run_info.iteration_timings));
end

%Termination conditions are based on max CEI value, number of iterations
%and time budget
while (iter==1) || ((max(gM.complete_expected_improvement) > user_input.max_cei_threshold) && ...
        (iter <= user_input.max_num_iterations) && (sum(run_info.iteration_timings) <= user_input.time_budget))
    
    
    
    %Begin rapid search iterations for fixed number of iterations
    for rapid_iter = 1:user_input.num_rapid_search_iterations
        iter = iter + 1;
        tic;
        
        %Find linear index of current optimal solution relative to ordering
        %in search set
        current_optimal_solution_search_set_index = find(rM.search_set == rM.current_optimal_solution);
        
        
            
        %Find linear indices of max CEI solutions relative to ordering in
        %search set
        max_cei_solutions_search_set_indices = find(ismember(rM.search_set, rM.max_cei_solutions));
        
        %Specify solutions to simualte and number of simulation
        %replications to be performed
        simulated_solutions = [current_optimal_solution_search_set_index; max_cei_solutions_search_set_indices];
        simulated_solutions_replications = user_input.replications(1) .* (rM.num_replications(simulated_solutions)==0) + ...
            user_input.replications(2) .* (rM.num_replications(simulated_solutions) > 0);
        simulated_solutions_coordinates = index_to_coordinate(gM.bound_constraints, [rM.current_optimal_solution; rM.max_cei_solutions]);
        
        
        %Perform simulation at simulated_solutions and store output
        if user_input.num_simulation_workers == 1

            simulation_output = cell(2, 1);
            simulation_output{1} = user_input.simulator(simulated_solutions_coordinates(1,:), simulated_solutions_replications(1));
            simulation_output{2} = user_input.simulator(simulated_solutions_coordinates(2,:), simulated_solutions_replications(2));
        else
            simulation_output = cell(user_input.num_simulation_workers, 1);
            parfor (i = 1:user_input.num_simulation_workers, user_input.num_simulation_workers)
                simulation_output{i} = user_input.simulator(simulated_solutions_coordinates(i,:),simulated_solutions_replications(i));
            end
        end
        
        %Update sample statistics for simulated solutions
        for i = 1:length(simulated_solutions)
            sim_out = simulation_output{i}';
            if rM.num_replications(simulated_solutions(i))>0
                running_cumulative_means = cumsum([rM.response_sum(simulated_solutions(i));sim_out])./cumsum([rM.num_replications(simulated_solutions(i)); ones(length(sim_out),1)]);
                rM.response_sum_of_square_differences(simulated_solutions(i)) = rM.response_sum_of_square_differences(simulated_solutions(i)) + ...
                sum((sim_out-running_cumulative_means(1:end-1)).*(sim_out-running_cumulative_means(2:end)));
            else
                running_cumulative_means = cumsum(sim_out)./linspace(1,length(sim_out),length(sim_out))';
                rM.response_sum_of_square_differences(simulated_solutions(i)) = sum((sim_out(2:end)-running_cumulative_means(1:end-1)).*(sim_out(2:end)-running_cumulative_means(2:end)));
            end
            
            rM.response_sum(simulated_solutions(i)) = rM.response_sum(simulated_solutions(i)) + sum(simulation_output{i});
            rM.num_replications(simulated_solutions(i)) = rM.num_replications(simulated_solutions(i)) + simulated_solutions_replications(i);
            
            
        end
        
            
        
               
        %Update all conditional precision matrix of for rapid metamodel
        %object
        rM.update_precision_matrices();
        
        %Find solution with lowest sample mean in rapid metamodel object
        rM.compute_current_optimal_solution();
        
        %Compute conditional statistics for solutions in rapid metamodel
        %object
        rM.compute_conditional_statistics();
        
        
        %Compute max CEI solution, or max q-CEI solutions if there is more
        %than 1 processor
        rM.compute_complete_expected_improvement();
        if user_input.num_simulation_workers == 1
            rM.compute_max_cei_solution();
        else
            rM.compute_max_q_cei_solutions();
        end
        
        
        
        %Record run information in run_info
        run_info.max_cei_solutions = [run_info.max_cei_solutions; rM.max_cei_solutions'];
        run_info.current_optimal_solutions = [run_info.current_optimal_solutions; rM.current_optimal_solution];
        run_info.max_ceis = [run_info.max_ceis; max(rM.complete_expected_improvement)];
        run_info.iteration_timings = [run_info.iteration_timings; toc];
        
        if (user_input.verbose)
            fprintf('Iteration: %d; Est Obj Fcn Val: %f; Max CEI: %f; Time Elapsed: %f\n',iter,full(rM.response_sample_means(current_optimal_solution_search_set_index)),max(rM.complete_expected_improvement),...
                sum(run_info.iteration_timings));
        end
        
    end
    
    %Beginning of global search iteration
    iter = iter + 1;

    
    %Find index of solution with lowest sample mean and max CEI solutions
    %relative to indexing of search set
    current_optimal_solution_search_set_index = find(rM.search_set == rM.current_optimal_solution);
    max_cei_solutions_search_set_indices = find(ismember(rM.search_set, rM.max_cei_solutions));
    
    %Specify solutions to simulate and number of simulation replications to
    %be performed
    simulated_solutions = [current_optimal_solution_search_set_index; max_cei_solutions_search_set_indices];
        simulated_solutions_replications = user_input.replications(1) .* (rM.num_replications(simulated_solutions)==0) + ...
            user_input.replications(2) .* (rM.num_replications(simulated_solutions) > 0);
        simulated_solutions_coordinates = index_to_coordinate(gM.bound_constraints, [rM.current_optimal_solution; rM.max_cei_solutions]);
        
        
    %Perform simulation at simulated_solutions and store output
    if user_input.num_simulation_workers == 1
        simulation_output = cell(2, 1);
        simulation_output{1} = user_input.simulator(simulated_solutions_coordinates(1,:), simulated_solutions_replications(1));
        simulation_output{2} = user_input.simulator(simulated_solutions_coordinates(2,:), simulated_solutions_replications(2));
    else
        simulation_output = cell(user_input.num_simulation_workers, 1);
        parfor (i = 1:user_input.num_simulation_workers, user_input.num_simulation_workers)

            simulation_output{i} = user_input.simulator(simulated_solutions_coordinates(i,:),simulated_solutions_replications(i));
        end
    end
    
    %Update sample statistics for simulated solutions
    for i = 1:length(simulated_solutions)
            sim_out = simulation_output{i}';
            if rM.num_replications(simulated_solutions(i))>0
                running_cumulative_means = cumsum([rM.response_sum(simulated_solutions(i));sim_out])./cumsum([rM.num_replications(simulated_solutions(i)); ones(length(sim_out),1)]);
                rM.response_sum_of_square_differences(simulated_solutions(i)) = rM.response_sum_of_square_differences(simulated_solutions(i)) + ...
                sum((sim_out-running_cumulative_means(1:end-1)).*(sim_out-running_cumulative_means(2:end)));
            else
                running_cumulative_means = cumsum(sim_out)./linspace(1,length(sim_out),length(sim_out))';
                rM.response_sum_of_square_differences(simulated_solutions(i)) = sum((sim_out(2:end)-running_cumulative_means(1:end-1)).*(sim_out(2:end)-running_cumulative_means(2:end)));
            end
            
            rM.response_sum(simulated_solutions(i)) = rM.response_sum(simulated_solutions(i)) + sum(simulation_output{i});
            rM.num_replications(simulated_solutions(i)) = rM.num_replications(simulated_solutions(i)) + simulated_solutions_replications(i);
            
    end
    
    


    
    %Update conditional precision matrices for solutions in rapid metamodel
    %object
    rM.update_precision_matrices();
    
    %Update global metamodel object with new simulation information gained
    %from rapid metamodel object
    update_global_sample_mean_precision_matrices_from_search_set(gM, rM);
    
    %Find solution with lowest sample mean across entire feasible region
    gM.compute_current_optimal_solution();
    
    %Compute conditional statistics for solutions in search set, given
    %partition
    compute_conditional_statistics_search_set(gM, rM);
    
    %Compute conditional statistics for solutions in fixed set, given
    %partition
    compute_conditional_statistics_fixed_set(gM, rM);
    
    %Compute CEI values for all solutions in feasible region jointly
    gM.compute_complete_expected_improvement();
    
    
    %Construct new partition of feasible region into search and fixed sets
    gM.construct_partition();

    %Using partiion, construct new rapid metamodel object
    rM = RapidMetamodel(gM);
    
    %Compute conditional statistics for solutions in rapid metamodel object
    rM.compute_conditional_statistics();
    
    %Find max CEI solution(s) out of solutions in rapid metamodel object
    if user_input.num_simulation_workers == 1
    
        rM.compute_max_cei_solution();
    else
        rM.compute_max_q_cei_solutions();
    end


    %Record run information in run_info
    if user_input.num_simulation_workers == 1
        run_info.max_cei_solutions = [run_info.max_cei_solutions; rM.max_cei_solutions];
    else
        run_info.max_cei_solutions = [run_info.max_cei_solutions; rM.max_cei_solutions'];
    end
    run_info.current_optimal_solutions = [run_info.current_optimal_solutions; gM.current_optimal_solution];
    run_info.max_ceis = [run_info.max_ceis; max(gM.complete_expected_improvement)];
    run_info.iteration_timings = [run_info.iteration_timings; toc];
    
    
    if (user_input.verbose)
            fprintf('Iteration: %d; Est Obj Fcn Val: %f; Max CEI: %f; Time Elapsed: %f\n',iter,full(gM.response_sample_means(gM.current_optimal_solution)),...
                max(gM.complete_expected_improvement), sum(run_info.iteration_timings));
    end
    
    
    
    

    
    
end

    
    
    
end
