function run_info = GMIA(user_input, eD, run_info)
%GMIA is the main function that runs the GMIA algorithm (see rGMIA for
%running rGMIA). It takes as input a struct, user_input, to specify run
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
%search_set_size: leave empty, [], for GMIA use
%num_simulation_workers: number of processors available for parallel use
%verbose: binary (1/0) flag whether to output iteration number, estimated
%objective function value, max CEI value at iteration and total time
%elapsed at each iteration
%num_rapid_search_iteration: leave empty, [], for GMIA use
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

%Termination conditions are based on max CEI value, number of iterations
%and time budget
    while (iter==1) || ((max(gM.complete_expected_improvement(gM.max_cei_solutions)) > user_input.max_cei_threshold) && ...
            (iter <= user_input.max_num_iterations) && (sum(run_info.iteration_timings) <= user_input.time_budget))

        tic;

        %Find current optimal solution
        gM.compute_current_optimal_solution();
        
        %Compute conditional statistics (means, variances, covariances) for
        %all solutions
        gM.compute_conditional_statistics_non_partitioned();

        %Compute CEI values for all solutions
        gM.compute_complete_expected_improvement();
        
        %Compute either max CEI solution if single processor is used or q
        %solutions that maximize q-CEI if q > 1 processors are used.
        if user_input.num_simulation_workers == 1
            gM.compute_max_cei_solution();
        else 
            gM.compute_max_q_cei_solutions();
        end
        
        %Specify solutions to simulate and number of simulation
        %replications to be performed
        simulated_solutions = [gM.current_optimal_solution; gM.max_cei_solutions];
        simulated_solutions_replications = user_input.replications(1) .* (gM.num_replications(simulated_solutions)==0) + ...
            user_input.replications(2) .* (gM.num_replications(simulated_solutions) > 0);
        simulated_solutions_coordinates = index_to_coordinate(gM.bound_constraints, simulated_solutions);
        
        
        
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
            if gM.num_replications(simulated_solutions(i))>0
                running_cumulative_means = cumsum([gM.response_sum(simulated_solutions(i));sim_out])./cumsum([gM.num_replications(simulated_solutions(i)); ones(length(sim_out),1)]);
                gM.response_sum_of_square_differences(simulated_solutions(i)) = gM.response_sum_of_square_differences(simulated_solutions(i)) + ...
                sum((sim_out-running_cumulative_means(1:end-1)).*(sim_out-running_cumulative_means(2:end)));
            else
                running_cumulative_means = cumsum(sim_out)./linspace(1,length(sim_out),length(sim_out))';
                gM.response_sum_of_square_differences(simulated_solutions(i)) = sum((sim_out(2:end)-running_cumulative_means(1:end-1)).*(sim_out(2:end)-running_cumulative_means(2:end)));
            end
            
            gM.response_sum(simulated_solutions(i)) = gM.response_sum(simulated_solutions(i)) + sum(simulation_output{i});
            gM.num_replications(simulated_solutions(i)) = gM.num_replications(simulated_solutions(i)) + simulated_solutions_replications(i);
            
        end
        
            
        
            


        %Update conditional precision matrix
        gM.update_conditional_precision_matrix_non_partitioned();



        %Record run information in run_info
        run_info.max_cei_solutions = [run_info.max_cei_solutions; gM.max_cei_solutions'];
        run_info.current_optimal_solutions = [run_info.current_optimal_solutions; gM.current_optimal_solution];
        run_info.max_ceis = [run_info.max_ceis; max(gM.complete_expected_improvement(gM.max_cei_solutions))];
        run_info.iteration_timings = [run_info.iteration_timings; toc];
        
        
        if (user_input.verbose)
            fprintf('Iteration: %d; Est Obj Fcn Val: %f; Max CEI: %f; Time Elapsed: %f\n',iter,full(gM.response_sample_means(gM.current_optimal_solution)),...
                max(gM.complete_expected_improvement(gM.max_cei_solutions)), sum(run_info.iteration_timings));
        end
        
        
        iter = iter + 1;



    end

end

    
    
    
    
    
    
    
    

