classdef GlobalMetamodel < Metamodel
    %GLOBALMETAMODEL subclass which inherits from Metamodel.
    %   The GlobalMetamodel class is a subclass of Metamodel that contains
    %all information pertaining to the GMRF of all solutions in GMIA and
    %rGMIA, assuming there are n solutions and d dimensions.
    properties (SetAccess = protected)
%         Properties only set accessible form class and subclass methods

%         Inherited properties from Metamodel
%         precision_matrix
%         intrinsic_precision_matrix
%         conditional_precision_matrix
%         conditional_means 
%         conditional_variances
%         conditional_covariances
%         complete_expected_improvement 
%         current_optimal_solution 
%         max_cei_solutions
%         num_solutions
%         search_set
%         MLE_mean
%         MLE_precision
%         num_simulation_workers
%         num_screening_solutions
%         pardiso_flag
%         response_sample_means (Dependent)
%         response_sample_variances (Dependent)
        
%         GlobalMetamodel specific properties
        bound_constraints %Dense dx2 matrix containing inclusive lower and upper-bounds of box constraints for each dimension.
        design_points %Dense kx1 vector containing linear index of design points, assuming column-major ordering.
        search_set_size %Scalar value indicating the number of solutions chosen to be in the search set. This is left empty in GMIA.       
        
    end
    
    methods (Access = private)

        function construct_precision_matrix(globalObj)
            %Given the MLE_mean and MLE_precision, construct the associated precision matrix and store in GlobalMetamodel.precision_matrix.
            grid_dimensions = globalObj.bound_constraints(:,2) - (globalObj.bound_constraints(:,1)-1);
            dim = length(grid_dimensions);
            
            row = [globalObj.MLE_precision(1) -globalObj.MLE_precision(1) * globalObj.MLE_precision(2)];
            for i = 3:length(globalObj.MLE_precision)
                row = [row zeros(1,prod(grid_dimensions(1:i-2)) - length(row)) -globalObj.MLE_precision(1)*globalObj.MLE_precision(i)];
            end
            row = [row zeros(1,globalObj.num_solutions-length(row))];
            prec_mat = sptoeplitz(row);
            
            offset = find(row~=0);
            offset = offset(2:end)-1;
            for i = 1:dim-1
                end_zeroed_rows = prod(grid_dimensions(1:i)):prod(grid_dimensions(1:i)):(size(prec_mat,1) - 1);
                if i == 1
                    zeroed_rows = end_zeroed_rows;
                else
                    begin_zeroed_rows = end_zeroed_rows - prod(grid_dimensions(1:i-1)) + 1;
                    zeroed_rows = [];
                    for j = 1:length(end_zeroed_rows)
                        zeroed_rows = [zeroed_rows, begin_zeroed_rows(j):end_zeroed_rows(j)];
                    end
                end
                prec_mat(zeroed_rows, zeroed_rows + offset(i)) = 0;
            end
            globalObj.precision_matrix = triu(prec_mat) + triu(prec_mat)' - spdiags(diag(prec_mat),0,prod(grid_dimensions),prod(grid_dimensions));
        end
    end
    
    
    methods
        
        function globalObj = GlobalMetamodel(experimentalDesignObj, user_input)
        %GlobalMetamodel constructor    
            globalObj.bound_constraints = experimentalDesignObj.bound_constraints;
            grid_dimensions = globalObj.bound_constraints(:,2) - (globalObj.bound_constraints(:,1)-1);
            globalObj.num_solutions = prod(grid_dimensions);
            globalObj.response_sum = experimentalDesignObj.response_sum;
            globalObj.response_sum_of_square_differences = experimentalDesignObj.response_sum_of_square_differences;
            globalObj.num_replications = experimentalDesignObj.num_replications;
            globalObj.design_points = experimentalDesignObj.design_points;
            globalObj.MLE_mean = experimentalDesignObj.MLE_mean;
            globalObj.MLE_precision = experimentalDesignObj.MLE_precision;
            globalObj.search_set_size = user_input.search_set_size;
            globalObj.pardiso_flag = user_input.pardiso_flag;
            globalObj.num_screening_solutions = user_input.num_screening_solutions;
            globalObj.num_simulation_workers = user_input.num_simulation_workers;
            
            
            globalObj.construct_precision_matrix();
            intrinsic_precision_matrix_diag = sparse(zeros(globalObj.num_solutions,1));
            intrinsic_precision_matrix_diag(globalObj.design_points) = globalObj.num_replications(globalObj.design_points)./globalObj.response_sample_variances(globalObj.design_points);
            globalObj.intrinsic_precision_matrix = spdiags(intrinsic_precision_matrix_diag,0,globalObj.num_solutions,globalObj.num_solutions);
            globalObj.conditional_precision_matrix = globalObj.precision_matrix + globalObj.intrinsic_precision_matrix;
            
          
        end
        
        
        function compute_max_q_cei_solutions(globalObj)
            %Compute the max_q_cei_solutions and store as vector in GlobalMetamodel.max_cei_solutions.
            global pardiso3 id3
            q = globalObj.num_simulation_workers - 1;
            screened_solutions = sortrows([globalObj.complete_expected_improvement linspace(1,globalObj.num_solutions,globalObj.num_solutions)'], -1);
            screened_solutions = [globalObj.current_optimal_solution; screened_solutions(1:globalObj.num_screening_solutions,2)];
            unscreened_solutions = setdiff(linspace(1,globalObj.num_solutions, globalObj.num_solutions)', screened_solutions);
            
            
            pardiso3.factorize(id3,tril(globalObj.conditional_precision_matrix(unscreened_solutions, unscreened_solutions)));
            
            screened_conditional_covariance_matrix = inv(globalObj.conditional_precision_matrix(screened_solutions, screened_solutions) -...
                globalObj.conditional_precision_matrix(screened_solutions, unscreened_solutions) *...
                pardiso3.solve(id3, tril(globalObj.conditional_precision_matrix(unscreened_solutions, unscreened_solutions)),...
                full(globalObj.conditional_precision_matrix(unscreened_solutions, screened_solutions))));
            
            
            max_q_cei_set = screened_solutions(2);
            max_q_cei_set_opt = [screened_solutions(1); max_q_cei_set];
            
            
            
            for i = 2:q
                candidate_solutions = setdiff(screened_solutions, max_q_cei_set_opt);
                q_cei_cand = NaN(length(candidate_solutions),1);
                parfor (i = 1:length(candidate_solutions), globalObj.num_simulation_workers)
                    cand_max_q_cei_set = [max_q_cei_set; candidate_solutions(i)];
                    cand_max_q_cei_set_opt = [max_q_cei_set_opt; candidate_solutions(i)];
                    cand_max_q_cei_set_opt_screened_ind = find(ismember(screened_solutions, cand_max_q_cei_set_opt));
                    q_cei_cand(i) = compute_q_cei(globalObj.conditional_means(cand_max_q_cei_set_opt), ...
                       screened_conditional_covariance_matrix(cand_max_q_cei_set_opt_screened_ind, cand_max_q_cei_set_opt_screened_ind));
                end
                [~, max_q_cei_cand_solution] = max(q_cei_cand);
                max_q_cei_set = [max_q_cei_set; candidate_solutions(max_q_cei_cand_solution)];
                max_q_cei_set_opt = [max_q_cei_set_opt; candidate_solutions(max_q_cei_cand_solution)];
            end
            
            globalObj.max_cei_solutions = max_q_cei_set;
            
            
            
     
        end
        
        
        
        
        
        
        function compute_current_optimal_solution(globalObj)
        %Find the solution with the lowest response sample mean (since we assume problem is a minimization problem) and store in GlobalMetamodel.max_cei_solutions.
            sample_means = [globalObj.response_sample_means linspace(1,length(globalObj.response_sample_means),length(globalObj.response_sample_means))'];
            sample_means = sample_means(sample_means(:,1)~=0,:);
            [~, index_current_optimal_solution] = min(sample_means(:,1));
            globalObj.current_optimal_solution = full(sample_means(index_current_optimal_solution,2));
        end
        
        function compute_max_cei_solution(globalObj)
        %Find single max CEI solution. This is called when there is 1 processor available.
            [~, globalObj.max_cei_solutions] = max(globalObj.complete_expected_improvement);
        end
        
        
        function compute_conditional_statistics_non_partitioned(globalObj)
        %Compute conditional statistics for solutions in search set and fixed set separately. This is only used in rGMIA and not in GMIA.
            if globalObj.pardiso_flag
                global id1 pardiso1
                pardiso1.factorize(id1,tril(globalObj.conditional_precision_matrix))
                globalObj.conditional_variances = full(diag(pardiso1.invert(id1,tril(globalObj.conditional_precision_matrix))));
                globalObj.conditional_means = globalObj.MLE_mean + pardiso1.solve(id1,tril(globalObj.conditional_precision_matrix),full(globalObj.intrinsic_precision_matrix*(globalObj.response_sample_means-globalObj.MLE_mean)));
                globalObj.conditional_covariances = pardiso1.solve(id1,tril(globalObj.conditional_precision_matrix),double(full((1:globalObj.num_solutions == globalObj.current_optimal_solution)')));
            else
                conditional_precision_matrix_cholesky_factor = chol(globalObj.conditional_precision_matrix)';
                globalObj.conditional_variances = full(diag(conditional_precision_matrix_cholesky_factor'\...
                    (conditional_precision_matrix_cholesky_factor\speye(globalObj.num_solutions,globalObj.num_solutions))));
                globalObj.conditional_means = full(globalObj.MLE_mean + conditional_precision_matrix_cholesky_factor'\...
                    (conditional_precision_matrix_cholesky_factor\(globalObj.intrinsic_precision_matrix*(globalObj.response_sample_means-globalObj.MLE_mean))));
                globalObj.conditional_covariances = conditional_precision_matrix_cholesky_factor'\...
                    (conditional_precision_matrix_cholesky_factor\(1:globalObj.num_solutions == globalObj.current_optimal_solution)');
            end
                
        end

        
        
        function update_conditional_precision_matrix_non_partitioned(globalObj)
        %Update conditional precision matrix that in a non-partitioned manner. This function is used in GMIA when there is no search set.
            globalObj.intrinsic_precision_matrix(globalObj.current_optimal_solution, globalObj.current_optimal_solution) = globalObj.num_replications(globalObj.current_optimal_solution)/...
                globalObj.response_sample_variances(globalObj.current_optimal_solution);
           globalObj.conditional_precision_matrix(globalObj.current_optimal_solution, globalObj.current_optimal_solution) = ...
                globalObj.precision_matrix(globalObj.current_optimal_solution, globalObj.current_optimal_solution) + ...
                globalObj.intrinsic_precision_matrix(globalObj.current_optimal_solution, globalObj.current_optimal_solution);
            
            
            max_cei_solution_mat_lin_index = sub2ind(size(globalObj.intrinsic_precision_matrix),globalObj.max_cei_solutions, globalObj.max_cei_solutions);
            
            globalObj.intrinsic_precision_matrix(max_cei_solution_mat_lin_index) = globalObj.num_replications(globalObj.max_cei_solutions)./...
                globalObj.response_sample_variances(globalObj.max_cei_solutions);
            globalObj.conditional_precision_matrix(max_cei_solution_mat_lin_index) = globalObj.precision_matrix(max_cei_solution_mat_lin_index) + ...
                globalObj.intrinsic_precision_matrix(max_cei_solution_mat_lin_index);
            
        end
        
        
        function update_global_sample_mean_precision_matrices_from_search_set(globalObj, rapidObj)
        %Update the vectors, matrices storing global sample statistics (for all solutions), after simulations of solutions in the search set have finished.
            
            globalObj.response_sum(globalObj.search_set) = rapidObj.response_sum;
            globalObj.response_sum_of_square_differences(globalObj.search_set) = rapidObj.response_sum_of_square_differences;
            globalObj.num_replications(globalObj.search_set) = rapidObj.num_replications;
            
            global_intrinsic_precision_matrix_diag = diag(globalObj.intrinsic_precision_matrix);
            rapid_intrinsic_precision_matrix_diag = diag(rapidObj.intrinsic_precision_matrix);
            global_intrinsic_precision_matrix_diag(globalObj.search_set) = rapid_intrinsic_precision_matrix_diag;
            
            globalObj.intrinsic_precision_matrix = diag(global_intrinsic_precision_matrix_diag);
            globalObj.conditional_precision_matrix = globalObj.precision_matrix + globalObj.intrinsic_precision_matrix;
            
        end
      
        
        function compute_conditional_statistics_search_set(globalObj, rapidObj)
            %Compute conditional statistics for solutions in the search set
            
            fixed_set = setdiff(linspace(1,globalObj.num_solutions,globalObj.num_solutions)',globalObj.search_set);
            
            
            rapidObj.conditional_covariance_matrix = inv(rapidObj.conditional_precision_matrix - rapidObj.intermediate_matrix_B);
            rapidObj.conditional_means = rapidObj.MLE_mean + rapidObj.conditional_covariance_matrix*(rapidObj.intrinsic_precision_matrix*(rapidObj.response_sample_means - ...
                rapidObj.MLE_mean) - rapidObj.intermediate_vector_a);
            rapidObj.conditional_variances = diag(rapidObj.conditional_covariance_matrix);
            
            if ismember(globalObj.current_optimal_solution, globalObj.search_set)
                current_optimal_solution_search_set_index = find(globalObj.search_set == globalObj.current_optimal_solution);
                rapidObj.conditional_covariances = rapidObj.conditional_covariance_matrix(:, current_optimal_solution_search_set_index);
            else
                current_optimal_solution_fixed_set_index = find(fixed_set == globalObj.current_optimal_solution);
                rapidObj.conditional_covariances = -rapidObj.conditional_covariance_matrix*rapidObj.intermediate_matrix_A(current_optimal_solution_fixed_set_index,:)';
            end
            
            globalObj.conditional_means(globalObj.search_set) = rapidObj.conditional_means;
            globalObj.conditional_variances(globalObj.search_set) = rapidObj.conditional_variances;
            globalObj.conditional_covariances(globalObj.search_set) = rapidObj.conditional_covariances;
            
        end
        
        function compute_conditional_statistics_fixed_set(globalObj, rapidObj)
            %Compute conditional statistics for solutions in the fixed set
            
            fixed_set = setdiff(linspace(1,globalObj.num_solutions,globalObj.num_solutions)',globalObj.search_set);
            current_optimal_solution_search_set_index = find(rapidObj.search_set == rapidObj.current_optimal_solution);
            
            
            conditional_covariance_matrix_upper_cholesky_factor = chol(rapidObj.conditional_covariance_matrix);
            product_of_A_and_cholesky_factor = rapidObj.intermediate_matrix_A * conditional_covariance_matrix_upper_cholesky_factor;
            A_conditional_covariance_matrix_A_diag = NaN(size(product_of_A_and_cholesky_factor,2),1);
            
            for i = 1:size(product_of_A_and_cholesky_factor,1)
                A_conditional_covariance_matrix_A_diag(i) = norm(product_of_A_and_cholesky_factor(i,:))^2;
            end
            
            if globalObj.pardiso_flag
                %Use this if using PARDISO to compute conditional means,
                %variances and covariances
                global id2 pardiso2
                globalObj.conditional_variances(fixed_set) = diag(pardiso2.invert(id2,tril(globalObj.conditional_precision_matrix(fixed_set, fixed_set)))) + A_conditional_covariance_matrix_A_diag;
                globalObj.conditional_means(fixed_set) = globalObj.MLE_mean + pardiso2.solve(id2,tril(globalObj.conditional_precision_matrix(fixed_set, fixed_set)),...
                    full(globalObj.intrinsic_precision_matrix(fixed_set,fixed_set) * (globalObj.response_sample_means(fixed_set) - globalObj.MLE_mean)))...
                    - rapidObj.intermediate_matrix_A * (globalObj.conditional_means(globalObj.search_set) - globalObj.MLE_mean);
                if ismember(globalObj.current_optimal_solution,globalObj.search_set)
                    globalObj.conditional_covariances(fixed_set) = -rapidObj.intermediate_matrix_A * rapidObj.conditional_covariance_matrix(:,current_optimal_solution_search_set_index);
                else
                    index_curr_opt_fixed_set = find(fixed_set==globalObj.current_optimal_solution);
                    globalObj.conditional_covariances(fixed_set) = pardiso2.solve(id2,tril(globalObj.conditional_precision_matrix(fixed_set, fixed_set)),...
                        double(1:length(fixed_set) == index_curr_opt_fixed_set)') + rapidObj.intermediate_matrix_A * rapidObj.conditional_covariance_matrix...
                        * rapidObj.intermediate_matrix_A(index_curr_opt_fixed_set,:)';
                end
            else
                %Use this if not using PARDISO to compute conditional
                %means, variances and covariances
                fixed_set_conditional_precision_matrix_cholesky_factor = chol(globalObj.conditional_precision_matrix(fixed_set, fixed_set))';
                globalObj.conditional_variances(fixed_set) = full(diag(fixed_set_conditional_precision_matrix_cholesky_factor'\...
                    (fixed_set_conditional_precision_matrix_cholesky_factor\speye(length(fixed_set),length(fixed_set)))));
                globalObj.conditional_means(fixed_set) = globalObj.MLE_mean + fixed_set_conditional_precision_matrix_cholesky_factor'\...
                    (fixed_set_conditional_precision_matrix_cholesky_factor\(globalObj.intrinsic_precision_matrix(fixed_set,fixed_set) *...
                    (globalObj.response_sample_means(fixed_set) - globalObj.MLE_mean))) - rapidObj.intermediate_matrix_A * (globalObj.conditional_means(globalObj.search_set) - globalObj.MLE_mean);
                if ismember(globalObj.current_optimal_solution,globalObj.search_set)
                    globalObj.conditional_covariances(fixed_set) = -rapidObj.intermediate_matrix_A * rapidObj.conditional_covariance_matrix(:,current_optimal_solution_search_set_index);
                else
                    index_curr_opt_fixed_set = find(fixed_set==globalObj.current_optimal_solution);
                    globalObj.conditional_covariances(fixed_set) = fixed_set_conditional_precision_matrix_cholesky_factor'\...
                    (fixed_set_conditional_precision_matrix_cholesky_factor\(1:length(fixed_set) == index_curr_opt_fixed_set)') + ...
                     rapidObj.intermediate_matrix_A * rapidObj.conditional_covariance_matrix * rapidObj.intermediate_matrix_A(index_curr_opt_fixed_set,:)';
                end
                    
                    
            end
                
                
                
        end          
        
        function construct_partition(globalObj)
        %Construct partiiton of the feasible region into search set and fixed set solutions, identified by their indices.
                        
            complete_expected_improvement = [globalObj.complete_expected_improvement linspace(1,globalObj.num_solutions,globalObj.num_solutions)'];
            complete_expected_improvement(globalObj.current_optimal_solution,:) = [];
            complete_expected_improvement = sortrows(complete_expected_improvement,-1);
            search_set = complete_expected_improvement(1:globalObj.search_set_size-1,2);
            search_set = [globalObj.current_optimal_solution;search_set];
            globalObj.search_set = search_set;
            
            fixed_set = setdiff(linspace(1,globalObj.num_solutions,globalObj.num_solutions)',globalObj.search_set);
            
            if globalObj.pardiso_flag
                global id2 pardiso2
                pardiso2.factorize(id2,tril(globalObj.conditional_precision_matrix(fixed_set, fixed_set)));
            end
        end
        
        
%         %Setter functions to perform error handling
%         
%         function set.bound_constraints(globalObj,value)
%             %Metamodel.bound_constraints needs to be a matrix of
%             %numerics, nonempty, integer-valued, finite, dense and
%             %appropriately sized
%             if ~isnumeric(value) || isempty(value) || sum(sum(value ~= floor(value))) || sum(sum(~isfinite(value))) || issparse(value) || size(value,2) ~= 2
%                 error('GlobalMetamodel.bound_constraints is invalid.');
%             else
%                 globalObj.bound_constraints = value;
%             end
%         end
%         
%         function set.design_points(globalObj,value)
%             %GlobalMetamodel.design_points needs to be a vector of
%             %numerics, nonempty, integer-valued, finite, dense and
%             %appropriately sized
%             if ~isnumeric(value) || isempty(value) || sum(value ~= floor(value)) || sum(~isfinite(value)) || issparse(value) || size(value,2) ~= 1
%                 error('GlobalMetamodel.bound_constraints is invalid.');
%             else
%                 globalObj.design_points = value;
%             end
%         end
%         
% 
%         
% 
%         
%         function set.search_set_size(globalObj,value)
%             %GlobalMetamodel.search_set_size needs to be an integer less
%             %than the number of feasible solutions if not empty
%             if ~isempty(value)
%                 if sum(value ~= floor(value)) || value > globalObj.num_solutions
%                     error('GLobalMetamodel.search_set_size is invalid.');
%                 else
%                     globalObj.search_set_size = value;
%                 end
%             end
%         end
        
        
      
        
        
        
    end

    
            
 
end

function q_cei =  compute_q_cei(conditional_means, conditional_covariance_matrix)
    %Helper function used to compute q-CEI, given conditional means and conditional covariance matrix.

            q = length(conditional_means)-1;

            sparse_id = speye(q,q);
            q_cei = 0;

            for k = 1:q
                D_k = [-sparse_id(:,1:k), sparse(ones(q,1)), -sparse_id(:,k+1:end)];
                M_k = D_k * conditional_means;
                Sigma_k = D_k * conditional_covariance_matrix * D_k';
                inner_sum = 0;


                for i = 1:q
                    set_excl_i = setdiff(linspace(1,q,q)',i);
                    C_ki = -M_k(set_excl_i) + (M_k(i)/Sigma_k(i,i))*Sigma_k(set_excl_i,i);
                    Cond_Sigma_ki = Sigma_k(set_excl_i,set_excl_i) - Sigma_k(set_excl_i,i)*Sigma_k(i,set_excl_i)/Sigma_k(i,i);

                    inner_sum = inner_sum + Sigma_k(i,1)*normpdf(0,M_k(i),sqrt(Sigma_k(i,i)))*mvncdf(C_ki,zeros(q-1,1),Cond_Sigma_ki);


                end
            q_cei = q_cei - (M_k(1) * mvncdf(zeros(q,1),M_k,Sigma_k) - inner_sum);

            end
        end
    
