classdef RapidMetamodel < Metamodel
    properties (SetAccess = ?GlobalMetamodel)
%         Properties only accessible from RapidMetamodel and
%         GlobalMetamodel classes

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
        
%         RapidMetamodel specific properties
        conditional_covariance_matrix %Dense nxn conditional covariance matrix of all solutions in search set.
        intermediate_matrix_A %Intermediate matrix A used to compute search set conditional statistics. Refer to paper for definition.
        intermediate_matrix_B %Intermediate matrix B used to compute search set conditional statistics. Refer to paper for definition.
        intermediate_vector_a %Intermediate vector a used to compute search set conditional statistics. Refer to paper for definition.

    end
    methods
        function rapidObj = RapidMetamodel(globalObj)
        %RapidMetamodel constructor   
            rapidObj.num_solutions = globalObj.search_set_size;
            rapidObj.search_set = globalObj.search_set;
            rapidObj.response_sum = full(globalObj.response_sum(globalObj.search_set));
            rapidObj.response_sum_of_square_differences = full(globalObj.response_sum_of_square_differences(globalObj.search_set));
            rapidObj.num_replications = full(globalObj.num_replications(globalObj.search_set));
            rapidObj.complete_expected_improvement = globalObj.complete_expected_improvement(globalObj.search_set);
            rapidObj.current_optimal_solution = globalObj.current_optimal_solution;
            rapidObj.max_cei_solutions = globalObj.max_cei_solutions;
            rapidObj.MLE_mean = globalObj.MLE_mean;
            rapidObj.MLE_precision = globalObj.MLE_precision;
            rapidObj.num_screening_solutions = globalObj.num_screening_solutions;
            rapidObj.num_simulation_workers = globalObj.num_simulation_workers;
            
            intr_prec_mat_diag = diag(globalObj.intrinsic_precision_matrix);
            rapidObj.precision_matrix = globalObj.precision_matrix(globalObj.search_set,globalObj.search_set);
            rapidObj.intrinsic_precision_matrix = diag(intr_prec_mat_diag(globalObj.search_set));
            rapidObj.conditional_precision_matrix = globalObj.conditional_precision_matrix(globalObj.search_set,globalObj.search_set);
            rapidObj.pardiso_flag = globalObj.pardiso_flag;
            
            
            fixed_set = setdiff(linspace(1,globalObj.num_solutions,globalObj.num_solutions)',globalObj.search_set);
            
            if globalObj.pardiso_flag
                global id2 pardiso2
                pardiso2.factorize(id2,tril(globalObj.conditional_precision_matrix(fixed_set,fixed_set)));
                rapidObj.intermediate_matrix_A = pardiso2.solve(id2,tril(globalObj.conditional_precision_matrix(fixed_set,fixed_set)),full(globalObj.conditional_precision_matrix(fixed_set,globalObj.search_set)));  
            else
                fixed_set_conditional_precision_matrix_cholesky_factor = chol(globalObj.conditional_precision_matrix(fixed_set, fixed_set))';
                rapidObj.intermediate_matrix_A = fixed_set_conditional_precision_matrix_cholesky_factor'\...
                    (fixed_set_conditional_precision_matrix_cholesky_factor\full(globalObj.conditional_precision_matrix(fixed_set,fixed_set)));
            end
            rapidObj.intermediate_matrix_B = globalObj.conditional_precision_matrix(fixed_set,globalObj.search_set)' * rapidObj.intermediate_matrix_A;
            rapidObj.intermediate_vector_a = rapidObj.intermediate_matrix_A' * (globalObj.intrinsic_precision_matrix(fixed_set,fixed_set) * (globalObj.response_sample_means(fixed_set) - globalObj.MLE_mean));
            
        end
        
        
        function compute_max_q_cei_solutions(rapidObj)
            %Find the q solutions that have the largest q-CEI
            global pardiso3 id3
            q = rapidObj.num_simulation_workers - 1;
            screened_solutions = sortrows([rapidObj.complete_expected_improvement linspace(1,rapidObj.num_solutions,rapidObj.num_solutions)'], -1);
            screened_solutions = [find(rapidObj.search_set==rapidObj.current_optimal_solution); screened_solutions(1:rapidObj.num_screening_solutions,2)];
            
            
            
            
            max_q_cei_set = screened_solutions(2);
            max_q_cei_set_opt = [screened_solutions(1); max_q_cei_set];
            
            
            
            for i = 2:q
                candidate_solutions = setdiff(screened_solutions, max_q_cei_set_opt);
                q_cei_cand = NaN(length(candidate_solutions),1);
                parfor (j = 1:length(candidate_solutions), rapidObj.num_simulation_workers)
                    cand_max_q_cei_set = [max_q_cei_set; candidate_solutions(j)];
                    cand_max_q_cei_set_opt = [max_q_cei_set_opt; candidate_solutions(j)];
                    q_cei_cand(j) = compute_q_cei(rapidObj.conditional_means(cand_max_q_cei_set_opt), ...
                       rapidObj.conditional_covariance_matrix(cand_max_q_cei_set_opt, cand_max_q_cei_set_opt));
                end
                [~, max_q_cei_cand_solution] = max(q_cei_cand);
                max_q_cei_set = [max_q_cei_set; candidate_solutions(max_q_cei_cand_solution)];
                max_q_cei_set_opt = [max_q_cei_set_opt; candidate_solutions(max_q_cei_cand_solution)];
            end
            
            rapidObj.max_cei_solutions = rapidObj.search_set(max_q_cei_set);
            
            
            
            
     
        end
        
        
        function compute_current_optimal_solution(rapidObj)
            %Find the current optimal solution (solution with the lowest
            %sample mean) within the search set
            sample_means = [rapidObj.response_sample_means linspace(1,length(rapidObj.response_sample_means),length(rapidObj.response_sample_means))'];
            sample_means = sample_means(sample_means(:,1)~=0,:);
            [~, index_current_optimal_solution] = min(sample_means(:,1));
            index_current_optimal_solution = sample_means(index_current_optimal_solution,2);
            rapidObj.current_optimal_solution = rapidObj.search_set(index_current_optimal_solution);
        end
        
        function compute_complete_expected_improvement(rapidObj)
            %Compue CEI values for all solutions in the search set
            current_optimal_solution_search_set_index = find(rapidObj.search_set == rapidObj.current_optimal_solution);
            
            variance_of_difference = rapidObj.conditional_variances + rapidObj.conditional_variances(current_optimal_solution_search_set_index) - 2*rapidObj.conditional_covariances;
            difference_of_means = rapidObj.conditional_means(current_optimal_solution_search_set_index) - rapidObj.conditional_means;
            critical_value = difference_of_means./sqrt(variance_of_difference);
            variance_of_difference = variance_of_difference .* (abs(variance_of_difference) >= 10e-10);
            complete_expected_improvement = difference_of_means.*normcdf(critical_value) + sqrt(variance_of_difference).*normpdf(critical_value);
            complete_expected_improvement(current_optimal_solution_search_set_index) = 0;
            rapidObj.complete_expected_improvement = complete_expected_improvement;
        end
        
        function compute_max_cei_solution(rapidObj)
            %Find the solution in the search set with the largest CEI value
            [~, max_cei_solution] = max(rapidObj.complete_expected_improvement);
            rapidObj.max_cei_solutions = rapidObj.search_set(max_cei_solution);            
        end
        
        function update_precision_matrices(rapidObj)
            %Update precision matrices with new simulation information
            current_optimal_solution_search_set_index = find(rapidObj.search_set == rapidObj.current_optimal_solution);
            max_cei_solutions_search_set_indices = find(ismember(rapidObj.search_set, rapidObj.max_cei_solutions));

            
            intr_prec_mat_diag = diag(rapidObj.intrinsic_precision_matrix);
            intr_prec_mat_diag(current_optimal_solution_search_set_index) = rapidObj.num_replications(current_optimal_solution_search_set_index)/...
                rapidObj.response_sample_variances(current_optimal_solution_search_set_index);
            
            intr_prec_mat_diag(max_cei_solutions_search_set_indices) = (rapidObj.num_replications(max_cei_solutions_search_set_indices)./...
                rapidObj.response_sample_variances(max_cei_solutions_search_set_indices)).*(rapidObj.num_replications(max_cei_solutions_search_set_indices) > 0);

            
            rapidObj.intrinsic_precision_matrix = diag(intr_prec_mat_diag);
            
            rapidObj.conditional_precision_matrix(current_optimal_solution_search_set_index,current_optimal_solution_search_set_index) = ...
                rapidObj.precision_matrix(current_optimal_solution_search_set_index,current_optimal_solution_search_set_index) + ...
                rapidObj.intrinsic_precision_matrix(current_optimal_solution_search_set_index,current_optimal_solution_search_set_index);
            
            max_cei_solution_mat_lin_index = sub2ind(size(rapidObj.intrinsic_precision_matrix),max_cei_solutions_search_set_indices, max_cei_solutions_search_set_indices);
            
            rapidObj.conditional_precision_matrix(max_cei_solution_mat_lin_index) = ...
                rapidObj.precision_matrix(max_cei_solution_mat_lin_index) + rapidObj.intrinsic_precision_matrix(max_cei_solution_mat_lin_index);
            
            
            
        end
        
        function compute_conditional_statistics(rapidObj)
           %Compute conditional statistics for all solutions in the search
           %set
           current_optimal_solution_search_set_index = find(rapidObj.search_set == rapidObj.current_optimal_solution); 
            
           rapidObj.conditional_covariance_matrix = inv(rapidObj.conditional_precision_matrix - rapidObj.intermediate_matrix_B);
           rapidObj.conditional_variances = diag(rapidObj.conditional_covariance_matrix);
           rapidObj.conditional_covariances = rapidObj.conditional_covariance_matrix(:,current_optimal_solution_search_set_index);
           rapidObj.conditional_means = rapidObj.MLE_mean + rapidObj.conditional_covariance_matrix*(rapidObj.intrinsic_precision_matrix*(rapidObj.response_sample_means - rapidObj.MLE_mean) - ...
               rapidObj.intermediate_vector_a);
        end        

                
        
        
        %Setter functions to perform error handling
        
%         function set.conditional_covariance_matrix(rapidObj,value)
%             %RapidMetamodel.conditional_covariance_matrix needs to be a
%             %matrix of numerics, nonempty, finite, dense and appropriately
%             %sized
%             if ~isnumeric(value) || isempty(value) || sum(sum(~isfinite(value))) || issparse(value) || sum(size(value) ~= [rapidObj.num_solutions, rapidObj.num_solutions])
%                 error('RapidMetamodel.conditional_covariance_matrix is invalid.');
%             else
%                 rapidObj.conditional_covariance_matrix = value;
%             end
%         end
%         
%         function set.intermediate_matrix_A(rapidObj,value)
%             %RapidMetamodel.intermediate_matrix_A needs to be a
%             %matrix of numerics, nonempty, finite, dense and appropriately
%             %sized
%             if ~isnumeric(value) || isempty(value) || sum(sum(~isfinite(value))) || issparse(value) || size(value,2) ~= rapidObj.num_solutions
%                 error('RapidMetamodel.conditional_covariance_matrix_A is invalid.');
%             else
%                 rapidObj.intermediate_matrix_A = value;
%             end
%         end
%         
%         function set.intermediate_matrix_B(rapidObj,value)
%             %RapidMetamodel.intermediate_matrix_B needs to be a
%             %matrix of numerics, nonempty, finite, dense and appropriately
%             %sized
%             if ~isnumeric(value) || isempty(value) || sum(sum(~isfinite(value))) || issparse(value) || sum(size(value) ~= [rapidObj.num_solutions, rapidObj.num_solutions])
%                 error('RapidMetamodel.conditional_covariance_matrix_B is invalid.');
%             else
%                 rapidObj.intermediate_matrix_B = value;
%             end
%         end        
%     
%         function set.intermediate_vector_a(rapidObj,value)
%             %RapidMetamodel.intermediate_vector_a needs to be a
%             %vector of numerics, nonempty, finite, dense and appropriately
%             %sized
%             if ~isnumeric(value) || isempty(value) || sum(sum(~isfinite(value))) || issparse(value) || size(value,2) ~= 1
%                 error('RapidMetamodel.conditional_covariance_vector_a is invalid.');
%             else
%                 rapidObj.intermediate_vector_a = value;
%             end
%         end             
%             
            
        

        

        
        
    end
        
        
    
  

                
end

function q_cei = compute_q_cei(conditional_means, conditional_covariance_matrix)
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
                    if min(eig(Cond_Sigma_ki)) <= 0
                        keyboard;
                    end
                    Cond_Sigma_ki = 0.5*(Cond_Sigma_ki + Cond_Sigma_ki');
                    inner_sum = inner_sum + Sigma_k(i,1)*normpdf(0,M_k(i),sqrt(Sigma_k(i,i)))*mvncdf(C_ki,zeros(q-1,1),Cond_Sigma_ki);


                end
	    Sigma_k = 0.5*(Sigma_k + Sigma_k');
            q_cei = q_cei - (M_k(1) * mvncdf(zeros(q,1),M_k,Sigma_k) - inner_sum);

            end
        end
    
