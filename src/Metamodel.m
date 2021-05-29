classdef (Abstract=true) Metamodel < handle
    %METAMODEL superclass from which GlobalMetamodel and RapidMetamodel
    %inherit common traits.
    %   The Metamodel class is an abstract superclass (no objects are
    %instantiated) that contains traits common to GlobalMetamodel and
    %RapidMetamodel objects, assuming there are n solutions and d
    %dimensions.
    properties
        response_sum %Sparse nx1 vector containing sum of response across replications for simulated solutions in Metamodel object, assuming column-major ordering.
        response_sum_of_square_differences %Sparse nx1 vector containing sum of square difference between response and mean across replications for simulated solutions in a Metamodel object, assuming column-major ordering.
        num_replications %Sparse nx1 vector containing the total number of replications performed for simulated solutions in a Metamodel object, assuming column-major ordering.
    end
    properties (SetAccess = protected)
        precision_matrix %Sparse nxn precision matrix of Metamodel object based on MLE estimation.
        intrinsic_precision_matrix %Sparse nxn intrinsic precision matrix of Metamodel object that is updated as simulation occurs.
        conditional_precision_matrix %Sparse nxn conditional precision matrix of Metamodel object that is updated as simulation occurs.
        conditional_means %Dense nx1 vector of conditional means for all solutions in Metamodel object.
        conditional_variances %Dense nx1 vector of conditional variances for all solutions in Metamodel object.
        conditional_covariances %Dense nx1 vector of conditional covariances relative to current optimal solution in Metamodel object.
        complete_expected_improvement %Dense nx1 vector of CEI values for all solutions in Metamodel object.
        current_optimal_solution %Scalar value containing the index of the current optimal solution, assuming column-major ordering.
        max_cei_solutions %Scalar/qx1 vector containing the index/indices of either the solution with the largest CEI value or the set of solutions estimated to have the largest q-CEI values.
        num_solutions %Scalar number of feasible solutions in the Metamodel object. For a GlobalMetamodel object, this is the number of solutions in the feasible region. For a RapidMetamodel object, this is the number of solutions in the search set.
        search_set %Dense vector of indices of the solutions chosen to be in the search set. This is left empty in GMIA.
        MLE_mean %Scalar value representing the MLE value of the mean of the GMRF.
        MLE_precision %Dense vector dx1 containing the MLE values of the conditional precision of solutions and conditional correlations between solutions.
        num_simulation_workers %Scalar number of available processors available to perform simulations in parallel.
        num_screening_solutions %Number of candidate solutions used to construct the q-CEI optimal set, based on their marginal CEI values.
        pardiso_flag %Binary (0/1) flag to determine whether PARDISO solver is to be used (1) or not (0).
    end
    properties (Dependent = true)
        response_sample_means %Sparse nx1 vector containing the response sample means of simulated solutions.
        response_sample_variances %Sparse nx1 vector containing the response sample variances of simulated solutions.
    end
       
    methods
        
        
        function compute_complete_expected_improvement(obj)
            %Compute CEI values for all solutions in GMRF and return as a vector, assuming column-major ordering.
            variance_of_difference = obj.conditional_variances + obj.conditional_variances(obj.current_optimal_solution) - 2*obj.conditional_covariances;
            difference_of_means = obj.conditional_means(obj.current_optimal_solution) - obj.conditional_means;
            critical_value = difference_of_means./sqrt(variance_of_difference);
            variance_of_difference = variance_of_difference .* (abs(variance_of_difference) >= 10e-10);
            complete_expected_improvement = difference_of_means.*normcdf(critical_value) + sqrt(variance_of_difference).*normpdf(critical_value);
            complete_expected_improvement(obj.current_optimal_solution) = 0;
            obj.complete_expected_improvement = complete_expected_improvement;
        end
        
        



            
            

        
        function value = get.response_sample_means(obj)
            %Getter method to return response sample means. If called on solutions have not been simulated, the value returned will be NaN.
            value = obj.response_sum./obj.num_replications;
        end
        
        function value = get.response_sample_variances(obj)
            %Getter method to return response sample variances. If called on solutions have not been simulated, the value returned will be NaN.
            value = obj.response_sum_of_square_differences./obj.num_replications;
        end
        
        
        %Setter methods to perform error handling
        
%         function set.response_sum(obj,value)
%             %Metamodel.response_sum needs to be a vector of
%             %numerics, nonempty, real-valued, finite and
%             %appropriately sized
%             if ~isnumeric(value) || isempty(value) || ~isreal(value) || sum(~isfinite(value)) || sum(size(value) ~= [obj.num_solutions, 1])
%                 error('Metamodel.response_sum is invalid.');
%             else
%                 obj.response_sum = value;
%             end
%         end
%         
%         function set.response_sum_of_square_differences(obj,value)
%             %Metamodel.response_sum_of_square_differences needs to be a vector of
%             %numerics, nonempty, real-valued, >= 0, finite and
%             %approprtiately sized
%             if ~isnumeric(value) || isempty(value) || ~isreal(value) || sum(~ge(value,0)) || sum(~isfinite(value)) || sum(size(value) ~= [obj.num_solutions, 1])
%                 error('Metamodel.response_sum_of_square_differences is invalid.');
%             else
%                 obj.response_sum_of_square_differences = value;
%             end
%         end
%         
%         function set.num_replications(obj,value)
%             %Metamodel.num_replications needs to be a vector of numerics,
%             %nonempty, integer, >=0, finite and appropriately
%             %sized
%             if ~isnumeric(value) || isempty(value) || sum(value ~= floor(value)) || sum(~ge(value,0)) || sum(~isfinite(value)) || sum(size(value) ~= [obj.num_solutions, 1])
%                 error('Metamodel.num_replications is invalid.');
%             else
%                 obj.num_replications = value;
%             end
%         end
%         
%         function set.precision_matrix(obj,value)
%             %Metamodel.precision_matrix needs to be a matrix of
%             %numerics, nonempty, finite, sparse and appropriately
%             %sized
%             if ~isnumeric(value) || isempty(value) || sum(sum(~isfinite(value))) || ~issparse(value) || sum(size(value) ~= [obj.num_solutions, obj.num_solutions])
%                 error('Metamodel.precision_matrix is invalid.');
%             else
%                 obj.precision_matrix = value;
%             end
%         end
%         
%         function set.intrinsic_precision_matrix(obj,value)
%             %Metamodel.intrinsic_precision_matrix needs to be a diagonal
%             %matrix of numerics, nonempty, finite, sparse and
%             %appropriately sized
%             if ~isnumeric(value) || isempty(value) || ~isdiag(value) || sum(sum(~isfinite(value))) || ~issparse(value) || sum(size(value) ~= [obj.num_solutions, obj.num_solutions])
%                 error('Metamodel.intrinsic_precision_matrix is invalid.');
%             else
%                 obj.intrinsic_precision_matrix = value;
%             end
%         end
%         
%         function set.conditional_precision_matrix(obj,value)
%             %Metamodel.conditional_precision_matrix needs to be a matrix of
%             %numerics, nonempty, finite, sparse and appropriately
%             %sized
%             if ~isnumeric(value) || isempty(value) || sum(sum(~isfinite(value))) || ~issparse(value) || sum(size(value) ~= [obj.num_solutions, obj.num_solutions])
%                 error('Metamodel.conditional_precision_matrix is invalid.');
%             else
%                 obj.conditional_precision_matrix = value;
%             end
%         end
%         
%          function set.conditional_means(obj,value)
%             %Metamodel.conditional_means needs to be a vector of
%             %numerics, nonempty, real-valued, finite, dense and
%             %appropriately sized
%             if ~isnumeric(value) || isempty(value) || ~isreal(value) || sum(~isfinite(value)) || issparse(value) || sum(size(value) ~= [obj.num_solutions, 1])            
%                 error('Metamodel.conditional_means is invalid.');
%             else
%                 obj.conditional_means = value;
%             end
%         end
%         
%         function set.conditional_variances(obj,value)
%             %Metamodel.conditional_variances needs to be a vector of
%             %numerics, nonempty, real-valued, >= 0, finite, dense and
%             %appropriately sized
%             if ~isnumeric(value) || isempty(value) || ~isreal(value) || sum(~ge(value,0)) || sum(~isfinite(value)) || issparse(value) || sum(size(value) ~= [obj.num_solutions, 1])
%                 error('Metamodel.conditional_variances is invalid.');
%             else
%                 obj.conditional_variances = value;
%             end
%         end                
% 
%         function set.conditional_covariances(obj,value)
%             %Metamodel.conditional_covariances needs to be a vector of
%             %numerics, nonempty, real-valued, finite, dense and
%             %appropriately sized
%             if ~isnumeric(value) || isempty(value) || ~isreal(value) || sum(~isfinite(value)) || issparse(value) || sum(size(value) ~= [obj.num_solutions, 1])
%                 error('Metamodel.conditional_covariances is invalid.');
%             else
%                 obj.conditional_covariances = value;
%             end
%         end        
% 
%         function set.complete_expected_improvement(obj,value)
%             %Metamodel.complete_expected_improvement needs to be a vector of
%             %numerics, nonempty, real-valued, >= 0, finite, dense and
%             %appropriately sized
%             if ~isnumeric(value) || isempty(value) || ~isreal(value) || sum(~ge(value,0)) || sum(~isfinite(value)) || issparse(value) || sum(size(value) ~= [obj.num_solutions, 1])
%                 error('Metamodel.complete_expected_improvement is invalid.');
%             else
%                 obj.complete_expected_improvement = value;
%             end
%         end
%         
%         function set.current_optimal_solution(obj,value)
%             %Metamodel.current_optimal_solution needs to be a numeric
%             %integer, nonempty, > 0, finite and appropriately sized
%             if ~isnumeric(value) || isempty(value) || value ~= floor(value) || ~gt(value,0) || ~isfinite(value) || sum(size(value) ~= [1, 1])
%                 error('Metamodel.current_optimal_solution is invalid.');
%             else
%                 obj.current_optimal_solution = value;
%             end
%         end
%         
%         function set.max_cei_solutions(obj,value)
%             %Metamodel.max_cei_solutions needs to be a numeric integer,
%             %nonempty, > 0, finite
%             if ~isnumeric(value) || isempty(value) || value ~= floor(value) || ~gt(value,0) || ~isfinite(value)
%                 error('Metamodel.max_cei_solutions is invalid.');
%             else
%                 obj.max_cei_solutions = value;
%             end
%         end
%         
%         function set.search_set(obj,value)
%             %Metamodel.search_set needs to be a vector of
%             %numerics, nonempty, integer-valued, finite and dense
%             if ~isnumeric(value) || isempty(value) || sum(value ~= floor(value)) || sum(~isfinite(value)) || issparse(value)
%                 error('Metamodel.search_set is invalid.');
%             else
%                 obj.search_set = value;
%             end
%         end        
%                 
%         function set.MLE_mean(obj,value)
%             %Metamodel.MLE_mean needs to be a nonempty numeric scalar
%             if ~isnumeric(value) || isempty(value) || sum(size(value) ~= [1, 1])
%                 error('Metamodel.MLE_mean is invalid.');
%             else
%                 obj.MLE_mean = value;
%             end
%         end
%         
%         function set.MLE_precision(obj,value)
%             %Metamodel.MLE_precision needs to be a vector of
%             %numerics, nonempty, finite and dense 
%             if ~isnumeric(value) || isempty(value) || sum(~isfinite(value)) || issparse(value)
%                 error('Metamodel.MLE_precision is invalid.');
%             else
%                 obj.MLE_precision = value;
%             end
%         end         
   

   
    
    
    end

end
            
        
        
