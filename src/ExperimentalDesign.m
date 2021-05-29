classdef ExperimentalDesign
    %EXPERIMENTALDESIGN contains the design points, associated simulation information and resulting GMRF MLEs to be used by GMIA and rGMIA.
    %   An ExperimentalDesign object stores summary statistics (e.g., MLE
    %estimate of mean). This object is used as an input to run GMIA/rGMIA.
    

    properties
        bound_constraints %Dense dx2 matrix containing inclusive lower and upper-bounds of box constraints for each dimension.
        response_sum %Sparse nx1 vector containing sum of response after num_replications. Vectorization is done in column-major ordering.
        response_sum_of_square_differences %Sparse nx1 vector containing sum of square difference between response and mean. Vectorization is done in column-major ordering.
        num_replications %Sparse nx1 number of replications simulated at each solution. Vectorization is done in column-major ordering.
        design_points %Dense kx1 vector containing linear index of design points, assuming column-major ordering.
        MLE_mean %Scalar mean of the GMRF.
        MLE_precision %Dense dx1 vector containing thetas: [\theta_0, \theta_1, ..., \theta_d]
    end
    
    methods
        function experimentalDesignObj = ExperimentalDesign(bound_constraints, response_sum, response_sum_of_square_differences, num_replications, design_points, MLE_mean, MLE_precision)
        %ExperimentalDesign Constructor Method
            experimentalDesignObj.bound_constraints = bound_constraints;
            experimentalDesignObj.response_sum = response_sum;
            experimentalDesignObj.response_sum_of_square_differences = response_sum_of_square_differences;
            experimentalDesignObj.num_replications = num_replications;
            experimentalDesignObj.design_points = design_points;
            experimentalDesignObj.MLE_mean = MLE_mean;
            experimentalDesignObj.MLE_precision = MLE_precision;
                
        end
    end
    
end

