function coordinates = index_to_coordinate(bound_constraints, index)
%index_to_coordinate takes in dense dx2 matrix containing inclusive lower
%and upper-bounds of box constraints for each dimension and linear index
%and outputs the appropriate coordinates of the solution.

[coordinates{1:size(bound_constraints,1)}] = ind2sub(bound_constraints(:,2) - (bound_constraints(:,1)-1), index);

coordinates = reshape(cell2mat(coordinates'),length(coordinates{1}),length(coordinates)) + repmat((bound_constraints(:,1)'-1),length(index),1);

    
end

