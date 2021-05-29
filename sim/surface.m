function output = surface(design_points, grid_dimensions, y, var_noise, rep)
%Surface simulator takes a deterministic surface as input, y, and adds white
%noise with variance var_noise to convert the surface to a simulation.
%design_points is the solution to simulate.
%grid_dimensions is a vector with the number of solutions along each
%dimension.
%y is the response surface
%var_noise is the variance of white noise added to the surface.
%rep is the number of replications to produce.



d = size(design_points,2);
design_points_index = zeros(size(design_points,1),1);
for i = d:-1:2
    design_points_index = design_points_index + (design_points(:,i)-1)*prod(grid_dimensions(1:i-1));
end
design_points_index = design_points_index + design_points(:,1);
    

var_noise = var_noise(design_points_index);
y = y(design_points_index);

output = repmat(y,1,rep) + repmat(sqrt(var_noise),1,rep).*randn(size(design_points,1),rep);


end

