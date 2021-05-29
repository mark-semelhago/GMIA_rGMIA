function Q = prec_mat_constr(theta, grid_dimensions)
%PREC_MAT_CONSTR 
%Creates sparse precision matrix


d = length(grid_dimensions);
row = [theta(1) -theta(1)*theta(2)];
for i = 3:length(theta)
    row = [row zeros(1,prod(grid_dimensions(1:i-2))-length(row)) -theta(1)*theta(i)];
end
row = [row zeros(1,prod(grid_dimensions)-length(row))];
Q = sptoeplitz(row);


%Take into account boundary effects
offset = find(row~=0);
offset = offset(2:end)-1;
for i = 1:d-1
    zeroed_rows = prod(grid_dimensions(1:i)):prod(grid_dimensions(1:i)):size(Q,1)-1;
    Q(zeroed_rows,zeroed_rows+offset(i)) = 0;
end
Q = triu(Q) + triu(Q)'- spdiags(diag(Q),0,prod(grid_dimensions),prod(grid_dimensions));



end

