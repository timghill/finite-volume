function z_node = interp_el_node(dmesh, z_el)
F = scatteredInterpolant(dmesh.tri.elements(:,1), dmesh.tri.elements(:, 2), z_el);
z_node = F(dmesh.tri.nodes(:, 1), dmesh.tri.nodes(:, 2));