figure
hold on

for ii=1:dmesh.tri.n_edges
    nn = dmesh.tri.connect_edge(ii,:);
    plot(dmesh.tri.nodes(nn,1), dmesh.tri.nodes(nn,2), 'k')
end

ii = 50;
iEdge = dmesh.tri.connect_el_edge(ii, 2);
plot(dmesh.tri.elements(ii, 1), dmesh.tri.elements(ii, 2), 'r*')

neighs = dmesh.tri.edge_stencil{iEdge};
plot(dmesh.tri.elements(neighs,1), dmesh.tri.elements(neighs,2), 'bo')