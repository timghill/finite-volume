function dmesh = supplement_dmesh(dmesh)
% Calculate additional fields that help optimize the FV methods

% Calculate normal components for edges
dmesh.tri.edge_nx = zeros(dmesh.tri.n_edges, 1);
dmesh.tri.edge_ny = zeros(dmesh.tri.n_edges, 1);

dmesh.tri.edge_area = zeros(dmesh.tri.n_edges, 2);

% figure
% hold on
% axis image

for kk=1:dmesh.tri.n_edges
    
%     nodes = dmesh.tri.connect_edge(kk,:);
%     plot(dmesh.tri.nodes(nodes,1), dmesh.tri.nodes(nodes, 2), 'k')
    
    
    neigh_els = dmesh.tri.connect_edge_el(kk,:);
    if neigh_els(1)>0
        iEdge = find(dmesh.tri.connect_el_edge(neigh_els(1),:)==kk);
        dmesh.tri.edge_nx(kk) = dmesh.tri.nx(neigh_els(1), iEdge);
        dmesh.tri.edge_ny(kk) = dmesh.tri.ny(neigh_els(1), iEdge);
    else
        iEdge = find(dmesh.tri.connect_el_edge(neigh_els(2),:)==kk);
        dmesh.tri.edge_nx(kk) = -dmesh.tri.nx(neigh_els(2), iEdge);
        dmesh.tri.edge_ny(kk) = -dmesh.tri.ny(neigh_els(2), iEdge);
    end

    if neigh_els(1)>0
        dmesh.tri.edge_area(kk,1) = dmesh.tri.area(neigh_els(1));
    else
        dmesh.tri.edge_area(kk,1) = dmesh.tri.area(neigh_els(2));
    end
    
    if neigh_els(2)>0
        dmesh.tri.edge_area(kk,2) = dmesh.tri.area(neigh_els(2));
    else
        dmesh.tri.edge_area(kk,2) = dmesh.tri.area(neigh_els(1));
    end
    
    
%     if mod(kk,25)==0 && min(neigh_els)>0
%         norm_scale = 0.033;
%         plot(dmesh.tri.edge_midpoints(kk,1), dmesh.tri.edge_midpoints(kk,2), 'r*')
%         plot(dmesh.tri.elements(neigh_els(1),1), dmesh.tri.elements(neigh_els(1),2), 'b*')
%         plot(dmesh.tri.elements(neigh_els(2),1), dmesh.tri.elements(neigh_els(2),2), 'g*')
%         plot([dmesh.tri.edge_midpoints(kk,1), dmesh.tri.edge_midpoints(kk,1) + norm_scale*dmesh.tri.edge_nx(kk)],...
%             [dmesh.tri.edge_midpoints(kk,2), dmesh.tri.edge_midpoints(kk,2) + norm_scale*dmesh.tri.edge_ny(kk)], 'b')
%        
%     end
   
end