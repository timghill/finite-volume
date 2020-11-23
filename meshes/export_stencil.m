close all; clear all

dmesh=load('dmesh_03.mat');

% Arbitrary interior element index
ii=600;

disp('Central element node coordinates:')
disp(dmesh.tri.nodes(dmesh.tri.connect(ii,:),:))

disp('Central element centroid coordinates')
disp(dmesh.tri.elements(ii,:))

% Edge index
edge=dmesh.tri.connect_el_edge(ii,2);

neigh_els=dmesh.tri.edge_stencil{edge};

nodes=dmesh.tri.connect(neigh_els,:);

disp('Edge midpoint:')
dmesh.tri.edge_midpoints(edge,:)


disp('Edge stencil node coordinates:')
for kk=1:length(nodes)
    disp(dmesh.tri.nodes(dmesh.tri.connect(neigh_els(kk),:),:))
end

disp('_______________________________________________________________')
% Now we look at the compact and extended element stencils
disp('Compact stencil nodes:')
cnodes=dmesh.tri.connect(dmesh.tri.node_stencil_compact{ii},:);
for kk=1:length(cnodes)
   disp(dmesh.tri.nodes(cnodes(kk,:),:))
end

disp('_______________________________________________________________')
% Now we look at the compact and extended element stencils
disp('Compact extended nodes:')
cnodes=dmesh.tri.connect(dmesh.tri.node_stencil_extended{ii},:);
for kk=1:length(cnodes)
   disp(dmesh.tri.nodes(cnodes(kk,:),:))
end