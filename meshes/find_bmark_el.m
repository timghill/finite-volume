bmark_el=zeros(dmesh.tri.n_elements,3);
% Find boundary mark for each edge of each element
for ii=1:dmesh.tri.n_elements
   nodes=dmesh.tri.connect(ii,:);
   nodes=[nodes,nodes(1)];
   
   edges=zeros(1,3);
   for jj=1:3
       edgei=find(dmesh.tri.connect_edge(:,1)==nodes(jj) & dmesh.tri.connect_edge(:,2)==nodes(jj+1));
       edgej=find(dmesh.tri.connect_edge(:,2)==nodes(jj) & dmesh.tri.connect_edge(:,1)==nodes(jj+1));
       
       if isempty(edgei)
          edgeix=edgej; 
       else
           edgeix=edgei;
       end
       
       edge_bmark=dmesh.tri.bmark_edge(edgeix);
       
       bmark_el(ii,jj)=edge_bmark;
   end
end