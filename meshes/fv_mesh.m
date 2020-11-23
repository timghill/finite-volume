% Script to make meshes for FV Shallow Water solver

%% 1 - Rectangular mesh
boundary_xy=[-1,-1;
             1,-1;
             1,1;
             -1,1];
         
bmark=[2;2;2;2];
bmark_edge=[2;2;2;2];

% max_area=0.2*0.01;
max_area=0.4*0.01;

rmesh=make_mesh_wrapper(boundary_xy,bmark,bmark_edge,max_area);

mesh_plot_tri(gca,rmesh.tri, 1, 1)

axis image

save('rect_mesh.mat','-struct','rmesh')

%% 2 - Circular mesh for circular dam break problem
ds=0.005;
s=0:ds:1-ds;
s=s';
bxy=[cos(s*2*pi),sin(s*2*pi)];
bmark=ones(size(s));
bmark_edge=bmark;
max_area=2/3*0.01;

cmesh=make_mesh_wrapper(bxy,bmark,bmark_edge,max_area);

figure
hold on
mesh_plot_tri(gca,cmesh.tri,1,1)
axis image
save('circ_mesh.mat','-struct','cmesh')

%% 3 - Hyperbolic constriction
dx=0.01;
r=0.5;
bx_upper=-1:dx:1-dx;
by_upper=cosh(bx_upper)-1;
by_upper=(r*by_upper./max(by_upper))+1-r;

dy=0.01;
by_right=(1):-dy:-(1)+dy;
bx_right=ones(size(by_right));

bx_lower=bx_upper(end:-1:1);
by_lower=-by_upper(end:-1:1);

bx_left=-ones(size(by_right));
by_left=-by_right(end:-1:1);

bx=[bx_upper,bx_right,bx_lower,bx_left];
by=[by_upper,by_right,by_lower,by_left];

bx=bx(end:-1:1);
by=by(end:-1:1);

plot(bx,by)

bxy=[bx',by'];

bmark=ones(size(bx));
bmark_edge=bmark;
max_area=2/3*0.01;

hmesh=make_mesh_wrapper(bxy,bmark,bmark_edge,max_area);

figure
hold on
mesh_plot_tri(gca,hmesh.tri,1,1)
axis image
save('hyp_mesh.mat','-struct','hmesh')

% plot(bx_upper,by_upper)

%% 4 - Coarse rectangular-ish mesh
dx=0.1;
% r=0.5;
bx_upper=-1+dx:dx:1-dx;
by_upper=ones(size(bx_upper));

dy=0.1;
by_right=1-dy:-dy:-(1)+dy;
bx_right=ones(size(by_right));

bx_lower=bx_upper(end:-1:1);
by_lower=-by_upper(end:-1:1);

bx_left=-ones(size(by_right));
% by_left=by_right(end:-1:1);
by_left=-1+dy:dy:1-dy;

bx=[bx_upper,bx_right,bx_lower,bx_left];
by=[by_upper,by_right,by_lower,by_left];

bx=bx(end:-1:1);
by=by(end:-1:1);

o=ones(size(bx_upper))';
bmark=[3*o;2*o;2*o;2*o];
bmark_edge=[3*o;2*o;2*o;2*o];
% bmark(end)=3;
max_area=0.01;

bxy=[bx',by'];
rmesh=make_mesh_wrapper(bxy,bmark,bmark_edge,max_area);

% Fix boundary
% rmesh.tri.bmark_edge(69)=3;

% mesh_plot_tri(gca,rmesh.tri, 1, 1)

figure
hold on
for ii=1:rmesh.tri.n_edges
    nodes=rmesh.tri.connect_edge(ii,:);
    nodesx=rmesh.tri.nodes(nodes,1);
    nodesy=rmesh.tri.nodes(nodes,2);
    co=get(gca,'ColorOrder');
    plot(nodesx,nodesy,'Color',co(rmesh.tri.bmark_edge(ii)+1,:))
end

axis image

save('rect_mesh2.mat','-struct','rmesh')

%% 5 - Scaled coarse rectangular-ish mesh
ymin=-500;
ymax=500;
xmin=0;
xmax=2.5e3;
% max_area=(ymax-ymin)*xmax/575;
max_area=(ymax-ymin)*xmax/400;


bx=(bx-min(bx))./(max(bx)-min(bx));
bx=bx*xmax;

by=2*(by-min(by))./(max(by)-min(by))-1;
by=by*ymax;

plot(bx,by)

bxy=[bx', by'];
rmesh=make_mesh_wrapper(bxy,bmark,bmark_edge,max_area);


rmesh.tri.bmark_el(237,2)=3;
rmesh.tri.bmark_edge(497)=2;
figure
hold on
mesh_plot_tri(gca,rmesh.tri, 1, 1)
% for ii=1:rmesh.tri.n_edges
%     text(rmesh.tri.edge_midpoints(ii,1),rmesh.tri.edge_midpoints(ii,2),sprintf('%d',ii))
% end

% figure
% hold on
co=get(gca,'ColorOrder');
for jj=1:rmesh.tri.n_elements
    text(rmesh.tri.elements(jj,1),rmesh.tri.elements(jj,2),sprintf('%d',jj),'color',co(max(rmesh.tri.bmark_el(jj,:))+1,:))
end

axis image

save('phys_mesh.mat','-struct','rmesh')

%% 6 - rectangular mesh at various resolutions
areas=[0.001,0.0025,0.005,0.0075,0.01,0.025,0.05];
for ii=1:length(areas)
    area=areas(ii)
    bxy=[-1,-1;
        1,-1;
        1,1;
        -1,1];

    bmark=[2;2;2;2];
    bmark_edge=[2;2;2;2];
    rmesh=make_mesh_wrapper(bxy,bmark,bmark_edge,area);

    save(sprintf('dmesh_%02d.mat',ii),'-struct','rmesh')
end