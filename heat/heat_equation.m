% heat_equation.m solves the two-dimensional linear heat equation on an
% unstructured triangular mesh using finite volume methods. The solver uses
% first-order FV discretization in space with built-in matlab
% time-stepping.

dmesh=load('../meshes/circ_mesh.mat');

%% Setup
% Time-stepping
tend=1;
t=0;
tspan=linspace(t,tend,51);

% Parameters
params.gamma=1e-1;

% The solver supports three types of boundary conditions, specified
% using params.bc:
%   (1) 'dirichlet':    Boundary has fixed value params.v_dirichlet
%   (2) 'neumann':      Boundary has no heat flux out of the domain
%                       (grad.v=0 on boundary) [These BC
%   (3) 'flux':         Allow flux out of domain by enforcing
%                       d(grad.v)/dn=0 on the boundary.
% Note that 'neumann' boundary conditions conserve thermal energy, while
% the other BCs allow heat to flow out of the domain through the boundary

params.bc='neumann';
params.v_dirichlet=0;   % Must specify value for Dirichlet condition

%% Initial conditions
u=zeros(size(dmesh.tri.elements(:,1)));
trix=dmesh.tri.elements(:,1);
triy=dmesh.tri.elements(:,2);
trinorm=sqrt(trix.^2+triy.^2);
% u(trix>-0.5 & trix<0.5 & triy>-0.5 & triy<0.5)=2;
u(trinorm<0.5)=1;

M0=sum(u.*dmesh.tri.area);

%% Solver
odefun=@(t,y) rhs_heat_unstructured(t,y,dmesh,params);
tic;
[tt,yout]=ode45(odefun,tspan,u);
toc;

%% Post-processing
u=yout';
% Animate the outputs
for ii=1:length(tspan)
    t=tspan(ii);
    figure(1)
    u_node=interp_el_node(dmesh,u(:,ii));
    if strcmp(params.bc,'dirichlet')
        u_node(dmesh.tri.bmark==1)=0;
    end
    trisurf(dmesh.tri.connect,dmesh.tri.nodes(:,1),dmesh.tri.nodes(:,2),u_node)
    axis([-1,1,-1,1,-0.25,1.25])
    colormap('winter')
    title(sprintf('t = %.3f', t))
    caxis([-0.25,1.25])
    drawnow
    pause(0.1)
end

kk=[10,30];
for ii=kk
    t=tspan(ii);
    figure
    u_node=interp_el_node(dmesh,u(:,ii));
    if strcmp(params.bc,'dirichlet')
        u_node(dmesh.tri.bmark==1)=0;
    end
    trisurf(dmesh.tri.connect,dmesh.tri.nodes(:,1),dmesh.tri.nodes(:,2),u_node)
    axis([-1,1,-1,1,0,1])
    cmocean('thermal')
    title(sprintf('t = %.3f', t))
    caxis([0,1])
    drawnow
    print(sprintf('heat_%s_%i',params.bc,ii),'-dpng','-r600')
end

% Compute total "Mass" in the domain. Depending on the boundary conditions
% this may be conserved.
M=sum(u.*dmesh.tri.area,1)/M0;
figure
plot(tt,M)