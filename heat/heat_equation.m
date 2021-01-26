% heat_equation.m solves the two-dimensional linear heat equation on an
% unstructured triangular mesh using finite volume methods with built-in
% matlab time-stepping.

dmesh=load('../meshes/circ_mesh.mat');
dmesh = supplement_dmesh(dmesh);

%% Setup
% Time-stepping
tend=1;
t=0;
tspan=linspace(t,tend,51);
dt = 0.01;

% Parameters
params.gamma=1e-1;      % Thermal conductivity constant

% alpha and beta control the flux parameterization. The flux is of the form
%       gamma * v^alpha * (grad.v)^beta.
% Therefore, alpha = 0 and beta = 1 corresponds to the linear heat
% equation, while alpha >0 and beta = 0.5 corresponds to a turbulent flow
% parameterization.
params.alpha = 0;
params.beta = 1;

% The solver supports three types of boundary conditions, specified
% using params.bc:
%   (1) 'dirichlet':    Boundary has fixed value params.v_dirichlet
%   (2) 'neumann':      Boundary has no heat flux out of the domain
%                       (grad.v=0 on boundary) [These BC
%   (3) 'flux':         Allow flux out of domain by enforcing
%                       d(grad.v)/dn=0 on the boundary.
% Note that 'neumann' boundary conditions conserve thermal energy, while
% the other BCs allow heat to flow out of the domain through the boundary
params.v_dirichlet=0;   % Must specify value for Dirichlet condition

% Switches
% The solver supports different methods to calculate gradients. The options
% are:
% (1) 'gg' to use simple Green-Gauss/Divergence theorem method. Boundary values
%      are computed as the average of neighbouring element values
% (2) 'gg-hybrid' to use Green-Gauss method with a least-squares approach
%     calculate boundary values.
% (3) 'lsq' to use weighted least-squares approach to calculate the
%      gradient value on each element.
gradient_solvers = {'gg', 'gg-hybrid', 'lsq'};
params.gradient = gradient_solvers{2};
rhsfunc = @rhs_heat_unstructured_optimized;
params.bc='dirichlet';

%% Initial conditions
u0=zeros(size(dmesh.tri.elements(:,1)));
trix=dmesh.tri.elements(:,1);
triy=dmesh.tri.elements(:,2);
trinorm=sqrt(trix.^2+triy.^2);
% u0(trinorm<0.5)=1;
u0 = exp(-trinorm.^2/0.25);

M0=sum(u0.*dmesh.tri.area);

params.derivs = true;
dx = rhsfunc(u0, dmesh, params);
params.derivs = false;
figure
trisurf(dmesh.tri.connect,dmesh.tri.nodes(:,1),dmesh.tri.nodes(:,2),interp_el_node(dmesh, dx.vx))

params.derivs = false;

%% Solver
odefun=@(t,y) rhsfunc(y,dmesh,params);
opts = odeset('Stats', 'on');
tic;
% [tt,yout]=odeRK(odefun,tspan,dt,u0,opts);
[tt,yout] = ode45(odefun,tspan,u0,opts);
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

kk=[20,51];
for ii=kk
    t=tspan(ii);
    figure
    u_node=interp_el_node(dmesh,u(:,ii));
    if strcmp(params.bc,'dirichlet')
        u_node(dmesh.tri.bmark==1)=0;
    end
    trisurf(dmesh.tri.connect,dmesh.tri.nodes(:,1),dmesh.tri.nodes(:,2),u_node, 'FaceColor', 'interp')
    axis([-1,1,-1,1,0,1])
    cmocean('thermal')
    title(sprintf('t = %.3f', t))
    caxis([0,1])
    drawnow
    print(sprintf('heat_%i_corrected',ii),'-dpng','-r600')
end

% Compute total "Mass" in the domain. Depending on the boundary conditions
% this may be conserved.
M=sum(u.*dmesh.tri.area,1)/M0;
figure
plot(tt,M)