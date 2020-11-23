% heat_equation.m solves the two-dimensional nonlinear heat equation on an
% unstructured triangular mesh using finite volume methods. The solver uses
% first-order FV discretization in space with built-in matlab
% time-stepping. The difference between this and ../heat/heat_equation.m is
% that here the flux is gamma*u*grad.u instead of just gamma*grad.u.

% dmesh=load('../meshes/circ_mesh.mat');
dmesh=load('../meshes/phys_mesh.mat');

%% Setup
% Time-stepping
tend=1e4;
t=0;
tspan=linspace(t,tend,51);

% Parameters
params.gamma=5e-1;
params.alpha=5/4;
params.beta=3/2;

% The solver supports three types of boundary conditions, specified
% using params.bc:
%   (1) 'dirichlet':    Boundary has fixed value params.v_dirichlet
%   (2) 'neumann':      Boundary has no heat flux out of the domain
%                       (grad.v=0 on boundary) [These BC
%   (3) 'flux':         Allow flux out of domain by enforcing
%                       d(grad.v)/dn=0 on the boundary.
% Note that 'neumann' boundary conditions conserve thermal energy, while
% the other BCs allow heat to flow out of the domain through the boundary

params.bc='flux';
params.v_dirichlet=1e-3;   % Must specify value for Dirichlet condition

%% Initial conditions
u=1e-3*ones(size(dmesh.tri.elements(:,1)));
trix=dmesh.tri.elements(:,1);
triy=dmesh.tri.elements(:,2);
trinorm=sqrt(trix.^2+triy.^2);
u(trix>1e3 & trix<2e3 & triy>-250 & triy<250)=2;
% u(trinorm<0.5)=1;


M0=sum(u.*dmesh.tri.area);

%% Solver
odefun=@(t,y) rhs_nonlinear_heat(t,y,dmesh,params,'solver','area');
opts=odeset('RelTol',1e-2,'AbsTol',1e-5);
tic;
[tt,yout]=ode45(odefun,tspan,u,opts);
toc;
u=yout';

%% Post-processing
% Animate the outputs
for ii=1:length(tspan)
    t=tspan(ii);
    figure(1)
    u_node=interp_el_node(dmesh,u(:,ii));
    if strcmp(params.bc,'dirichlet')
        u_node(dmesh.tri.bmark==1)=0;
    end
    trisurf(dmesh.tri.connect,dmesh.tri.nodes(:,1),dmesh.tri.nodes(:,2),u_node)
%     axis([-1,1,-1,1,-0.25,1.25])
    colormap('winter')
    title(sprintf('t = %.3f', t))
    caxis([-0.25,1.25])
    drawnow
    pause(0.1)
end


kk=[4,18];
for ii=kk
    t=tspan(ii);
    figure
    u_node=interp_el_node(dmesh,u(:,ii));
    if strcmp(params.bc,'dirichlet')
        u_node(dmesh.tri.bmark==1)=0;
    end
    trisurf(dmesh.tri.connect,dmesh.tri.nodes(:,1),dmesh.tri.nodes(:,2),u_node)
%     axis([-1,1,-1,1,0,1])
    cmocean('thermal')
    title(sprintf('t = %.3f', t))
    caxis([0,1])
    drawnow
    print(sprintf('nonlinear_heat_%s_%i',params.bc,ii),'-dpng','-r600')
end

M=sum(u.*dmesh.tri.area,1)/M0;
figure
plot(tt,M)

% figure
% vprime=rhs_nonlinear_heat(t,yout(end,:)',dmesh,params);
% u_node=interp_el_node(dmesh,vprime);
% if strcmp(params.bc,'dirichlet')
%     u_node(dmesh.tri.bmark==1)=0;
% end
% trisurf(dmesh.tri.connect,dmesh.tri.nodes(:,1),dmesh.tri.nodes(:,2),u_node)
% %     axis([-1,1,-1,1,0,1])
% cmocean('balance')
% title(sprintf('t = %.3f', t))
% caxis([0,1])

figure
subplot(2,1,1)
hold on
% Make some fake data
uu=1*0.03*trix;
vprime=rhs_nonlinear_heat(0,uu,dmesh,params,'test','first-order');
u_node=interp_el_node(dmesh,vprime.vx);
trisurf(dmesh.tri.connect,dmesh.tri.nodes(:,1),dmesh.tri.nodes(:,2),u_node)
%     axis([-1,1,-1,1,0,1])
cmocean('balance')
caxis([0.03-0.02,0.03+0.02])
colorbar
axis image
title(sprintf('t = %.3f', t))

subplot(2,1,2)
hold on
% Make some fake data
uu=1*0.03*trix;
vprime=rhs_nonlinear_heat(0,uu,dmesh,params,'test','area');
u_node=interp_el_node(dmesh,vprime.vx);
trisurf(dmesh.tri.connect,dmesh.tri.nodes(:,1),dmesh.tri.nodes(:,2),u_node)
%     axis([-1,1,-1,1,0,1])
cmocean('balance')
caxis([0.03-0.02,0.03+0.02])
colorbar
axis image
title(sprintf('t = %.3f', t))
