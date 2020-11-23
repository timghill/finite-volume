% shallow_water_unstructured.m solves the shallow water equations on an
% unstructured triangular mesh using FV methods. The solver uses first
% order Lax-Friedrichs type flux spatially, the built-in matlab time stepper
% ode 45. The code has three cases set up:
%  1: Square domain, initial "dam" in center of domain
%  2: Circular domain with initial circular dam
%  3: Wave propagating through domain with constriction in the middle

casenum=1;

%% Setup
% Time-stepping
tend=3;
t=0;
tspan=linspace(t,tend,51);

%% Initial conditions, in rectangular coordinates, then interp onto FV mesh

if casenum==1 % Square bulge
    dmesh=load('../meshes/rect_mesh.mat');

    h=ones(size(dmesh.tri.elements(:,1)));
    trix=dmesh.tri.elements(:,1);
    triy=dmesh.tri.elements(:,2);
    h(trix>-0.5 & trix<0.5 & triy>-0.5 & triy<0.5)=2;
    
elseif casenum==2 % Circular dam break
    dmesh=load('../meshes/circ_mesh.mat');

    h=ones(size(dmesh.tri.elements(:,1)));
    trix=dmesh.tri.elements(:,1);
    triy=dmesh.tri.elements(:,2);
    trinorm=sqrt(trix.^2 + triy.^2);
    h(trinorm<0.5)=2;
    
elseif casenum==3 % Hyperbolic constriction
    dmesh=load('../meshes/hyp_mesh.mat');

    h=ones(size(dmesh.tri.elements(:,1)));
    trix=dmesh.tri.elements(:,1);
    h(trix<-0.75)=2;
end

% mx and my are easy -- zeros
mx=zeros(dmesh.tri.n_elements,1);
my=zeros(dmesh.tri.n_elements,1);

% Initial state vector
v0=[h;mx;my];

% Compute initial total energy and mass for later analysis
E0=sum(h.^2.*dmesh.tri.area);
M0=sum(h.*dmesh.tri.area);

%% Solver
params=struct;  % We don't have parameters for SW equations since we set
                % g=0, but we will have parameters for other models
odefun=@(t,v) rhs_sw_unstructured(t,v,dmesh,params);

tic;
[tt,yout]=ode45(odefun,tspan,v0);
toc;

%% Post-processing

% Unpack model outputs
N=dmesh.tri.n_elements;
h=yout(:,1:N)';
mx=yout(:,N+1:2*N)';
my=yout(:,2*N+1:end)';

%% Animate the outputs
for ii=1:length(tspan)
    t=tspan(ii);
    figure(1)
    h_node=interp_el_node(dmesh,h(:,ii));
    trisurf(dmesh.tri.connect,dmesh.tri.nodes(:,1),dmesh.tri.nodes(:,2),h_node)
    axis([-1,1,-1,1,0.5,2.5])
    colormap('winter')
    title(sprintf('t = %.3f', t))
    caxis([0.5,2.5])
    drawnow
    pause(0.1)
end

kk=[4,18,28];
for ii=kk
    t=tspan(ii);
    figure
    h_node=interp_el_node(dmesh,h(:,ii));
    trisurf(dmesh.tri.connect,dmesh.tri.nodes(:,1),dmesh.tri.nodes(:,2),h_node,'EdgeColor','none');%,'FaceColor','interp')
    axis([-1,1,-1,1,0.75,2])
    colormap(palettes('blue-4'))
    title(sprintf('t = %.3f', t))
    caxis([0.75,2])
    drawnow
    print(sprintf('sw_solution_%i',ii),'-dpng','-r600')
%     pause(0.1)
end

%% Analysis - Compute total energy and mass
M=sum(h.*dmesh.tri.area,1)'./M0;
E=sum( h.^2.*dmesh.tri.area + 0.5*(mx.^2./h + my.^2./h).*dmesh.tri.area,1)'./E0;

figure
hold on
plot(tt,M)
plot(tt,E)
xlabel('Time (s)')

legend({'Mass','Energy'},'box','off','Location','east')

print('sw_mass_energy_conservation','-dpng','-r600')