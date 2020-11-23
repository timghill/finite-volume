% shallow_water_2d.m simulates 2D shallow water equations for a Riemann
% problem where there is initially a bulge of fluid in the center of the
% domain. Uses Lax-Friedrichs type flux.

% Spatial parameters
xmin=-1;
xmax=1;
ymin=-1;
ymax=1;

m=25;   % Number of cells in each direction

dx=(xmax-xmin)/m;
dy=(ymax-ymin)/m;
area=dx*dy;

x=xmin:dx:xmax-dx;
y=ymin:dy:ymax-dy;

[xx,yy]=meshgrid(x,y);  % 2D coordinate grids for plotting later on

% Temporal parameters
tend=3;
t=0;
c=0.8;      % CFL safety constant (should be <1)

% Plotting parameters
% T is an array of times to make a plot at. The time stepping loop will
% adjust time step to simulate system at exactly these times, and will show
% a plot of h at these times after the simulation runs
T=[1.35, 3];

% Script can optionally display an animation of simulation
play_animation=true;

% Initial conditions, in terms of conserved quantitites, with ghost cells
h = ones(m,m);
% h(m/4+2:3*m/4+1, m/4+2:3*m/4+1)=2;
h(xx>yy)=2;
mx=zeros(m,m);
my=zeros(m,m);

h_avg=mean(h(:));

h_new=h;
mx_new=mx;
my_new=my;

% h_out will be shape (n_times, nx, ny), storing h values where you have
% specified to make plots
h_out=[];

% tt is the actual times used in the simulation, and h_total is the height
% integrated over the domain (which equals the total volume of fluid) at
% each time step
tt=[];
h_total=[];

tj=1;   % Index for keeping track of outputs at times given in T

% Fixed cell normal vector components
nx=[1;0;-1;0];
ny=[0;1;0;-1];

for q=1:4
    nvec=[nx(q);ny(q)]
    tvec=[-ny(q);nx(q)]
end

while t<tend
    h_new=h;
    mx_new=mx;
    my_new=my;
% for i=[1]
    save_output=false;
    
    % Unpack physical variables
    u=mx./h;
    v=my./h;
    
    % Compute eigenvalues and timestep
    lambda_x=abs(u)+sqrt(h);
    lambda_y=abs(v)+sqrt(h);
    dt=0.5*c*min( min(dx./lambda_x(:)), min(dy./lambda_y(:)));
    
    % This conditional checks if we are near one of the times we need to
    % save for plotting
    if t+dt>T(tj)
        dt = T(tj)-t;
        tj=tj+1;
        save_output=true;
    end
    
    % This conditional checks if we are near the end of the simulation time
    if t+dt>tend
        dt=tend-t;
    end
    
    t=t+dt;
    
%     lambda=abs(mx./h)+sqrt(h);
%     lambda_max_plus=max(lambda(2:end-1,3:end), lambda(2:end-1,2:end-1));
%     lambda_max_plus=max(lambda_max_plus(:));
%     lambda_max_minus=max(lambda(2:end-1,1:end-2), lambda(2:end-1,2:end-1));
%     lambda_max_minus=max(lambda_max_minus(:));

    for ii=1:m
        for jj=1:m
            for q=1:4
                % Compute normal and tangent vectors, using right handed
                % coordinate system such that nvec cross tvec = +1.
                nvec=[nx(q);ny(q)];
                tvec=[-ny(q);nx(q)];
                
                % Compute indices of adjacent element
                adj_ii=ii+nvec(2);
                adj_jj=jj+nvec(1);
                
                % Compute center and adjacent normal and tangent momentum
                % components
                h1=h(ii,jj);
                mn1=nvec(1)*mx(ii,jj) + nvec(2)*my(ii,jj);
                mt1=tvec(1)*mx(ii,jj) + tvec(2)*my(ii,jj);
                
                if adj_ii<1 || adj_ii>m || adj_jj<1 || adj_jj>m
                   mn2=-mn1;
                   mt2=mt1;
                   h2=h1;
                else
                    h2=h(adj_ii,adj_jj);
                    mn2=nvec(1)*mx(adj_ii,adj_jj) + nvec(2)*my(adj_ii,adj_jj);
                    mt2=tvec(1)*mx(adj_ii,adj_jj) + tvec(2)*my(adj_ii,adj_jj);
                end
                % Construct center and adjacent state vectors
                u1=[h(ii,jj);mn1;mt1];
                u2=[h2;mn2;mt2];
                
%                 f=dt/dx;
%                 f=dx/dt;
%                 f=dx/dt;
                % FV discretization
                [f1,fn,ft]=fv_sw_single_bndry(u1,u2);
                
                % Compute x and y components of flux from normal and
                % tangent components
                fx=nvec(1)*fn+tvec(1)*ft;
                fy=nvec(2)*fn+tvec(2)*ft;
                
                % Update the state -- could use a built-in matlab solver
                % instead of this time loop
                h_new(ii,jj)=h_new(ii,jj) - dt*f1*dx/area;
                mx_new(ii,jj)=mx_new(ii,jj) - dt*fx*dx/area;
                my_new(ii,jj)=my_new(ii,jj) - dt*fy*dx/area;
                
            end
        end
    end
    
    h=h_new;
    mx=mx_new;
    my=my_new;
    
%     my(end,12)
    
    %% Boundary conditions
    % we use ghost cells and reverse the normal component of momentum in
    % the ghost cells to approximate no-normal-flux boundary conditions
    
    % Boundaries with normal in x-direction
%     mx(:,1)=-mx(:,2);
%     mx(:,end)=-mx(:,end-1);
%     
%     my(:,1)=my(:,2);
%     my(:,end)=my(:,end-1);
%     
%     h(:,1)=h(:,2);
%     h(:,end)=h(:,end-1);
%     
%     % Boundaries with normal in y-direction
%     my(1,:)=-my(2,:);
%     my(end,:)=-my(end-1,:);
%     mx(1,:)=mx(2,:);
%     mx(end,:)=mx(end-1,:);
%     
%     h(1,:)=h(2,:);
%     h(end,:)=h(end-1,:);
    
    % If we determined we are near a specified output time, save the
    % output height
    if save_output
        h_out=cat(3,h_out,h);
    end
    
    % Make an animation if we specified we want one.
    if play_animation
        figure(1)
        mesh(xx,yy,h)
        title(sprintf('t = %.3f', t))
        axis([-1, 1, -1, 1, h_avg-1, h_avg+1])
        caxis([h_avg-1, h_avg+1])
        colormap(palettes('blue-green-sat-on-ends'))
        drawnow
    end
    
    % Compute total volume of liquid to check conservation
    h_total=[h_total sum(h(:)*dx*dy)];
    tt=[tt t];
end
% 
% Make the figures
for i=1:length(T)
    figure
    mesh(xx,yy,h_out(:,:,i))
    title(sprintf('h (t = %.3f)', T(i)))
    xlabel('x')
    ylabel('y')
    colormap(palettes('blue-green-sat-on-ends'))
    drawnow
end

figure
plot(tt,h_total)
grid on
xlabel('Time')
ylabel('Volume of water')

% plot(mx(30,:))
