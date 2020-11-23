function vprime=rhs_heat_unstructured(t,v,dmesh,params)
% function vprime=rhs_heat_unstructured computes the spatial discretization
% of the heat equation using first-order finite volume methods on an
% unstructured triangular mesh. The solver supports three types of boundary
% conditions, specified using params.bc:
%   (1) 'dirichlet':    Boundary has fixed value params.v_dirichlet
%   (2) 'neumann':      Boundary has no heat flux out of the domain
%                       (grad.v=0 on boundary)
%   (3) 'flux':         Allow flux out of domain by enforcing
%                       d(grad.v)/dn=0 on the boundary.
% This function is written to pass to the matlab solvers as
% odefun=@(t,y) rhs_heat_unstructured(t,y,dmesh,params)
% [tout,yout]=odeXXX(odefun,tspan,v0)

vprime=zeros(size(v));  % Time derivative of state variable v
vx=zeros(size(v));      % x-component of gradient of v
vy=zeros(size(v));      % y-component of gradient of v

% First compute the gradient of v by integrating over triangle edges (e.g.
% by applying divergence theorem)
for ii=1:dmesh.tri.n_elements
    v_gradx=0;
    v_grady=0;
    for q=1:3
        nvec=[dmesh.tri.nx(ii,q),dmesh.tri.ny(ii,q)];
        r=dmesh.tri.ds(ii,q)/dmesh.tri.area(ii);
        
        v1=v(ii);
        adj_i=dmesh.tri.connect_el_el(ii,q);
        
        if adj_i==-1    % Apply boundary conditions
            if strcmp(params.bc,'dirichlet')
                v2=params.v_dirichlet;
            elseif strcmp(params.bc,'neumann')
                v2=v1;
            elseif strcmp(params.bc,'flux')
                v2=v1;
            end
            
            v_bndry=0.5*(v1+v2);
            v_gradx=v_gradx+v_bndry*nvec(1)*r;
            v_grady=v_grady+v_bndry*nvec(2)*r;
        else
            v2=v(adj_i);
            v_bndry=0.5*(v1+v2);
            v_gradx=v_gradx+v_bndry*nvec(1)*r;
            v_grady=v_grady+v_bndry*nvec(2)*r;
        end
        
    end
    vx(ii)=v_gradx;
    vy(ii)=v_grady;
end

% Now that we have gradient, compute the numerical flux across the
% boundaries
for ii=1:dmesh.tri.n_elements
    bndry_el=false;
    
    f=0;
    for q=1:3
        nvec=[dmesh.tri.nx(ii,q),dmesh.tri.ny(ii,q)];
        r=dmesh.tri.ds(ii,q)/dmesh.tri.area(ii);
        
        vx1=vx(ii);
        vy1=vy(ii);
        adj_i=dmesh.tri.connect_el_el(ii,q);
        if adj_i==-1    % Apply boundary conditions
            if strcmp(params.bc,'dirichlet')
                vx2=+vx1;
                vy2=+vy1;
            elseif strcmp(params.bc,'neumann')
                vx2=-vx1;
                vy2=-vy1;
            elseif strcmp(params.bc,'flux')
                vx2=vx1;
                vy2=vy1;
            end
                
            bndry_el=true;
        else
            vx2=vx(adj_i);
            vy2=vy(adj_i);
        end
        
        f_x1=params.gamma*vx1;
        f_y1=params.gamma*vy1;
        
        f_x2=params.gamma*vx2;
        f_y2=params.gamma*vy2;
        
        fx=0.5*(f_x1 + f_x2);
        fy=0.5*(f_y1 + f_y2);
        
        f=f+(fx*nvec(1) + fy*nvec(2))*r;
        
        if strcmp(params.bc,'dirichlet') && bndry_el
            f=0;
        end
    end
    vprime(ii)=f;
end