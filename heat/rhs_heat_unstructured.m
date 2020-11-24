function vprime=rhs_heat_unstructured(v,dmesh,params)
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
vx=zeros(size(v));      % x-component of gradient of v, defined on elements
vy=zeros(size(v));      % y-component of gradient of v, defined on elements

% First compute the gradient of v by integrating over triangle edges (e.g.
% by applying divergence theorem)
for ii=1:dmesh.tri.n_elements
    for kk=1:3
        nvec=[dmesh.tri.nx(ii,kk),dmesh.tri.ny(ii,kk)];
        r=dmesh.tri.ds(ii,kk)/dmesh.tri.area(ii);
        iEdge=dmesh.tri.connect_el_edge(ii,kk);

        if strcmp(params.gradient, 'gg-hybrid')
            neigh_els=dmesh.tri.edge_stencil{iEdge};

            edgex=dmesh.tri.edge_midpoints(iEdge,1);
            edgey=dmesh.tri.edge_midpoints(iEdge,2);

            dx=dmesh.tri.elements(neigh_els,1)-edgex;
            dy=dmesh.tri.elements(neigh_els,2)-edgey;

            A = [ones(size(dx)), dx, dy];
            b = v(neigh_els);
            W=diag(1./sqrt(dx.^2 + dy.^2));
            W=W/sum(W(:));
            x_lsq = (W*A)\(W*b);
            v_bndry = x_lsq(1);
            vx(ii) = vx(ii) + r*v_bndry*nvec(1);
            vy(ii) = vy(ii) + r*v_bndry*nvec(2);

        elseif strcmp(params.gradient, 'gg')
            v1 = v(ii);
            adj_i = dmesh.tri.connect_el_el(ii,kk);

            if adj_i>0
                v2 = v(adj_i);

            else        % Apply boundary conditions
                if strcmp(params.bc,'dirichlet')
                    v2=params.v_dirichlet;
                elseif strcmp(params.bc,'neumann')
                    v2=v1;
                elseif strcmp(params.bc,'flux')
                    v2=v1;
                end
            end

            v_bndry=0.5*(v1+v2);
            vx(ii) = vx(ii) + r*v_bndry*nvec(1);
            vy(ii) = vy(ii) + r*v_bndry*nvec(2);
        end
           
    end
    
    if strcmp(params.gradient, 'lsq')
        neigh_els = dmesh.tri.node_stencil_extended{ii};

        elx = dmesh.tri.elements(ii, 1);
        ely = dmesh.tri.elements(ii, 2);

        neighx = dmesh.tri.elements(neigh_els, 1);
        neighy = dmesh.tri.elements(neigh_els, 2);

        dx = neighx - elx;
        dy = neighy - ely;

        A = [dx, dy];
        b = v(neigh_els) - v(ii);
        W = diag(1./sqrt(dx.^2 + dy.^2));
        W = W/sum(W(:));
        x_lsq = (W*A)\(W*b);

        vx(ii) = x_lsq(1);
        vy(ii) = x_lsq(2);
    end
end

% Now that we have the gradient, compute the numerical flux across the
% boundaries
for ii=1:dmesh.tri.n_elements
    bndry_el=false;
    f=0;
    for q=1:3
        nvec=[dmesh.tri.nx(ii,q),dmesh.tri.ny(ii,q)];
        r=dmesh.tri.ds(ii,q)/dmesh.tri.area(ii);
        
        vx1=vx(ii);
        vy1=vy(ii);
        v1 = v(ii);
        adj_i=dmesh.tri.connect_el_el(ii,q);
        if adj_i==-1    % Apply boundary conditions
            if strcmp(params.bc,'dirichlet')
                vx2=+vx1;
                vy2=+vy1;
                v2 = v1;
            elseif strcmp(params.bc,'neumann')
                vx2=-vx1;
                vy2=-vy1;
                v2 = v1;
            elseif strcmp(params.bc,'flux')
                vx2=vx1;
                vy2=vy1;
                v2 = v1;
            end
                
            bndry_el=true;
        else
            v2 = v(adj_i);
            vx2=vx(adj_i);
            vy2=vy(adj_i);
        end
        
        
        f_x1=-params.gamma*(abs(v1))^params.alpha*(abs(vx1))^params.beta*sign(vx1);
        f_y1=-params.gamma*(abs(v1))^params.alpha*(abs(vy1))^params.beta*sign(vy1);
        
        f_x2=-params.gamma*(abs(v2))^params.alpha*(abs(vx2))^params.beta*sign(vx2);
        f_y2=-params.gamma*(abs(v2))^params.alpha*(abs(vy2))^params.beta*sign(vy2);
        
        fx=0.5*(f_x1 + f_x2);
        fy=0.5*(f_y1 + f_y2);
        
        f = f - (fx*nvec(1) + fy*nvec(2))*r;
        
        if strcmp(params.bc,'dirichlet') && bndry_el
            f = 0;
        end
    end
    vprime(ii) = f;
end

if params.derivs
    vprime = struct;
    vprime.vx = vx;
    vprime.vy = vy;
end