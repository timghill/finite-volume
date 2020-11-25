function vprime=rhs_heat_unstructured_optimized(v,dmesh,params)
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

vx=zeros(size(v));      % x-component of gradient of v, defined on elements
vy=zeros(size(v));      % y-component of gradient of v, defined on elements

% First compute the gradient of v by integrating over triangle edges (e.g.
% by applying divergence theorem)
if strcmp(params.gradient, 'gg-hybrid') || strcmp(params.gradient, 'gg')
    for kk=1:dmesh.tri.n_edges
        nvec = [dmesh.tri.edge_nx(kk), dmesh.tri.edge_ny(kk)];
        r1 = dmesh.tri.edge_length(kk)/dmesh.tri.edge_area(kk,1);
        r2 = dmesh.tri.edge_length(kk)/dmesh.tri.edge_area(kk,2);
        if dmesh.tri.bmark_edge(kk)==0
            % For interior elements, calculate gradient according to
            % specified method
            if strcmp(params.gradient, 'gg-hybrid')

                stencil_els = dmesh.tri.edge_stencil{kk};
                edgex = dmesh.tri.edge_midpoints(kk, 1);
                edgey = dmesh.tri.edge_midpoints(kk, 2);

                neigh_els = dmesh.tri.connect_edge_el(kk,:);

                dx = dmesh.tri.elements(stencil_els, 1) - edgex;
                dy = dmesh.tri.elements(stencil_els, 2) - edgey;

                A = [ones(size(dx)), dx, dy];
                b = v(stencil_els);
                W=diag(1./sqrt(dx.^2 + dy.^2));
                W=W/sum(W(:));
                x_lsq = (W*A)\(W*b);
                v_bndry = x_lsq(1);

                if neigh_els(1)>0
                    vx(neigh_els(1)) = vx(neigh_els(1)) + r1*v_bndry*nvec(1);
                    vy(neigh_els(1)) = vy(neigh_els(1)) + r1*v_bndry*nvec(2);
                end

                if neigh_els(2)>0
                    vx(neigh_els(2)) = vx(neigh_els(2)) - r2*v_bndry*nvec(1);
                    vy(neigh_els(2)) = vy(neigh_els(2)) - r2*v_bndry*nvec(2);
                end


            elseif strcmp(params.gradient, 'gg')
                neigh_els = dmesh.tri.connect_edge_el(kk,:);
                    v1 = v(neigh_els(1));
                    v2 = v(neigh_els(2));
                    v_bndry = 0.5*(v1 + v2);

                    vx(neigh_els(1)) = vx(neigh_els(1)) + r1*v_bndry*nvec(1);
                    vy(neigh_els(1)) = vy(neigh_els(1)) + r1*v_bndry*nvec(2);

                    vx(neigh_els(2)) = vx(neigh_els(2)) - r2*v_bndry*nvec(1);
                    vy(neigh_els(2)) = vy(neigh_els(2)) - r2*v_bndry*nvec(2);
            end
        else
            neigh_els = dmesh.tri.connect_edge_el(kk,:);
            % On boundary elements, use simple GG method
            if neigh_els(1)>0 && neigh_els(2)<0
                v1 = v(neigh_els(1));
            elseif neigh_els(1)<0 && neigh_els(2)>0
                v1 = v(neigh_els(2));
                nvec = -nvec;
            end

            if strcmp(params.bc,'dirichlet')
                v2 = params.v_dirichlet;
            elseif strcmp(params.bc,'neumann')
                v2 = v1;
            elseif strcmp(params.bc,'flux')
                v2 = v1;

            end

            v_bndry = 0.5*(v1 + v2);

            if neigh_els(1)>0
                vx(neigh_els(1)) = vx(neigh_els(1)) + r1*v_bndry*nvec(1);
                vy(neigh_els(1)) = vy(neigh_els(1)) + r1*v_bndry*nvec(2);
            else
                vx(neigh_els(2)) = vx(neigh_els(2)) - r2*v_bndry*nvec(1);
                vy(neigh_els(2)) = vy(neigh_els(2)) - r2*v_bndry*nvec(2);
            end
        end
    end
elseif strcmp(params.gradient, 'lsq')
    for ii=1:dmesh.tri.n_elements
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

% Efficiently compute flux on elements
f_x = -params.gamma*(abs(v)).^params.alpha.*(abs(vx)).^params.beta.*sign(vx);
f_y = -params.gamma*(abs(v)).^params.alpha.*(abs(vy)).^params.beta.*sign(vy);
vprime= zeros(dmesh.tri.n_elements, 1);
for kk=1:dmesh.tri.n_edges
    nvec = [dmesh.tri.edge_nx(kk), dmesh.tri.edge_ny(kk)];
    r1 = dmesh.tri.edge_length(kk)/dmesh.tri.edge_area(kk,1);
    r2 = dmesh.tri.edge_length(kk)/dmesh.tri.edge_area(kk,2);
    neigh_els = dmesh.tri.connect_edge_el(kk,:);

    if dmesh.tri.bmark_edge(kk)==0
        fx1 = f_x(neigh_els(1));
        fy1 = f_y(neigh_els(1));
        fx2 = f_x(neigh_els(2));
        fy2 = f_y(neigh_els(2));
    else
        if neigh_els(1)>0
            fx1 = f_x(neigh_els(1));
            fy1 = f_y(neigh_els(1));
            
            if strcmp(params.bc,'dirichlet') % Want rate of change to be zero
                fx2=fx1;
                fy2=fy1;
            elseif strcmp(params.bc,'neumann') % Generally, also want rate of change to be zero
                fx2=-fx1;
                fy2=-fy1;
            elseif strcmp(params.bc,'flux')
                fx2 = fx1;
                fy2 = fy1;
            else
                disp('BCs not found')
            end
            
        else
            fx2 = f_x(neigh_els(2));
            fy2 = f_y(neigh_els(2));
            
            if strcmp(params.bc,'dirichlet') % Want rate of change to be zero
                fx1=fx2;
                fy1=fy2;
            elseif strcmp(params.bc,'neumann') % Generally, also want rate of change to be zero
                fx1=-fx2;
                fy1=-fy2;
            elseif strcmp(params.bc,'flux')
                fx1 = fx2;
                fy1 = fy2;
            else
                disp('BCs not found')
            end
        end
    end
    fx = 0.5*(fx1 + fx2);
    fy = 0.5*(fy1 + fy2);
    f = (fx*nvec(1) + fy*nvec(2));
    
    
    if neigh_els(1)>0
        vprime(neigh_els(1)) = vprime(neigh_els(1)) - r1*f;
    end
    if neigh_els(2)>0
        vprime(neigh_els(2)) = vprime(neigh_els(2)) + r2*f;
    end
    
end

if strcmp(params.bc, 'dirichlet')
    vprime(sum(dmesh.tri.bmark_el,2)>0)=0;
end

if params.derivs
    vderiv = vprime;
    vprime = struct;
    vprime.vx = vx;
    vprime.vy = vy;
    
    vprime.fx = f_x;
    vprime.fy = f_y;
    
    vprime.vprime = vderiv;
end