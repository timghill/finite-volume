function vprime=rhs_sw_unstructured(t,v,dmesh,params)
% function vprime=rhs_sw_unstructured computes the spatial derivatives of
% v=[h;mx;my] on unstructured triangular mesh dmesh using first-order
% Lax-Friedrichs flux function. Note that params is not used, but this type
% of function is what we will use later for the full model that has
% numerous parameters to keep track of

% Unpack the state vector
N=dmesh.tri.n_elements;
h=v(1:N);
mx=v(N+1:2*N);
my=v(2*N+1:end);

hprime =zeros(N,1);
mxprime=zeros(N,1);
myprime=zeros(N,1);

% Loop over the elements to compute spatial derivatives
for ii=1:N
       for q=1:3
           nvec=[dmesh.tri.nx(ii,q),dmesh.tri.ny(ii,q)];
           tvec=-[-nvec(2),nvec(1)];

            adj_i=dmesh.tri.connect_el_el(ii,q);
            
            h1=h(ii);
            mn1=nvec(1)*mx(ii) + nvec(2)*my(ii);
            mt1=tvec(1)*mx(ii) + tvec(2)*my(ii);

            if adj_i==-1
                h2=h1;
                mn2=-mn1; % This is the ideal wall condition!
                mt2=mt1;
            else
                h2=h(adj_i);
                mn2=nvec(1)*mx(adj_i) + nvec(2)*my(adj_i);
                mt2=tvec(1)*mx(adj_i) + tvec(2)*my(adj_i);
            end

            % Construct center and adjacent state vectors
            u1=[h1;mn1;mt1];
            u2=[h2;mn2;mt2];
            [f1,fn,ft]=fv_sw_single_bndry(u1,u2);

            area=dmesh.tri.area(ii);
            l=dmesh.tri.ds(ii,q);

            fx=nvec(1)*fn + tvec(1)*ft;
            fy=nvec(2)*fn + tvec(2)*ft;

            hprime(ii )=hprime(ii ) - f1*l/area;
            mxprime(ii)=mxprime(ii) - fx*l/area;
            myprime(ii)=myprime(ii) - fy*l/area;
       end
end

vprime=[hprime;mxprime;myprime];