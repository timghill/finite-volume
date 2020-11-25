function [T,Y] = odeRK(odefun,tspan,dt,v0,varargin)
% A simple RK4 timestepping scheme built to use the same syntax as the
% built-in matlab ode integrators (e.g. ode45), but with an additional
% argument for the timestep dt. Ignores any additional arguments (e.g.
% 'Stats', opts, etc.)

% Definition of runge-kutta constants
Y=zeros(length(tspan),length(v0));
T=zeros(length(tspan),1);

Y(1,:)=v0';

% Take first step
t=tspan(1);
v=v0 + dt*odefun(t,v0);
t=t+dt;

jj=2; % Index for stepping through tspan array
while t < tspan(end)
    % Determine next time step
    if t+dt >= tspan(jj)
        dt_ii=tspan(jj)-t;
        save_state=true;
    else
        dt_ii=dt;
        save_state=false;
    end
    
    % This is the advancement of the state    
    k1 = odefun(t, v);
    k2 = odefun(t + dt_ii/2, v + dt_ii*k1/2);
    k3 = odefun(t + dt_ii/2, v + dt_ii*k2/2);
    k4 = odefun(t + dt_ii, v + dt_ii*k3);
    
    v = v + dt_ii*(k1 + 2*k2 + 2*k3 + k4)/6;
    t=t+dt_ii;
    
    if save_state
        Y(jj,:) = v';
        T(jj)=t;
        jj=jj+1;
    end
    
end