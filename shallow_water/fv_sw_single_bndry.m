function [f1,fn,ft]=fv_sw_single_bndry(u1,u2)
% Compute numerical fluxes for 2D shallow water equations, where
% ui=[h,m_normal,m_tangent], where normal and tangent are in reference to
% the outward unit normal of the cell.

h1=u1(1);
mn1=u1(2);
mt1=u1(3);

h2=u2(1);
mn2=u2(2);
mt2=u2(3);

% This defines the equations
f_h1=mn1;
f_n1=mn1^2/h1 + 0.5*h1^2;
f_t1=mn1*mt1/h1;

f_n1(h1==0)=0;
f_t1(h1==0)=0;

f_h2=mn2;
f_n2=mn2^2/h2 + 0.5*h2^2;
f_t2=mn2*mt2/h2;

f_n2(h2==0)=0;
f_t2(h2==0)=0;

% Compute eigenvalues
lambda_max=max(abs(mn1/h1)+sqrt((h1)),abs(mn2/h2)+sqrt((h2)));
% lambda_

% lambda_max=0;

% Compute the numerical fluxes
f1=0.5*(f_h1+f_h2) - 0.5*lambda_max*(h2-h1);
fn=0.5*(f_n1+f_n2) - 0.5*lambda_max*(mn2-mn1);
ft=0.5*(f_t1+f_t2) - 0.5*lambda_max*(mt2-mt1);