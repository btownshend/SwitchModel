function xdot = binding(t,x);
%
% Solves a set of differential equations describing RNA turnover including
% ribozyme cleavage, etc, and carried through to protein translation, etc.
%
% [t,x] = ode45(?cellkin?,[0 150],x0)
% integrates from t = 0 to t = 150 sec, with initial conditions
% Rn0 = x0(1), Rc0 = x0(2), and P0=x0(3) and x0 is a column vector
% Nov 6 2012 Leo d'Espaux
% based on http://homepages.rpi.edu/~bequeb/courses/cpc/modules/ode01.pdf
%
% since the states are passed to this routine in the x vector,
% convert to natural notation
    R = x(1);
    P = x(2);
    C = x(3);
  
 
% rates in units of mol sec etc
    ka=0;
    kd=1.9e-3;
    kcl=0.5/60;
    
    
% the modeling equations are:
    dRdt = -kcl*R - ka*R*P + kd*C;
    dPdt = -ka*P*R + kd*C;
    dCdt = ka*P*R - kd*C;
  
 
% now, create the column vector of state derivatives
    xdot = [dRdt;dPdt;dCdt];
    
% end of file
