function xdot = cellkin(t,x);
%
% Solves a set of differential equations describing RNA turnover including
% ribozyme cleavage, etc, and carried through to protein translation.
%
% [t,x] = ode45(‘cellkin’,[0 150],x0)
% integrates from t = 0 to t = 150 sec, with initial conditions
% Rn0 = x0(1), Rc0 = x0(2), and P0=x0(3) and x0 is a column vector
% which needs to be declared in matlab prompt, initially as [0;0;0]
% Nov 6 2012 Leo d'Espaux


% convert variables to x(i) notation
    Rn = x(1);
    Rc = x(2);
    P = x(3);
  

% now all the constants

% Initial nuclear RNA. Made up 
    Rn0 = 25;
% rate of elongation in nucleotides/sec div by length of 3’UTR
    kelong = 25*(1/2000);      
% rate of rz cleavage in /sec 
    kcl=0.5/60; 
% rate of nuclear RNA decay in /sec. Made up
    kndecay=0.01;          
% rate of RNA processing, transport in /sec. Made up  
    kproc=0.000001;                      
% rate of RNA cytoplasmic RNA decay in /sec. Made up 
    kcdecay=0.01;           
% rate of protein translation in /sec
% 8AA/sec, 233AA/protein, 0.64ribosomes/100nt, give .153 prot/sec
% slow down to allow for ribosome loading delay
    ktrans=(0.0001)*0.153;
% rate of protein decay in /min. made up
    kpdecay=0.01;

    
% the modeling equations are:
  dRndt = kelong*Rn0 - kcl*Rn - kndecay*Rn - kproc*Rn;
  dRcdt = kproc*Rn - kcl*Rc - kcdecay*Rc;
  dPdt = ktrans*Rc - kpdecay*P;
  
 
% now, create the column vector of state derivatives
  xdot = [dRndt;dRcdt;dPdt];
% end of file cellkin.m
