x0=[0;0;0];
[t,x]=ode45('cellkin', [0,500], x0);
figure(1); clf;
subplot(3,1,1),plot(t,x(:,1),'k'), xlabel('t'), ylabel('Rn');
subplot(3,1,2),plot(t,x(:,2),'k'), xlabel('t'), ylabel('Rc');
subplot(3,1,3),plot(t,x(:,3),'k'), xlabel('t'), ylabel('P');

% Parameter fitting
% Need to find realistic constants to allow for reasonable RNA and protein half lives, expression profiles at steady state (comparing switches of different kcl, time to reach steady state, etc.
% We know that sTRSV has ~2% expression. So can assume that only have ~2% Rc left by the time you start making proteins

% Protein binding
x0=[];
x0(1)=1.25e-16; %5pM assuming 40kD
x0(2)=100e-9;
x0(3)=0;

[t,x]=ode45('binding', [0,500], x0);
figure(2); clf;
subplot(3,1,1),plot(t,x(:,1),'k'), xlabel('t'), ylabel('RNA');
subplot(3,1,2),plot(t,x(:,2),'k'), xlabel('t'), ylabel('Protein');
subplot(3,1,3),plot(t,x(:,3),'k'), xlabel('t'), ylabel('Complex');