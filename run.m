x0=[0;0;0];   % Initial conditions
timerange=[0,500];   % Time (in minutes)

% Parameters of model
% Each is defined with a 3-element vector consisting of lowest possible, best estimate, highest possible
rng=[];
% Initial nuclear RNA. Made up 
rng.Rn0 = [25,25,25];
% rate of elongation in nucleotides/sec div by length of 3’UTR
rng.kelong = [25,25,25]*(1/2000);      
% rate of rz cleavage in /sec 
rng.kcl=[0.5,0.5,0.5]/60; 
% rate of nuclear RNA decay in /sec. Made up
rng.kndecay=[0.01,.01,.01];          
% rate of RNA processing, transport in /sec. Made up  
rng.kproc=[1e-6,1e-6,1e-6];
% rate of RNA cytoplasmic RNA decay in /sec. Made up 
rng.kcdecay=[1e-2,1e-2,1e-2];           
% rate of protein translation in /sec
% 8AA/sec, 233AA/protein, 0.64ribosomes/100nt, give .153 prot/sec
% slow down to allow for ribosome loading delay
rng.ktrans=[(0.0001)*0.153,(0.0001)*0.153,(0.0001)*0.153];
% rate of protein decay in /min. made up
rng.kpdecay=[0.001,0.01,0.1]/60;

% Form matrix for each of the parameter points
% Vectors are [Rn,Rc,P]
ntrials=10;
t50=[];
odeoptions=odeset('MaxStep',1);   % maximum 1 second steps
for trial=1:ntrials
  pick=struct('Rn0',rpick(rng.Rn0),'kelong',rpick(rng.kelong),'kcl',rpick(rng.kcl),'kndecay',rpick(rng.kndecay),...
              'kproc',rpick(rng.kproc),'kcdecay',rpick(rng.kcdecay),'ktrans',rpick(rng.ktrans),'kpdecay',rpick(rng.kpdecay));
  A=[-pick.kcl-pick.kndecay-pick.kproc, 0, 0
     pick.kproc, -pick.kcl-pick.kcdecay, 0
     0, pick.ktrans, -pick.kpdecay];
  c=[pick.kelong*pick.Rn0,0,0]';
  [t,x]=ode45(@(t,x) A*x + c, timerange, x0,odeoptions);

  % Plot results
  figure; clf;
  col='rgbcym';
  variables={'[RNA_{nuclear}]','[RNA_{cyto}]','[Protein]'};
  for i=1:size(x,2)
    semilogy(t,x(:,i),col(i));
    hold on;
    xlabel('Time (sec)');
    ylabel('Concentration');
  end
  legend(variables);
  title(sprintf('Trial %d',trial));
  % Find t50
  t50(trial)=t(find(x(:,3)>x(end,3)/2,1));
  fprintf('Trial %d: Rn_0=%g, k_{elong}=%g, k_{cl}=%g, k_{ndecay}=%g, k_{proc}=%g, k_{cdecay}=%g, k_{trans}=%g, k_{pdecay}=%g, t_{50}=%.1f\n', trial, pick.Rn0,pick.kelong,pick.kcl,pick.kndecay,pick.kproc,pick.kcdecay,pick.ktrans,pick.kpdecay,t50(trial));
  picks{trial}=pick;
end
figure;clf;
for i=1:ntrials
  plot(picks{i}.kpdecay,t50(i),'.');
  hold on;
end
xlabel('k_{pdecay}');
ylabel('t_{50%}');

% Parameter fitting
% Need to find realistic constants to allow for reasonable RNA and protein half lives, expression profiles at steady state (comparing switches of different kcl, time to reach steady state, etc.
% We know that sTRSV has ~2% expression. So can assume that only have ~2% Rc left by the time you start making proteins

