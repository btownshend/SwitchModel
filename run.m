x0=[0;0;0];   % Initial conditions
timerange=[0,4]*60*60;   % Time (in seconds)

% Parameters of model
% Each is defined with a 3-element vector consisting of lowest possible, best estimate, highest possible
rng=[];
% Rate of nuclear RNA production in molecules/second (dominated by initiation)
% Since system is linear, this just scales everything else, so not much point in letting it vary
rng.Rn0 = [25,25,25];
% rate of elongation in nucleotides/sec div by length of 3’UTR
% Not relevant to model since Rn0 controls production rate
rng.kelong = [25,25,25]*(1/2000);      
% rate of rz cleavage in /sec 
rng.kcl=[0.1,0.5,1.0]/60; 
% rate of nuclear RNA decay in /sec. Made up
rng.kndecay=[0.001,.01,.05];          
% rate of RNA processing, transport in /sec. Made up  
rng.kproc=[0.5e-6,1e-6,2e-6];
% rate of RNA cytoplasmic RNA decay in /sec. Made up 
rng.kcdecay=[0.5e-2,1e-2,2e-2];           
% rate of protein translation in /sec
% 8AA/sec, 233AA/protein, 0.64ribosomes/100nt, give .153 prot/sec
% slow down to allow for ribosome loading delay
rng.ktrans=(0.0001)*0.153*[0.5,1,2];
% rate of protein decay in /min. made up
rng.kpdecay=[0.005,0.01,0.02]/60;

% Form matrix for each of the parameter points
% Vectors are [Rn,Rc,P]
ntrials=50;
t50=[];
finalP=[];
odeoptions=odeset('MaxStep',60);   % maximum 1 minute steps
col='rgbcym';
variables={'RNA_{nuclear}','RNA_{cyto}','Protein'};

for trial=1:ntrials
  pick=struct('Rn0',rpick(rng.Rn0),'kelong',rpick(rng.kelong),'kcl',rpick(rng.kcl),'kndecay',rpick(rng.kndecay),...
              'kproc',rpick(rng.kproc),'kcdecay',rpick(rng.kcdecay),'ktrans',rpick(rng.ktrans),'kpdecay',rpick(rng.kpdecay));
  A=[-pick.kcl-pick.kndecay-pick.kproc, 0, 0
     pick.kproc, -pick.kcl-pick.kcdecay, 0
     0, pick.ktrans, -pick.kpdecay];
  c=[pick.Rn0,0,0]';
  [t,x]=ode45(@(t,x) A*x + c, timerange, x0,odeoptions);

  if trial<3
    % Plot results
    setfig(sprintf('trial %d',trial));
    clf;
    for i=1:size(x,2)
      semilogy(t/60,x(:,i),col(i));
      hold on;
      xlabel('Time (min)');
      ylabel('Concentration');
    end
    legend(variables);
    title(sprintf('Trial %d',trial));
    c=axis;
    c(3)=1e-10;
    axis(c);
  end
  
  % Find t50
  t50(trial)=t(find(x(:,3)>x(end,3)/2,1));
  finalP(trial)=x(end,3);
  fprintf('Trial %d: Rn_0=%g, k_{elong}=%g, k_{cl}=%g, k_{ndecay}=%g, k_{proc}=%g, k_{cdecay}=%g, k_{trans}=%g, k_{pdecay}=%g, t_{50}=%.1f,[P]_{final}=%g\n', trial, pick.Rn0,pick.kelong,pick.kcl,pick.kndecay,pick.kproc,pick.kcdecay,pick.ktrans,pick.kpdecay,t50(trial),finalP(trial));
  picks{trial}=pick;
end
fn=fieldnames(rng);
for j=1:length(fn)
  field=fn{j};
  if rng.(field)(3)-rng.(field)(1)==0
    % Not modified
    continue;
  end
  setfig(field);
  clf;
  subplot(211);
  for i=1:ntrials
    plot(picks{i}.(field),t50(i)/60,'.');
    hold on;
  end
  xlabel(field);
  ylabel('t_{50%} (min)');
  subplot(212);
  for i=1:ntrials
    semilogy(picks{i}.(field),finalP(i),'.');
    hold on;
  end
  xlabel(field);
  ylabel(sprintf('[Protein] at t=%d',timerange(2)));
end
% Parameter fitting
% Need to find realistic constants to allow for reasonable RNA and protein half lives, expression profiles at steady state (comparing switches of different kcl, time to reach steady state, etc.
% We know that sTRSV has ~2% expression. So can assume that only have ~2% Rc left by the time you start making proteins

