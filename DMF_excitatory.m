function [fc, frN]=DMF_excitatory(scPath,simTime,dt,G,noiseAmp)

%scPath     : structural connectivity matrix path
%simTime    : simulation time
%dt         : Step size of Euler's method
%G          : Global coupling strength
%noiseAmp   : noise amplitude

%Load Structural Connectivity
SC = h5read(scPath,'/C');
nAreas = size(SC,1); % number of brain areas
SC(1:nAreas+1:nAreas*nAreas) = 0;
ds   = 2000;    % BOLD downsampling rate
dtt = 1e-3;


% loading DMF model cortical pool parameters
tmax = simTime; % 12 mins(total time to generate synaptic activity)
% tspan = 0:dt:tmax; 
taon = 100; % time constant for NMDA synapses
gamma = 0.641;
sigma = noiseAmp;
JN = 0.2609; % synaptic coupling for NMDA synapses
I0 = 0.3;  % external input current  
wEE = 0.9; % wEE or w+ synpatics weight for intra area excitatory-excitatory
neuro_act = zeros(tmax,nAreas);
frN = zeros(tmax,nAreas);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulating synaptic activity

disp('simulating global brain dynamics...');
sn=0.001*ones(nAreas,1); % excitatory(NMDA) synaptic gating variable
j = 0;
nn = 0;
nIters = length(1:dt:tmax) - 1;
for i=1:nIters
  xn = wEE*JN*sn + G*JN*SC*sn + I0; % input current to excitatory(NMDA) population
  rn = phie_excitatory(xn); % excitatory(NMDA) population firing rate
  sn = sn + dt*(-sn/taon+(gamma/1000.0)*(1-sn).*rn) + sqrt(dt)*sigma*randn(nAreas,1);
  sn(sn>1) = 1;
  sn(sn<0) = 0;
  j = j + 1;
  if(j == 1/dt)
    nn = nn + 1;
    j = 0;
    neuro_act(nn,:) = sn';
    frN(nn,:) = rn';
  end  
end
disp('simulation completed and synaptic activity computed');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOLD empirical calculation using Hemodynamic model : Buxton et al, BALLOON MODEL
disp('computing BOLD signals using synaptic activity...')
%T = tmax*dt; % Total time in seconds
T = nn*dtt;
% B = Generate_Bold(T,neuro_act(1:nn,1)',dtt); 
B = Bold_ideal(T,neuro_act(1:nn,1)',dtt); 
BOLD_act = zeros(length(B),nAreas);
BOLD_act(:,1) = B;

for i=2:nAreas
%    B = Generate_Bold(T,neuro_act(1:nn,i)',dtt);
   B = Bold_ideal(T,neuro_act(1:nn,i)',dtt);
   BOLD_act(:,i) = B;
end
disp('BOLD signals computed')

% Downsampling and reordering removing the first 500ms
bds=BOLD_act(500:ds:end,:);
% % BOLD correlation matrix = Simulated Functional Connectivity
fc  = corrcoef(bds);
end
