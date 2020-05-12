function [fc, maxfrN, curr2, neuro_act, wEI]=DMF_lesioned(scPath,scName,fic,ffi,simTime,dt,G,noiseAmp,lesionAreas,ficWeights,setFICWeights,Iext)
%%
%scPath          : structural connectivity(SC) matrix path
%fic             : true or false
%ffi             : true or false
%simTime         : simulation time
%dt              : Step size of Euler's method
%G               : Global coupling strength
%noiseAmp        : noise amplitude
%lesionAreas     : Index of the areas which are lesioned in the SC
%ficWeights      : local inhibitory synaptic coupling initialization values
%setFICWeights   : [True/False] if set true ficWeights will be used for
%                  initialization
%%

%Load Structural Connectivity
SC = h5read(scPath,scName);
nAreas = size(SC,1); % number of brain areas
SC(1:nAreas+1:nAreas*nAreas) = 0;
ds   = 2000;    % BOLD downsampling rate
dtt = 1e-3;
nonLesionAreas = setdiff(1:nAreas,lesionAreas);

% loading DMF model cortical pool parameters
tmax = simTime; % 12 mins(total time to generate synaptic activity)
% tspan = 0:dt:tmax; 
tmax1 = 10000;    % 10 secs(simulation time for FIC adjustment)
% tspan1 = 0:dt:tmax1;
taon = 100; % time constant for NMDA synapses
taog = 10;  % time constant for GABA synapses
gamma = 0.641;
sigma = noiseAmp;
JN = 0.15; % synaptic coupling for NMDA synapses
if(setFICWeights)
    wEI = ficWeights;
elseif(fic)
    wEI = 0.9*ones(nAreas,1);
else
    wEI = ones(nAreas,1);
end
I0 = 0.382;  % external input current  
wextE = 1.0; % scaling factor for I0 for excitatory populations
wextI = 0.7; % scaling factor for I0 for excitatory populations
wEE = 1.4; % wEE or w+ synpatics weight for intra area excitatory-excitatory
curr1 = zeros(tmax1,nAreas);
curr2 = zeros(simTime,nAreas);
neuro_act = zeros(tmax,nAreas);
frN = zeros(tmax,nAreas);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Balanced feedback Inhibiton control commented
% Iext = zeros(nAreas,1);    % external current due to stimulus
% Iext([10 44]) = 0.02; % give external stimulus to right and left LOCC
delta = 0.02*ones(nAreas,1);
flag = 0;
if(fic)
    disp('FIC started')
    while(flag < length(nonLesionAreas))
     sn = 0.001*ones(nAreas,1); % excitatory(NMDA) synaptic gating variable
     sg = 0.001*ones(nAreas,1); % inhibitory(GABA) synaptic gating variable
     sn(lesionAreas) = 0;
     sg(lesionAreas) = 0;
     nIters = length(1:dt:tmax1) - 1;
     j = 0;
     nn = 0;
     for i = 1:nIters
      xn = I0*wextE + wEE*JN*sn + G*JN*SC*sn - wEI.*sg + Iext; % input current to excitatory(NMDA) population
      if(ffi)
        xg = I0*wextI + JN*sn - sg + G*JN*SC*sn; % input current to inhibitory(GABA) population with long range feed-forward inhibition(FFI)
      else
        xg = I0*wextI + JN*sn - sg;  % input current to inhibitory(GABA) population without FFI
      end
      xn(lesionAreas) = 0;
      xg(lesionAreas) = 0;
      rn = phie(xn); % excitatory(NMDA) population firing rate
      rg = phii(xg); % inhibitory(GABA) population firing rate
      rn(lesionAreas) = 0;
      rg(lesionAreas) = 0;
      sn = sn + dt*(-sn/taon+(gamma/1000.0)*(1-sn).*rn) + sqrt(dt)*sigma*randn(nAreas,1);
      sn(lesionAreas) = 0;
      sn(sn>1) = 1;
      sn(sn<0) = 0;
      sg = sg + dt*(-sg/taog+rg./1000.0) + sqrt(dt)*sigma*randn(nAreas,1);
      sg(lesionAreas) = 0;
      sg(sg>1) = 1;
      sg(sg<0) = 0;
      j = j + 1;
      if(j == 1/dt)
        nn = nn + 1;
        j = 0;
        curr1(nn,:) = xn' - 125.0/310.0;
      end
     end
     currm = mean(curr1(1000:end,:),1);
     flag = 0;
     for i = 1:length(nonLesionAreas)
      n = nonLesionAreas(i);
      if abs(currm(n)+0.026) > 0.005
       if currm(n) < -0.026
        wEI(n) = wEI(n) - delta(n);
        delta(n) = delta(n) - 0.001;
        if delta(n) < 0.001
           delta(n) = 0.001;
        end
       else
        wEI(n) = wEI(n) + delta(n);
%         delta(n) = delta(n) - 0.001;
%         if delta(n) < 0.001;
%            delta(n) = 0.001;
%         end
       end
      else
       flag = flag+1;
      end
     end
    disp(strcat('No.of areas satisfied FIC: ',int2str(flag)))
    end
    disp('FIC completed')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulating synaptic activity

disp('simulating global brain dynamics...');
sn=0.001*ones(nAreas,1); % excitatory(NMDA) synaptic gating variable
sg=0.001*ones(nAreas,1); % inhibitory(GABA) synaptic gating variable
sn(lesionAreas) = 0;
sg(lesionAreas) = 0;
j = 0;
nn = 0;
nIters = length(1:dt:tmax) - 1;
for i=1:nIters
  %noise = randn(nAreas,1);
  xn = I0*wextE + wEE*JN*sn + G*JN*SC*sn - wEI.*sg + Iext; % input current to excitatory(NMDA) population
  if(ffi)
    xg = I0*wextI + JN*sn - sg + G*JN*SC*sn; % input current to inhibitory(GABA) population with long range feed-forward inhibition(FFI)
  else
    xg = I0*wextI + JN*sn - sg;  % input current to inhibitory(GABA) population without FFI
  end
  xn(lesionAreas) = 0;
  xg(lesionAreas) = 0;
  rn = phie(xn); % excitatory(NMDA) population firing rate
  rg = phii(xg); % inhibitory(GABA) population firing rate
  rn(lesionAreas) = 0;
  rg(lesionAreas) = 0;
  sn = sn + dt*(-sn/taon+(gamma/1000.0)*(1-sn).*rn) + sqrt(dt)*sigma*randn(nAreas,1);
  sn(lesionAreas) = 0;
  sn(sn>1) = 1;
  sn(sn<0) = 0;
  sg = sg + dt*(-sg/taog+rg./1000.0) + sqrt(dt)*sigma*randn(nAreas,1);
  sg(lesionAreas) = 0;
  sg(sg>1) = 1;
  sg(sg<0) = 0;
  j = j + 1;
  if(j == 1/dt)
    nn = nn + 1;
    j = 0;
    neuro_act(nn,:) = sn';
    frN(nn,:) = rn';
    curr2(nn,:) = xn' - 125.0/310.0;
  end
end
maxfrN = max(frN);
disp('simulation completed and synaptic activity computed');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOLD empirical calculation using Hemodynamic model : Buxton et al, BALLOON MODEL
disp('computing BOLD signals using synaptic activity...')
%T = tmax*dt; % Total time in seconds
T = nn*dtt;
boldFlag = 1;
for j = 1:length(nonLesionAreas)
    i = nonLesionAreas(j);
    if(boldFlag)
%         B = Bold_ideal(T,neuro_act(1:nn,i)',dtt);
        B = BOLD(T,neuro_act(1:nn,i)',dtt);
        BOLD_act = zeros(length(B),nAreas);
%         BOLD_act(:,i) = atanh(B);
        BOLD_act(:,i) = B;
        boldFlag = 0;
    else
%         B = Bold_ideal(T,neuro_act(1:nn,i)',dtt);
        B = BOLD(T,neuro_act(1:nn,i)',dtt);
%         BOLD_act(:,i) = atanh(B);
        BOLD_act(:,i) = B;
    end
end
disp('BOLD signals computed')

% Downsampling and reordering removing the first 500ms
bds=BOLD_act(500:ds:end,:);
% % BOLD correlation matrix = Simulated Functional Connectivity
fc  = corrcoef(bds);
fc(lesionAreas,:) = 0;
fc(:,lesionAreas) = 0;
end
