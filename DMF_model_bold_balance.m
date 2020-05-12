function [Cb_FI, neuro_act, frN]=DMF_model_bold_balance(scPath,fic,simTime,dt)

% DMF model simulation based on Gustavo's code
% Model needs a SC matrix to simulate neuronal resting activity of about 12 minutes.
% Generated neuronal signal is passed to a HEMODYNAMIC model to simulate BOLD signals.
% Script also include recently proposed model by Deco et al. (2013,2014) with Feedback Inhibition control
% to maintain balance ex-inh over entire cortex composed of multiple cortical regions.

if(nargin<4)
   disp('Using default values for FIC:false,simTime:12 mins,dt:0.001')
   fic = false;
   simTime = 12*60*1000; % simulation time in milli seconds
   dt = 0.1; % Sampling rate of simulated neuronal activity (seconds)
end

%Load Structural Connectivity
SC = h5read(scPath,'/C');
Nnew = size(SC,1);
Isubdiag = find(tril(ones(Nnew),-1)); % indices of lower diagonal elements of SC matrix
scSubDiag = SC(Isubdiag); % lower diagonal elements of SC matrix
% dtt   = 1e-3;   % Sampling rate of simulated neuronal activity (seconds)
ds   = 100;    % BOLD downsampling rate


% loading DMF model cortical pool parameters
tmax = simTime; % 12 mins(total time to generate synaptic activity)
% tspan = 0:dt:tmax; 
tmax1 = 1000;    % 10 secs(simulation time for FIC adjustment)
% tspan1 = 0:dt:tmax1;
dtt = 0.001;
taon = 100; % time constant for NMDA synapses
taog = 10;  % time constant for GABA synapses
gamma = 0.641;
sigma = 0.01;
JN = 0.15; % synaptic coupling for NMDA synapses
wEI = 0.9*ones(Nnew,1);
I0 = 0.382;  % external input current  
wextE = 1.0; % scaling factor for I0 for excitatory populations
wextI = 0.7; % scaling factor for I0 for excitatory populations
wEE = 1.4; % wEE or w+ synpatics weight for intra area excitatory-excitatory
curr = zeros(tmax1,Nnew);
% curr_new = zeros(tmax,Nnew);
neuro_act = zeros(tmax,Nnew);
frN = zeros(tmax,Nnew);
%frG = zeros(tmax,Nnew);
% neuro_act = zeros(tmax/10,Nnew);
% neuro_act2 = zeros(tmax1,Nnew);
G = 1.5;
%G=0.1:0.025:3.0; %relevant Global scaling parameter to simulate global cortical model
%max_firing_rate = zeros(1,66);%firing rate of all cortical pools Intialized
%for we=G;
%we
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Balanced feedback Inhibiton control commented
Iext = zeros(Nnew,1);    % external current due to stimulus
delta = 0.02*ones(Nnew,1);
flag = 0;
if(fic)
    disp('FIC started')
    while(flag < Nnew)
     sn = 0.001*ones(Nnew,1); % excitatory(NMDA) synaptic gating variable
     sg = 0.001*ones(Nnew,1); % inhibitory(GABA) synaptic gating variable
     j = 0;
     nn = 0;
     for i = 0:dt:tmax1
      noiseN = randn(Nnew,1);
      noiseG = randn(Nnew,1);
      xn = I0*wextE + wEE*JN*sn + G*JN*SC*sn - wEI.*sg + Iext; % input current to excitatory(NMDA) population
      xg = I0*wextI + JN*sn - sg; % input current to inhibitory(GABA) population
      rn = phie(xn); % excitatory(NMDA) population firing rate
      rg = phii(xg); % inhibitory(GABA) population firing rate
      sn = sn + dt*(-sn/taon+(1-sn)*gamma.*rn./1000.) + sqrt(dt)*sigma*noiseN;
      sn(sn>1) = 1;
      sn(sn<0) = 0;
      sg = sg + dt*(-sg/taog+rg./1000.) + sqrt(dt)*sigma*noiseG;
      sg(sg>1) = 1;
      sg(sg<0) = 0;
      j = j + 1;
      if(j == 1/dt)
        nn = nn + 1;
        curr(nn,:) = xn' - 125/310;
        j = 0;
      end
     end
     currm = mean(curr(1000:end,:),1);
     flag = 0;
     for n = 1:1:Nnew
      %if (n==2&3) && abs(currm(n)+0.026)>0.005 
      if abs(currm(n)+0.026) > 0.005
    %    disp(currm(n))
       if currm(n) < -0.026
        wEI(n) = wEI(n) - delta(n);
        delta(n) = delta(n) - 0.001;
        if delta(n) < 0.001;
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
%      disp(strcat('No.of areas satisfied FIC: ',int2str(flag)))
    end
    disp('FIC completed')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulating synaptic activity

disp('simulating global brain dynamics...');
sn=0.001*ones(Nnew,1); % excitatory(NMDA) synaptic gating variable
sg=0.001*ones(Nnew,1); % inhibitory(GABA) synaptic gating variable
j = 0;
nn = 0;
for i=0:dt:tmax
  noiseN = randn(Nnew,1);
  noiseG = randn(Nnew,1);
  xn = I0*wextE + wEE*JN*sn + G*JN*SC*sn - wEI.*sg + Iext; % input current to excitatory(NMDA) population
  xg = I0*wextI + JN*sn - sg; % input current to inhibitory(GABA) population
  rn = phie(xn); % excitatory(NMDA) population firing rate
  rg = phii(xg); % inhibitory(GABA) population firing rate
  sn = sn + dt*(-sn/taon+(1-sn)*gamma.*rn./1000.) + sqrt(dt)*sigma*noiseN;
  sn(sn>1) = 1;
  sn(sn<0) = 0;
  sg = sg + dt*(-sg/taog+rg./1000.) + sqrt(dt)*sigma*noiseG;
  sg(sg>1) = 1;
  sg(sg<0) = 0;
  j = j + 1;
  if(j == 1/dt)
      nn = nn + 1;
      neuro_act(nn,:) = sn';
      frN(nn,:) = rn';
      j = 0;
  end
end
disp('simulation completed and synaptic activity computed');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% computes max excitatory firing rate from Each pool
% to compute bifurcation diagram as a function of Global
% scaling parameter 'G'

%max_firing_rate(ix) = max(rn);
%ix=ix+1;
%end

% 
% %%% BOLD empirical calculation using Hemodynamic model
% 
% %Friston BALLOON MODEL
disp('computing BOLD signals using synaptic activity...')
% T = tmax*dt; % Total time in seconds
T = nn*dtt;

B = BOLD(T,neuro_act(1:nn,1)',dtt); % B=BOLD activity, bf=Foutrier transform, f=frequency range)
% disp(size(B))
% disp(length(B))
BOLD_act = zeros(length(B),Nnew);
BOLD_act(:,1) = B;  
% 
for nnew=2:Nnew
   B = BOLD(T,neuro_act(1:nn,nnew),dtt);
   BOLD_act(:,nnew) = B;
end
disp('BOLD sginals computed')
% 
% 
% % Downsampling and reordering removing the first 500ms
bds_FI=BOLD_act(500:ds:end,:);
% 
% % Global regression
% %[b]=Global_regression(bds');
% 
% % clear bds BOLD_filt
% 
% % BOLD correlation matrix = Simulated Functional Connectivity
Cb_FI  = corrcoef(bds_FI);

% cb  = atanh(Cb_FI(Isubdiag)); % Vector containing all the FC values below the diagonal 
%imagesc(Cb_FI);colormap(darkb2r(-0.5,1));colorbar();set(gca,'Fontsize', 20);
% % 
% % %%%%%%%
% % 
%Coef_restFC    = corrcoef(Cb,FC_emp);
% Coef_rest    = corrcoef(cb,sc)
% % % % %fittcorr(ix)=Coef(2)
%fittcorr_rest(ix)=Coef_restFC(2)
%ix=ix+1
%end
%figure
%plot(wee,fittcorr);
% plot(wee,fittcorr2);
% 
% 
% load Human_66.mat C Order
% figure;
% subplot(2,2,1)
% imagesc(C);
% colorbar;
% subplot(2,2,2)
% imagesc(Cnew);
% colorbar;
% subplot(2,1,2)
% plot(wee,fittcorr);

end