%% DMF model simulation based on Gustavo's code
%% Model needs a SC matrix to simulate neuronal resting activity of about 15 seconds.
%% Generated neuronal signal is passed to a HEMODYNAMIC model to simulate BOLD signals.
%% Script also include recently proposed model by Deco et al. (2013,2014) with Feedback Inhibition control
%% to maintain balance ex-inh over entire cortex composed of multiple cortical regions.
fcs = cell(length(0.1:0.1:4));
iters = 0;
for we = 0.1:0.1:4
iters = iters + 1
close all;
clearvars -except fcs we iters;

%Load Structural Connectivity

connectionmatrixpath = '../../data/connectivity_66/';
weights = load([connectionmatrixpath,'Human_66.mat']);
%weights_FC = load([connectionmatrixpath, 'CN_20120927_fMRI_new']);
FC_emp = weights.FC_emp(1:66,1:66);
%weights = load([connectionmatrixpath, 'Haggman_SC.mat']);
%connectionmatrixpath = '/Users/dipanjanroy/Documents/DRCodes2013/Node_model_RS/BrainNetworkModels_3.0.2/';
%weights = load([connectionmatrixpath, 'SC_Haggman.mat']);
%weights = load([connectionmatrixpath, 'Human_66.mat']);
% connectionmatrixpath='/Users/dipanjanroy/Documents/DRCodes2013/Node_model_RS/BOLD_BALANCE_MODEL_REST/';
% weights = load([connectionmatrixpath, '']);
%weights = load([connectionmatrixpath, 'CN_20120927_SC.mat']);

%C=weights.SC_cap_agg_bwflav2_norm(1:10,1:10);
%C=C.*10^-3;
%load 'FC_Haggman_new.mat' FC_emp %load empirical correlation matrix
%load SC_Haggman.mat C % load structural connectivity matrix obtained from DSI/DTI
C=weights.C(1:66,1:66);
%FC_new=FC_emp;
Nnew = size(C,1);
Isubdiag = find(tril(ones(Nnew),-1));
fc = atanh(FC_emp(Isubdiag));
sc = C(Isubdiag);

%weights = load([connectionmatrixpath, 'XB_20120831_SC.mat']);
%SC_cap_agg_counts=(weights.SC_cap_agg_bwflav2); %Normalization factor across all-SC 

%C=Cnew(Order,Order);
%C = weights.C(1:66,1:66);
%C=SC_cap_agg_counts(1:68,1:68);

%Isubdiag = find(tril(ones(68),-1)); % Indexes of all the values below the diagonal.

%fc       = atanh(FC_emp(Isubdiag)); % Vector containing all the FC values below the
%Nnew = size(C,1);
dtt   = 1e-3;   % Sampling rate of simulated neuronal activity (seconds)
ds   = 100;    % BOLD downsampling rate
%%

% loading DMF model cortical pool parameters

%%%%%%%%%

dt=1;
tmax=8*60*1000;
tspan=0:dt:tmax;
tmax1=4000;
tspan1=0:dt:tmax1;

taon=100;
% taog=10;
gamma=0.641;
sigma=0.001;
JN=0.2609;
% J=ones(Nnew,1);
I0=0.3;
Jexte=1.0;
% Jexti=0.7;
w=0.9;
%Iext=0.02;
curr=zeros(tmax1,Nnew);
curr_new=zeros(tmax,Nnew);
neuro_act=zeros(tmax,Nnew);
frN = zeros(tmax,Nnew);
% neuro_act2=zeros(tmax1,Nnew);
ix=1;
% we=1.4;
%G=0.1:0.025:3.0; %relevant Global scaling parameter to simulate global cortical model
%max_firing_rate = zeros(1,66);%firing rate of all cortical pools Intialized
%for we=G;
%we
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Balanced feedback Inhibiton control commented
%k=1;
% sn=0.001*ones(Nnew,1);
% sg=0.001*ones(Nnew,1);
% %Iext=zeros(Nnew,1);
% delta=0.02*ones(Nnew,1);
% 
% for k=1:5000
%  sn=0.001*ones(Nnew,1);
%  sg=0.001*ones(Nnew,1);
%  nn=1;
%  j=0;
% for i=2:1:length(tspan1)
%   xn=I0*Jexte+w*JN*sn+we*JN*C*sn-J.*sg;
%  % xg=I0*Jexti+JN*sn-sg;
%   rn=phie(xn);
%  % rg=phii(xg);
%   sn=sn+dt*(-sn/taon+(1-sn)*gamma.*rn./1000.)+sqrt(dt)*sigma*randn(Nnew,1);
%   sn(sn>1) = 1;  
%   sn(sn<0) = 0;      
%  % sg=sg+dt*(-sg/taog+rg./1000.)+sqrt(dt)*sigma*randn(Nnew,1);
%  % sg(sg>1) = 1;        
%  % sg(sg<0) = 0;
%   j=j+1;
%   if j==10
%    curr(nn,:)=xn'-125/310;
%    neuro_act2(nn,:)=rn';
%    nn=nn+1;
%    j=0;
%   end
%  end
% 
%  currm=mean(curr(1000:end,:),1);
%  flag=0;
%  for n=1:1:Nnew
%   if abs(currm(n)+0.026)>0.005 
%    if currm(n)<-0.026 
%     J(n)=J(n)-delta(n);
%     delta(n)=delta(n)-0.001;
%     if delta(n)<0.001;
%        delta(n)=0.001;
%     end
%    else 
%     J(n)=J(n)+delta(n);
%    end
%   else
%    flag=flag+1;
%   end
%  end
% 
%  if flag==Nnew
%   break;
%  end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 sn=0.001*ones(Nnew,1);
% sg=0.001*ones(Nnew,1);
%Iext = [0.0;0.0;0.0;0.01;0.0;0.0;0.0;0.0;0.0;0.0];
 
 nn=1;
 j=0;
    
 for i=2:1:length(tspan)
     
  xn=I0*Jexte+w*JN*sn+we*JN*C*sn;
%  xg=I0*Jexti+JN*sn-sg;
  rn=phie_single(xn);
%  rg=phii(xg);
  sn=sn+dt*(-sn/taon+(1-sn)*gamma.*rn./1000.)+sqrt(dt)*sigma*randn(Nnew,1);
  sn(sn>1) = 1; 
  sn(sn<0) = 0;             
%  sg=sg+dt*(-sg/taog+rg./1000.)+sqrt(dt)*sigma*randn(Nnew,1);
%  sg(sg>1) = 1;        
%  sg(sg<0) = 0;
  j=j+1;
  if j == (1/dt)
%   curr_new(nn,:)=xn'-125/310;   
   frN(nn,:) = rn';
   neuro_act(nn,:)=sn';
   nn=nn+1;
   j=0;
  end
 end
% computes max excitatory firing rate from Each pool
% to compute bifurcation diagram as a function of Global
% scaling parameter 'G'

%max_firing_rate(ix) = max(rn);
%ix=ix+1;
%end
nn=nn-1;
% 
% %%% BOLD empirical calculation using Hemodynamic model
% 
% %Friston BALLOON MODEL
T = nn*dtt; % Total time in seconds

B = BOLD(T,neuro_act(1:nn,1)'); % B=BOLD activity, bf=Foutrier transform, f=frequency range)
BOLD_act = zeros(length(B),Nnew);
BOLD_act(:,1) = B;  
% 
for nnew=2:Nnew

   B = BOLD(T,neuro_act(1:nn,nnew));
   BOLD_act(:,nnew) = B;
end
% 
% 
% % Downsampling and reordering removing the first 500ms
bds=BOLD_act(500:ds:end,:);
% 
% % Global regression
% %[b]=Global_regression(bds');
% 
% % clear bds BOLD_filt
% 
% % BOLD correlation matrix = Simulated Functional Connectivity
fc_sim  = corrcoef(bds);
fcs{iters} = fc_sim;
% cb  = atanh(Cb(Isubdiag)); % Vector containing all the FC values below the diagonal 
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