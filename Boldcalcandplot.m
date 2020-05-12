%%% Pearson correlation coefficient computed from BOLD simulated signals

Friston BALLOON MODEL
T = nn*dtt; % Total time in seconds

B = BOLD(T,neuro_act(1:nn,1)'); % B=BOLD activity, bf=Foutrier transform, f=frequency range)
BOLD_act = zeros(length(B),Nnew);
BOLD_act(:,1) = B;  

for nnew=2:Nnew

    B = BOLD(T,neuro_act(1:nn,nnew));
    BOLD_act(:,nnew) = B;
end


% Downsampling and reordering removing the first 500ms
bds=BOLD_act(500:ds:end,:);

% Global regression
%[b]=Global_regression(bds');

% clear bds BOLD_filt

% BOLD correlation matrix = Simulated Functional Connectivity
Cb  = corrcoef(bds);
cb  = atanh(Cb(Isubdiag)); % Vector containing all the FC values below the diagonal

imagesc(Cb);colormap(darkb2r(-0.5,1));colorbar();set(gca,'Fontsize', 20);