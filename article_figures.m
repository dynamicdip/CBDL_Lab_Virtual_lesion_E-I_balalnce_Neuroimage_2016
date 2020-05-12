%%%%%% Empirical Healthy SC %%%%%%
% sc_hc = h5read('../../data/connectivity_68/Human_68.hdf5','/C');
% load('SC_color_map.mat');
% imagesc(sc_hc)
% colormap(cmp)
% c = colorbar('YTick',[0 0.2 0.4 0.6 0.8]);
% axpos = get(gca,'Position');
% cpos = get(c,'Position');
% cpos(3) = 0.5*cpos(3);
% set(c,'Position',cpos);
% set(gca,'Position',axpos);
% set(gca,'XTick',[1 10 20 30 40 50 60 68],'YTick',[1 10 20 30 40 50 60 68])
% xlabel('Area ID','FontWeight','bold')
% ylabel('Area ID','FontWeight','bold')
% savePath = '../../reports/figures/sc_hc';
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Virtually Lesioned SC %%%%%%
% load('LesionColorMap.mat')
% c = plot_lesioned_sc('../../data/connectivity_68/lesioned/sc_24_PCUN_nn1.hdf5','LesionColorMapnn1.mat','../../data/connectivity_68/fs_region_abbrev_sort.txt')
% savePath = '../../reports/figures/sc_lpcun';
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% optimal G %%%%%%
% load('../../results/68/optimalG_68_fic_dt1_sigma0.001_8mins_EI.mat','fcCorrs')
% area(0.025:0.025:1,fcCorrs,'FaceColor',[0.5 0.5 0.5])
% hold on;
% plot(fcCorrs,'k','LineWidth',2)
% plot([0.6],[0.6],'gx','MarkerSize',10)
% % text(0.62,0.6,'\leftarrowr=0.6,G=0.6','FontSize',20,'Color',[.0 .0 .0])
% xlabel('Global Coupling Strength(G)','FontWeight','bold')
% ylabel('FC Correlation Fit','FontWeight','bold')
% savePath = '../../reports/figures/optimalG';
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% optimal G Firing rates%%%%%%
% load('../../results/68/optimalG_68_fic_dt1_sigma0.001_8mins_EI.mat','maxfrNs')
% barh(maxfrNs{24},'FaceColor',[0.5 0.5 0.5])
% box off;
% set(gca,'XTick',[1 5 10],'YTick',[1 10 20 30 40 50 60 68])
% xlabel('Firing Rate','FontWeight','bold')
% ylabel('Area ID','FontWeight','bold')
% savePath = '../../reports/figures/optimalG_fr';
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Empirical FC %%%%%%%
% fc_hc_emp = h5read('../../data/connectivity_68/Human_68.hdf5','/CC');
% figure();
% imagesc(fc_hc_emp)
% colormap('hot')
% caxis([-0.117 1.0])
% c = colorbar('YTick',[-0.1 0.2 0.5 0.8 1],'FontSize',15);
% axpos = get(gca,'Position');
% cpos = get(c,'Position');
% cpos(3) = 0.5*cpos(3);
% set(c,'Position',cpos);
% set(gca,'Position',axpos);
% set(gca,'XTick',[1 10 20 30 40 50 60 68],'YTick',[1 10 20 30 40 50 60 68],'FontSize',15)
% xlabel('Area ID','FontWeight','bold')
% ylabel('Area ID','FontWeight','bold')
% savePath = '../../reports/figures/fc_hc_emp';
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Simulated FC healthy %%%%%%%
% fc_hc_sim = load('sim_hc.mat','fc');
% figure();
% imagesc(fc_hc_sim.fc)
% colormap('hot')
% c = colorbar('YTick',[-0.1 0.2 0.5 0.8 1],'FontSize',15);
% axpos = get(gca,'Position');
% cpos = get(c,'Position');
% cpos(3) = 0.5*cpos(3);
% set(c,'Position',cpos);
% set(gca,'Position',axpos);
% set(gca,'XTick',[1 10 20 30 40 50 60 68],'YTick',[1 10 20 30 40 50 60 68],'FontSize',15)
% xlabel('Area ID','FontWeight','bold')
% ylabel('Area ID','FontWeight','bold')
% savePath = '../../reports/figures/fc_hc_sim';
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Loss of balance upon lesioning PCUN Prior%%%%%
% load('sim_hc.mat','wEI')
% simTime = 10*1000;
% [fc1,frN1,cur1] = DMF_lesioned('../../data/connectivity_68/lesioned/sc_24_PCUN_nn1.hdf5',false,false,simTime,1,0.6,0.001,[24],wEI,true);
% curm1 = mean(cur1(1000:end,:),1);
% abscurm1 = abs(curm1 + 0.026);
% abscurm1(24) = 0;
% barh(abscurm1);
% hold on;
% xlim([0 0.025])
% ylim([1 69])
% plot([0.005 0.005],[0 70],'Color','r')
% ylabel('Area ID','FontWeight','bold')
% xlabel('E-I balance test condition*','FontWeight','bold')
% t = [0:5:60];
% set(gca,'YTick',[1 t(2:end) 64 68],'XTick',[0 0.01 0.02],'FontSize',13)
% box off;
% savePath = '../../reports/figures/balance_loss1_1';
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Loss of balance upon lesioning PCUN Post %%%%%
% load('sim_hc.mat','wEI')
% simTime = 10*1000;
% [fc2,frN2,cur2] = DMF_lesioned('../../data/connectivity_68/lesioned/sc_24_PCUN_nn1.hdf5',true,false,simTime,1,0.6,0.001,[24],wEI,true);
% curm2 = mean(cur2(1000:end,:),1);
% abscurm2 = abs(curm2 + 0.026);
% abscurm2(24) = 0;
% barh(abscurm2);
% hold on;
% xlim([0 0.025])
% ylim([1 69])
% plot([0.005 0.005],[0 70],'Color','r')
% ylabel('Area ID','FontWeight','bold')
% xlabel('E-I balance test condition*','FontWeight','bold')
% t = [0:5:60];
% set(gca,'YTick',[1 t(2:end) 64 68],'XTick',[0 0.01 0.02],'FontSize',13)
% box off;
% savePath = '../../reports/figures/balance_loss1_2';
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% No.of areas that lost balance for each lesion center%%%%%%%%%%%
% scPaths = {'../../data/connectivity_68/lesioned/sc_1_BSTS_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_2_CAC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_3_CMF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_4_CUN_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_5_ENT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_6_FUS_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_7_IP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_8_IT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_9_ISTH_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_10_LOCC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_11_LOF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_12_LING_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_13_MOF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_14_MT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_15_PARH_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_16_PARC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_17_POPE_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_18_PORB_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_19_PTRI_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_20_PCAL_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_21_PCNT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_22_PC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_23_PREC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_24_PCUN_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_25_RAC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_26_RMF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_27_SF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_28_SP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_29_ST_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_30_SMAR_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_31_FP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_32_TP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_33_TT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_34_INS_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_35_BSTS_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_36_CAC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_37_CMF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_38_CUN_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_39_ENT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_40_FUS_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_41_IP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_42_IT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_43_ISTH_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_44_LOCC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_45_LOF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_46_LING_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_47_MOF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_48_MT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_49_PARH_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_50_PARC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_51_POPE_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_52_PORB_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_53_PTRI_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_54_PCAL_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_55_PCNT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_56_PC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_57_PREC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_58_PCUN_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_59_RAC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_60_RMF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_61_SF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_62_SP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_63_ST_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_64_SMAR_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_65_FP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_66_TP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_67_TT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_68_INS_nn1.hdf5'};
% load('sim_hc.mat','wEI');
% simTime = 10*1000;
% balLoss = zeros(length(scPaths),1);
% for i = 1:length(scPaths)
%     nbrs = h5read(scPaths{i},'/Neighbours');
%     [fc1,frN1,cur1] = DMF_lesioned(scPaths{i},false,false,simTime,1,0.6,0.001,nbrs,wEI,true);
%     curm1 = mean(cur1(1000:end,:),1);
%     abscurm1 = abs(curm1 + 0.026);
%     abscurm1(nbrs) = 0;
%     balLoss(i) = length(find(abscurm1 > 0.005));
% end
% regionNames = strsplit(fileread('../../data/connectivity_68/fs_region_abbrev_sort.txt'),'\n');
% regionNames = regionNames(1:68);
% barh(balLoss);
% ylabel('Lesion Center','FontWeight','bold')
% xlabel('No. of Areas lost E-I balance','FontWeight','bold')
% box off;
% set(gca,'YTick',1:68,'YTickLabel',regionNames,'FontSize',5,'XTick',[0:5:60],'FontSize',10)
% set(gca,'YTick',[1 t(2:end) 64 68],'FontSize',13,'XTick',[0:10:60],'FontSize',13)
% savePath = '../../reports/figures/balance_loss_2';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% No.of areas that lost balance for each lesion center%%%%%%%%%%%
% scPaths = {'../../data/connectivity_68/lesioned/sc_1_BSTS_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_2_CAC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_3_CMF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_4_CUN_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_5_ENT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_6_FUS_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_7_IP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_8_IT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_9_ISTH_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_10_LOCC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_11_LOF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_12_LING_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_13_MOF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_14_MT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_15_PARH_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_16_PARC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_17_POPE_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_18_PORB_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_19_PTRI_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_20_PCAL_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_21_PCNT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_22_PC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_23_PREC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_24_PCUN_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_25_RAC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_26_RMF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_27_SF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_28_SP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_29_ST_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_30_SMAR_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_31_FP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_32_TP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_33_TT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_34_INS_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_35_BSTS_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_36_CAC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_37_CMF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_38_CUN_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_39_ENT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_40_FUS_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_41_IP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_42_IT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_43_ISTH_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_44_LOCC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_45_LOF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_46_LING_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_47_MOF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_48_MT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_49_PARH_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_50_PARC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_51_POPE_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_52_PORB_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_53_PTRI_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_54_PCAL_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_55_PCNT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_56_PC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_57_PREC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_58_PCUN_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_59_RAC_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_60_RMF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_61_SF_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_62_SP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_63_ST_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_64_SMAR_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_65_FP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_66_TP_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_67_TT_nn1.hdf5',...
% '../../data/connectivity_68/lesioned/sc_68_INS_nn1.hdf5'};
% % load('sim_hc.mat','wEI');
% % scEmp = h5read('../../data/connectivity_68/Human_68.hdf5','/C');
% % [commStrucSCEmp,modulSCEmp] = community_louvain(scEmp,1);
% % partCoefSCEmp = participation_coef(scEmp,commStrucSCEmp,1);
% % simTime = 10*1000;
% % for i = 1:length(scPaths)
% %     nbrs = h5read(scPaths{i},'/Neighbours');
% %     [fc1,frN1,cur1] = DMF_lesioned(scPaths{i},false,false,simTime,1,0.6,0.001,nbrs,wEI,true);
% %     curm1 = mean(cur1(1000:end,:),1);
% %     abscurm1 = abs(curm1 + 0.026);
% %     abscurm1(nbrs) = 0;
% %     balLoss(i) = length(find(abscurm1 > 0.005));
% % end
% [rho ,pVal] = corrcoef(partCoefSCEmp,balLoss);
% rho = rho(1,2);
% pVal = pVal(1,2);
% figure(2);
% plot(partCoefSCEmp,balLoss,'ko','MarkerFaceColor','k','MarkerSize',5);
% ylabel('No. of Areas lost E-I balance','FontWeight','bold');
% xlabel('Participation Coefficient','FontWeight','bold');
% hold on;
% box off;
% p=polyfit(partCoefSCEmp,balLoss',1);
% pred = polyval(p,partCoefSCEmp);
% plot(partCoefSCEmp,pred,'k');
% text(0.025,60,['r=',sprintf('%.2f',rho),', p<0.001'],'FontSize',15,'BackgroundColor',[0.3 0.7 0.7]);
% xlim([0 0.7])
% ylim([1 60])
% t = 0:10:60;
% set(gca,'XTick',0:0.2:0.7,'YTick',[1 t(2:end)],'FontSize',13)
% savePath = '../../reports/figures/balance_loss_3';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% FC Simulated lPCUN prior %%%%%%%%%%%%%%
% fc_lpcun_prior = h5read('../../results/68/fc_24_PCUN_G0.6_sigma0.001_ficprior_nn1_8mins.hdf5','/avgCC');
% imagesc(fc_lpcun_prior)
% colormap(hot)
% c = colorbar('YTick',[0 0.2 0.4 0.6 0.8 1],'FontSize',15);
% axpos = get(gca,'Position');
% cpos = get(c,'Position');
% cpos(3) = 0.5*cpos(3);
% set(c,'Position',cpos);
% set(gca,'Position',axpos);
% set(gca,'XTick',[1 10 20 30 40 50 60 68],'YTick',[1 10 20 30 40 50 60 68],'FontSize',15)
% xlabel('Area ID','FontWeight','bold')
% ylabel('Area ID','FontWeight','bold')
% savePath = '../../reports/figures/fc_lpcun_prior';
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% FC Simulated lPCUN post %%%%%%%%%%%%%%
% fc_lpcun_post = h5read('../../results/68/fc_24_PCUN_G0.6_sigma0.001_ficpost_nn1_8mins.hdf5','/avgCC');
% imagesc(fc_lpcun_post)
% colormap('hot')
% c = colorbar('YTick',[0 0.2 0.4 0.6 0.8 1],'FontSize',15);
% axpos = get(gca,'Position');
% cpos = get(c,'Position');
% cpos(3) = 0.5*cpos(3);
% set(c,'Position',cpos);
% set(gca,'Position',axpos);
% set(gca,'XTick',[1 10 20 30 40 50 60 68],'YTick',[1 10 20 30 40 50 60 68],'FontSize',15)
% xlabel('Area ID','FontWeight','bold')
% ylabel('Area ID','FontWeight','bold')
% savePath = '../../reports/figures/fc_lpcun_post';
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Graph Metrics:Immediate FC effects %%%%%%%%%%%%%%%%%%%%
% graph_metrics
% savePath = '../../reports/figures/gm1';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Graph Metrics:FC effects after Re-establishing E-I balance %%%%%%%%%%%%%%%%%%%%
% graph_metrics
% savePath = '../../reports/figures/gm2';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%% Balance loss vs particpation Coefficient for SC downsampled
% %%%%%%%%%%%%% from 998 areas SC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2);
% plot(partCoefSCEmp,balLoss,'ko','MarkerFaceColor','k','MarkerSize',5);
% ylabel('No. of Areas lost E-I balance','FontWeight','bold');
% xlabel('Participation Coefficient','FontWeight','bold');
% hold on;
% box off;
% p=polyfit(partCoefSCEmp,balLoss,1);
% pred = polyval(p,partCoefSCEmp);
% plot(partCoefSCEmp,pred,'k');
% text(3,18,['r=',sprintf('%.2f',rho),', p<0.001'],'FontSize',15,'BackgroundColor',[0.3 0.7 0.7]);
% xlim([2 22])
% ylim([2 20])
% set(gca,'XTick',2:4:22,'YTick',2:4:20,'FontSize',13)
% savePath = '../../reports/figures/balance_loss_vs_partCoef_resampled_roi_avg';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% FC recovery across subjects and location%%%%%%%%%%%%%%%
figure(1);
% individual_graph_metrics;
% boxplot(fcDistancePrior(:,idxs));
% box off;
% ylim([1,25]);
% set(gca,'XTick',1:length(idxs),'XTickLabel',{'16','21','22','23','24','28','30','50','55','56','57','58','62','64'},'FontSize',8);
% text(2,-4,'Lesion Location Index','FontWeight','bold','FontSize',12)
% ylabel('FC Distance','FontWeight','bold','FontSize',12)
% savePath = '../../reports/figures/fcrecovery_subjects_prior'

boxplot(fcDistancePost(:,idxs));
box off;
ylim([1,25]);
set(gca,'XTick',1:length(idxs),'XTickLabel',{'16','21','22','23','24','28','30','50','55','56','57','58','62','64'},'FontSize',8);
text(2,-4,'Lesion Location Index','FontWeight','bold','FontSize',12)
ylabel('FC Distance','FontWeight','bold','FontSize',12)
savePath = '../../reports/figures/fcrecovery_subjects_post'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%