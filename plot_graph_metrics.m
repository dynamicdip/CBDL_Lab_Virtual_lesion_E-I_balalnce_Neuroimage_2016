%%%%%%% Graph Metric plots WOFIC%%%%%%%
% subplot(2,3,1);
% hold on;
% xlabel('Clustering Coefficient')
% ylabel('FC\_Emp vs FC\_Sim Correlation')
% p=polyfit(clstrcoef,fccorr_wofic,1);
% r = corrcoef(clstrcoef,fccorr_wofic);
% r = r(1,2);
% pred = polyval(p,clstrcoef);
% plot(clstrcoef,fccorr_wofic,'ko');
% plot(clstrcoef,pred);
% gtext(['r=',num2str(r)],'FontSize',20)
% 
% subplot(2,3,2);
% hold on;
% xlabel('Node Betweenness Centrality')
% ylabel('FC\_Emp vs FC\_Sim Correlation')
% r = corrcoef(nbc,fccorr_wofic);
% r = r(1,2);
% p=polyfit(nbc,fccorr_wofic,1);
% pred = polyval(p,nbc);
% plot(nbc,fccorr_wofic,'ko');
% plot(nbc,pred);
% gtext(['r=',num2str(r)],'FontSize',20)
% 
% subplot(2,3,3);
% hold on;
% xlabel('Path Length')
% ylabel('FC\_Emp vs FC\_Sim Correlation')
% r = corrcoef(pathlength,fccorr_wofic);
% r = r(1,2);
% p=polyfit(pathlength,fccorr_wofic,1);
% pred = polyval(p,pathlength);
% plot(pathlength,fccorr_wofic,'ko');
% plot(pathlength,pred);
% gtext(['r=',num2str(r)],'FontSize',20)
% 
% subplot(2,3,4);
% hold on;
% xlabel('Modularity')
% ylabel('FC\_Emp vs FC\_Sim Correlation')
% r = corrcoef(modularity,fccorr_wofic);
% r = r(1,2);
% p=polyfit(modularity,fccorr_wofic,1);
% pred = polyval(p,modularity);
% plot(modularity,fccorr_wofic,'ko');
% plot(modularity,pred);
% gtext(['r=',num2str(r)],'FontSize',20)
% 
% subplot(2,3,5);
% hold on;
% xlabel('Global Efficiency')
% ylabel('FC\_Emp vs FC\_Sim Correlation')
% r = corrcoef(globalefficiency,fccorr_wofic);
% r = r(1,2);
% p=polyfit(globalefficiency,fccorr_wofic,1);
% pred = polyval(p,globalefficiency);
% plot(globalefficiency,fccorr_wofic,'ko');
% plot(globalefficiency,pred);
% gtext(['r=',num2str(r)],'FontSize',20)

%%%%%%% Graph Metric plots FIC%%%%%%%
subplot(2,3,1);
hold on;
xlabel('Clustering Coefficient')
ylabel('FC\_Emp vs FC\_Sim Correlation')
r = corrcoef(clstrcoef,fccorr_fic);
r = r(1,2);
p=polyfit(clstrcoef,fccorr_fic,1);
pred = polyval(p,clstrcoef);
plot(clstrcoef,fccorr_fic,'ko');
plot(clstrcoef,pred);
gtext(['r=',num2str(r)],'FontSize',20)

subplot(2,3,2);
hold on;
xlabel('Node Betweenness Centrality')
ylabel('FC\_Emp vs FC\_Sim Correlation')
r = corrcoef(nbc,fccorr_fic);
r = r(1,2);
p=polyfit(nbc,fccorr_fic,1);
pred = polyval(p,nbc);
plot(nbc,fccorr_fic,'ko');
plot(nbc,pred);
gtext(['r=',num2str(r)],'FontSize',20)

subplot(2,3,3);
hold on;
xlabel('Path Length')
ylabel('FC\_Emp vs FC\_Sim Correlation')
r = corrcoef(pathlength,fccorr_fic);
r = r(1,2);
p=polyfit(pathlength,fccorr_fic,1);
pred = polyval(p,pathlength);
plot(pathlength,fccorr_fic,'ko');
plot(pathlength,pred);
gtext(['r=',num2str(r)],'FontSize',20)

subplot(2,3,4);
hold on;
xlabel('Modularity')
ylabel('FC\_Emp vs FC\_Sim Correlation')
r = corrcoef(modularity,fccorr_fic);
r = r(1,2);
p=polyfit(modularity,fccorr_fic,1);
pred = polyval(p,modularity);
plot(modularity,fccorr_fic,'ko');
plot(modularity,pred);
gtext(['r=',num2str(r)],'FontSize',20)

subplot(2,3,5);
hold on;
xlabel('Global Efficiency')
ylabel('FC\_Emp vs FC\_Sim Correlation')
r = corrcoef(globalefficiency,fccorr_fic);
r = r(1,2);
p=polyfit(globalefficiency,fccorr_fic,1);
pred = polyval(p,globalefficiency);
plot(globalefficiency,fccorr_fic,'ko');
plot(globalefficiency,pred);
gtext(['r=',num2str(r)],'FontSize',20)