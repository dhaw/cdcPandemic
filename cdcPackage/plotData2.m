figure;
fs=10; lw=2;
cmap=parula(10);
colormap parula
hold on;
plot(xus,nansum(yusx,2)/sum(NNbarUS)*1000,'k-','linewidth',3);
plot(xca,nansum(ycax,2)/sum(NNbarCA)*1000,'color',cmap(1,:),'linewidth',lw);
plot(xco,nansum(ycox,2)/sum(NNbarCO)*1000,'color',cmap(2,:),'linewidth',lw);
plot(xct,nansum(yctx,2)/sum(NNbarCT)*1000,'color',cmap(3,:),'linewidth',lw);
plot(xga,nansum(ygax,2)/sum(NNbarGA)*1000,'color',cmap(4,:),'linewidth',lw);
plot(xmd,nansum(ymdx,2)/sum(NNbarMD)*1000,'color',cmap(5,:),'linewidth',lw);
plot(xmn,nansum(ymnx,2)/sum(NNbarMN)*1000,'color',cmap(6,:),'linewidth',lw);
plot(xnmNew,nansum(ynmNewx,2)/sum(NNbarNM)*1000,'color',cmap(7,:),'linewidth',lw);
plot(xny,nansum(ynyx,2)/sum(NNbarNY)*1000,'color',cmap(8,:),'linewidth',lw);
plot(xor,nansum(yorx,2)/sum(NNbarOR)*1000,'color',cmap(9,:),'linewidth',lw);
plot(xtn,nansum(ytnx,2)/sum(NNbarTN)*1000,'color',cmap(10,:),'linewidth',lw);
plot(xus,nansum(yusx,2)/sum(NNbarUS)*1000,'k-','linewidth',2);
legend('US','CA','CO','CT','GA','MD','MN','NM','NY','OR','TN','location','NW')
box on; grid on; grid minor;
xlabel('Time (weeks)')
ylabel('Hosp./1000')
axis([18,52,0,.12])%max(nansum(yusx,2))]);
set(gca,'fontsize',fs)

numStates=10;
numStacks=2;
numAge=5;
labels={'CA','CO','CT','GA','MD','MN','NM','NY','OR','TN'};
X=zeros(numStates,numStacks,numAge);
X(1,1,:)=nansum(ycax(xca<35,:),1)/sum(NNbarCA)*1000;
X(1,2,:)=nansum(ycax(xca>34,:),1)/sum(NNbarCA)*1000;
X(2,1,:)=nansum(ycox(xco<35,:),1)/sum(NNbarCO)*1000;
X(2,2,:)=nansum(ycox(xco>34,:),1)/sum(NNbarCO)*1000;
X(3,1,:)=nansum(yctx(xct<35,:),1)/sum(NNbarCT)*1000;
X(3,2,:)=nansum(yctx(xct>34,:),1)/sum(NNbarCT)*1000;
X(4,1,:)=nansum(ygax(xga<35,:),1)/sum(NNbarGA)*1000;
X(4,2,:)=nansum(ygax(xga>34,:),1)/sum(NNbarGA)*1000;
X(5,1,:)=nansum(ymdx(xmd<35,:),1)/sum(NNbarMD)*1000;
X(5,2,:)=nansum(ymdx(xmd>34,:),1)/sum(NNbarMD)*1000;
X(6,1,:)=nansum(ynyx(xny<35,:),1)/sum(NNbarMN)*1000;
X(6,2,:)=nansum(ynyx(xny>34,:),1)/sum(NNbarMN)*1000;
X(7,1,:)=nansum(ynmNewx(xnmNew<35,:),1)/sum(NNbarNM)*1000;
X(7,2,:)=nansum(ynmNewx(xnmNew>34,:),1)/sum(NNbarNM)*1000;
X(8,1,:)=nansum(yorx(xor<35,:),1)/sum(NNbarNY)*1000;
X(8,2,:)=nansum(yorx(xor>34,:),1)/sum(NNbarNY)*1000;
X(9,1,:)=nansum(ycax(xca<35,:),1)/sum(NNbarOR)*1000;
X(9,2,:)=nansum(ycax(xca>34,:),1)/sum(NNbarOR)*1000;
X(10,1,:)=nansum(ytnx(xtn<35,:),1)/sum(NNbarTN)*1000;
X(10,2,:)=nansum(ytnx(xtn>34,:),1)/sum(NNbarTN)*1000;
plotBarStackGroups(X,labels);
