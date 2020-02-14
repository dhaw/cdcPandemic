function f=cdcUncertaintyPlot(xsto,ydata)%,simTheta)
lx=size(xsto,1);
burn=5000;
int=200;
xdata=1:36;
xbounds=[5,9,21,28];
xcut=16;%Match data - removed from start of year
thresh=9;%Match with cdcLhoodsW5
pands=zeros(length(xdata),floor((lx-burn)/int));
for i=burn+1:int:lx
    pandsimi=cdcPandemicSimulationW5(xsto(i,:),xdata,0,0,ydata,243);
    pands(:,i)=sum(pandsimi,2);
end
fs=12; lw=.5;
maxy=max(max(max([pands])),max(sum(ydata,2)));%,ydata(1:thresh,:)]));
col=lines(7);
figure
hold on
%plot([thresh+xcut,thresh+xcut],[0,maxy],'--','linewidth',2,'color',col(1,:))
plot([xbounds(1),xbounds(1)],[0,maxy],'--','linewidth',2,'color',col(1,:))
plot([xbounds(2),xbounds(2)],[0,maxy],'--','linewidth',2,'color',col(1,:))
plot(xdata(1:thresh)+xcut,pands(1:thresh,:),'linewidth',lw,'color',[.5,.5,.5]);
plot(xdata(thresh:end)+xcut,pands(thresh:end,:),'linewidth',lw,'color','k');
%simTheta=sum(simTheta(1:thresh,:),2);
%plot((1:thresh)+xcut,simTheta,'k-','linewidth',2)
plotData=sum(ydata,2);
plot((1:thresh)+xcut,plotData(1:thresh),'k-','linewidth',2)
plot((thresh:xdata(end))+xcut,plotData(thresh:end),'k--','linewidth',2)
set(gca,'fontsize',fs)
xlabel('Time (weeks)')
ylabel('Hospitalisations')
axis([xdata(1)+xcut,xdata(end)+xcut,0,maxy])
grid on
grid minor
box on