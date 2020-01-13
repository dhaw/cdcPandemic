function f=cdcUncertaintyPlot(xsto,ydata)
lx=size(xsto,1);
burn=1000;
int=100;
xdata=1:35;
thresh=17;%Match with cdcLhoodsW5
pands=zeros(length(xdata),floor((lx-burn)/int));
for i=burn+1:int:lx
    pandsimi=cdcPandemicSimulationW5(xsto(i,:),xdata,0,0,ydata,243);
    pands(:,i)=sum(pandsimi,2);
end
fs=12; lw=.5;
maxy=max(max(pands));
col=lines(7);
figure
hold on
plot([thresh,thresh],[0,maxy],'--','linewidth',2,'color',col(1,:))
plot(xdata,pands,'linewidth',lw,'color','k');
set(gca,'fontsize',fs)
xlabel('Time (weeks)')
ylabel('Incidence')
axis([xdata(1),xdata(end),0,maxy])
grid on
grid minor
box on