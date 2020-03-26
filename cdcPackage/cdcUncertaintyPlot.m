<<<<<<< HEAD
function f=cdcUncertaintyPlot(xsto,ydata)%,simTheta)
=======
function f=cdcUncertaintyPlot(xsto,ydata,theta)%,simTheta)
>>>>>>> 7109d6dbccee4bd48cd4f166634736f1a390a7be
lx=size(xsto,1);
burn=5000;
int=200;
xdata=1:36;
<<<<<<< HEAD
xbounds=[5,9,21,28];
xcut=16;%Match data - removed from start of year
thresh=9;%Match with cdcLhoodsW5
=======
%xbounds=15;%[5,9,21,28];
xcut=16;%Match data - removed from start of year
thresh=15;%Match with cdcLhoodsW5
xbounds=[1,15]+xcut;%[5,9,21,28];

>>>>>>> 7109d6dbccee4bd48cd4f166634736f1a390a7be
pands=zeros(length(xdata),floor((lx-burn)/int));
j=1;
for i=burn+1:int:lx
    pandsimi=cdcPandemicSimulationW5(xsto(i,:),xdata,0,0,ydata,243);
    pands(:,j)=sum(pandsimi,2);
    j=j+1;
end
%}
pandsq=prctile(pands,[25,50,75],2);
y1=pandsq(:,1);
y2=pandsq(:,2);
y3=pandsq(:,3);

fs=12; lw=.5;
<<<<<<< HEAD
maxy=max(max(max([pands])),max(sum(ydata,2)));%,ydata(1:thresh,:)]));
=======
maxy=max(y3);%max(max(pands)),max(sum(ydata,2)));%,ydata(1:thresh,:)]));
>>>>>>> 7109d6dbccee4bd48cd4f166634736f1a390a7be
col=lines(7);
col1=col(1,:);
figure
hold on
%plot([thresh+xcut,thresh+xcut],[0,maxy],'--','linewidth',2,'color',col(1,:))
<<<<<<< HEAD
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
=======
%plot([xbounds(1),xbounds(1)],[0,maxy],'--','linewidth',2,'color',col(1,:))
plot([xbounds(2),xbounds(2)],[0,maxy],'--','linewidth',2,'color',col(1,:))
%{
plot(xdata(1:thresh)+xcut,pands(1:thresh,:),'linewidth',lw,'color',[.5,.5,.5]);
plot(xdata(thresh:end)+xcut,pands(thresh:end,:),'linewidth',lw,'color','k');
%}
%
tvec=xdata+xcut;
tvec2=[tvec,fliplr(tvec)];
inBetween=[y1;flipud(y3)];
fill(tvec2,inBetween,col1,'facealpha',.2);
plot(tvec,y1,'-','linewidth',1,'color',col1);
plot(tvec,y3,'-','linewidth',1,'color',col1);
h1=plot(tvec,y2,'-','linewidth',2,'color',col1);
%}
%
plotData=sum(ydata,2);
h2=plot((1:thresh)+xcut,plotData(1:thresh),'k-','linewidth',2,'color','k');%col(2,:));
plot((thresh:xdata(end))+xcut,plotData(thresh:end),'k--','linewidth',2,'color',[0,0,0])%col(2,:))
%}
%Plot a simulation:
%{
mle=cdcPandemicSimulationW5(theta(4:end),xdata,0,0,0,243);
h3=plot(xdata+xcut,sum(mle,2),'-','linewidth',2,'color',.5*[1,1,1]);%col(3,:));
%}
set(gca,'fontsize',fs)
xlabel('Time (weeks)')
ylabel('Incidence')
axis([xdata(1)+xcut,xdata(end)+xcut,0,maxy])
%legend([h2,h1,h3],'Data','MCMC','MLE','location','NW')
>>>>>>> 7109d6dbccee4bd48cd4f166634736f1a390a7be
grid on
grid minor
box on

f=pands;