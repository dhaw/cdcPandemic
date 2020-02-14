function [f,h2mean]=subUncertaintyPlot(NNbar,xsto,xdata,ydata,theta)%,simTheta)
lx=size(xsto,1);
burn=1000;
int=100;
%xdata=1:36;
%xbounds=15;%[5,9,21,28];
thresh=15;%Match with cdcLhoodsW5
xbounds=[1,15];%[5,9,21,28];

mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
nbar=size(mu,2);

toRun=floor((lx-burn)/int);
pands=zeros(length(xdata),toRun);
Z2=zeros(toRun,nbar);
j=1;
for i=burn+1:int:lx
    [pandsimi,~,zi]=subPandemicSimulation(NNbar,xsto(i,:),xdata,0,0,ydata,243);
    pands(:,j)=sum(pandsimi,2);
    Z2(j,:)=zi;
    j=j+1;
end
h2mean=sum(mean(Z2,1)./mu);
%}
pandsq=prctile(pands,[25,50,75],2);
y1=pandsq(:,1);
y2=pandsq(:,2);
y3=pandsq(:,3);

fs=12; lw=.5;
maxy=max(y3);%max(max(pands)),max(sum(ydata,2)));%,ydata(1:thresh,:)]));
col=lines(7);
col1=col(1,:);
figure
hold on
%plot([thresh+xcut,thresh+xcut],[0,maxy],'--','linewidth',2,'color',col(1,:))
%plot([xbounds(1),xbounds(1)],[0,maxy],'--','linewidth',2,'color',col(1,:))
plot([xbounds(2),xbounds(2)],[0,maxy],'--','linewidth',2,'color',col(1,:))
%{
plot(xdata(1:thresh)+xcut,pands(1:thresh,:),'linewidth',lw,'color',[.5,.5,.5]);
plot(xdata(thresh:end)+xcut,pands(thresh:end,:),'linewidth',lw,'color','k');
%}
%
tvec=xdata;
tvec2=[tvec;flipud(tvec)];
inBetween=[y1;flipud(y3)];
fill(tvec2,inBetween,col1,'facealpha',.2);
plot(tvec,y1,'-','linewidth',1,'color',col1);
plot(tvec,y3,'-','linewidth',1,'color',col1);
h1=plot(tvec,y2,'-','linewidth',2,'color',col1);
%}
%
plotData=sum(ydata,2);
h2=plot(xdata,plotData,'k-','linewidth',2,'color','k');%col(2,:));
%plot((thresh:xdata(end))+xcut,plotData(thresh:end),'k--','linewidth',2,'color',[0,0,0])%col(2,:))
%}
%Plot a simulation:
%
[mle,~,~]=subPandemicSimulation(NNbar,theta(4:end),xdata,0,0,0,243);
h3=plot(xdata,sum(mle,2),'-','linewidth',2,'color',.5*[1,1,1]);%col(3,:));
%}
set(gca,'fontsize',fs)
xlabel('Time (weeks)')
ylabel('Incidence')
axis([xdata(1),xdata(end),0,maxy])
legend([h2,h1,h3],'Data','MCMC','MLE','location','NW')
grid on
grid minor
box on

f=pands;