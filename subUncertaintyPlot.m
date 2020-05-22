function [f,h2mean]=subUncertaintyPlot(NNbar,xsto,xdata,ydata,theta)%,simTheta)
thisage=0;%=0 for total
col=[0.2422    0.1504    0.6603;
    0.2504    0.1650    0.7076;
    0.2578    0.1818    0.7511;
    0.2647    0.1978    0.7952;
    0.2706    0.2147    0.8364];%lines(7);
if thisage>0
    col1=col(thisage,:);
else
    col1=.5*[1,1,1];
end
vlinecol=0*[1,1,1];%col(1,:);

if thisage==0
    plotData=sum(ydata,2);
else
    plotData=ydata(:,thisage);
end

lx=size(xsto,1);
burn=1000;
int=100;
%xdata=1:36;
%xbounds=15;%[5,9,21,28];
thresh=15;%Match with cdcLhoodsW5
xbounds=[1,35];%[5,9,21,28];

mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
nbar=size(mu,2);

toRun=floor((lx-burn)/int);
pands=zeros(length(xdata),toRun);
Z2=zeros(toRun,nbar);
j=1;
for i=burn+1:int:lx
    [pandsimi,~,zi]=subPandemicSimulation(NNbar,xsto(i,:),xdata,0,0,ydata,243);
    %Sum or specific age group
    if thisage==0
        pands(:,j)=sum(pandsimi,2);
    else
        pands(:,j)=pandsimi(:,thisage);
    end
    Z2(j,:)=zi;
    j=j+1;
end
h2mean=sum(mean(Z2,1)./mu);
%}
pandsq=prctile(pands,[25,50,75],2);
y1=pandsq(:,1);
y2=pandsq(:,2);
y3=pandsq(:,3);

fs=10; lw=.5;
maxy=max(max(y3),max(plotData));%max(max(pands)),max(sum(ydata,2)));%,ydata(1:thresh,:)]));
figure
hold on
%plot([thresh+xcut,thresh+xcut],[0,maxy],'--','linewidth',2,'color',col(1,:))
%plot([xbounds(1),xbounds(1)],[0,maxy],'--','linewidth',2,'color',col(1,:))
plot([xbounds(2),xbounds(2)],[0,maxy],'--','linewidth',2,'color',vlinecol)
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
h2=plot(xdata,plotData,'k-','linewidth',2,'color','k');%col(2,:));
%plot((thresh:xdata(end))+xcut,plotData(thresh:end),'k--','linewidth',2,'color',[0,0,0])%col(2,:))
%}
maxy=max(maxy,max(plotData));
%Plot a simulation:
%
[mle,~,~]=subPandemicSimulation(NNbar,theta,xdata,0,0,0,243);
if thisage==0
    mleSum=sum(mle,2);
else
    mleSum=mle(:,thisage);
end
h3=plot(xdata,mleSum,'--','linewidth',2,'color',.5*[1,1,1]);%col(3,:));
maxy=max(maxy,max(mleSum));
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