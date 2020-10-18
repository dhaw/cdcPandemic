<<<<<<< HEAD
function f=cdcPlotDataUncertainty(y)%,mu,sig)
mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
sig=[129.6984 329.7708 231.8866 134.0024 58.4746];
nbar=5;
xdata=3:36;
xcut=16;
ly=size(y,1);
murep=repmat(mu,ly,1);
sigrep=repmat(sig,ly,1);
ymu=y.*mu;
yminus=y.*(murep-sigrep);
yplus=y.*(murep+sigrep);
fs=12; lw=1;
cmap=lines(nbar);
figure
hold on
for i=1:nbar
    cmapi=cmap(i,:);
    plot(xdata+xcut,ymu(xdata,i),'-','linewidth',2,'color',cmapi)
    plot(xdata+xcut,yminus(xdata,i),'-','linewidth',lw,'color',cmapi)
    plot(xdata+xcut,yplus(xdata,i),'-','linewidth',lw,'color',cmapi)
end
xlabel('Time (weeks)')
ylabel('Incidence')
set(gca,'fontsize',fs)
axis([xdata(1)+xcut,xdata(end)+xcut,0,max(max(yplus))])
grid on
grid minor
box on
=======
function f=cdcPlotDataUncertainty(y)%,mu,sig)
mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
sig=[129.6984 329.7708 231.8866 134.0024 58.4746];
nbar=5;
xdata=3:36;
xcut=0;%16;
ly=size(y,1);
murep=repmat(mu,ly,1);
sigrep=repmat(sig,ly,1);
y2=y.*murep;
y1=y.*(murep-sigrep);
y3=y.*(murep+sigrep);
fs=12; lw=1;
cmap=lines(nbar);
figure
maxy=max(max(y3));
tvec=xdata;%+xcut
tvec2= [tvec,fliplr(tvec)];
hold on
for i=1:nbar
    cmapi=cmap(i,:);
    
    inBetween=[y1(xdata,i);flipud(y3(xdata,i))];
    fill(tvec2,inBetween,cmapi,'facealpha',.2);
    
    plot(tvec,y1(xdata,i),'-','linewidth',1,'color',cmapi);
    plot(tvec,y3(xdata,i),'-','linewidth',1,'color',cmapi);
    
    plot(tvec,y2(xdata,i),'-','linewidth',2,'color',cmapi)
    %plot(xdata+xcut,yminus(xdata,i),'-','linewidth',lw,'color',cmapi)
    %plot(xdata+xcut,yplus(xdata,i),'-','linewidth',lw,'color',cmapi)
end

xlabel('Time (weeks)')
ylabel('Incidence')
set(gca,'fontsize',fs)
%axis([xdata(1)+xcut,xdata(end)+xcut,0,max(max(y3))])
axis([xdata(1),xdata(end),0,maxy])
grid on
grid minor
box on
>>>>>>> 7109d6dbccee4bd48cd4f166634736f1a390a7be
