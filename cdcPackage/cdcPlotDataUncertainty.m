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
