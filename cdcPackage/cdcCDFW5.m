function [f,g]=cdcCDFW5(xsto,factor)
plotStuff=0;
lx=length(xsto);
burn=200;
qmean=1./[393.0256 999.3054 702.6867 406.0680 177.1958]*factor;%Row - as z2
qsd=[129.6984 329.7708 231.8866 134.0024 58.4746];
%Nvec=(160000:1000:260000)';
Nvec=(90000:200:240000)';%(550000:100:558000)';%
Nmat=repmat(Nvec,1,1);%Last arg=1 or 4
lN=length(Nvec);
x=zeros(lN,1);%Last arg=1 or 4
for i=burn+1:10:lx
    [~,~,z2]=cdcPandemicSimulationW5([xsto(i,1),xsto(i,2)],4:7,0,0,0,243);
    mu=repmat(z2.*qmean,lN,1);%lN times 4
    sig=mu.*(1-repmat(qmean,lN,1));
    sigsum=sum(sig,2);
    sigsum=10*sqrt(sigsum);%Match sd in likelihoods %W5*
    %xadd=normcdf(2*mu-Nmat,mu,sig);
    musum=sum(mu,2);
    xadd=normcdf(2*musum-Nmat,musum,sigsum);
    x=x+xadd;
end
lx2=floor((lx-burn)/10);
x=x/lx2;
f=x;
g=Nvec;%flipud(cumsum(flipud(x)));

if plotStuff==1
    fs=12; lw=2;
    cmap=lines(7);
    figure
    hold on
    plot(Nvec,x,'-','linewidth',lw,'color',cmap(1,:));
    xlabel('N')
    ylabel('p(H>N)')
    set(gca,'fontsize',fs)
    axis([Nvec(1),Nvec(end),0,1])
    grid on
    grid minor
    box on
end