function [f,g]=cdcCDF(xsto,factor)
plotStuff=0;
lx=length(xsto);
burn=200;
qmean=1./[143.4400,364.7100,148.2000,64.6700]*factor;%Row - as z2
%Nvec=(160000:1000:260000)';
Nvec=(90000:200:240000)';%(550000:100:558000)';%
Nmat=repmat(Nvec,1,1);%Last arg=1 or 4
lN=length(Nvec);
x=zeros(lN,1);%Last arg=1 or 4
for i=burn+1:10:lx
    [~,~,z2]=cdcPandemicSimulation([xsto(i,1),xsto(i,2)],4:7,0,0,0,243);
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
    %{
    %For schematic with shaded area:
    here=11;
    colFill=.9*[1,1,1];
    Nstar=Nvec(here);
    xstar=x(here);
    set(gca,'xtick',[Nstar],'xticklabel',{'N^*'})
    x2=[Nvec(here:end);flipud(Nvec(here:end))];
    inBetween=[x(here:end);zeros(length(x)-here+1,1)];
    fill(x2,inBetween,colFill)
    plot([Nstar,Nstar],[0,xstar],'k-','linewidth',lw)
    text(110650,.2,'Y^*','fontsize',15)
    %}
    plot(Nvec,x,'-','linewidth',lw,'color',cmap(1,:));
    xlabel('N')
    ylabel('p(H>N)')
    %legend('0-4','5-19','2--64','65+')
    set(gca,'fontsize',fs)
    axis([Nvec(1),Nvec(end),0,1])
    grid on
    grid minor
    box on
end