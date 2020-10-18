%%
%{
function [Z2,realZ2]=subWave2Z1(NNbar,xdata,ydataNX,xsto,vaxparams)
burn=10000;
int=500;
NNtot=sum(NNbar);
realZ2=nansum(ydataNX(xdata>34,:),1);
realZ2=realZ2'/NNtot*1000;
mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
sig=[129.6984 329.7708 231.8866 134.0024 58.4746];
prior=xsto(burn+1:int:end,:);
lp=size(prior,1);
ly=size(ydataNX,1);
nbar=length(mu);
%numNorm=100;%Number of multiplier samples per MCMC
%mults=normrnd(mu',sig',[nbar,numNorm]);
Z2=zeros(nbar,lp);
t0=243;
for i=1:lp
    ri=normrnd(mu,sig);
    ydata=ydataNX.*repmat(ri,ly,1);
    y0=nansum(ydata(xdata<35,:),1)';
    I0=ydata(xdata==34,:)'/30;
    I0(isnan(I0))=0;
    if isempty(I0)==1
        I0=zeros(nbar,1);
    end
    y0=y0-I0;
    [~,~,z2]=subPandemicSimulationVaxZ1(NNbar,prior(i,:),xdata,0,0,ydata,243,vaxparams,t0,y0,I0);%mcmc=1
    Z2(:,i)=z2./ri/NNtot*1000;
end
%}
%%
%
function [Z2,realZ2]=subWave2Z1(NNbar,xdata,ydata,ydataNX,thetac,vaxparams)
NNtot=sum(NNbar);
realZ2=nansum(ydataNX(xdata>34,:),1);
realZ2=realZ2'/NNtot*1000;
[xsto,~,~,~,~]=subMCMC(NNbar,xdata,ydataNX,thetac,vaxparams);
[~,z2]=subUncertaintyPlot(NNbar,xsto,xdata,ydata,thetac(3:end),vaxparams);
lx=size(z2,1);
Z2=z2/NNtot*1000;%/normrnd(repmat(mu',1,lx),repmat(sig',1,lx))*1000;
%}
%{
y0=ydataNX(xdata<18,:);
y0=sum(y0,1)';
y0mu=y0.*mu;
y0sig=y0.*sig;
%}
