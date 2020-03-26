function [Z2,realZ2]=subWave2(NNbar,xdata,ydataNX,xsto)
NNtot=sum(NNbar);
realZ2=nansum(ydataNX(xdata>34,:),1);
realZ2=realZ2'/NNtot*1000;

burn=10000;
int=500;
prior=xsto(burn+1:int:end,:);
lp=size(prior,1);
ly=size(ydataNX,1);

mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
sig=[129.6984 329.7708 231.8866 134.0024 58.4746];
nbar=length(mu);
%numNorm=100;%Number of multiplier samples per MCMC
%mults=normrnd(mu',sig',[nbar,numNorm]);

Z2=zeros(nbar,lp);
for i=1:lp
    ri=normrnd(mu,sig);
    [~,~,z2]=subPandemicSimulation(NNbar,prior(i,:),xdata,0,0,ydataNX.*repmat(ri,ly,1),243);%mcmc=1
    Z2(:,i)=z2'./ri'/NNtot*1000;
end
%{
y0=ydataNX(xdata<18,:);
y0=sum(y0,1)';
y0mu=y0.*mu;
y0sig=y0.*sig;
%}
