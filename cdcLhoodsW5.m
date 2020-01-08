function f=cdcLhoodsW5(params,ydataNX)
%Input data - hospitalisations (i.e. without multipliers)
xdata=1:35;%Match MLE
lx=length(xdata);
NNbar=[1464566;3790730;10236474;3984200;2412129];
nbar=length(NNbar);
ydataNX=ydataNX(xdata,:);
hosp=reshape(ydataNX,nbar*lx,1);
%
%Scale incidence back up to absolute:
%ydata=ydata.*NNvec;
%}
%Input data without multipliers

mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
sig=[129.6984 329.7708 231.8866 134.0024 58.4746];

muvec=reshape(repmat(mu,lx,1),nbar*lx,1);
sigvec=reshape(repmat(sig,lx,1),nbar*lx,1);

ysim=cdcPandemicSimulationW5(params,xdata,0,0,0,243);
%ysim=ysim(xdata,:);%simCut in cdcPandemicSimulationW5
ysim=reshape(ysim,nbar*lx,1);

y=ysim./hosp;
y(hosp==0)=0;

%Normal likelihood:
L=-log(normpdf(y,muvec,sigvec));
f=sum(L);
end