function f=cdcLhoodsW5old(params,ydata)
%Input data - hospitalisations (i.e. without multipliers)
xdata=16:32;%Match MLE
lx=length(xdata);
NNbar=[1464566;3790730;10236474;3984200;2412129];
NNmat=repmat(NNbar',lx,1);
NNvec=kron(NNbar,ones(lx,1));
nbar=length(NNbar);
ydata=ydata(xdata,:);
ydata=reshape(ydata,nbar*lx,1);
%
%Scale incidence back up to absolute:
%ydata=ydata.*NNvec;
%}
%Input data without multipliers
qmean=1./[393.0256 999.3054 702.6867 406.0680 177.1958];
qsd=1./[129.6984 329.7708 231.8866 134.0024 58.4746];
qvec=repmat(qmean,lx,1);
qvec=reshape(qvec,nbar*lx,1);

%ydata=ydata./qvec;

ysim=cdcPandemicSimulationW5(params,xdata,0,0,0,243).*NNmat;
ymean=reshape(ysim,nbar*lx,1).*qvec;
ysd=ymean.*(1-qvec);
ysd=10*sqrt(ysd);%W5*

%Normal likelihood:
L=-log(normpdf(ydata,ymean,ysd));
f=sum(L);
end