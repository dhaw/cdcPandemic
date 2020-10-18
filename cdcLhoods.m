function f=cdcLhoods(params,ydata)
%Input data - hospitalisations (i.e. without multipliers)
xdata=4:8;%Match MLE
lx=length(xdata);
ydata=ydata(xdata,:);
ydata=reshape(ydata,4*lx,1);
%
%Scale incidence back up to absolute:
NNbar=[19169690;62121035;184015269;39570590];
NNmat=repmat(NNbar',lx,1);
NNvec=kron(NNbar,ones(lx,1));
ydata=ydata.*NNvec;
%}
%Input data without multipliers
na=4;%Number of age groups output - 1 or 4 - ****match with byAge in pandemic1DallV****
qmean=1./[143.4400,364.7100,148.2000,64.6700];
qvec=repmat(qmean,lx,1);
qvec=reshape(qvec,4*lx,1);

%ydata=ydata./qvec;

ysim=cdcPandemicSimulation(params,xdata,0,0,0,243).*NNmat;
ymean=reshape(ysim,na*lx,1).*qvec;
ysd=ymean.*(1-qvec);
ysd=10*sqrt(ysd);

%Normal likelihood:
L=-log(normpdf(ydata,ymean,ysd));
f=sum(L);
end