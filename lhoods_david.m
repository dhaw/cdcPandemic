function f=lhoods_david(params,ydata)
xdata=4:12;%Match MLE

%Input data without multipliers
na=4;%Number of age groups output - 1 or 4 - ****match with byAge****
qmean=1./[143.4400,364.7100,148.2000,64.6700];
lx=length(xdata);
qvec=repmat(qmean,length(xdata),1);
qvec=reshape(qvec,4*lx,1);

ydata=ydata./qvec;

ysim=pandemic1DallV(params,xdata,0,0,0);
ymean=reshape(ysim,na*lx,1);
ysd=sqrt(ymean.*(1-qvec))/100;
%Normal likelihood:
L=-sum(log(normpdf(ydata,ymean,ysd)));
f=L;
end
%{
function f=Lqtheta(qvec,params,xdata,ydata,lx)
y=pandemic1DallV(params,xdata,0,1,0);
ymean=reshape(y,4*lx,1);
ysd=sqrt(ymean.*(1-qvec));
f=normpdf(ydata,ymean,ysd);
end
%}