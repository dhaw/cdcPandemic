function f=lhoods(data,factor)
%Input data without multipliers
na=4;%Number of age groups output - 1 or 4 - ****match with byAge****
if na==1
    data=sum(data,2);
    NNbar=[19169690;62121035;184015269;39570590];
    %To do: ****define global q****
else
    qmean=1./[143.4400,364.7100,148.2000,64.6700];
end
xdata=4:7;%4:7 for 1st wave
lx=length(xdata);

ydata=data(xdata,:);
ydata=reshape(ydata,na*lx,1);


qvec=repmat(qmean,length(xdata),1);
qvec=reshape(qvec,na*lx,1);

%params=[factor,.9,.77];
%params=[3.2740,0.9284,0.7383];
params=[3.2740,0.9284,factor];

ydata=ydata./qvec;
%f=Lqtheta(qvec,params,xdata,ydata,lx);

ysim=pandemic1DallV(params,xdata,0,0,0);
ymean=reshape(ysim,na*lx,1);
ysd=sqrt(ymean.*(1-qvec));
%Normal likelihood:
L=log(normpdf(ydata,ymean,ysd));
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