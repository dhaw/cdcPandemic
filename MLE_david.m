function f=MLE_david(ydata)

%ydata=data(xdata,:);
%ydata=reshape(ydata,lx,1);
%params0=[3,.58,243,.95];
params0=[3,.58,243,.95];
%params0=[2.3,.95,.61];
%[phat,pci]=mle(ydata,'pdf',@(params)pandemic1DallV(params,xdata,0,0,0),'start',params0);
%[phat,pci]=mle(yvec,'logpdf',@(params)lhoods_david(yvec,xdata,params),'start',params0);

options=optimset('Display','off','MaxIter',10000,'TolX',10^-30,'TolFun',10^-30);
[theta,fval,exitflag,output,grad,hessian]=fminunc('lhoods_david',params0,options,ydata);
%exp(theta);
f=theta;
end