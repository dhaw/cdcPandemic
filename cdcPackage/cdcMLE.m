function f=cdcMLE(ydata)

params0=[3.1134    0.5924];

options=optimset('Display','off','MaxIter',10000,'TolX',10^-30,'TolFun',10^-30);
[theta,fval,exitflag,output,grad,hessian]=fminunc('cdcLhoods',params0,options,ydata);
f=theta;
end