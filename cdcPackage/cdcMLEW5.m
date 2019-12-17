function f=cdcMLEW5(ydataNX)

%params0=[3.1134    0.5924];
params0=[3.0114    0.5778    1.4065    0.3662];

options=optimset('Display','off','MaxIter',10000,'TolX',10^-30,'TolFun',10^-30);
[theta,fval,exitflag,output,grad,hessian]=fminunc('cdcLhoodsW5',params0,options,ydataNX);
f=theta;
end