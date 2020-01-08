function f=cdcMLEW5(ydataNX)

params0=[3.4449    0.3646    0.6082    0.0263    1.8760    0.3932];

options=optimset('Display','off','MaxIter',10000,'TolX',10^-30,'TolFun',10^-30);
%[theta,fval,exitflag,output,grad,hessian]=fminunc('cdcLhoodsW5',params0,options,ydataNX);
theta=fminunc('cdcLhoodsW5',params0,options,ydataNX);
f=theta;
end