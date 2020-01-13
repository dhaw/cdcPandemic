function f=cdcMLEW5(ydataNX)

params0=[2.6741    0.4705    0.8489    0.0000    1.3119    0.3362];

options=optimset('Display','off','MaxIter',10000,'TolX',10^-30,'TolFun',10^-30);
%[theta,fval,exitflag,output,grad,hessian]=fminunc('cdcLhoodsW5',params0,options,ydataNX);
theta=fminunc('cdcLhoodsW5',params0,options,ydataNX);
f=theta;
end