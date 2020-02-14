function f=cdcMLEW5(ydataNX,theta)

Cc=[1.9200    0.4268    0.5260    0.2554    0.1665;
    1.7600    8.7522    2.2855    1.0876    1.2190;
    4.0700    4.5939    6.6160    4.5939    2.9494;
    0.9000    0.8885    1.6180    2.3847    1.6919;
    0.2300    0.2975    0.5712    0.8756    1.8930];
Cvec=reshape(Cc,25,1)';

params0=theta;%[-3.4997    0.5707    0.9578    1.2509    0.3657];

options=optimset('Display','off');%,'MaxIter',10000,'TolX',10^-30,'TolFun',10^-30);
%[theta,fval,exitflag,output,grad,hessian]=fminunc('cdcLhoodsW5',params0,options,ydataNX);
theta=fminunc('cdcLhoodsW5',params0,options,ydataNX);
f=theta;
end