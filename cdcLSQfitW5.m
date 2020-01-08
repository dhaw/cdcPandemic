function f=cdcLSQfitW5(data)

Cc=[1.9200    0.4268    0.4600    0.4434    0.1665;
        1.7600    8.7522    2.4600    1.8630    1.2190;
        0.4500    0.4911    2.5900    0.8132    0.3412;
        4.5200    4.9913    5.6500    6.9995    4.3001;
        0.2300    0.2975    0.4300    0.6903    1.8930]';
Cparam=reshape(Cc,25,1)';
Cm=.9*Cparam;
Cp=1.1*Cparam;

xdata=1:35;
ydata=data(xdata,:);
fun=@(x,xdata)cdcPandemicSimulationW5(x,xdata,0,0,0,243);
%x0=[3,.58,1.4,.38];
x0=[3.5,.4,.6,Cparam,1.46/.775,.38];%[3.5,.45,.2,1.34,.38]; [7,.8,.2];
lb=[2,0,0,Cm,1,.25];%1.2,.25];
ub=[10,1,1,Cp,3,.5];%,2.5,.5];%1.7,1];
options=optimset('MaxFunEvals',10000,'MaxIter',10000);
x=lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options);
f=x;
end

function f=fun(x,xdata)
f=cdcPandemicSimulationW5(x,xdata,0,0,0,243);
end