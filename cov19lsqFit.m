function f=cov19lsqFit(data,thetaIn)
xdata=[5:9,21:28];
ydata=data(xdata,:);
fun=@(x,xdata)cdcPandemicSimulationW5(x,xdata,0,0,0,243);
x0=thetaIn;
lb=[.2,.5,60,1,.25];
ub=[1,1,140,2,1.25];

%options=optimset('MaxFunEvals',1000000,'MaxIter',1000000);
x=lsqcurvefit(fun,x0,xdata,ydata,lb,ub);%,options);
f=x;
end

function f=fun(x,xdata)
pandsim=cdcPandemicSimulationW5(x,xdata,0,0,0,243);
f=pandsim(:,1:2);
end