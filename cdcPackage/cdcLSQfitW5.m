function f=cdcLSQfitW5(data)

xdata=1:35;
ydata=data(xdata,:);
fun=@(x,xdata)cdcPandemicSimulationW5(x,xdata,0,0,0,243);
%x0=[3,.58,1.4,.38];
x0=[3.5,.45,.2,1.34,.38];
lb=[3.5,.1,0,1.3,.25];%1.2,.25];
ub=[6,1,.5,1.7,.5];%1.7,1];
x=lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
f=x;
end

function f=fun(x,xdata)
f=cdcPandemicSimulationW5(x,xdata,0,0,0,243);
end