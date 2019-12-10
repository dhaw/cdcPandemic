function f=cdcLSQfit(data)

xdata=4:7;%4:12 for 2 peaks, 4:7 for first only
ydata=data(xdata,:);
fun=@(x,xdata)cdcPandemicSimulation(x,xdata,0,0,0,243);
x0=[1.4,3.8];
x=lsqcurvefit(fun,x0,xdata,ydata);
f=x;
end

function f=fun(x,xdata)
f=cdcPandemicSimulation(x,xdata,0,0,0,243);
end