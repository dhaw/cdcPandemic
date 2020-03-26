function f=cdcLSQfit(data)

xdata=4:12;%4:12 for 2 peaks, 4:7 for first only
ydata=data(xdata,:);
fun=@(x,xdata)cdcPandemicSimulation(x,xdata,0,0,0,243);
%x0=[3,.58,1.4,.38];
x0=[3,.7,.98];
lb=[2.5,.2,.9];
ub=[8,1,1];
x=lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
f=x;
end

function f=fun(x,xdata)
f=cdcPandemicSimulation(x,xdata,0,0,0,243);
end
%{
xdata=4:12;%4:12 for 2 peaks, 4:7 for first only
ydata=data(xdata,:);

NNbar=[19169690;62121035;184015269;39570590];
NNrep=repmat(NNbar',length(xdata),1);
ydata=ydata.*NNrep;

fun=@(x,xdata)cdcPandemicSimulation(x,xdata,0,0,0,243);
x0=[3,.7,.95];
x=lsqcurvefit(fun,x0,xdata,ydata);
f=x;
end

function f=fun(x,xdata,NNrep)
pand=cdcPandemicSimulation(x,xdata,0,0,0,243);
f=sum(pand.*NNrep,2);
end
%}