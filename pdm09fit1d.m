function f=pdm09fit1d(data)

xdata=4:7;%Not from 1 - in case truncate 1st month, also to avoid previous season
ydata=data(xdata,:);
fun=@(x,xdata)pandemic1DallV(x,xdata,0,0,0);
x0=[3.2740    0.9284    0.7383];
x=lsqcurvefit(fun,x0,xdata,ydata);
f=x;
end

function f=fun(x,xdata)
f=pandemic1Dall2(x,xdata,0,0);
end

%ydata=reshape(ydata,4*length(xdata),1);
%lb=[2,200,0,0,1,-1]; ub=[8,250,1,1,2,0];