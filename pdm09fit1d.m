function f=pdm09fit1d(data)

xdata=5:12;%Not from 1 - in case truncate 1st month, also to avoid previous season
ydata=data(xdata,:);
fun=@(x,xdata)pandemic1Dall(x,xdata,0,0);
x0=[2.4,220,.77,.7,1.4,-.8];
x=lsqcurvefit(fun,x0,xdata,ydata);
f=x;