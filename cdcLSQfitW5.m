function f=cdcLSQfitW5(data,thetaIn)
%{
Cc=[1.9200    0.4268    0.5260    0.2554    0.1665;
    1.7600    8.7522    2.2855    1.0876    1.2190;
    4.0700    4.5939    6.6160    4.5939    2.9494;
    0.9000    0.8885    1.6180    2.3847    1.6919;
    0.2300    0.2975    0.5712    0.8756    1.8930];
Cvec=reshape(Cc,25,1)';
%}
%{
%Cvec=thetaIn(3:27);
Cm=.9*Cvec;
Cp=1.1*Cvec;
%}
xdata=[5:9,21:28];
ydata=data(xdata,:);
fun=@(x,xdata)cdcPandemicSimulationW5(x,xdata,0,0,0,243);
%
%x0=[thetaIn(1:3),Cvec,thetaIn(end-2:end)];
x0=thetaIn;
lb=[.2,.5,60,1,.25];
ub=[1,1,140,2,1.25];
%}
%{
x0=[Cvec    -3.4997    1.2509    0.3657];
lb=[Cm,-30,1,.25];
ub=[Cp,60,3,.5];
%}
%{
x0=[0,.6,.85,Cparam,1.2,.35];
lb=[-30,0,0,Cm,1,.25];
ub=[60,1,1,Cp,3,.5];
%}
%options=optimset('MaxFunEvals',1000000,'MaxIter',1000000);
x=lsqcurvefit(fun,x0,xdata,ydata,lb,ub);%,options);
f=x;
end

function f=fun(x,xdata)
pandsim=cdcPandemicSimulationW5(x,xdata,0,0,0,243);
f=pandsim(:,1:2);
end