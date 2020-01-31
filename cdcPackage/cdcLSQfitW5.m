function f=cdcLSQfitW5(data,thetaIn)
%{
Cc=[1.9200    0.4268    0.5260    0.2554    0.1665;
    1.7600    8.7522    2.2855    1.0876    1.2190;
    4.0700    4.5939    6.6160    4.5939    2.9494;
    0.9000    0.8885    1.6180    2.3847    1.6919;
    0.2300    0.2975    0.5712    0.8756    1.8930]';
Cvec=reshape(Cc,25,1)';
Cm=.9*Cvec;
Cp=1.1*Cvec;
%}
%
Cvec=thetaIn(3:27);
Cm=.9*Cvec;
Cp=1.1*Cvec;
%}
xdata=[1:36];%,21:31];%thetamle:[6:10,21:30];
ydata=data(xdata,1:end);%Match agesOut
fun=@(x,xdata)cdcPandemicSimulationW5(x,xdata,0,0,0,243);
%{
x0=thetaIn;
lb=[0,.5,0,0,1.1,.25];
ub=[1,1,1,150,1.3175,1.25];
%}
%{
x0=[thetaIn(1:3),thetaIn(end-2:end)]; 
%x0=[Cvec    -3.4997    1.2509    0.3657];
lb=[.2,.5,.5,0,1.1,.25];
ub=[1,1,1,140,1.3175,1.25];
%}
%
x0=[thetaIn(1:2),Cvec,thetaIn(end-2:end)];
lb=[.5,0,Cm,60,1.05,.25];
ub=[1,1,Cp,140,1.32,1.2];
%}
%options=optimset('MaxFunEvals',100000,'MaxIter',100000);
x=lsqcurvefit(fun,x0,xdata,ydata,lb,ub);%,options);
f=x;
end

function f=fun(x,xdata)
pandsim=cdcPandemicSimulationW5(x,xdata,0,0,0,243);
f=pandsim;
end