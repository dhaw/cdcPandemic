function f=subLsqFit(NNbar,xdata,ydata,thetaIn)
age2mats=0;
nbar=5;

threshl=17;%17;%Lower (included in fit)
threshu=52;%Upper

ydata(xdata<threshl,:)=[];
xdata(xdata<threshl)=[];
ydata(xdata>threshu,:)=[];
xdata(xdata>threshu)=[];

plim=[1,0;
    1,0;
    .5,0;
    1,0;
    1,0;
    220,120;%240,0;
    2,1;%max,min
    2.4,.6]';
if age2mats==1
    Cc1=reshape(thetaIn(2:2+nbar^2-1),nbar^2,1);
    Cc2=reshape(thetaIn(2+nbar^2:2+2*nbar^2-1),nbar^2,1);
    plim=[plim(:,2)';[1.2*Cc1,.8*Cc1];[1.2*Cc2,.8*Cc2];plim(:,end-5:end)'];
else
    Cvec=reshape(thetaIn(3:3+nbar^2-1),nbar^2,1);
    plim=[plim(:,1:2)';[2*Cvec,0*Cvec];plim(:,end-5:end)'];
end

x0=thetaIn;
ub=plim(:,1);
lb=plim(:,2);

options=optimset('MaxFunEvals',100000,'MaxIter',100000);
fun=@(params,xdata)subPandemicSimulation(NNbar,params,xdata,0,0,0,243);
theta=lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options);
f=theta;
end

function f=funx(x,xdata)
pandsim=cdcPandemicSimulationW5(x,xdata,0,0,0,243);
f=pandsim(:,1:2);
end