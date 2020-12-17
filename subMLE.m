function f=subMLE(NNbar,xdata,ydataNX,thetaIn,vaxparams)
age2mats=0;
%{
Cc=[1.9200    0.4268    0.5260    0.2554    0.1665;
    1.7600    8.7522    2.2855    1.0876    1.2190;
    4.0700    4.5939    6.6160    4.5939    2.9494;
    0.9000    0.8885    1.6180    2.3847    1.6919;
    0.2300    0.2975    0.5712    0.8756    1.8930]';
Cvec=reshape(Cc,25,1);
%}
nbar=5;

params0=thetaIn;
%%
%Bounds:
%
plim=[.3,0;%1,.4;
    31,-31;%.2,0;%35,-5;%60,-30;%1,0;%XX
    .4,0;
    %plim=[.5,0;
    .605,.495;%.66,.44;%.605,.495;
    .6,.4;%.55,.45;
    100,-100;%30,-5;%180,-60%XX
    2,1;%1.595,1.305;%1.595,1.305;%max,min
    2/3,1/4]';%1/2.7
%}
%{
plim=[1,0;
    60,-60;
    .5,0;
    .605,.495;
    .55,.45;
    60,-60;%240,0;
    3,1;%max,min
    1,1/3.3]';%1/2.7
%}
%{
plim=[1,0;
    .5,0;
    1,0;
    1,0;
    1,0;
    100,-60;%240,0;
    2,1;%max,min
    1,1/3.3]';%1/2.7
%}
if age2mats==1
    Cc1=reshape(thetaIn(3:3+nbar^2-1),nbar^2,1);
    Cc2=reshape(thetaIn(3+nbar^2:3+2*nbar^2-1),nbar^2,1);
    %Cc2=reshape(thetaIn(2+nbar^2:2+2*nbar^2-1),nbar^2,1);
    plim=[plim(:,1:2)';[1.2*Cc1,.8*Cc1];[1.2*Cc2,.8*Cc2];plim(:,end-5:end)'];
    %plim=[plim(:,2)';[1.1*Cc1,.9*Cc1];plim(:,end-5:end)'];
    %plim=[plim(:,1:2)';[1.1*Cc1,.9*Cc1];plim(:,end-5:end)'];%plim(:,1:2)';
else
    Cvec=reshape(thetaIn(3:3+nbar^2-1),nbar^2,1);
    plim=[plim(:,1:2)';[1.2*Cvec,.8*Cvec];plim(:,end-5:end)'];
    %Cvec=reshape(thetaIn(1:nbar^2),nbar^2,1);
    %plim=[plim(:,1:2)';[1.1*Cvec,.9*Cvec];plim(:,end-5:end)'];
    %plim=[plim(:,2)';[1.2*Cvec,.8*Cvec];plim(:,end-5:end)'];%W2 only
    %plim=[1.2*thetaIn',.8*thetaIn'];
    %plim(1,:)=[1,0];
    %plim(2,:)=[1,0];
    %plim(28,1)=.2;
end
ub=plim(:,1);
lb=plim(:,2);
%%
options=optimset('MaxFunEvals',100000,'algorithm','sqp');
fun=@(params)subLhoods(NNbar,params,xdata,ydataNX,vaxparams);
%theta=fmincon(fun,params0,[],[],[],[],lb',ub',[],options);%options);

rng default%For reproducibility
options=optimoptions(@fmincon,'MaxFunctionEvaluations',100000,'algorithm','interior-point','UseParallel',true);
problem=createOptimProblem('fmincon','x0',params0,'objective',fun,'lb',lb,'ub',ub,'options',options);
ms=MultiStart;
[xoptim,~]=run(ms,problem,20);
%save('optimOut2.mat','xoptim','yoptim','exitflag')

%theta=fminunc(fun,params0,options);%options);
f=xoptim;
end
%{
nbar=size(ydataNX,2);
xdata=1:36;
lx=length(xdata);
ydataNX=ydataNX(xdata,:);
hosp=reshape(ydataNX,nbar*lx,1);
[theta,~]=mle(hosp,'pdf',@(params,hosp)lhoodFct(params,hosp),'start',params0);
%}