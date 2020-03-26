function f=subMLE(NNbar,xdata,ydataNX,thetaIn)
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
plim=[1,0;
    1,0;
    1,0;
    1,0;
    1,0;
    240,0;
    2,1;%max,min
    2.4,.6]';
if age2mats==1
    Cc1=reshape(thetaIn(2:2+nbar^2-1),nbar^2,1);
    Cc2=reshape(thetaIn(2+nbar^2:2+2*nbar^2-1),nbar^2,1);
    plim=[plim(:,2)';[1.2*Cc1,.8*Cc1];[1.2*Cc2,.8*Cc2];plim(:,end-5:end)'];
else
    Cvec=reshape(thetaIn(3:3+nbar^2-1),nbar^2,1);
    plim=[plim(:,1:2)';[100*Cvec,0*Cvec];plim(:,end-5:end)'];
    %plim=[plim(:,2)';[100*Cvec,.01*Cvec];plim(:,end-5:end)'];%W2 only
end
ub=plim(:,1);
lb=plim(:,2);
%%
options=optimset('MaxFunEvals',100000);
fun=@(params)subLhoods(NNbar,params,xdata,ydataNX);
theta=fmincon(fun,params0,[],[],[],[],lb',ub',[],options);%options);
f=theta;
end
%{
nbar=size(ydataNX,2);
xdata=1:36;
lx=length(xdata);
ydataNX=ydataNX(xdata,:);
hosp=reshape(ydataNX,nbar*lx,1);
[theta,~]=mle(hosp,'pdf',@(params,hosp)lhoodFct(params,hosp),'start',params0);
%}