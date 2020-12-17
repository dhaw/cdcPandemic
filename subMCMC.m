function [xsto, outsto, history, accept_rate,covmat]=subMCMC(NNbar,xdata,ydataNX,thetac,vaxparams)%,xpriors,xbinEdges)%,x0
% F:          Function giving log-posterior density for a parameter set x
% x0:         Initial value of parameter set x
% n:          Number of iterations
% cov0:       Initial covariance matrix
% fac:        Scaling factor for covariance matrix. Set fac = 1 for default
% fixinds:    Elements of x that should be held fixed. Set fixinds = [] for full MCMC
% blockinds:  Number of 'epi parameters' (e.g. beta, X2) in x, if we want to vary epi and non-epi parameters as independent 'blocks'. Set blockinds = [] if we want a full covariance matrix
% displ:      Structure with display options. Set displ = true to show progress

%MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ);

C1=[6.3137    1.9526    1.1026    0.2633    0.0770    4.9328   15.7934    2.4278    0.5151    0.1393    5.6393    5.6127    8.5182  2.8388    0.1215    1.6730    1.6056    2.7618    1.9235    0.1121    0.8147    1.5306    0.4504    0.3602    0.1609]';
C2=[7.5103    2.0601    1.1629    0.2796    0.0773    5.2350   20.8923    2.6155    0.5974    0.1406    6.0040    6.7120    8.8204  2.9078    0.1250    1.9520    2.5667    2.8596    1.9930    0.1165    0.8168    1.5605    0.4641    0.3906    0.1830]';
C=[1.9200    0.4268    0.5260    0.2554    0.1665;
    1.7600    8.7522    2.2855    1.0876    1.2190;
    4.0700    4.5939    6.6160    4.5939    2.9494;
    0.9000    0.8885    1.6180    2.3847    1.6919;
    0.2300    0.2975    0.5712    0.8756    1.8930]';
C=reshape(C,25,1);
age2mats=0;

cut=34;%34;%Cut from here +1
ydataNX(xdata>cut,:)=[];
xdata(xdata>cut)=[];
%{
cut=20;%Lower (included in fit)
ydataNX(xdata<cut,:)=[];
xdata(xdata<cut)=[];
%}
n=100000;
sigma=1;
fixinds=[];
blockind=[];
displ=false;

%
plim=[.3,0;%1,.4;
    31,-31;%.2,0;%35,-5;%60,-30;%1,0;%XX
    .4,0;
    %plim=[.5,0;
    .605,.495;%.66,.44;%.605,.495;
    .6,.4;%.55,.45;
    100,-100;%30,-5;%180,-60%XX
    1.595,1.305;%1.595,1.305;%max,min
    2/3,1/4]';%1/2.7
%}
%{
plim=[.4,0;%1,.4;
    100,-100;%.2,0;%35,-5;%60,-30;%1,0;%XX
    .5,0;
    %plim=[.5,0;
    .66,.44;%.605,.495;
    .6,.4;%.55,.45;
    100,-100;%30,-5;%180,-60%XX
    1.7,1.3;%1.595,1.305;%max,min
    1,1/4]';%1/2.7
%}
%{
plim=[1,.4;
    .2,0;%35,-5;%60,-30;%1,0;%XX
    .5,0;
    %plim=[.5,0;
    1,0%.8,.44;%.66,.44;%.605,.495;%XX
    .6,0;%.6,.4;%.55,.45;%XX
    100,0;%30,-5;%180,-60%XX
    2,1;%1.595,1.305;%max,min
    2,1/4]';%1/2.7
%}
%{
plim=[1,0;
    100,-100;
    1,0;
    1,0;
    1,0;
    120,00;%240,0;
    3,1;%max,min
    3,1/6]';%1/2.7
%}
%plim=[2*thetac;0*thetac];
%%
%{
Cc=[1.9200    0.4268    0.5260    0.2554    0.1665;
    1.7600    8.7522    2.2855    1.0876    1.2190;
    4.0700    4.5939    6.6160    4.5939    2.9494;
    0.9000    0.8885    1.6180    2.3847    1.6919;
    0.2300    0.2975    0.5712    0.8756    1.8930]';
Cvec=reshape(Cc,25,1);
%}
if age2mats==1
    Cvec1=thetac(3:27)';
    Cvec2=thetac(28:52)';
    plim=[plim(:,1:2)';[1.2*C1,.8*C1];[1.2*C2,.8*C2];plim(:,end-5:end)']';
    %plim=[[1.2*Cvec1,.8*Cvec1];plim']';
    x0=[thetac(1:2),Cvec1',Cvec2',thetac(end-5:end)];%+.2*(rand(1,length(thetac)-1)-.5).*thetac(2:end);
else
    Cvec=thetac(3:27)';
    %{
    %From C - set cut=34
    plim=[[1.1*Cvec,.9*Cvec];plim']';
    x0=thetac(3:end);%+.2*(rand(1,length(thetac)-2)-.5).*thetac(3:end);
    %}
    plim=[plim(:,1:2)';[1.2*Cvec,.8*Cvec];plim(:,end-5:end)']';
    %plim=[1.2*thetac;.8*thetac];
    x0=thetac;%+.2*(rand(1,length(thetac)-2)-.5).*thetac(3:end);
    %x0=prior;
    %x0=thetac(4:end);%.*(1+.02*rand(1,length(thetac)-3));
    %plim=kron([1.2;.8],thetac(3:end));
    %{
    %Including x22/phi2 etc - set cut=52
    plim=[plim(:,1:2)';[1.1*Cvec,.9*Cvec];plim(:,3:end)']';
    x0=thetac;%+.2*(rand(1,length(thetac))-.5).*thetac;
    %}
end

F=@(params)fcn(NNbar,params,xdata,ydataNX,plim,vaxparams);
%F=@(params)fcnPrior(NNbar,params,xdata,ydataNX,plim,xpriors,xbinEdges);
[xsto, outsto, history, accept_rate,covmat] = MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ);
end

function f=unif(x,plim)
val=1./(plim(1,:)-plim(2,:));
in=(x-plim(1,:)).*(x-plim(2,:));
in(in>0)=0;
in(in<0)=1;
f=val.*in;
end

function f=randic(plim)
f=plim(1,:)+rand(1,size(plim,2)).*(plim(2,:)-plim(1,:));
end

function f=fcn(NNbar,params,xdata,ydataNX,plim,vaxparams)
%if params(end)>.1
p1=params((plim(2,:)-params)>0);
p2=params((plim(1,:)-params)<0);
if isempty(p1)==1 && isempty(p2)==1
    f=-subLhoods(NNbar,params,xdata,ydataNX,vaxparams)+sum(log(unif(params,plim)));
    %+log(unif(params(1),plim1))+log(unif(params(2),plim2));%+log(unif(params(3),plim3));%+log(unif(params(4),plim4)));
else
    f=-inf;
end
end

function f=fcnPrior(NNbar,params,xdata,ydataNX,plim,xpriors,xbinEdges)%Fix
p1=params((plim(2,:)-params)>0);
p2=params((plim(1,:)-params)<0);
lp=length(params);%Feed in
if isempty(p1)==1 && isempty(p2)==1
    %Gamma:
    f=-subLhoods(NNbar,params,xdata,ydataNX)+sum(log(gampdf(params',xpriors(1,:)',xpriors(2,:)')));
    %{
    %Histogram:
    p0=0;
    for i=1:lp
        index=discretize(params(i),xbinEdges(i,:));
        p0=p0+log(xpriors(i,index));
    end
    f=-subLhoods(NNbar,params,xdata,ydataNX)+p0;
    %}
else
    f=0;
end
end