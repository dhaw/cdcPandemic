function [xsto, outsto, history, accept_rate,covmat]=subMCMCxphi2(NNbar,xdata,ydataNX,thetac)%,xpriors,xbinEdges)%,x0age2mats=0;
%mcmc=0 in simulation

cut=17;%Lower (included in fit)
ydataNX(xdata>cut,:)=[];
xdata(xdata>cut)=[];

n=100000;
sigma=1;
fixinds=[];
blockind=[];
displ=false;

%plim=[.85,.75;%max, min
    %.04,0]';

plim=[.81,.79;%max, min
    .001,0;
    1,0;
    1,0;
    1,0;
    240,0;
    2,1;%max,min
    2.4,.6]';
Cvec1=thetac(3:27)';
plim=[plim(:,1:2)';[1.1*Cvec1,.9*Cvec1];plim(:,3:end)']';

%x0=randic(plim);
x0=thetac;%(1:2);%+.1*(rand(1,2)-.5).*thetac(1:2);
%x0=prior;
%x0=thetac(4:end);%.*(1+.02*rand(1,length(thetac)-3));

params2=thetac(3:end);

%F=@(params1)fcn(NNbar,params1,params2,xdata,ydataNX,plim);
F=@(params)fcn2(NNbar,params,xdata,ydataNX,plim);
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

function f=fcn(NNbar,params1,params2,xdata,ydataNX,plim)
%if params(end)>.1
p1=params1((plim(2,:)-params1)>0);
p2=params1((plim(1,:)-params1)<0);
if isempty(p1)==1 && isempty(p2)==1
    f=-subLhoods(NNbar,[params1,params2],xdata,ydataNX)+sum(log(unif(params1,plim)));
    %+log(unif(params(1),plim1))+log(unif(params(2),plim2));%+log(unif(params(3),plim3));%+log(unif(params(4),plim4)));
else
    f=-inf;
end
end

function f=fcn2(NNbar,params,xdata,ydataNX,plim)
%if params(end)>.1
p1=params((plim(2,:)-params)>0);
p2=params((plim(1,:)-params)<0);
if isempty(p1)==1 && isempty(p2)==1
    f=-subLhoods(NNbar,params,xdata,ydataNX)+sum(log(unif(params,plim)));
    %+log(unif(params(1),plim1))+log(unif(params(2),plim2));%+log(unif(params(3),plim3));%+log(unif(params(4),plim4)));
else
    f=-inf;
end
end