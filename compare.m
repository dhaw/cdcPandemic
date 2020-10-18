function [xsto, outsto, history, accept_rate,covmat]=compare(NNbar,xdata,ydataNX,thetac)%,xpriors,xbinEdges)%,x0
% F:          Function giving log-posterior density for a parameter set x
% x0:         Initial value of parameter set x
% n:          Number of iterations
% cov0:       Initial covariance matrix
% fac:        Scaling factor for covariance matrix. Set fac = 1 for default
% fixinds:    Elements of x that should be held fixed. Set fixinds = [] for full MCMC
% blockinds:  Number of 'epi parameters' (e.g. beta, X2) in x, if we want to vary epi and non-epi parameters as independent 'blocks'. Set blockinds = [] if we want a full covariance matrix
% displ:      Structure with display options. Set displ = true to show progress

%MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ);

age2mats=0;

cut=35;%Lower (included in fit)
ydataNX(xdata>cut,:)=[];
xdata(xdata>cut)=[];
%{
cut=20;%Lower (included in fit)
ydataNX(xdata<cut,:)=[];
xdata(xdata<cut)=[];
%}
n=10000;
sigma=1;
fixinds=[];
blockind=[];
displ=false;

%plim=[1,.0;
    %1,0;
    %1,0;
    plim=[1,0;
    1,0;
    1,0;
    240,0;
    2,1;%max,min
    2.4,.6]';
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
    Cvec1=thetac(2:26)';
    %Cvec2=thetac(27:51)';
    %plim=[[10*Cvec1,0*Cvec1];[10*Cvec2,0*Cvec2];plim']';
    plim=[[1.2*Cvec1,.8*Cvec1];plim']';
    x0=[Cvec1',thetac(end-5:end)];%+.2*(rand(1,length(thetac)-1)-.5).*thetac(2:end);
else
    Cvec=thetac(3:27)';
    plim=[[1.2*Cvec,.8*Cvec];plim']';
    %plim=[plim(:,1:3)';[1.1*Cvec,.9*Cvec];plim(:,end-4:end)']';
    %
    %x0=randic(plim);
    x0=thetac(3:end)+.2*(rand(1,length(thetac)-2)-.5).*thetac(3:end);
    %x0=prior;
    %x0=thetac(4:end);%.*(1+.02*rand(1,length(thetac)-3));
end

F=@(params)fcn(NNbar,params,xdata,ydataNX,plim);
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

function f=fcn(NNbar,params,xdata,ydataNX,plim)
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
© 2020 GitHub, Inc.