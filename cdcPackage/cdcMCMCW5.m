function [xsto, outsto, history, accept_rate]=cdcMCMCW5(ydataNX,thetac)
% F:          Function giving log-posterior density for a parameter set x
% x0:         Initial value of parameter set x
% n:          Number of iterations
% cov0:       Initial covariance matrix
% fac:        Scaling factor for covariance matrix. Set fac = 1 for default
% fixinds:    Elements of x that should be held fixed. Set fixinds = [] for full MCMC
% blockinds:  Number of 'epi parameters' (e.g. beta, X2) in x, if we want to vary epi and non-epi parameters as independent 'blocks'. Set blockinds = [] if we want a full covariance matrix
% displ:      Structure with display options. Set displ = true to show progress

%MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ);

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
    140,60;
    1.7,1.1;%max,min
    2.2,.25]';
%{
Cc=[1.9200    0.4268    0.5260    0.2554    0.1665;
    1.7600    8.7522    2.2855    1.0876    1.2190;
    4.0700    4.5939    6.6160    4.5939    2.9494;
    0.9000    0.8885    1.6180    2.3847    1.6919;
    0.2300    0.2975    0.5712    0.8756    1.8930]';
Cvec=reshape(Cc,25,1);
%}
Cvec=thetac(4:28)';
plim=[[1.1*Cvec,.9*Cvec];plim']';
%plim=[plim(:,1:3)';[1.1*Cvec,.9*Cvec];plim(:,end-4:end)']';
%}
    
%x0=randic(plim);
x0=thetac(4:end);

F=@(params)fcn(params,ydataNX,plim);
[xsto, outsto, history, accept_rate] = MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ);
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

function f=fcn(params,ydataNX,plim)
f=-cdcLhoodsTest(params,ydataNX)+sum(log(unif(params,plim)));
%+log(unif(params(1),plim1))+log(unif(params(2),plim2));%+log(unif(params(3),plim3));%+log(unif(params(4),plim4)));
end