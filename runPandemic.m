function [f,g,D]=runPandemic(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,isdual,solvetype,numseed,eps,randic,tauend)
%ZfinalSizeAllMulti2
%isdual: 0=SM, 1=DM, 2=IM
%solvetype: 1=FSC, 2=ODE, 3=SCM
if solvetype==1
    error('Final size calculations not supported')
end
ICdepAge=1;
if isdual==0
    beta=betaS;
elseif isdual==1
    beta=betaD;
elseif isdual==2
    beta=betaI;
else
    beta=beta3;
end
randFact=1;
time=(1:tauend);
t0=0; tend=1000;
phi1=1; phi2=0;
NN0=NNrep; NN0(NNrep==0)=1;
%alpha=1;%TSIR/sub-exp parameter
%vover=[1/5,1/14,1/45,1/16];%[1/5,1/14,1/45,0];
%toAge=round(NNbar(1:n)*vover(1));
%%
Ni=repmat(NNrep,1,nbar); Nj=Ni';
Niover=1./Ni; Niover(Ni==0)=1;
Mj=(Kbar')*NNbar;%C in denom?? .*Cbar %(Kbar')*
Mj(Mj==0)=1;
Mjover=1./Mj;
Mjover=repmat(Mjover',nbar,1);
if isdual==0
    D=Kbar.*Mjover.*Cbar.*Nj;
elseif isdual==1
    D=(K1.*Mjover)*((Kbar').*Cbar);
    D=D.*Nj;
elseif isdual==2
    D=Kbar'.*Niover.*Cbar.*Nj;
elseif isdual==3
    D=Kbar.*Cbar;
end
%%
%INITIAL CONDITION:
if randic==1
    if solvetype==3
        if ICdepAge==1
            rec0=round(randFact*rand*NNbar(nn));
        else
            error('This bit of code isnt written yet')
        end
    else
        if ICdepAge==1
            rec0=randFact*rand*NNbar(nn);
        else
            error('This bit of code isnt written yet')
        end
    end
else
    rec0=zeros(nbar,lp);
end
S0=NNbar-rec0;
%%
zn=zeros(nbar,1);
seed=10^-(numseed);%*NNprob;
seedvec=zeros(nbar,1); seedvec(2*n+1:3*n)=seed*ones(n,1);
%seedvec=seed*ones(nbar,1);
thresh=0;%Remove from ODE solver
%%
%SIMULATE (UP TO ATTACK RATES):   
Zsol=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,minNind,maxNind,D,S0,beta,t,t0,tend,zn,phi1,phi2,seed,solvetype,thresh,alpha,seedvec,plotTau);
%%
f=Zsol;
end

function f=solveZi(Zi,Z0,beta,gamma,D,Nages,addbit)
lZi=log(Nages-Zi); lZi(Zi>1)=NaN;
f=lZi-log(Nages-Z0)+beta/gamma*D*(Zi-Z0)+addbit;
end
%%
function f=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,minNind,maxNind,D,S0,beta,tau,t0,tend,zn,phi1,phi2,seed,solvetype,thresh,alpha,seedvec,plotTau)
y0=[S0;zn;NNbar-S0];
if solvetype==2
    [tout,yout]=ode23(@(t,y)integr8all(t,y,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,alpha,seedvec),[t0,tend],y0);
    %Incidence curve in here:
    %
    figure
    fs=12; lw=2;
    Y=yout(:,nbar+1:2*nbar); %Z=yout(:,2*nbar+1:3*nbar);
    %Y=sum(Y,2);
    %Y1=Y(:,minNind)+Y(:,minNind+n)+Y(:,minNind+2*n)+Y(:,minNind+3*n); Y1=Y1/NN(minNind);
    %Y2=Y(:,maxNind)+Y(:,maxNind+n)+Y(:,maxNind+2*n)+Y(:,maxNind+3*n); Y2=Y2/NN(maxNind);
    Ysum=sum(Y,2);
    Yall=Y(:,1:n)+Y(:,n+1:2*n)+Y(:,2*n+1:3*n)+Y(:,3*n+1:end);
    Yall=abs(Yall);%****cheat ;)
    %
    %Unlogged plots:
    hold on
    %plot(tout,Y1,'--','linewidth',lw,'color',[.165,.31,.431]);%[.165,.31,.431][.447,.553,.647]
    %plot(tout,Y2,'-','linewidth',lw,'color',[.165,.31,.431]);
    plot(tout,Ysum,'k','linewidth',lw);
    plot(tout,Yall);
    %}
    %{
    %Logged plots:
    semilogy(tout,Ysum,'k','linewidth',lw);
    hold on
    semilogy(tout,Yall);
    %}
    xlabel('Time (days)','FontSize',fs);
    ylabel('Prevalence','FontSize',fs);
    set(gca,'FontSize',fs);
    maxY=max(Ysum);
    %axis([0,tend,0,maxY])
    axis ([0,tend,0,maxY])
    %legend('Min','Max','location','NE')
    grid on
    grid minor
    hold off
    %}
    rec1=yout(end,2*nbar+1:end)+yout(end,nbar+1:2*nbar);%******** Truncation - add infectious to rercovered?
    f=rec1';
elseif solvetype==3
    f=stochSim(y0,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,tau,alpha);
else
    error('Invalid solvetype')
end
end
%%
function f=integr8all(t,y,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,alpha,seedvec)
mu=1/1800;%5*360=1800
%phi=phi1-phi2*cos(pi*t/180);
phi=1;%+sin(pi*t/180);
S=y(1:nbar);
I=y(nbar+1:2*nbar);
R=y(2*nbar+1:end);
%seed1=seed*S./NN0;%*exp(-t);
seed1=seedvec.*S./NN0;%*exp(-t);
Sfoi=phi*(beta*S.*(D*(I./NN0))+seed1);%seed1);
Sdot=-Sfoi;%+mu*R;
Idot=Sfoi-gamma*I;
Rdot=gamma*I;%-mu*R;
f=[Sdot;Idot;Rdot];
end
%%
function f=stochSim(y,beta,gamma,n,nbar,NN,N0,D,seed,phi1,phi2,tau,alpha)
%Feed in mu if required
factor=6;
tend=360*factor; beta=beta/factor; gamma=gamma/factor;
Vec=zeros(nbar,tend);
%
S=y(1:nbar);
I=y(nbar+1:2*nbar);
R=y(2*nbar+1:end);

%Test seed:
%[Nmax,maxInd]=max(NN);
%ii=4*(maxInd-1)+3;
%seed=0; I(ii)=10; S(ii)=S(ii)-10;
%
%I=floor(N0/10000); S=S-I;

%Different from here:
i=1;
threshold=30;%Number of time-steps for necessary simulation/seed
while i<tend && (i<threshold || sum(I)>0)%At least 30 time-steps
phi=1;%phi1-phi2*cos(pi*i*f1/180);%Just seed here
Sout=1-exp(-phi*(beta*(D*(I./N0).^alpha)+seed*heaviside(threshold-i)));%+mu*R;.^alpha
Sout(Sout>1)=1;
Sout=binornd(S,Sout); Sout(S==0)=0;
S=S-Sout; I=I+Sout;
%
Iout=1-exp(-gamma);
Iout=binornd(I,Iout);
I=I-Iout; R=R+Iout;
Vec(:,i)=I;
i=i+1;
end
%{
Vsum=sum(Vec,1);%/sum(NN);
maxV=max(Vsum);
fs=12; lw=2;
Vtot=Vsum(factor:factor:end);
Vall=Vec(1:n,factor:factor:end)+Vec(n+1:2*n,factor:factor:end)+Vec(2*n+1:3*n,factor:factor:end)+Vec(3*n+1:end,factor:factor:end);
%
%Unlogged plots:
figure
plot(1:tend/factor,Vtot,'-','linewidth',2,'color',[.447,.553,.647]);
plot(1:tend/factor,Vall);
axis([0,tend/factor,0,maxV]);
xlabel('Time (days)'); ylabel('Prevalence'); set(gca,'fontsize',fs)
grid on
grid minor
%}
%
%{
%Logged plots:
figure
semilogy(1:tend/factor,Vtot,'k','linewidth',lw)
hold on
semilogy(1:tend/factor,Vall);
axis([0,tend/factor,0,maxV]);
xlabel('Time (days)'); ylabel('Prevalence'); set(gca,'fontsize',fs)
grid on
grid minor
%}
f=R;
if sum(isnan(R))>0
    fuck=1;
    print('Somenting is NaN')
end
end