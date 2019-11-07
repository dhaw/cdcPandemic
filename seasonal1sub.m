function [f,g,D]=seasonal1sub(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,isdual,solvetype,numseed,eps,randic,tauend)
%ZfinalSizeAllMulti2
%isdual: 0=SM, 1=DM, 2=IM
%solvetype: 1=FSC, 2=ODE, 3=SCM
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
demog=1;
plotTau=0;
time=(1:tauend);
lt=length(time);
t0=0; tend=1000;
mu=1/80;%In ODE code
phi1=1; phi2=0;
NN0=NNrep; NN0(NNrep==0)=1;
Nages=NNbar./NN0;
alpha=1;%TSIR/sub-exp parameter
thetaGeom=1;
rate=eps;
%
%Theta:
Prow=[0,0,0,1-eps,eps];
%Prow=[0,0,0,0,1-eps,eps];
if min(Prow)<0
    error('Check drift parameters')
end
%%
%DO NOT TOUCH PARAMETERS:
lp=length(Prow);
Prow=Prow/sum(Prow);
%
vover=[1/5,1/14,1/45,1/16];%[1/5,1/14,1/45,0];
toAge=round(NNbar(1:n)*vover(1));
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
B=zeros(nbar,lp);
if randic==1
    if solvetype==3
        if ICdepAge==1
            b=zeros(nbar,lp);
            b(:,1)=NNbar;
            for nn=1:nbar
                if NNbar(nn)>0
                    r1=round((lp-1)*rand);%How many numbers of years non-zero
                    r1=randsample(2:lp,r1);%Which years non-zero
                    toMove=round(randFact*rand*NNbar(nn));
                    if toMove>0 && isempty(r1)==0
                        toHere=datasample(r1,toMove);
                        v=ones(size(toHere));
                        toAdd=accumarray(toHere(:),v(:));
                        b(nn,2:length(toAdd))=toAdd(2:end);
                        b(nn,1)=b(nn,1)-sum(toAdd);
                    end
                end 
            end
        else
                error('This bit of code isnt written yet')
        end
    else
        if ICdepAge==1
            b=zeros(nbar,lp);
            b(:,1)=NNbar;
            for nn=1:nbar
                if NNbar(nn)>0
                    r1=round((lp-1)*rand);%How many numbers of years non-zero
                    r1=randsample(2:lp,r1);%Which years non-zero
                    toMove=randFact*rand*NNbar(nn);
                    if toMove>0 && isempty(r1)==0
                        nnvec=rand(1,length(r1));
                        nnvec=nnvec/sum(nnvec)*toMove;
                        b(nn,r1)=nnvec;
                        b(nn,1)=b(nn,1)-toMove;
                    end
                end 
            end
        else
            error('This bit of code needs checking')
            %Need to rewrite in terms of B
            maxYears=lp-1;
            Prows=zeros(n,maxYears);
            numNonZero=floor((maxYears+1)*rand(n,1));%0-6 years - pick whow many are non zero
            for i=1:n
                rsi=randsample(maxYears,numNonZero(i));
                Prows(i,rsi)=1;
            end
            findNZ=find(Prows);
            lNZ=length(findNZ);
            randNZ=rand(lNZ,1);
            Prows(findNZ)=randNZ;
            rand2=rand(n,1); repRand=repmat(rand2,1,maxYears);
            sumP=sum(Prows,2); repP=repmat(sumP,1,maxYears);
            B=Prows./repP.*repRand;
            B(isinf(B)==1)=0; B(isnan(B)==1)=0;
            NNratio=repmat(NNbar./NNrep,1,maxYears);
            B=[zeros(nbar,1),repmat(B,4,1).*NNratio];
        end
    end
else
    b=zeros(nbar,lp);
    b(:,1)=NNbar;
end
S0=b(:,1);
Z0=sum(b(:,2:end),2);%./NNbar;%*(1-mu);
%%
zn=zeros(nbar,1);
A1=zeros(n,lt);
A2=A1;
seed=10^-(numseed);%*NNprob;
seedvec=zeros(nbar,1); seedvec(2*n+1:3*n)=seed*ones(n,1);
%seedvec=seed*ones(nbar,1);
thresh=0;%Remove from ODE solver
%%
%SIMULATE (UP TO ATTACK RATES):
for t=1:lt    
if solvetype==1
    %Final size:
    addbit=0;%seed;%tend*
    %If seed=0, no epidemic happens - because use ODE solver for initial condition
    IC=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,minNind,maxNind,D,Z0,beta,t,t0,tend,zn,phi1,phi2,seed,2,thresh,alpha,seedvec); IC=IC./NN0;%New IC here (if ever necessary)
    options=optimset('Display','off');
    funt=@(Zi)solveZi(Zi,Z0,beta,gamma,D,Nages,addbit);
    Zsol=fsolve(funt,IC,options);
else%solvetype=2/3
    Zsol=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,minNind,maxNind,D,S0,beta,t,t0,tend,zn,phi1,phi2,seed,solvetype,thresh,alpha,seedvec,plotTau);
    %Zsol=Zsol./NN0;
end
nu=Zsol-Z0;
A1(:,t)=sum(reshape(Zsol./NN0,n,na),2);%Prop immune for spatial cell (before antigenic drift)
A2(:,t)=sum(reshape(nu./NN0,n,na),2);%AR for spatial cell
%%
%IF THETA GEOMETRIC:
if thetaGeom==1
if solvetype==3
    error('This bit of code isnt written yet')
else
    if demog==1
        Z0=(1-rate)*Zsol;
        %Necessary for randic==1:
        btake1=Z0(1:n,:)*vover(1);
        btake2=Z0(n+1:2*n,:)*vover(2);
        btake3=Z0(2*n+1:3*n,:)*vover(3);
        btake4=Z0(3*n+1:end,:)*vover(4);
        Z0=Z0-[btake1;btake2;btake3;btake4]+[zeros(n,1);btake1;btake2;btake3];
        b(1:n,1)=b(1:n,1)+sum(btake4,2);
        S0=NNbar-Z0;
    else
        Z0=(1-rate)*Zsol;
        S0=NNbar-Z0;
    end
end
%if demog==1
%    error('This bit of code isnt written yet')
%end
else
%%
%ASSIGNING IMMUNITY:
if solvetype==3
    badd=zeros(nbar,lp);
    b(:,1)=b(:,1)-nu;
    for ii=1:lp-1%Short loop - check repeated use of binomial
        baddii=binornd(nu,Prow(ii));%Were sus to both, now sus to s2/immune to just s1
        badd(:,ii)=badd(:,ii)+baddii;
        nu=nu-baddii;
    end
    badd(:,end)=nu;
    %b=b+badd;
else
    b(:,1)=b(:,1)-nu;
    %b=b+repmat(nu,1,lp).*Prow;
end
%%
%DRIFT:
b=[b(:,1)+b(:,2),b(:,3:end),zeros(nbar,1)];
%%
if solvetype==3
    b=b+badd;
else
    b=b+repmat(nu,1,lp).*repmat(Prow,nbar,1);
end
%%
%AGEING:
if demog==1
    if solvetype==3
        b1=b(1:n,:);
        b2=b(n+1:2*n,:);
        b3=b(2*n+1:3*n,:);
        b4=b(3*n+1:end,:);
        btake1=zeros(n,lp);
        btake2=btake1;
        btake3=btake1;
        btake4=btake1;
        for nn=1:n
            movenn=toAge(nn);
            %
            btakenn=zeros(1,lp);
            bnn=b1(nn,:);
            findnn=find(bnn>0);
            popsnn=bnn(findnn);
            fromnn=randsample(1:NNbar(nn),movenn);
            [histVec,~]=histcounts(fromnn,[0,cumsum(popsnn)+.5]);
            btakenn(findnn)=histVec;
            btake1(nn,:)=btakenn;
            %
            btakenn=zeros(1,lp);
            bnn=b2(nn,:);
            findnn=find(bnn>0);
            popsnn=bnn(findnn);
            fromnn=randsample(1:NNbar(n+nn),movenn);
            [histVec,~]=histcounts(fromnn,[0,cumsum(popsnn)+.5]);
            btakenn(findnn)=histVec;
            btake2(nn,:)=btakenn;
            %
            btakenn=zeros(1,lp);
            bnn=b3(nn,:);
            findnn=find(bnn>0);
            popsnn=bnn(findnn);
            fromnn=randsample(1:NNbar(2*n+nn),movenn);
            [histVec,~]=histcounts(fromnn,[0,cumsum(popsnn)+.5]);
            btakenn(findnn)=histVec;
            btake3(nn,:)=btakenn;
            %
            btakenn=zeros(1,lp);
            bnn=b4(nn,:);
            findnn=find(bnn>0);
            popsnn=bnn(findnn);
            fromnn=randsample(1:NNbar(3*n+nn),movenn);
            [histVec,~]=histcounts(fromnn,[0,cumsum(popsnn)+.5]);
            btakenn(findnn)=histVec;
            btake4(nn,:)=btakenn;
        end
    else
        btake1=b(1:n,:)*vover(1);
        btake2=b(n+1:2*n,:)*vover(2);
        btake3=b(2*n+1:3*n,:)*vover(3);
        btake4=b(3*n+1:end,:)*vover(4);
    end
    b=b-[btake1;btake2;btake3;btake4]+[zeros(n,lp);btake1;btake2;btake3];
    b(1:n,1)=b(1:n,1)+sum(btake4,2);
    if min(min(min(b)))<0
        error('Shit the bed')
    end
    
end
S0=b(:,1);
Z0=sum(b(:,2:end),2);%*(1-mu);
end
end
%%
%}
f=A1;
g=A2;
end

function f=solveZi(Zi,Z0,beta,gamma,D,Nages,addbit)
lZi=log(Nages-Zi); lZi(Zi>1)=NaN;
f=lZi-log(Nages-Z0)+beta/gamma*D*(Zi-Z0)+addbit;
end
%%
function f=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,minNind,maxNind,D,S0,beta,tau,t0,tend,zn,phi1,phi2,seed,solvetype,thresh,alpha,seedvec,plotTau)
%icR=ic.*NNrep;
%y0=[NNbar-icR;zn;icR];
y0=[S0;zn;NNbar-S0];
%{
cond=sum(icR<thresh);
if cond==0
    seed=0;
end
%}
if solvetype==2
    [tout,yout]=ode23(@(t,y)integr8all(t,y,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,alpha,seedvec),[t0,tend],y0);
    %Incidence curve in here:
    %
    if tau==plotTau
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
    end
    %}
    rec0=yout(end,2*nbar+1:end)+yout(end,nbar+1:2*nbar);%******** Truncation - add infectious to rercovered?
    f=rec0';
elseif solvetype==3
    y1=y0;%+[-ic1;ic1;zn];
    f=stochSim(y1,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,tau,alpha);
%}
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