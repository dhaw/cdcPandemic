function f=pandemic1Dall3(params,xdata)%(R0,phi1,phi2,tlag,seednum,tswitch,closureFactor,betacModifier)
byAge=1;
%%
%Fixed parameters:
R0modifier=1;%/.775;
R0=1.4012*R0modifier;
tswitch=220;
betacModifier=.7737;%.9;
closureFactor=.7001;
monthShift=-.8013;
phi1=1;
phi2=0;
tlag=30;%Days
stoch=0;%Stochasticity - not included yet
pop=3e10;%Population
gamma=1/2.6;
seednum=3;
seasonality=1;
ftimes=1;
tclose=10^4;
%%
%Input parameters:
seednum=params(1);
%tswitch=params(2);
%betacModifier=params(2);
%phi1=params(2);
%phi2=params(3);
%closureFactor=params(3);
%R0=params(4);
%monthShift=params(5);
%%
n=1; nbar=4;
NNbar=[19169690;62121035;184015269;39570590];
NN=sum(NNbar);
NNrep=repmat(NN,4,1);
%Age mixing:
Cc=[27.57920413,8.051767033,4.975736133,0.850626995;9.165259795,43.43045174,8.195858852,2.158756533;5.941537452,5.863025518,14.20166331,5.533694466;0.600583289,0.807369258,1.39444674,7.848296781];
Co=Cc;
Cc(2,2)=closureFactor*Co(2,2);
%Calculate betas:
Sstart=repmat(NNbar,1,nbar);
Mj=NNbar; Mj(Mj==0)=1; Mjover=1./Mj;
Mjover=repmat(Mjover',nbar,1);
Dc=(Sstart.*Mjover).*Cc;
%Dc=Sstart*Cc/NN;
Gc=1/gamma*Dc;
d=eigs(Gc,1); R0c=max(d); betac=R0/R0c*betacModifier;
Do=(Sstart.*Mjover).*Co;
%Do=Sstart*Co/NN;
Go=1/gamma*Do;
d=eigs(Go,1); R0o=max(d); betao=R0/R0o;
%}
%%
%For simulation:
t0=0; tend=720;
NN0=NNrep; NN0(NNrep==0)=1;
Ni=repmat(NNrep,1,nbar); Nj=Ni';
Dc=Cc.*Mjover*NN;
Do=Co.*Mjover*NN;
zn=zeros(nbar,1);
y0=[NNbar;zn;zn];
seed=10^(-seednum);
%%
%Simulate:
    [tout,yout]=ode45(@(t,y)integr8all(t,y,betac,betao,gamma,nbar,NN0,Dc,Do,seasonality,phi1,phi2,seed,tlag,tswitch,tclose),[t0,tend],y0);
    Y=yout(:,1:nbar);
    Y=-diff(Y,1,1);
    tdiff=diff(tout);
    Y=Y./repmat(tdiff,1,nbar);
    tout=tout(2:end);
    if byAge==1
        NNdiv=repmat(NNbar',size(Y,1),1);
        Y=Y./NNdiv;
    else
        Y=sum(Y,2)/NN;
    end
    tmonth=1+ceil(tout/28+monthShift);
    tmonth(tmonth<1)=1;
    if byAge==1
        f1=accumarray(tmonth,Y(:,1));
        f2=accumarray(tmonth,Y(:,2));
        f3=accumarray(tmonth,Y(:,3));
        f4=accumarray(tmonth,Y(:,4));
        fall=ftimes*[f1,f2,f3,f4];
    else
        fall=ftimes*accumarray(tmonth,Y);
    end
    f=fall(xdata,:);
end
%%
function f=integr8all(t,y,betac,betao,gamma,nbar,NN0,Dc,Do,seasonality,phi1,phi2,seed,tlag,tswitch,tclose)
if seasonality==1
    phi=phi1-phi2*cos(pi*(t-tlag)/180);
else
    phi=1;
end
if t<tswitch || t>tclose
    XX=Dc;
    beta=betac;
else
    XX=Do;
    beta=betao;
end
S=y(1:nbar);
I=y(nbar+1:2*nbar);
seed1=seed.*S./NN0;
Sfoi=phi*(beta*S.*(XX*(I./NN0)+seed1));
Sdot=-Sfoi;
Idot=Sfoi-gamma*I;
Rdot=gamma*I;
f=[Sdot;Idot;Rdot];
end