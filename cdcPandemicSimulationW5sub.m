function [f,g,z2]=cdcPandemicSimulationW5sub(params,xdata,plotComp,plotEpis,ydata,tswitch)%(R0,phi1,phi2,tlag,seednum,tswitch,closureFactor,betacModifier)
%V: antiViral treatment
%plotComp: plot comparison (with data)
%plotEpis: plot incidence (non-aggregated)
foi=1;%1 for sum(I_j/N_j), 2 for (sum I_j)/N
ages=5;
if ages==4
    %death=[.005,.0072,.04,.079]'/100;
    hosp=[.042,.016,.029,.166]';
    NNbar=[19169690;62121035;184015269;39570590];
    nbar=length(NNbar);
    Cc=[27.57920413,8.051767033,4.975736133,0.850626995;
        9.165259795,43.43045174,8.195858852,2.158756533;
        5.941537452,5.863025518,14.20166331,5.533694466;
        0.600583289,0.807369258,1.39444674,7.848296781]';
    tdays=30;%Days per month
    simCut=3;%Cut this many months from start of year
    legString={'0-4','5-19','20-64','65+'};
    from57=13/45;
    t1=1;%Plot data from month t1
elseif ages==5
    NNbar=[1464566;3790730;10236474;3984200;2412129];
    nbar=length(NNbar);
    hosp=1./[393.0256 999.3054 702.6867 406.0680 177.1958]';
    Cc=[1.9200    0.4268    0.4600    0.4434    0.1665;
        1.7600    8.7522    2.4600    1.8630    1.2190;
        0.4500    0.4911    2.5900    0.8132    0.3412;
        4.5200    4.9913    5.6500    6.9995    4.3001;
        0.2300    0.2975    0.4300    0.6903    1.8930]';
    tdays=7;%Days per week
    simCut=16;%Cut this many weeks from start of year
    legString={'0-4','5-17','18-24','25-64','65+'};
    from57=8/30;
    t1=18;%Plot data from week t1    
end
%
tswitchAttack=243;%Always calculate attck rate from here
%Changes for closure:
logPlots=0;
byAge=1;%=0 for global incidence plot, =1 to stratify by age - f
ageInc=1;%Total or age-specific incidence out (before aggregated) - g
relInc=0;%Relative incidence - fraction of age group population - both
immuneFactor=0;%In over 52s (born pre-1957) %W5*
%%
%Fixed parameters:
R0=1.6/.775;%1.4065;  1.4082;%
seednum=4;%2.9756;
%tswitch=243;
betacModifier=1;
closureFactor=0.5844;
tshift=-1;
seasonality=1;%On
ftimes=1;
tclose=10^4;
phi1=1;
phi2=0;
tlag=0;%Days
gamma=.38;%.3362;
tau=0;
tv=0;
propSym=.55;%Data
relInf=.5;%Data
mu=[.005,.0072,.04,.079,1.57]'/100;%Data
mu=repmat(mu,2,1);
sub=.1;%Relative size of sub-population
eps=0;%School closure coupling
%%
%Input parameters:
%
seednum=params(1);
closureFactor=params(2);
immuneFactor=params(3);%[4,.42,.2,1.38,.38]
%sub=params(4);
%phi2=params(4);
Cc=reshape(params(4:4+nbar^2-1),5,5);
%}
%
%propSymp=params(end-2);
R0=params(end-1);
gamma=params(end);
%}
phi1=1-phi2;
%%
%Calculate remaining/dependent parameters:
%gamma=1/TR;
gammabar=1/(1/gamma-1/tau-1);
NN=sum(NNbar);
NNrep=repmat(NN,nbar,1);
%Age mixing:
Co=Cc;
%Cc(2:3,2:3)=closureFactor*Co(2:3,2:3);
Cc(2,2)=closureFactor*Co(2,2);
%
%Calculate betas (sub population only):
Sstart=repmat(NNbar*sub,1,nbar);%Ni
Mj=NNbar'*sub; Mj(Mj==0)=1; Mjover=1./Mj;
Mjover=repmat(Mjover,nbar,1);%1/Nj
if foi==1
    Dc=(Sstart.*Mjover).*Cc;
    %Dc=Sstart*Cc/NN;
    Gc=1/gamma*Dc;
    d=eigs(Gc,1); R0c=max(d); betac=R0/R0c*betacModifier;
    Do=(Sstart.*Mjover).*Co;
    Go=1/gamma*Do;
    d=eigs(Go,1); R0o=max(d); betao=R0/R0o;
%}
elseif foi==2
    error('NN/Mjover need to be modified for sub')
    Dc=Sstart.*Cc/NN/sub;
    Gc=1/gamma*Dc;
    d=eigs(Gc,1); R0c=max(d); betac=R0/R0c*betacModifier;
    Do=(Sstart.*Mjover).*Co;
    Go=1/gamma*Do;
    d=eigs(Go,1); R0o=max(d); betao=R0/R0o;
end
%Initial condition (sub-population only):
if ages==4
    zn=zeros(nbar,1);
    y0=[NNbar;zn;zn;zn;zn];
    a3out=y0(nbar-1)*immuneFactor*from57;
    y0(4)=y0(4)-a3out; y0(end-nbar-1)=a3out;
    a4out=y0(nbar)*immuneFactor;
    y0(5)=y0(5)-a4out; y0(end-nbar)=a4out;
elseif ages==5
    zn=zeros(nbar,1);
    y0=[NNbar*sub;zn;zn;zn;zn];
    a4out=y0(nbar-1)*immuneFactor*from57;
    y0(4)=y0(4)-a4out; y0(end-nbar-1)=a4out;
    a5out=y0(nbar)*immuneFactor;
    y0(5)=y0(5)-a5out; y0(end-nbar)=a5out;
end
betac=betao*betacModifier;%Assume R0 fitted to Oct wave
%betao=betac;
%}
%%
%For simulation:
t0=0; tend=450;
%Full FOI:
if foi==1
    Dc=kron([1,eps;eps,1],Cc);%Tuned for eps=0 (calculating betas above)
    Do=kron([sub,1-sub;1-sub,sub],Co);%repmat(Co,2,2)/2;
    NN0=NNbar;
    NN0(NN0==1)=1;
elseif foi==2
    error('NN/Mjover need to be modified for sub')
    Dc=kron([1,eps;eps,1],Cc.*Mjover*NN);
    Do=repmat(Co.*Mjover*NN,2,2)/2;
    NN0=NNrep;
    NN0(NNrep==0)=1;
end
seed=10^(-seednum);
%%
%Do=kron(ones(2),Do);
%Dc=kron([1,eps;eps,1],Dc);
%NNbar=[sub*NNbar;(1-sub)*NNbar];
NN0=[sub*NN0;(1-sub)*NN0];
y0=[sub*y0;(1-sub)*y0];
%%
%Simulate:
%ODE solver:
[tout,yout]=ode45(@(t,y)integr8all(t,y,betac,betao,gamma,tau,gammabar,nbar,NN0,Dc,Do,seasonality,phi1,phi2,seed,hosp,tlag,tswitch,tclose,tv,propSym,relInf,mu),[t0,tend],y0);
%Incidence curve in here:
Y=yout(:,1:nbar);
Y=-diff(Y,1,1);
tdiff=diff(tout);
Y=Y./repmat(tdiff,1,nbar);
tout=tout(2:end);
if ageInc==1
    g=[tout,Y];
else
    g=[tout,sum(Y,2)];
end
z2=sum(Y(tout>tswitchAttack,:),1);
if relInc==1
if byAge==1
    NNdiv=repmat(NNbar',size(Y,1),1);
    Y=Y./NNdiv;
else
    Y=sum(Y,2)/NN;
end
elseif byAge==0
    Y=sum(Y,2);
end
%
if plotEpis==1
    figure
    fs=12; lw=2;
    colormap lines
    if logPlots==0
        plot(tout,Y,'linewidth',lw);
    elseif logPlots==1
        semilogy(tout,Y,'linewidth',lw);
    end
    xlabel('Time (days)','FontSize',fs);
    ylabel('Incidence','FontSize',fs);
    set(gca,'FontSize',fs);
    maxY=max(max(Y));
    axis([0,tend,0,maxY]);
    if byAge==1
        legend(legString,'location','NW')
    end
    grid on
    grid minor
    box on
    hold off
end
tmonth=ceil(tout/tdays);%+tshift);
tmonth(tmonth==0)=1;
if byAge==1
    f1=accumarray(tmonth,Y(:,1));
    f2=accumarray(tmonth,Y(:,2));
    f3=accumarray(tmonth,Y(:,3));
    f4=accumarray(tmonth,Y(:,4));
    if ages==4
        fall=ftimes*[f1,f2,f3,f4];
    else
        f5=accumarray(tmonth,Y(:,5));
        fall=ftimes*[f1,f2,f3,f4,f5];
    end
    fall(1:simCut,:)=[];
    if plotComp==1
        figure
        fs=12; lw=2;
        maxf=max(max([fall;ydata]));
        cmap=lines(nbar);
        hold on
        simVec=simCut+1:simCut+size(fall,1);
        h1=plot(simVec,ftimes*fall(:,1),'linewidth',lw,'color',cmap(1,:));
        h2=plot(simVec,ftimes*fall(:,2),'linewidth',lw,'color',cmap(2,:));
        h3=plot(simVec,ftimes*fall(:,3),'linewidth',lw,'color',cmap(3,:));
        h4=plot(simVec,ftimes*fall(:,4),'linewidth',lw,'color',cmap(4,:));
        if ages==5
            h5=plot(simVec,ftimes*fall(:,5),'linewidth',lw,'color',cmap(5,:));
        end
        plotVec=t1:t1+size(ydata,1)-1;
        plot(plotVec,ydata(:,1),'--','linewidth',lw,'color',cmap(1,:))
        plot(plotVec,ydata(:,2),'--','linewidth',lw,'color',cmap(2,:))
        plot(plotVec,ydata(:,3),'--','linewidth',lw,'color',cmap(3,:))
        plot(plotVec,ydata(:,4),'--','linewidth',lw,'color',cmap(4,:))
        if ages==5
            plot(plotVec,ydata(:,5),'--','linewidth',lw,'color',cmap(5,:))
        end
        xlabel('Time (weeks)','FontSize',fs);
        ylabel('Incidence','FontSize',fs);
        set(gca,'FontSize',fs);
        axis([simCut+1,size(fall,1),0,maxf]);%tend
        if ages==4
            legend([h1,h2,h3,h4],legString,'location','NW')
        else
            legend([h1,h2,h3,h4,h5],legString,'location','NW')
        end
        %legend([h1,h2,h3,h4,h5],legString,'location','NW')
        grid on
        grid minor
        box on
        hold off
    end
else
    fall=ftimes*accumarray(tmonth,Y);
    if relInc==1
        ydata=sum(ydata,2)/NN;
    else
        ydata=sum(ydata,2);
    end
    if plotComp==1 %&& 
        figure
        fs=12; lw=2;
        fall=accumarray(tmonth,Y);
        fall(1:simCut)=[];
        maxf=max([fall;ydata]);
        hold on
        %plot([1,1],[0,maxf],'k--','linewidth',lw)
        %plot([12,12],[0,maxf],'k--','linewidth',lw)
        h1=plot(1:length(fall),ftimes*fall,'linewidth',lw);
        h2=plot(t1:size(ydata,1),ydata(t1:end,1),'--','linewidth',lw);
        xlabel('Time (months)','FontSize',fs);
        ylabel('Incidence','FontSize',fs);
        set(gca,'FontSize',fs);
        axis([0,size(ydata,1),0,maxf]);%tend
        legend([h1,h2],'Sim','Data','location','NW')
        grid on
        grid minor
        box on
        hold off
    end
end
f=fall(xdata,:);
end
%%
function f=integr8all(t,y,betac,betao,gamma,tau,gammabar,nbar,NNin,Dc,Do,seasonality,phi1,phi2,seed,hosp,tlag,tswitch,tclose,tv,propSym,relInf,mu)
%propSym=.55;%Data
%relInf=.5;%Data
%mu=[.005,.0072,.04,.079,1.57]'/100;%Data
%mu=repmat(mu,2,1);
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
S=[y(1:nbar);y(5*nbar+1:6*nbar)];
%I=y(nbar+1:2*nbar);
%IV=y(2*nbar+1:3*nbar);
IS=[y(nbar+1:2*nbar);y(6*nbar+1:7*nbar)];
IA=[y(2*nbar+1:3*nbar);y(7*nbar+1:8*nbar)];
if t<30
    seed1=seed.*S./NNin;
    seed1(5*nbar+1:end)=0;%Seed onbly sub population
else
    seed1=0;
end
if t<tv
    taux=0;
else
    taux=tau;
end
Sfoi=phi*(beta*S.*(XX*((IS+relInf*IA)./NNin)+seed1));
%{
Sdot=-Sfoi;
Idot=Sfoi-gamma*I-.55*taux*I.*hosp;
IVdot=.55*taux*I.*hosp-gammabar*IV;
Rdot=gamma*I+gammabar*IV;
f=[Sdot;Idot;IVdot;Rdot];
%}
Sdot=-Sfoi;
ISdot=propSym*Sfoi-gamma*IS-mu.*IS;%-.55*taux*I.*hosp;
IAdot=(1-propSym)*Sfoi-gamma*IA;
%IVdot=.55*taux*I.*hosp-gammabar*IV;
Rdot=gamma*(IS+IA);%+gammabar*IV;
Ddot=mu.*IS;
f=[Sdot(1:nbar);ISdot(1:nbar);IAdot(1:nbar);Rdot(1:nbar);Ddot(1:nbar);
    Sdot(nbar+1:end);ISdot(nbar+1:end);IAdot(nbar+1:end);Rdot(nbar+1:end);Ddot(nbar+1:end)];
end