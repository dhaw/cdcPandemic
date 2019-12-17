function [f,g,z2]=cdcPandemicSimulationW5(params,xdata,plotComp,plotEpis,ydata,tswitch)%(R0,phi1,phi2,tlag,seednum,tswitch,closureFactor,betacModifier)
%V: antiViral treatment
%plotComp: plot comparison (with data)
%plotEpis: plot incidence (non-aggregated)
ages=5;
if ages==4
    %death=[.005,.0072,.04,.079]'/100;%W5*
    hosp=[.042,.016,.029,.166]';
    NNbar=[19169690;62121035;184015269;39570590];
    Cc=[27.57920413,8.051767033,4.975736133,0.850626995;9.165259795,43.43045174,8.195858852,2.158756533;5.941537452,5.863025518,14.20166331,5.533694466;0.600583289,0.807369258,1.39444674,7.848296781];
    tdays=30;%Days per timestep
    simCut=3;%Cut this many months from start of year
    legString={'0-4','5-19','20-64','65+'};
    from57=13/45;
    plotTo=12;%Months
    t1=1;%Plot data from month t1
else
    NNbar=[1464566;3790730;10236474;3984200;2412129];
    hosp=1./[393.0256 999.3054 702.6867 406.0680 177.1958]';
    Cc=[1.9200    0.4268    0.4600    0.4434    0.1665;1.7600    8.7522    2.4600    1.8630    1.2190;0.4500    0.4911    2.5900    0.8132    0.3412;4.5200    4.9913    5.6500    6.9995    4.3001;0.2300    0.2975    0.4300    0.6903    1.8930]';
    tdays=7;%Days per timestep
    simCut=17;%Cut this many weeks from start of year
    legString={'0-4','5-17','18-24','25-64','65+'};
    from57=8/30;
    plotTo=52;%Weeks
    t1=18;%Plot data from month t1
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
R0=1.4;%1.4065;  1.4082;%
seednum=4;%2.9756;
%tswitch=243;
betacModifier=1;%.9;
closureFactor=0.5844;
tshift=-1;
seasonality=1;%On
ftimes=1;
tclose=10^4;
phi1=1;
%phi2=0;
tlag=0;%Days
gamma=.3362;
tau=0;
tv=0;
%%
%Input parameters:
%
seednum=params(1);
closureFactor=params(2);
immuneFactor=params(3);%[4,.42,.2,1.38,.38]
%monthShift=params(3);
%tswitch=params(3);
%phi1=params(3);
%}
%
R0=params(end-1);
gamma=params(end);
%}
phi2=1-phi1;
%%
%gamma=1/TR;
gammabar=1/(1/gamma-1/tau-1);
nbar=length(NNbar);
NN=sum(NNbar);
NNrep=repmat(NN,nbar,1);
%Age mixing:
%USA:
Co=Cc;
Cc(2:3,2:3)=closureFactor*Co(2:3,2:3);
%Cc(2,2)=closureFactor*Co(2,2);
%Calculate betas:
%
Sstart=repmat(NNbar,1,nbar);
Mj=NNbar; Mj(Mj==0)=1; Mjover=1./Mj;
Mjover=repmat(Mjover',nbar,1);
Dc=(Sstart.*Mjover).*Cc;
%Dc=Sstart*Cc/NN;
Gc=1/gamma*Dc;
d=eigs(Gc,1); R0c=max(d); betac=R0/R0c*betacModifier;
Do=(Sstart.*Mjover).*Co;
Go=1/gamma*Do;
d=eigs(Go,1); R0o=max(d); betao=R0/R0o;
betac=betao*betacModifier;%Assume R0 fitted to Oct wave
%betao=betac;
%}
%%
%For simulation:
t0=0; tend=450;
NN0=NNrep; NN0(NNrep==0)=1;
Ni=repmat(NNrep,1,nbar); Nj=Ni';
Dc=Cc.*Mjover*NN;
Do=Co.*Mjover*NN;
zn=zeros(nbar,1);
%
y0=[NNbar;zn;zn;zn];
%
a3out=y0(nbar-1)*immuneFactor*from57;
a4out=y0(nbar)*immuneFactor;
y0(4)=y0(4)-a3out; y0(14)=a3out;
y0(5)=y0(5)-a4out; y0(15)=a4out;
%}
seed=10^(-seednum);
%%
%Simulate:
%ODE solver:
%if solvetype==2
    [tout,yout]=ode45(@(t,y)integr8all(t,y,betac,betao,gamma,tau,gammabar,nbar,NN0,Dc,Do,seasonality,phi1,phi2,seed,hosp,tlag,tswitch,tclose,tv),[t0,tend],y0);
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
            semilogy(tout/7,Y,'linewidth',lw);
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
            axis([0,size(fall,1),0,maxf]);%tend
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
            fall(1:plotTo)=[];
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
function f=integr8all(t,y,betac,betao,gamma,tau,gammabar,nbar,NN0,Dc,Do,seasonality,phi1,phi2,seed,hosp,tlag,tswitch,tclose,tv)
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
IV=y(2*nbar+1:3*nbar);
if t<30
    seed1=seed.*S./NN0;
else
    seed1=0;
end
if t<tv
    taux=0;
else
    taux=tau;
end
Sfoi=phi*(beta*S.*(XX*(I./NN0)+seed1));
Sdot=-Sfoi;
Idot=Sfoi-gamma*I-.55*taux*I.*hosp;
IVdot=.55*taux*I.*hosp-gammabar*IV;
Rdot=gamma*I+gammabar*IV;
f=[Sdot;Idot;IVdot;Rdot];
end