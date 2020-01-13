function [f,g,z2]=cdcPandemicSimulation2(params,xdata,plotComp,plotEpis,ydata,tswitch,seednum,closureFactor)%(R0,phi1,phi2,tlag,seednum,tswitch,closureFactor,betacModifier)
tswitchAttack=243;
%V: antiViral treatment
%plotComp: plot comparison (with data)
%plotEpis: plot incidence (non-aggregated)
%
%Changes for closure:
logPlots=0;
byAge=1;%=0 for global incidence plot, =1 to stratify by age
ageInc=1;%Total or age-specific incidence out
relInc=1;%Relative incidence - fraction of age group population
t1=1;%Plot data from month t1
immuneFactor=0;%In over 52s (born pre-1957)
%%
%Fixed parameters:
R0=1.4065;
%seednum=2.9756;%3.0114;%2.8394;%
%tswitch=243;%220;
betacModifier=1;%.9;
%closureFactor=0.5844;%.5778;%.5924;%0.5672;%
monthShift=-1;
seasonality=1;%On
ftimes=1;
tclose=10^4;
phi1=1;
%phi2=0;
tlag=0;%Days
gamma=.3362;%1/2.6;
tau=0;
tv=0;
%%
%Input parameters:
%{
seednum=params(1);
closureFactor=params(2);
%immuneFactor=params(3);
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
hosp=[.042,.016,.029,.166]';
death=[.005,.0072,.04,.079]'/100;
n=1; nbar=4;
NNbar=[19169690;62121035;184015269;39570590];
%NNbar=[19037307;62045041;182377351;38799891];
NN=sum(NNbar);
NNrep=repmat(NN,4,1);
%{
NN=pop;
NNbar=pop*[1/5;1/14;1/45;1/16];
if stoch==1
    NNbar(NNbar<5)=5;
    NNbar=round(NNbar);
    NN=reshape(NNbar,n,4);
    NN=sum(NN,2);
end
%}
%Age mixing:
%USA:
Cc=[27.57920413,8.051767033,4.975736133,0.850626995;9.165259795,43.43045174,8.195858852,2.158756533;5.941537452,5.863025518,14.20166331,5.533694466;0.600583289,0.807369258,1.39444674,7.848296781];
Co=Cc;
Cc(2,2)=closureFactor*Co(2,2);
%Calculate betas:
%{
Sstart=repmat(NNbar,1,nbar);
Mj=Kbar'*NNbar;
Mj(Mj==0)=1; Mjover=1./Mj;
Mjover=repmat(Mjover',nbar,1);
DD=(Sstart.*K1.*Mjover)*(Kbar'.*Cbar);
GD=1/gamma*DD;
D=(K1.*Mjover)*((Kbar').*Cbar);
D=D.*Nj;
Db=(K1.*Mjover)*((Kbar').*CbarB);
Db=Db.*Nj*.7;
%}
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
t0=0; tend=720;
NN0=NNrep; NN0(NNrep==0)=1;
Ni=repmat(NNrep,1,nbar); Nj=Ni';
Dc=Cc.*Mjover*NN;
Do=Co.*Mjover*NN;
zn=zeros(nbar,1);
%
y0=[NNbar;zn;zn;zn];
a3out=y0(3)*immuneFactor*13/45;
a4out=y0(4)*immuneFactor;
y0(3)=y0(3)-a3out; y0(11)=a3out;
y0(4)=y0(4)-a4out; y0(12)=a4out;
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
        axis([0,400,0,maxY]);
        if byAge==1
            legend('0-4','5-19','20-64','65+')
        end
        grid on
        grid minor
        box on
        hold off
    end
    tmonth=1+ceil(tout/30+monthShift);
    if byAge==1
        f1=accumarray(tmonth,Y(:,1));
        f2=accumarray(tmonth,Y(:,2));
        f3=accumarray(tmonth,Y(:,3));
        f4=accumarray(tmonth,Y(:,4));
        fall=ftimes*[f1,f2,f3,f4];
        if plotComp==1
            figure
            fs=12; lw=2;
            maxf=max(max([fall;ydata]));
            cmap=lines(4);
            hold on
            plot([1,1],[0,maxf],'k--','linewidth',lw)
            plot([12,12],[0,maxf],'k--','linewidth',lw)
            h1=plot(t1:size(fall,1),ftimes*fall(t1:end,1),'linewidth',lw,'color',cmap(1,:));
            h2=plot(t1:size(fall,1),ftimes*fall(t1:end,2),'linewidth',lw,'color',cmap(2,:));
            h3=plot(t1:size(fall,1),ftimes*fall(t1:end,3),'linewidth',lw,'color',cmap(3,:));
            h4=plot(t1:size(fall,1),ftimes*fall(t1:end,4),'linewidth',lw,'color',cmap(4,:));
            plot(t1:size(ydata,1),ydata(t1:end,1),'--','linewidth',lw,'color',cmap(1,:))
            plot(t1:size(ydata,1),ydata(t1:end,2),'--','linewidth',lw,'color',cmap(2,:))
            plot(t1:size(ydata,1),ydata(t1:end,3),'--','linewidth',lw,'color',cmap(3,:))
            plot(t1:size(ydata,1),ydata(t1:end,4),'--','linewidth',lw,'color',cmap(4,:))
            xlabel('Time (months)','FontSize',fs);
            ylabel('Incidence','FontSize',fs);
            set(gca,'FontSize',fs);
            axis([0,13,0,maxf]);%tend
            legend([h1,h2,h3,h4],'0-4','5-19','20-64','65+','location','NW')
            grid on
            grid minor
            box on
            hold off
        end
    else
        fall=ftimes*accumarray(tmonth,Y);
        ydata=sum(ydata.*repmat(NNbar',size(ydata,1),1),2)/NN;
        if plotComp==1 %&& 
            figure
            fs=12; lw=2;
            fall=accumarray(tmonth,Y);
            maxf=max([fall;ydata]);
            hold on
            plot([1,1],[0,maxf],'k--','linewidth',lw)
            plot([12,12],[0,maxf],'k--','linewidth',lw)
            h1=plot(1:length(fall),ftimes*fall,'linewidth',lw);
            h2=plot(5:length(ydata),ydata(5:end),'--','linewidth',lw);
            xlabel('Time (months)','FontSize',fs);
            ylabel('Incidence','FontSize',fs);
            set(gca,'FontSize',fs);
            axis([0,13,0,maxf]);%tend
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