function [f,gg,z2]=subPandemicSimulationVaxScenarios(NNbar,params,xdata,plotComp,plotEpis,ydata,tswitch,scenarios)%(R0,phi1,phi2,tlag,seednum,tswitch,closureFactor,betacModifier)
gg=0;
numDays=52*7;

%hosp=scenarios.hosp;%./mean multiplier;
trigger=scenarios.trigger;
tau=scenarios.tau;
tv=scenarios.tv;
vaxparams=scenarios.vaxparams;
R0=scenarios.R0;

t0=243;%params(end-2);%181;%243; %scenarios.t0; %Change for 2009
tcount=243;%181;%1st July %243 %Change for 2009

vaxon=scenarios.vaxon;%+t0;%Change for 2009 - +t0 for hyp
tclose=scenarios.dayx+243;%t0;% Change for 2009 - +243 for 2009, +t0 or 243 for hyp
av=scenarios.av;

mcmc=0;
y0in=100*NNbar/sum(NNbar);%%Change for 2009 - 0 for 2009
beta2=0;
age2mats=1;
cf=.6078;%params(1);%.6078;%Change for 2009
ph2=.0009;%params(2);%.0009;%Change for 2009
t0shift=0;
dMonth=[-inf,273,304,334,365,396,424,10^4];%455];%10^4 in case go past April ,243
%dMonth replaced by vaxon
%%
%V: antiViral treatment
%plotComp: plot comparison (with data)
%plotEpis: plot incidence (non-aggregated)
foi=1;%1 for (sum I_j)/N, 2 for sum(I_j/N_j)
ages=5;
agesOut=(1:5);
hospOut=0;
if ages==4
    %death=[.005,.0072,.04,.079]'/100;
    hosp=[.042,.016,.029,.166]';
    repmu=repmat(1./hosp,length(xdata),1);
    %NNbar=[19169690;62121035;184015269;39570590];
    nbar=length(NNbar);
    Cc=[27.57920413,8.051767033,4.975736133,0.850626995;
        9.165259795,43.43045174,8.195858852,2.158756533;
        5.941537452,5.863025518,14.20166331,5.533694466;
        0.600583289,0.807369258,1.39444674,7.848296781];
    mu=zeros(4,1);
    tdays=30;%Days per month
    simCut=3;%Cut this many months from start of year
    legString={'0-4','5-19','20-64','65+'};
    from57=13/45;
    t1=1;%Plot data from month t1
elseif ages==5
    %NNbar=[1464566;3790730;10236474;3984200;2412129];
    nbar=length(NNbar);
    
    meanHosp=[393.0256  999.3054  702.6867  406.0680  177.1958]./scenarios.muTimes';
    hosp=1./meanHosp';%Rate
    repmu=repmat(meanHosp,length(xdata),1);
    Cc=[1.9200    0.4268    0.5260    0.2554    0.1665;
    1.7600    8.7522    2.2855    1.0876    1.2190;
    4.0700    4.5939    6.6160    4.5939    2.9494;
    0.9000    0.8885    1.6180    2.3847    1.6919;
    0.2300    0.2975    0.5712    0.8756    1.8930]';
    %{
    Cc=[1.9200    1.3000    3.1500    0.7600    0.3900
    0.5762    8.7413    4.5240    1.0857    0.8492
    0.6872    2.3496    6.6612    2.3094    1.1139
    0.3025    0.8968    3.2809    2.3960    1.2319
    0.1260    0.4516    1.9047    1.5044    1.9544];
    %}
    %Cc=ones(5);
    %mu=[.005,.0072,.04,.079,1.57]'./scenarios.muTimes/100;%Data
    %mu=hosp.*[.44,.16,3.4,4,7.6]'/100;
    mu=[.005,.0072,.04,.079,.27]'./scenarios.muTimes/100;
    
    tdays=7;%Days per week
    %simCut=16;%Cut this many weeks from start of 'year
    legString={'0-4','5-17','18-49','50-64','65+'};
    from57=24/30;
    t1=17;%Plot data from week t1
end
%
tswitchAttack=243;%Always calculate attck rate from here
%Changes for closure:
logPlots=0;
byAge=1;%=0 for global incidence plot, =1 to stratify by age - f
ageInc=1;%Total or age-specific incidence out (before aggregated) - g
relInc=0;%Relative incidence - fraction of age group population - both
%%
%Fixed parameters:
seednum=6;
%tswitch=243;
closureFactor=cf;
betacModifier=1;
phi1=1;
phi2=ph2;
seasonality=1;
ftimes=1;
%tclose=10^4;
tlag=31;%Days
%tau=0;
%tv=0;
tend=720;%End of April=484
%%
%Input parameters:
Cc=reshape(params(1:1+nbar^2-1),nbar,nbar);
%Exclude for MCMC:
if mcmc==0
    %Presented other way round to main simulation script
    if age2mats==1%So far only relevant for MLE fit
        phi2=params(1);
        tlag=params(2);
        Cc1=reshape(params(3:3+nbar^2-1),nbar,nbar);%Cc2;%Open=Cc1
        Cc2=Cc1;
        Cc2(1)=Cc1(1)-1.4149;
        Cc2(7)=Cc1(7)-6.1702;
    else
        %
        closureFactor=params(1);
        %betacModifier=params(2);
        %adultsDown=params(3);
        phi2=params(2);
        Cc=reshape(params(3:3+nbar^2-1),nbar,nbar);
        %}
        %{
        %W2 only:
        phi2=params(1);
        Cc=reshape(params(2:2+nbar^2-1),nbar,nbar);
        %}
    end
else
    if age2mats==1%So far only relevant for MLE fit
        %{
        Cc1=reshape(params(1:nbar^2),nbar,nbar);
        Cc2=reshape(params(nbar^2+1:2*nbar^2),nbar,nbar);
        %}
        %
        phi2=ph2;
        Cc1=reshape(params(1:nbar^2),nbar,nbar);
        Cc2=reshape(params(nbar^2+1:2*nbar^2),nbar,nbar);
        %
        %Better 1st wave - p2mn2c:
        %{
            Cc2=Cc1.*[4.0291    0.7121    0.5461    0.2756    0.5610
                    40.6147    0.6680    6.8128    1.7091    0.8669
                    1.3654    0.1002    0.6735    0.6426    1.0441
                    6.5303    0.4719    1.4841    0.6083    0.8570
                    2.0527    0.6075    0.8821    0.3976    3.1918];
        %}
        %{
            Cc2=Cc1.*[4.0131    0.7130    0.5473    0.2757    0.5507;
                    40.9427    0.6652    7.0092    1.7160    0.8702;
                    1.3737    0.1002    0.6736    0.6437    1.0445;
                    6.5318    0.4683    1.5042    0.6019    0.8633;
                    2.1013    0.6041    0.8915    0.4044    3.1849];
        %}
        %{
         Cc2=Cc1.*[3.0748    1.1730    0.8862    0.0583    4.7613;
                19.1443    0.9805    0.5643  351.9906  184.0789;
                0.9829    0.7847    0.4540    0.0071    0.5806;
                0.4814    0.6616    0.8005    1.8496    0.8138;
                140.7858    3.0448    1.4042    0.0097  132.3989];
        %}
    else
        Cc=reshape(params(1:nbar^2),nbar,nbar);
    end
end
%
%propSymp=params(end-2);
immuneFactor=0;%params(end-5); %Change for 2009
propSym=.55;%params(end-4); %Change for 2009
relInf=.5;%params(end-3); %Change for 2009
%t0=params(end-2)-120+t0shift;
%R0=1.46; %params(end-1);
gamma=1/3;%params(end); %1/3; %Change for 2009
%}
%Requiring imput parameters:
%phi1=1-phi2;
seedOn=t0+30+t0shift;
%R0=R0*(propSym+(1-propSym)*relInf);*Account for asym?
%Cc(3,:)=Cc(3,:)*adultsDown;
%Cc(:,3)=Cc(:,3)*adultsDown;
%Cc(3,3)=Cc(3,3)*adultsDown;
%%
%gamma=1/TR;
gammabar=1/(1/gamma-1/tau-1);
%nbar=length(NNbar);
NN=sum(NNbar);
NNrep=repmat(NN,nbar,1);
%Age mixing:
%USA:
if age2mats==1
    Co=Cc1;
    Cc=Cc2;
else
    Co=Cc;
end
if age2mats==0
    %Cc(2:3,2:3)=closureFactor/betacModifier*Co(2:3,2:3);
    if beta2==0
        Cc(2,2)=closureFactor/betacModifier*Co(2,2);
    else
        Cc(2,:)=closureFactor/betacModifier*Co(2,:);
    end
    %Cc(:,2)=closureFactor/betacModifier*Co(:,2);
end
%Calculate betas:
%
Sstart=repmat(NNbar,1,nbar);%Ni
Mj=NNbar'; Mj(Mj==0)=1; Mjover=1./Mj;
Mjover=repmat(Mjover,nbar,1);%1/Nj
if foi==1
    %{
    Dc=(Sstart.*Mjover).*Cc;
    
    Dc=[propSym*Dc,propSym*relInf*Dc;(1-propSym)*Dc,(1-propSym)*relInf*Dc];
    %Dc=Sstart*Cc/NN;
    Gc=Dc/gamma;%1/gamma*Dc;
    d=eigs(Gc,1); R0c=max(d); %betac=R0/R0c*betacModifier;
    %}
    Do=(Sstart.*Mjover).*Co;
    Do=[propSym*Do,propSym*relInf*Do;(1-propSym)*Do,(1-propSym)*relInf*Do];
    Go=1/gamma*(1+phi2)*Do;
    d=eigs(Go,1); R0o=max(d); betao=R0/R0o;
    %}
else
    Dc=Sstart.*Cc/NN;
    Gc=1/gamma*Dc;
    d=eigs(Gc,1); R0c=max(d); %betac=R0/R0c*betacModifier;
    Do=(Sstart.*Mjover).*Co;
    Go=1/gamma*Do;
    d=eigs(Go,1); R0o=max(d); betao=R0/R0o;
end
if ages==4
    zn=zeros(nbar,1);
    y0=[NNbar-y0in;zn;zn;zn;zn;zn;zn;zn+y0in];
    a3out=y0(nbar-1)*immuneFactor*from57;
    y0(4)=y0(4)-a3out; y0(end-nbar-1)=a3out;
    a4out=y0(nbar)*immuneFactor;
    y0(5)=y0(5)-a4out; y0(end-nbar)=a4out;
elseif ages==5
    zn=zeros(nbar,1);
    y0=[NNbar-y0in;zn;zn;propSym*y0in;zn;(1-propSym)*y0in;zn;zn];%+y0in];
    a4out=y0(nbar-1)*immuneFactor*from57;
    y0(4)=y0(4)-a4out; y0(end-nbar-1)=a4out;
    a5out=y0(nbar)*immuneFactor;
    y0(5)=y0(5)-a5out; y0(end-nbar)=a5out;
end
%betao=betao/(propSym+(1-propSym)*relInf);%*Account for asym?
betac=betao*betacModifier;%Assume R0 fitted to Oct wave

%betao=betao*betaModifier2;

%betao=betac;
%}
%%
%For simulation:
%Ni=repmat(NN0,1,nbar); Nj=Ni';
if foi==1
    %
    Dc=Cc;
    Do=Co;
    NN0=NNbar;
    NN0(NN0==1)=1;
    %}
    %{
    %As 4 age-group only, wrong:
    Dc=Cc.*Mjover*NN;
    Do=Co.*Mjover*NN;
    NN0=NNrep;
    NN0(NNrep==0)=1;
    %}
elseif foi==2
    Dc=Cc;%.*Mjover*NN;
    Do=Co;%.*Mjover*NN;
    NN0=NNrep;
    NN0(NNrep==0)=1;
end
seed=0;%10^(-seednum);%Change for 2009 - 0 for hyp
%%
%Simulate:
%ODE solver:
%if solvetype==2
global trig
trig=0;
    [tout,yout]=ode45(@(t,y)integr8all(t,y,betac,betao,gamma,tau,gammabar,nbar,NN0,Dc,Do,seasonality,phi1,phi2,seed,seedOn,hosp,tlag,tswitch,tclose,trigger,tv,propSym,relInf,mu,vaxparams,vaxon,av,dMonth),(t0:1:tend),y0);
    %Loads cut here (at bottom of script)
    yout=yout((tout>tcount-1),:);
    yout=yout(1:numDays+1,:);
    IS=-diff(yout(:,1:nbar)+yout(:,nbar+1:2*nbar)+yout(:,2*nbar+1:3*nbar),1)*propSym;%yout(1:numDays,2*nbar+1:3*nbar);
    %{
    if tau==0
        AV=zeros(numDays,5);
    else
        gammaS=(1-av)+av*1/(1/gamma-1);
        AV=IS.*repmat(av',numDays,1).*repmat(gammaS',numDays,1);
    end
    %}
    AV=diff(yout(:,4*nbar+1:5*nbar));
    H=IS.*repmat(hosp',numDays,1);%*gammaS;
    D=yout(:,end-nbar+1:end);

    D=diff(D,1);
    f=[IS,H,D,AV];
    %gg=0;
    z2=0;
end
%%
function f=integr8all(t,y,betac,betao,gamma,tau,gammabar,nbar,NNin,Dc,Do,seasonality,phi1,phi2,seed,seedOn,hosp,tlag,tswitch,tclose,trigger,tv,propSym,relInf,mu,vaxparams,vaxon,av,dMonth)
%propSym=.55;%Data
%relInf=.5;%Data
%mu=[.005,.0072,.04,.079,1.57]'/100;%Data
na=5;
if seasonality==1 %&& t>60
    phi=phi1+phi2*cos(2*pi*(t-tlag)/365);
else
    phi=phi1-phi2;
end
S=y(1:nbar);
SVH=y(nbar+1:2*nbar);
SV=y(2*nbar+1:3*nbar);
IS=y(3*nbar+1:4*nbar);
IV=y(4*nbar+1:5*nbar);
IA=y(5*nbar+1:6*nbar);
%g=max(g,sum([IS;IV],2)/sum(NNin));%g=max(g,max(sum([IS;IV],2))/sum(NNin));
global trig
%trig=max(trig,sum(IS+IV)/sum(NNin));
if trig>trigger/7
    trig=2;
end
if t<tswitch 
    XX=Dc;
    beta=betac;
elseif t>tclose && t<243+vaxon+35
    XX=Dc;
    beta=betac;
elseif trig>1 && t<243+vaxon+35
    XX=Dc;
    beta=betac;
else
    XX=Do;
    beta=betao;
end
if t<seedOn%t>seedOn && t<seedOn+14
    seed1=seed.*S./NNin;
    %seed1([1,2,4:nbar])=0;
else
    seed1=0;
end
if t<tv || tau==0
    taux=0;
    gammaS=gamma*ones(na,1);
    av=zeros(na,1);
else
    %taux=tau;
    gammaS=1./((1-av)*1./gamma+av*1./gamma*2/3);%1/(1/gamma-1);%1+hosp*(1/(1/gamma-1)-gamma);
    taux=0;
end
%
if t<243+vaxon%Change - this is hyp
    v1=0;
    v2=0;
elseif t<243+vaxon+35%sum(vaxparams(:,2).*NNin*(t-vaxon))>1.25e8
    v1=vaxparams(:,2);
    v2=vaxparams(:,3);
else
    v1=vaxparams(:,1);
    v2=vaxparams(:,3);
end
%}
Sfoi=phi*(beta*(XX*((IS+relInf*IA)./NNin)+seed1));
SVfoi=phi*(beta*(1-v2).*(XX*((IS+relInf*IA)./NNin)+seed1));
Sdot=-Sfoi.*S-v1.*S;
SVHdot=-SVfoi.*SVH+v1.*S-SVH/14;%14 days to work
SVdot=SVH/14-SVfoi.*SV;
ISdot=propSym*(Sfoi.*(S+SVH)+SVfoi.*SV)-gammaS.*IS-mu.*IS;%-.55*taux*I.*hosp;

IVdot=av.*IS;

IAdot=(1-propSym)*(Sfoi.*(S+SVH)+SVfoi.*SV)-gamma*IA;
Rdot=gammaS.*IS+gamma*IA;
Ddot=mu.*IS.*gammaS;
f=[Sdot;SVHdot;SVdot;ISdot;IVdot;IAdot;Rdot;Ddot];

trig=max(trig,propSym*sum(Sfoi.*(S+SVH)+SVfoi.*SV)/sum(NNin));
%{
Sfoi=phi*(beta*S.*(XX*((IS+relInf*IA)./NNin)+seed1));
SVfoi=phi*(beta*(1-v2).*SV.*(XX*((IS+relInf*IA)./NNin)+seed1));
Sdot=-Sfoi-v1.*S;
SVdot=-SVfoi+v1.*S;
%ISdot=propSym*(Sfoi+SVfoi)-.55*taux*IS.*hosp-gammaS.*IS-mu.*IS;%-.55*taux*I.*hosp;
ISdot=propSym*(Sfoi+SVfoi)-gammaS.*IS-mu.*IS;

IVdot=av.*IS;%.55*taux*IS.*hosp-gammabar*IV;

IAdot=(1-propSym)*(Sfoi+SVfoi)-gamma*IA;
Rdot=gammaS.*IS+gamma*IA;%+gammabar*IV;
Ddot=mu.*IS;
f=[Sdot;SVdot;ISdot;IVdot;IAdot;Rdot;Ddot];
%}
end
%%
%{
%Vax params:
vstart=243;
vrate=.0005;%Per capita doses administered/day
veff=.6;%beta -> (1-veff)*beta
vaxparams=[vstart,vrate,veff];
%}
%{
%W2 only:
t0shift=tswitch-params(end-2)+120;
y0in=ydata(xdata<18,:);
y0in=sum(y0in,1)';
%beforeShift=120;
%t0shift=tswitch-params(end-2)+120-beforeShift;%****
%}
%%
%{
    %Incidence curve in here:
    Y=yout(:,1:nbar)+yout(:,nbar+1:2*nbar);
    Y=-diff(Y,1,1);
    tdiff=diff(tout);
    Y=Y./repmat(tdiff,1,nbar);
    tout=tout(2:end);
    if ageInc==1
        gg=[tout,Y];
    else
        gg=[tout,sum(Y,2)];
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
        fs=10; lw=2;
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
    %tmonth=ceil(tout/tdays);%+tshift);
    %tmonth(tmonth<=0)=1;
    tunif=(1:floor(tout(end)))';
    tmonth=ceil(tunif/tdays);
    if byAge==1
        %{
        f1=accumarray(tmonth,Y(:,1));
        f2=accumarray(tmonth,Y(:,2));
        f3=accumarray(tmonth,Y(:,3));
        f4=accumarray(tmonth,Y(:,4));
        %}
        f1=accumarray(tmonth,interp1(tout,Y(:,1),tunif));
        f2=accumarray(tmonth,interp1(tout,Y(:,2),tunif));
        f3=accumarray(tmonth,interp1(tout,Y(:,3),tunif));
        f4=accumarray(tmonth,interp1(tout,Y(:,4),tunif));
        if ages==4
            fall=ftimes*[f1,f2,f3,f4];
        else
            %f5=accumarray(tmonth,Y(:,5));
            f5=accumarray(tmonth,interp1(tout,Y(:,5),tunif));
            fall=ftimes*[f1,f2,f3,f4,f5];
        end
        %fall(1:simCut,:)=[];
        if hospOut==1
            f=fall(1:90,:)./repmu;%xdata
        else
            f=fall(1:90,:);%xdata
        end
        if plotComp==1
            maxf=max(max([f;ydata]));
            figure
            fs=10; lw=2;
            %maxf=max(max([fall;ydata]));
            cmap=lines(nbar);
            hold on
            simVec=xdata(1):90;%xdata(end);%1:size(fall,1);%simCut+1:simCut+size(fall,1);
            if hospOut==1
                h1=plot(simVec,ftimes*fall(simVec,1)*hosp(1),'linewidth',lw,'color',cmap(1,:));
                h2=plot(simVec,ftimes*fall(simVec,2)*hosp(2),'linewidth',lw,'color',cmap(2,:));
                h3=plot(simVec,ftimes*fall(simVec,3)*hosp(3),'linewidth',lw,'color',cmap(3,:));
                h4=plot(simVec,ftimes*fall(simVec,4)*hosp(4),'linewidth',lw,'color',cmap(4,:));
                if ages==5
                    h5=plot(simVec,ftimes*fall(simVec,5)*hosp(5),'linewidth',lw,'color',cmap(5,:));
                end
                xlabel('Time (weeks)','FontSize',fs);
                ylabel('Hospitalisations','FontSize',fs);
            else
                h1=plot(simVec,ftimes*fall(simVec,1),'linewidth',lw,'color',cmap(1,:));
                h2=plot(simVec,ftimes*fall(simVec,2),'linewidth',lw,'color',cmap(2,:));
                h3=plot(simVec,ftimes*fall(simVec,3),'linewidth',lw,'color',cmap(3,:));
                h4=plot(simVec,ftimes*fall(simVec,4),'linewidth',lw,'color',cmap(4,:));
                if ages==5
                    h5=plot(simVec,ftimes*fall(simVec,5),'linewidth',lw,'color',cmap(5,:));
                end
                xlabel('Time (weeks)','FontSize',fs);
                ylabel('Incidence','FontSize',fs);
            end
            plotVec=xdata;%1:size(ydata,1);%t1:t1+size(ydata,1)-1;
            plot(plotVec,ydata(:,1),'--','linewidth',lw,'color',cmap(1,:))
            plot(plotVec,ydata(:,2),'--','linewidth',lw,'color',cmap(2,:))
            plot(plotVec,ydata(:,3),'--','linewidth',lw,'color',cmap(3,:))
            plot(plotVec,ydata(:,4),'--','linewidth',lw,'color',cmap(4,:))
            if ages==5
                plot(plotVec,ydata(:,5),'--','linewidth',lw,'color',cmap(5,:))
            end
            
            set(gca,'FontSize',fs);
            axis tight%([1,52-simCut,0,maxf]);%([simCut+1,size(fall,1),0,maxf]);%tend
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
        error('Adjust for hospitalisations as input')
        fall=ftimes*accumarray(tmonth,Y);
        if relInc==1
            ydata=sum(ydata,2)/NN;
        else
            ydata=sum(ydata,2);
        end
        if plotComp==1 %&& 
            figure
            fs=10; lw=2;
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
    f=f(:,agesOut);
    f(isnan(f)==1)=0;
    %{
    try
        f=fall(xdata,:)./repmu;
    catch
        ME
    end
    %}
    %}