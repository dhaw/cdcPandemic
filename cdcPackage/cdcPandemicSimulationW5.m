function [f,g,z2]=cdcPandemicSimulationW5(params,xdata,plotComp,plotEpis,ydata,tswitch)%(R0,phi1,phi2,tlag,seednum,tswitch,closureFactor,betacModifier)
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
    NNbar=[19169690;62121035;184015269;39570590];
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
    NNbar=[1464566;3790730;10236474;3984200;2412129];
    nbar=length(NNbar);
    
    meanHosp=[393.0256  999.3054  702.6867  406.0680  177.1958];
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
    mu=[.005,.0072,.04,.079,1.57]'/100;%Data
    tdays=7;%Days per week
    simCut=16;%Cut this many weeks from start of 'year
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
closureFactor=.712;%
betacModifier=1;
adultsDown=1;
immuneFactor=0;
tshift=-1;
seasonality=1;
ftimes=1;
tclose=10^4;
phi1=1;
phi2=0.0058;
tlag=30;%Days
tau=0;
tv=0;
propSym=.55;%Data
relInf=.5;%Data
R0=1.1804;%1.46/.775;
gamma=2.0987;%1/2.6
t0=79.9982;%0;
tend=450;
%%
%Input parameters:
%
closureFactor=params(1);
%betacModifier=params(2);
%adultsDown=params(3);
phi2=params(2);
immuneFactor=params(3);
Cc=reshape(params(4:4+nbar^2-1),nbar,nbar);
%betaModifier2=params(3);
%}
%Cc=reshape(params(1:1+nbar^2-1),nbar,nbar);
%
%propSymp=params(end-2);
propSym=params(end-4);
relInf=params(end-3);
t0=params(end-2);
R0=params(end-1);
gamma=params(end);
%}
%Requiring imput parameters:
%phi1=1-phi2;
seedOn=t0+30;
%R0=R0*(propSym+(1-propSym)*relInf);*Account for asym?
%Cc(3,:)=Cc(3,:)*adultsDown;
%Cc(:,3)=Cc(:,3)*adultsDown;
Cc(3,3)=Cc(3,3)*adultsDown;
%%
%gamma=1/TR;
gammabar=1/(1/gamma-1/tau-1);
%nbar=length(NNbar);
NN=sum(NNbar);
NNrep=repmat(NN,nbar,1);
%Age mixing:
%USA:
Co=Cc;
%Cc(2:3,2:3)=closureFactor/betacModifier*Co(2:3,2:3);
%Cc(2,2)=closureFactor/betacModifier*Co(2,2);
Cc(2,:)=closureFactor/betacModifier*Co(2,:);
%Cc(:,2)=closureFactor/betacModifier*Co(:,2);
%Calculate betas:
%
Sstart=repmat(NNbar,1,nbar);%Ni
Mj=NNbar'; Mj(Mj==0)=1; Mjover=1./Mj;
Mjover=repmat(Mjover,nbar,1);%1/Nj
if foi==1
    Dc=(Sstart.*Mjover).*Cc;
    %Dc=Sstart*Cc/NN;
    Gc=1/gamma*Dc;
    d=eigs(Gc,1); R0c=max(d); %betac=R0/R0c*betacModifier;
    Do=(Sstart.*Mjover).*Co;
    Go=1/gamma*Do;
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
    y0=[NNbar;zn;zn;zn;zn];
    a3out=y0(nbar-1)*immuneFactor*from57;
    y0(4)=y0(4)-a3out; y0(end-nbar-1)=a3out;
    a4out=y0(nbar)*immuneFactor;
    y0(5)=y0(5)-a4out; y0(end-nbar)=a4out;
elseif ages==5
    zn=zeros(nbar,1);
    y0=[NNbar;zn;zn;zn;zn];
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
seed=10^(-seednum);
%%
%Simulate:
%ODE solver:
%if solvetype==2
    [tout,yout]=ode45(@(t,y)integr8all(t,y,betac,betao,gamma,tau,gammabar,nbar,NN0,Dc,Do,seasonality,phi1,phi2,seed,seedOn,hosp,tlag,tswitch,tclose,tv,propSym,relInf,mu),[t0,tend],y0);
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
        fall(1:simCut,:)=[];
        if hospOut==1
            f=fall(xdata,:)./repmu;
        else
            f=fall(xdata,:);
        end
        if plotComp==1
            maxf=max(max([f;ydata]));
            figure
            fs=12; lw=2;
            %maxf=max(max([fall;ydata]));
            cmap=lines(nbar);
            hold on
            simVec=1:size(fall,1);%simCut+1:simCut+size(fall,1);
            if hospOut==1
                h1=plot(simVec,ftimes*fall(:,1)*hosp(1),'linewidth',lw,'color',cmap(1,:));
                h2=plot(simVec,ftimes*fall(:,2)*hosp(2),'linewidth',lw,'color',cmap(2,:));
                h3=plot(simVec,ftimes*fall(:,3)*hosp(3),'linewidth',lw,'color',cmap(3,:));
                h4=plot(simVec,ftimes*fall(:,4)*hosp(4),'linewidth',lw,'color',cmap(4,:));
                if ages==5
                    h5=plot(simVec,ftimes*fall(:,5)*hosp(5),'linewidth',lw,'color',cmap(5,:));
                end
                xlabel('Time (weeks)','FontSize',fs);
                ylabel('Hospitalisations','FontSize',fs);
            else
                h1=plot(simVec,ftimes*fall(:,1),'linewidth',lw,'color',cmap(1,:));
                h2=plot(simVec,ftimes*fall(:,2),'linewidth',lw,'color',cmap(2,:));
                h3=plot(simVec,ftimes*fall(:,3),'linewidth',lw,'color',cmap(3,:));
                h4=plot(simVec,ftimes*fall(:,4),'linewidth',lw,'color',cmap(4,:));
                if ages==5
                    h5=plot(simVec,ftimes*fall(:,5),'linewidth',lw,'color',cmap(5,:));
                end
                xlabel('Time (weeks)','FontSize',fs);
                ylabel('Incidence','FontSize',fs);
            end
            plotVec=1:size(ydata,1);%t1:t1+size(ydata,1)-1;
            plot(plotVec,ydata(:,1),'--','linewidth',lw,'color',cmap(1,:))
            plot(plotVec,ydata(:,2),'--','linewidth',lw,'color',cmap(2,:))
            plot(plotVec,ydata(:,3),'--','linewidth',lw,'color',cmap(3,:))
            plot(plotVec,ydata(:,4),'--','linewidth',lw,'color',cmap(4,:))
            if ages==5
                plot(plotVec,ydata(:,5),'--','linewidth',lw,'color',cmap(5,:))
            end
            
            set(gca,'FontSize',fs);
            axis([1,52-simCut,0,maxf]);%([simCut+1,size(fall,1),0,maxf]);%tend
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
    f=f(:,agesOut);
    f(isnan(f)==1)=0;
    %{
    try
        f=fall(xdata,:)./repmu;
    catch
        ME
    end
    %}
end
%%
function f=integr8all(t,y,betac,betao,gamma,tau,gammabar,nbar,NNin,Dc,Do,seasonality,phi1,phi2,seed,seedOn,hosp,tlag,tswitch,tclose,tv,propSym,relInf,mu)
%propSym=.55;%Data
%relInf=.5;%Data
%mu=[.005,.0072,.04,.079,1.57]'/100;%Data
if seasonality==1 && t>60
    phi=phi1-phi2*cos(pi*(t-tlag)/180);
else
    phi=phi1-phi2;
end
if t<tswitch || t>tclose% && t>200%
    XX=Dc;
    beta=betac;
else
    XX=Do;
    beta=betao;
end
S=y(1:nbar);
%I=y(nbar+1:2*nbar);
%IV=y(2*nbar+1:3*nbar);
IS=y(nbar+1:2*nbar);
IA=y(2*nbar+1:3*nbar);
if t<seedOn%t>seedOn && t<seedOn+14
    seed1=seed.*S./NNin;
    seed1([1,3:nbar])=0;
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
f=[Sdot;ISdot;IAdot;Rdot;Ddot];
end