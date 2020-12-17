function [f,g,z2]=subPandemicSimulationVax(NNbar,params,xdata,plotComp,plotEpis,ydata,tswitch,vaxparams)%(R0,phi1,phi2,tlag,seednum,tswitch,closureFactor,betacModifier)
mcmc=0;%Flag for omitting some parameters from fit
y0in=0; 
beta2=0;
age2mats=1;
cf=.6079;%.6074;
ph2=.0009;
t0shift=0;
tlag=-45;%In sinusoid %30 days
dMonth=[-inf,273,304,334,365,396,424,10^4];%455];%10^4 in case go past April ,243
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
    Cc=[27.57920413,8.051767033,4.975736133,0.850626995
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
adultsDown=1;
phi1=1;
phi2=ph2;
immuneFactor=0;%.3691;%us .3881;%mn .3691;%us 
tshift=-1;
seasonality=1;
ftimes=1;
tclose=10^4;
%tlag=30;%Days %30
tau=0;
tv=0;
propSym=.55;%Data
relInf=.5;%Data
R0=1.1804;%1.46/.775;
gamma=2.0987;%1/2.6
t0=79.9982;%0;
tend=365;%End of April=484
%%
%Input parameters:
Cc=reshape(params(1:1+nbar^2-1),nbar,nbar);
%Exclude for MCMC:
if mcmc==0
    if age2mats==1
        phi2=params(1);
        tlag=params(2);
        
        %Cc1=reshape(params(3:3+nbar^2-1),nbar,nbar);%****Cc1 or Cc2
        %Cc1=reshape(params(3+nbar^2:3+2*nbar^2-1),nbar,nbar);
        %{
        Cc1=[7.5103    5.2350    6.0040    1.9520    0.8168
        2.0601   20.8923    6.7120    2.5667    1.5605;
        1.1629    2.6155    8.8204    2.8596    0.4641;
        0.2796    0.5974    2.9078    1.9930    0.3906;
        0.0773    0.1406    0.1250    0.1165    0.1830];
        %}
        %{
        Cc2=Cc1-[1.1966    0.3022    0.3647    0.2790    0.0021;
                0.1075    5.0989    1.0993    0.9611    0.0299;
                0.0602    0.1877    0.3022    0.0978    0.0136;
                0.0163    0.0823    0.0691    0.0695    0.0305;
                0.0003    0.0013    0.0035    0.0045    0.0220];
        %}
        %{
        Cc1=Cc2+[1.1966    0.3022    0.3647    0.2790    0.0021;
                0.1075    5.0989    1.0993    0.9611    0.0299;
                0.0602    0.1877    0.3022    0.0978    0.0136;
                0.0163    0.0823    0.0691    0.0695    0.0305;
                0.0003    0.0013    0.0035    0.0045    0.0220];
        %}
        %{
        Cc1=Cc2;
        Cc1(1)=Cc1(1)+1.4149;
        Cc1(7)=Cc1(7)+6.1702;
        %}
        Cc1=reshape(params(3:3+nbar^2-1),nbar,nbar);%Cc2;%Open=Cc1
        Cc2=Cc1;
        Cc2(1)=Cc1(1)-1.4149;
        Cc2(7)=Cc1(7)-6.1702;
        %{
        Cc1=Cc2.*[1.1895    1.0613    1.0647    1.1667    1.0026
                1.0551    1.3228    1.1959    1.5986    1.0196
                1.0546    1.0773    1.0355    1.0354    1.0302
                1.0617    1.1598    1.0243    1.0361    1.0846
                1.0042    1.0095    1.0292    1.0397    1.1370];
        %}
    else
        %
        %closureFactor=params(1);
        %betacModifier=params(2);
        %adultsDown=params(3);
        closureFactor=params(1);
        phi2=params(2);
        tlag=31;%params(2);
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
        Cc2=reshape(params(1:nbar^2),nbar,nbar);
        %Cc2=reshape(params(nbar^2+1:2*nbar^2),nbar,nbar);
        %{
        Cc1=[7.5103    5.2350    6.0040    1.9520    0.8168
        2.0601   20.8923    6.7120    2.5667    1.5605;
        1.1629    2.6155    8.8204    2.8596    0.4641;
        0.2796    0.5974    2.9078    1.9930    0.3906;
        0.0773    0.1406    0.1250    0.1165    0.1830];
        Cc2=Cc1-[1.1966    0.3022    0.3647    0.2790    0.0021;
                0.1075    5.0989    1.0993    0.9611    0.0299;
                0.0602    0.1877    0.3022    0.0978    0.0136;
                0.0163    0.0823    0.0691    0.0695    0.0305;
                0.0003    0.0013    0.0035    0.0045    0.0220];
        %}
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
        %
        Cc1=Cc2+[1.1966    0.3022    0.3647    0.2790    0.0021;
                0.1075    5.0989    1.0993    0.9611    0.0299;
                0.0602    0.1877    0.3022    0.0978    0.0136;
                0.0163    0.0823    0.0691    0.0695    0.0305;
                0.0003    0.0013    0.0035    0.0045    0.0220];
        %}
        %{
        Cc1=Cc2.*[1.1895    1.0613    1.0647    1.1667    1.0026
                1.0551    1.3228    1.1959    1.5986    1.0196
                1.0546    1.0773    1.0355    1.0354    1.0302
                1.0617    1.1598    1.0243    1.0361    1.0846
                1.0042    1.0095    1.0292    1.0397    1.1370];
        %}
    else
        Cc=reshape(params(1:nbar^2),nbar,nbar);
    end
end
%
%propSymp=params(end-2);
immuneFactor=params(end-5);
propSym=params(end-4);
relInf=params(end-3);
t0=params(end-2);%-120+t0shift;
R0=params(end-1);
gamma=params(end);
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
    y0=[NNbar-y0in;zn;zn;zn;zn;zn;zn+y0in];
    a3out=y0(nbar-1)*immuneFactor*from57;
    y0(4)=y0(4)-a3out; y0(end-nbar-1)=a3out;
    a4out=y0(nbar)*immuneFactor;
    y0(5)=y0(5)-a4out; y0(end-nbar)=a4out;
elseif ages==5
    zn=zeros(nbar,1);
    y0=[NNbar-y0in;zn;zn;zn;zn;zn;zn+y0in];
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
    [tout,yout]=ode45(@(t,y)integr8all(t,y,betac,betao,gamma,tau,gammabar,nbar,NN0,Dc,Do,seasonality,phi1,phi2,seed,seedOn,hosp,tlag,tswitch,tclose,tv,propSym,relInf,mu,vaxparams,dMonth),[t0,tend],y0);
    %Incidence curve in here:
    Y=yout(:,1:nbar)+yout(:,nbar+1:2*nbar)+yout(:,2*nbar+1:3*nbar);
    Y=-diff(Y,1,1)*propSym;
    tdiff=diff(tout);
    Y=Y./repmat(tdiff,1,nbar);
    tout=tout(2:end);
    if ageInc==1
        g=[tout,Y];
    else
        g=[tout,sum(Y,2)];
    end
    %z2=sum(Y(tout>tswitchAttack,:),1);
    tsw=find(tout>tswitchAttack);
    tsw=tsw(1);
    z2=yout(end,5*nbar+1:6*nbar)+yout(end,6*nbar+1:end)-yout(tsw,5*nbar+1:6*nbar)-yout(tsw,6*nbar+1:end);
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
            f=fall(xdata,:)./repmu;
        else
            f=fall(xdata,:);
        end
        if plotComp==1
            maxf=max(max([f;ydata]));
            figure
            fs=10; lw=2;
            %maxf=max(max([fall;ydata]));
            cmap=lines(nbar);
            hold on
            simVec=xdata(1):xdata(end);%1:size(fall,1);%simCut+1:simCut+size(fall,1);
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
        
        %z2=ftimes*sum(fall(simVec>34,:),1);
        
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
end
%%
function f=integr8all(t,y,betac,betao,gamma,tau,gammabar,nbar,NNin,Dc,Do,seasonality,phi1,phi2,seed,seedOn,hosp,tlag,tswitch,tclose,tv,propSym,relInf,mu,vaxparams,dMonth)
%propSym=.55;%Data
%relInf=.5;%Data
%mu=[.005,.0072,.04,.079,1.57]'/100;%Data
if seasonality==1 %&& t>60
    phi=phi1+phi2*cos(2*pi*(t-tlag)/365);
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
SVH=y(nbar+1:2*nbar);
SV=y(2*nbar+1:3*nbar);
IS=y(3*nbar+1:4*nbar);
IA=y(4*nbar+1:5*nbar);
if t<seedOn%t>seedOn && t<seedOn+14
    seed1=seed.*S./NNin;
    %seed1([1,2,4:nbar])=0;
else
    seed1=0;
end
if t<tv
    taux=0;
else
    taux=tau;
end
%{
if t>vaxparams(1)
    vrate=vaxparams(2);
else
    vrate=0;
end
%}
thisMonth=discretize(t-14,dMonth);%-14 - vax delay
v1=vaxparams(:,thisMonth,1)/30;
v2=vaxparams(:,thisMonth,2);
Sfoi=phi*(beta*(XX*((IS+relInf*IA)./NNin)+seed1));
%SVfoi=phi*(beta*(1-vaxparams(3))*SV.*(XX*((IS+relInf*IA)./NNin)+seed1));
SVfoi=phi*(beta*(1-v2).*(XX*((IS+relInf*IA)./NNin)+seed1));
%{
Sdot=-Sfoi;
Idot=Sfoi-gamma*I-.55*taux*I.*hosp;
IVdot=.55*taux*I.*hosp-gammabar*IV;
Rdot=gamma*I+gammabar*IV;
f=[Sdot;Idot;IVdot;Rdot];
%}
Sdot=-Sfoi.*S-v1.*S;
SVHdot=-Sfoi.*SVH+v1.*S-SVH/14;%14 days to work
SVdot=SVH/14-SVfoi.*SV;
ISdot=propSym*(Sfoi.*(S+SVH)+SVfoi.*SV)-gamma*IS-mu.*IS;%-.55*taux*I.*hosp;
IAdot=(1-propSym)*(Sfoi.*(S+SVH)+SVfoi.*SV)-gamma*IA;
Rdot=gamma*(IS+IA);
Ddot=mu.*IS;
f=[Sdot;SVHdot;SVdot;ISdot;IAdot;Rdot;Ddot];
end
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