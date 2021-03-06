function f=pandemic1DallVraw(params,xdata,plotComp,plotEpis,ydata)%(R0,phi1,phi2,tlag,seednum,tswitch,closureFactor,betacModifier)
%V: antiViral treatment
%plotComp: plot comparison (with data)
%plotEpis: plot incidence (non-aggregated)
%
%Changes for closure:
%tclose as input, output f
%plotEpis=1;
logPlots=0;
byAge=1;
t1=1;%Plot data from month t1
%%
%Fixed parameters:
R0modifier=1;%/.775;
R0=1.4*R0modifier;
seednum=3;
tswitch=243;%220;
betacModifier=1;%.9;
closureFactor=1;
monthShift=-1;
seasonality=1;
ftimes=1;
tclose=10^4;
phi1=1;
phi2=0;
tlag=0;%Days
gamma=1/2.6;
tau=0;
tv=243;%400;
%%
%Input parameters:
seednum=params(1);
%tswitch=params(2);
phi1=1;%params(2);
phi2=1-phi1;%params(3);
%betacModifier=params(2);
closureFactor=params(2);
%R0=params(5);
%monthShift=params(4);
%ftimes=params(3);
%tau=params(4);
%%
gammabar=1/(1/gamma-1/tau-1);
hosp=[.042,.016,.029,.166]';
death=[.005,.0072,.04,.079]'/100;
n=1; nbar=4;
NNbar=[19169690;62121035;184015269;39570590];
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
%UK:
%Cc=[1.7112,0.4927,0.4049,0.0293;0.3919,2.5527,0.4753,0.0348;0.3234,0.5548,0.8996,0.0728;0.0528,0.1904,0.3744,0.3830];
%USA:
Cc=[27.57920413,8.051767033,4.975736133,0.850626995;9.165259795,43.43045174,8.195858852,2.158756533;5.941537452,5.863025518,14.20166331,5.533694466;0.600583289,0.807369258,1.39444674,7.848296781];
%Cc=Cc';
%China:
%{
Cnum=[6.92,.25,.77,.45;.19,3.51,.57,.2;.42,.38,1.4,.17;.36,.44,1.03,1.83];
Cdur=[3.88,.28,1.04,.49;.53,2.51,.75,.5;1.31,.8,1.14,.47;1,.85,.88,1.73];
Cc=Cnum.*Cdur;
%}
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
%Do=Sstart*Co/NN;
Go=1/gamma*Do;
d=eigs(Go,1); R0o=max(d); betao=R0/R0o;
betac=betao*betacModifier;%Assume R0 fitted to Oct wave
%}
%d=eigs(Cc,1)/gamma; R0c=max(d); betac=R0/R0c*betacModifier;
%d=eigs(Co,1)/gamma; R0o=max(d); betao=R0/R0o;
%%
%For simulation:
t0=0; tend=720;
NN0=NNrep; NN0(NNrep==0)=1;
Ni=repmat(NNrep,1,nbar); Nj=Ni';
%Dc=(Mjover*Cc).*Nj;
Dc=Cc.*Mjover*NN;
%Do=(Mjover*Co).*Nj;
Do=Co.*Mjover*NN;
zn=zeros(nbar,1);
y0=[NNbar;zn;zn;zn];
seed=10^(-seednum);
%%
%Simulate:
%ODE solver:
%if solvetype==2
    [tout,yout]=ode45(@(t,y)integr8all(t,y,betac,betao,gamma,tau,gammabar,nbar,NN0,Dc,Do,seasonality,phi1,phi2,seed,hosp,tlag,tswitch,tclose,tv),[t0,tend],y0);
    %Incidence curve in here:
    %Y=yout(:,nbar+1:2*nbar);
    Y=yout(:,1:nbar);
    Y=-diff(Y,1,1);
    tdiff=diff(tout);
    Y=Y./repmat(tdiff,1,nbar);
    tout=tout(2:end);
    %Ysum=sum(Y,2)/sum(NN);
    if byAge==1
        NNdiv=repmat(NNbar',size(Y,1),1);
        Y=Y./NNdiv;
    else
        Y=sum(Y,2)/NN;
    end
    %
    if plotEpis==1
        figure
        fs=12; lw=2;
        colormap lines
        if logPlots==0
            %plot(tout/7,Ysum,'k','linewidth',lw);
            plot(tout,Y,'linewidth',lw);
        elseif logPlots==1
            %semilogy(tout/7,Ysum,'k','linewidth',lw);
            semilogy(tout/7,Y,'linewidth',lw);
        end
        xlabel('Time (days)','FontSize',fs);
        ylabel('Incidence','FontSize',fs);
        set(gca,'FontSize',fs);
        maxY=max(max(Y));%,.05]);
        axis([0,400,0,maxY]);%tend
        %title(strcat('\tau=',num2str(tau)))
        if byAge==1
            legend('0-4','5-19','20-64','65+')
        end
        grid on
        grid minor
        box on
        hold off
    end
    %rec0=yout(end,2*nbar+1:end)+yout(end,nbar+1:2*nbar);%Truncation - add infectious to rercovered
    %f=Y;%rec0';
    tmonth=1+ceil(tout/30+monthShift);
    %tmonth(tmonth<1)=1;%Cheating
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
        if plotComp==1
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
    %f=reshape(fall(xdata,:),4*length(xdata),1);%For MLE
    %f=[tout,Y];
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