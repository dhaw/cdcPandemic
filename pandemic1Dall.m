function f=pandemic1Dall(params,xdata,plotComp,ydata)%(R0,phi1,phi2,tlag,seednum,tswitch,closureFactor,betacModifier)
solvetype=3;
%Changes for closure:
%tclose as input, output f
plotEpis=1; 
logPlots=0;
byAge=1;
%%
%Input parameters:
seednum=params(1);
tswitch=params(2);
betacModifier=params(3);
closureFactor=params(4);
R0=params(5);
monthShift=params(6);
ftimes=1;
tclose=10^4;
%%
%Fixed parameters:
%R0=1.7;
phi1=.5;
phi2=.5;
tlag=0;
%%
%Population:
pop=10000; NN=pop;
n=1; nbar=4;
NNrep=repmat(NN,4,1);
NNbar=pop*[1/5;1/14;1/45;1/16];
if solvetype==3
    NNbar(NNbar<5)=5;
    NNbar=round(NNbar);
    NN=reshape(NNbar,n,4);
    NN=sum(NN,2);
end
%Natural history:
gamma=1/2.6;
%seednum=4;
seasonality=0;
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
y0=[NNbar;zn;zn];
seed=10^(-seednum);
%%
%Simulate:
%ODE solver:
if solvetype==2
    [tout,yout]=ode45(@(t,y)integr8all(t,y,betac,betao,gamma,nbar,NN0,Dc,Do,seasonality,phi1,phi2,seed,tlag,tswitch,tclose),[t0,tend],y0);
elseif solvetype==3
    [tout,yout]=stochSim(y0,betac,betao,gamma,nbar,NN0,Dc,Do,seasonality,phi1,phi2,seed,tlag,tswitch,tclose);
end
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
        axis([0,100,0,maxY]);%tend
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
    tmonth=1+ceil(tout/28+monthShift);
    tmonth(tmonth<1)=1;%Cheating
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
            h1=plot(1:size(fall,1),ftimes*fall(:,1),'linewidth',lw,'color',cmap(1,:));
            h2=plot(1:size(fall,1),ftimes*fall(:,2),'linewidth',lw,'color',cmap(2,:));
            h3=plot(1:size(fall,1),ftimes*fall(:,3),'linewidth',lw,'color',cmap(3,:));
            h4=plot(1:size(fall,1),ftimes*fall(:,4),'linewidth',lw,'color',cmap(4,:));
            plot(5:size(ydata,1),ydata(5:end,1),'--','linewidth',lw,'color',cmap(1,:))
            plot(5:size(ydata,1),ydata(5:end,2),'--','linewidth',lw,'color',cmap(2,:))
            plot(5:size(ydata,1),ydata(5:end,3),'--','linewidth',lw,'color',cmap(3,:))
            plot(5:size(ydata,1),ydata(5:end,4),'--','linewidth',lw,'color',cmap(4,:))
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
    %f=[tout,Y];
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
Sfoi=phi*(beta*S.*(XX*(I./NN0))+seed1);
Sdot=-Sfoi;
Idot=Sfoi-gamma*I;
Rdot=gamma*I;
f=[Sdot;Idot;Rdot];
end

function [tout,yout]=stochSim(y,betac,betao,gamma,nbar,NN0,Dc,Do,seasonality,phi1,phi2,seed,tlag,tswitch,tclose)
%Feed in mu if required
factor=6;
tend=360*factor; betac=betac/factor; betao=betao/factor; gamma=gamma/factor;
Vec=zeros(3*nbar,tend);
S=y(1:nbar);
I=y(nbar+1:2*nbar);
R=y(2*nbar+1:end);
i=1;
threshold=30*factor;%Number of time-steps for necessary simulation/seed
while i<tend && (i<threshold || sum(I)>0)%At least 30 time-steps
    if seasonality==1
        phi=phi1-phi2*cos(pi*(t-tlag)/180);
    else
        phi=1;
    end
    if i<tswitch*factor || i>tclose*factor
        XX=Dc;
        beta=betac;
    else
        XX=Do;
        beta=betao;
    end
    seed1=seed.*S./NN0/factor;
    Sfoi=phi*(beta*S.*(XX*(I./NN0))+seed1);
    Sout=1-exp(-Sfoi);
    %Sout(Sout>1)=1;
    Sout=binornd(S,Sout); Sout(S==0)=0;
    S=S-Sout; I=I+Sout;
    %
    Iout=1-exp(-gamma);
    Iout=binornd(I,Iout);
    I=I-Iout; R=R+Iout;
    Vec(:,i)=[S;I;R];
    i=i+1;
end
tout=(1:i-1)'/factor;
yout=Vec(:,1:i-1)';
end