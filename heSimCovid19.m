function [f,g]=heSimCovid19(pr,beta,tvec,Dvec,n,nbar,NNvec,phi1,phi2,seedvec,S0,tau,plotTau)

%{
%Feed in to function - from prep
NNages=NNvec(:,1);
NNages=reshape(NNages,n,na);
NNages=sum(NNages,1);
%}
%DEout,Rout
%sigma,omega,gamma,hvec,muvec,pvec,qvec,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3
solvetype=2;
zn=zeros(nbar,1);
NNbar=NNvec(:,1);
if solvetype==2
    lt=length(tvec);
    t0=tvec(1);%tvec(1)=start time for while simulation
    %y0=[S0;repmat(zn,11,1);NNbar-S0];%HE
    y0=[S0;repmat(zn,10,1);NNbar-S0];%IC
    toutAll=[];
    Sout=[];
    Hout=Sout;
    Dout=Sout;
    Iout=[];
    DEout=zeros(nbar,lt);
    Rout=DEout;
    Rt=zeros(lt-1,1);
    for i=1:lt-1
        D=Dvec(:,:,i);
        tend=tvec(i+1);
        NNfeed=NNvec(:,i);
        NNfeed(NNfeed==0)=1;
        %
        %pr.z=0;
        %{
        if i>=2%3
            pr.q1=pr.qnew;%Won't turn off****
            pr.q2=pr.qnew;
            %pr.z=pr.znew;
            pr.odds=pr.oddsnew;
        end
        %}
        %topen=1;%tvec(3);
        [tout,Sclass,Hclass,Dclass,DEcum,Rcum,Itot,y0]=integr8(pr,beta,nbar,NNfeed,D,phi1,phi2,seedvec,t0,tend,y0,i);
        toutAll=[toutAll;tout(2:end)];
        Sout=[Sout;Sclass(2:end,:)];
        Hout=[Hout;Hclass(2:end,:)];
        Dout=[Dout;Dclass(2:end,:)];
        DEout(:,i)=DEcum;
        %Rout(:,i)=Rcum;
        Iout=[Iout;Itot(2:end,:)];
        t0=tend;
        %%
        %HE:
        %Calculate in advance?
        %Need to add age structure in here - only uses n, not na****
        Rt(i)=heComputeEigs(pr,beta,D,NNfeed,nbar,Sclass(end,:)');
        if i<lt-1
            Xh2w=NNvec(1:n-1,i+1)-NNvec(1:n-1,i);%Addition to each wp next intervention step
            Xw2h=-Xh2w; Xw2h(Xw2h<0)=0; 
            Xw2h=Xw2h./NNvec(1:n-1,i);
            Xw2h(NNvec(1:n-1,i)==0)=0;
            if NNvec(n,i)>0
                Xh2w(Xh2w<0)=0;
                Xh2w=Xh2w/NNvec(n,i);
            else
                Xh2w=0;
            end
            %Move all infection statuses:
            %y0=reshape(y0,[n,13]);%HE
            y0=reshape(y0,[n,12]);%IC
            %y0w2h=y0(1:n-1,:).*repmat(Xw2h,1,13);%HE
            y0w2h=y0(1:n-1,:).*repmat(Xw2h,1,12);%IC
            y0w2h=[-y0w2h;sum(y0w2h,1)];
            y0h2w=y0(n,:);
            y0h2w=kron(y0h2w,Xh2w);
            y0h2w=[y0h2w;-sum(y0h2w,1)];
            y0=y0+y0w2h+y0h2w;
            %y0=reshape(y0,13*nbar,1);%HE
            y0=reshape(y0,12*nbar,1);%IC
        end
        %%
    end
    if tau==plotTau
            dodiff=1;
            plotEpi(toutAll,Sout,Hout,n,dodiff,tvec);%-tvec(2)+38
    end
elseif solvetype==1
    error('Final size calculations not possible')
elseif solvetype==3
    error('Code not written yet')
    %f=stochSim(y0,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,tau,alpha);
end
%For plots:
f=Rt;%toutAll;%(2:end);
%g=max(sum(Hout(toutAll>tvec(3),:),2));
g=[max(sum(Hout(toutAll>tvec(3),:),2)),Rt(end)];%sum(Rout(end,:))];%Rout+DEout;%sum(Hout,2);%-diff(sum(Sout,2))./diff(toutAll);
%g=cumsum(Iout);
%g=sum(Dout,2);
end

function [tout,Sclass,Hclass,Dclass,DEcum,Rcum,Itot,y0new]=integr8(pr,beta,nbar,NN0,D,phi1,phi2,seedvec,t0,tend,y0,topen)
%ncomps=13;%Number of compartments
    [tout,yout]=ode45(@(t,y)integr8covid(t,y,pr,beta,nbar,NN0,D,phi1,phi2,seedvec,topen),[t0,tend],y0);
    Sclass=yout(:,1:nbar);
    %{
    %HE:
    Hclass=yout(:,10*nbar+1:11*nbar);
    Dclass=yout(:,11*nbar+1:12*nbar);
    DEcum=yout(end,11*nbar+1:12*nbar);%Deaths
    Rcum=yout(end,12*nbar+1:end);
    Itot=sum(yout(:,4*nbar+1:8*nbar),2);
    %}
    %
    %IC:
    Hclass=yout(:,9*nbar+1:10*nbar);
    Dclass=yout(:,10*nbar+1:11*nbar);
    DEcum=yout(end,10*nbar+1:11*nbar);%Deaths
    Rcum=yout(end,11*nbar+1:end);
    Itot=sum(yout(:,4*nbar+1:7*nbar),2);
    y0new=yout(end,:)';
    %}
end


function f=plotEpi(tout,Y,H,n,dodiff,tvec)
yvar='Inc./hosp.';%'Susceptibles'; 'Hospitalisations';
solvetype=2;
tend=tvec(end);%720;%For plot only
na=size(Y,2)/n;
cmap=lines(7);
if solvetype==2
    if dodiff==1
        Y=-diff(Y,1);
        tdiff=diff(tout,1);
        Y=Y./repmat(tdiff,1,n*na);%repmat - older version?
        tout(1)=[];
        H(1,:)=[];
    end
    figure
    fs=10; lw=2;
    if na==4
        Yall=[sum(Y(:,1:n),2),sum(Y(:,n+1:2*n),2),sum(Y(:,2*n+1:3*n),2),sum(Y(:,3*n+1:end),2)];
        Hall=[sum(H(:,1:n),2),sum(H(:,n+1:2*n),2),sum(H(:,2*n+1:3*n),2),sum(H(:,3*n+1:end),2)];
        %Yall=Y(:,1:n)+Y(:,n+1:2*n)+Y(:,2*n+1:3*n)+Y(:,3*n+1:end);
    elseif na==5
        Yall=[sum(Y(:,1:n),2),sum(Y(:,n+1:2*n),2),sum(Y(:,2*n+1:3*n),2),sum(Y(:,3*n+1:4*n),2),sum(Y(:,4*n+1:end),2)];
    elseif na==3
        Yall=[sum(Y(:,1:n),2),sum(Y(:,n+1:2*n),2),sum(Y(:,2*n+1:3*n),2)];
        Hall=[sum(H(:,1:n),2),sum(H(:,n+1:2*n),2),sum(H(:,2*n+1:3*n),2)];
    elseif na==1
        Yall=sum(Y,2);
        Hall=sum(H,2);
    else
        error('Number of age groups not recognised for plotting')
    end
    %maxY=max(max(Yall));
    maxY=max(max(max(Yall)),max(max(Hall)));
    %
    %Unlogged plots:
    h=zeros(1,na);
    hold on
    for i=2:length(tvec)-1
        plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
    end
    for i=1:na
        %h(i)=plot(tout,Yall(:,i),'linewidth',lw,'color',cmap(i,:));
        h1=plot(tout,Yall(:,i),'linewidth',lw,'color',cmap(1,:));%(i,:)
        h2=plot(tout,Hall(:,i),'-','linewidth',lw,'color',cmap(2,:));%'--', (i,:)
    end
    %plot(tout,Yall,'linewidth',.5,'color',cmap(2,:));
    %plot(tout,sum(Yall,2),'linewidth',lw,'color','k');
    %}
    %{
    %Logged plots:
    hold on
    semilogy(tout,Yall);
    %}
    

    lt=length(tvec);
    points=tvec+10;
    pointsy=.93*maxY;
    txt={'1','2','3','4','5','6'};
    for i=3:lt-1
        text(points(i),pointsy,txt{i-2},'fontsize',20)
    end
    
    xlabel('Time (days since 1st Jan)','FontSize',fs);
    ylabel('Population','FontSize',fs);%yvar
    set(gca,'FontSize',fs);
    axis ([0,tend,0,maxY])
    legend([h1,h2],'Incidence','Hospitalisations','location','NW')
    %
    if na==4
        legend(h,{'0-4','5-19','20-64','65+'},'location','NE')
    elseif na==5
        legend(h,{'0-4','5-17','18-49','50-64','65+'},'location','NE')
    elseif na==3
        legend(h,{'0-15','16-64','65+'},'location','NE')
    end
    %}
    grid on
    grid minor
    box on
    hold off
%}
end
end

function f=integr8covid(t,y,pr,betaIn,nbar,NN0,D,phi1,phi2,seedvec,itime)
%phi=phi1-phi2*cos(pi*t/180);%Seasonality****
phi=phi1;
%%
S=y(1:nbar);
E=y(nbar+1:2*nbar);
Ia=y(2*nbar+1:3*nbar);
%{
Ip=y(3*nbar+1:4*nbar);
Inm=y(4*nbar+1:5*nbar);
Ism=y(5*nbar+1:6*nbar);
Ins=y(6*nbar+1:7*nbar);
Iss=y(7*nbar+1:8*nbar);
Qm=y(8*nbar+1:9*nbar);
Qs=y(9*nbar+1:10*nbar);
H=y(10*nbar+1:11*nbar);
%H2=y(11*nbar+1:12*nbar);
I=2/3*Ia+2/3*Ip+Inm+Ism+Ins+Iss;%All infectious
%DE=y(10*nbar+1:11*nbar);
%R=y(11*nbar+1:end);
%}
%
%IC:
Inm=y(3*nbar+1:4*nbar);
Ism=y(4*nbar+1:5*nbar);
Ins=y(5*nbar+1:6*nbar);
Iss=y(6*nbar+1:7*nbar);
Qm=y(7*nbar+1:8*nbar);
Qs=y(8*nbar+1:9*nbar);
H=y(9*nbar+1:10*nbar);
%H2=y(11*nbar+1:12*nbar);
I=pr.red*Ia+Inm+Ism+Ins+Iss;%All infectious
%DE=y(10*nbar+1:11*nbar);
%R=y(11*nbar+1:end);
%}
%%
seed1=seedvec.*S./NN0;
beta=betaIn;
%{
if t>topen && t<topen+30
    beta=betaIn*(.25+.75*(t-topen)/30);
else
    beta=betaIn;
end
%}
Sfoi=phi*(beta*S.*(D*(I./NN0))+seed1);
%%
Sdot=-Sfoi;
Edot=Sfoi-pr.sigma*E;
%{
%HE model:
Iadot=(1-pr.p1)*pr.sigma*E-(1+pr.odds)*pr.g1*Ia;%2**
Ipdot=pr.p1*pr.sigma*E-(1+pr.odds)*pr.omega*Ip;%2**
Inmdot=(1-pr.p2)*(1-pr.p3)*pr.omega.*Ip-pr.g2*Inm;
Ismdot=(1-pr.p2)*pr.p3*pr.omega.*Ip-(pr.g2+pr.q1)*Ism;
Insdot=pr.p2*(1-pr.p4)*pr.omega.*Ip-pr.h*Ins;%*(1-h).*Ins;
Issdot=pr.p2*pr.p4*pr.omega.*Ip-(pr.h+pr.q2)*Iss;%(1-h).*Iss;
Qmdot=pr.g1*pr.odds*Ia+(1-pr.p2)*pr.omega*pr.odds.*Ip+pr.q1*Ism-pr.g4*Qm;
Qsdot=pr.p2*pr.omega*pr.odds.*Ip+pr.q2*Iss-pr.h*Qs;
%}
%
%HE model IC:
Iadot=(1-pr.p1)*pr.sigma*E-(1+pr.odds)*pr.g1*Ia;%2**
%Ipdot=pr.p1*pr.sigma*E-(1+pr.odds)*pr.omega*Ip;%2**
Inmdot=pr.p1*pr.sigma*(1-pr.p2)*(1-pr.p3)*E-pr.g2*Inm;
Ismdot=pr.p1*pr.sigma*(1-pr.p2)*pr.p3*E-(pr.g2+pr.q1)*Ism;
Insdot=pr.p1*pr.sigma*pr.p2*(1-pr.p4)*E-(pr.h+pr.gX)*Ins;%*(1-h).*Ins;
Issdot=pr.p1*pr.sigma*pr.p2*pr.p4*E-(pr.h+pr.q2+pr.gX)*Iss;%(1-h).*Iss;
Qmdot=pr.g1*pr.odds*Ia+pr.q1*Ism-pr.g4*Qm;
Qsdot=pr.q2*Iss-pr.h*Qs;
%}
Hdot=pr.h*(Ins+Iss+Qs)-pr.g3*H;
DEdot=pr.g3*pr.mu.*H;
Rdot=pr.g1*Ia+pr.g2*(Inm+Ism)+pr.g3*(1-pr.mu).*H+pr.g4*Qm+pr.gX*(Ins+Iss);
%f=[Sdot;Edot;Iadot;Ipdot;Inmdot;Ismdot;Insdot;Issdot;Qmdot;Qsdot;Hdot;DEdot;Rdot];
f=[Sdot;Edot;Iadot;Inmdot;Ismdot;Insdot;Issdot;Qmdot;Qsdot;Hdot;DEdot;Rdot];
end

%Stochastic variant - needs update to C19 flowchart
function f=stochSim(y,beta,gamma,n,nbar,NN,N0,D,seed,phi1,phi2,tau,alpha)
%Still flu-like/SIR****
%Feed in mu if required
factor=6;
tend=360*factor; beta=beta/factor; gamma=gamma/factor;
Vec=zeros(nbar,tend);
%
S=y(1:nbar);
I=y(nbar+1:2*nbar);
R=y(2*nbar+1:end);
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
f=R;
if sum(isnan(R))>0
    fuck=1;
    print('Somenting is NaN')
end
end

%{
Iadot=(1-pr.p1)*pr.sigma*E-pr.g1*Ia;
Ipdot=pr.p1*pr.sigma*E-pr.omega*Ip;
Inmdot=(1-pr.p2)*(1-pr.p3)*pr.omega.*Ip-pr.g2*Inm;
Ismdot=(1-pr.p2)*pr.p3*pr.omega.*Ip-(pr.g2+pr.q1)*Ism;
Insdot=pr.p2*(1-pr.p4)*pr.omega.*Ip-pr.h*Ins;%*(1-h).*Ins;
Issdot=pr.p2*pr.p4*pr.omega.*Ip-(pr.h+pr.q2)*Iss;%(1-h).*Iss;
Qmdot=pr.q1*Ism-pr.g4*Qm;
Qsdot=pr.q2*Iss-pr.h*Qs;
%}
%{
%With quarantine of all infs:
Iadot=(1-pr.p1)*pr.sigma*E-(pr.g1+pr.q1*pr.z)*Ia;
Ipdot=pr.p1*pr.sigma*E-(pr.omega+pr.q1*pr.z)*Ip;%q1=q2****
Inmdot=(1-pr.p2)*(1-pr.p3)*pr.omega.*Ip-(pr.g2+pr.q1*pr.z)*Inm;
Ismdot=(1-pr.p2)*pr.p3*pr.omega.*Ip-(pr.g2+pr.q1)*Ism;
Insdot=pr.p2*(1-pr.p4)*pr.omega.*Ip-(pr.h+pr.q2*pr.z)*Ins;%*(1-h).*Ins;
Issdot=pr.p2*pr.p4*pr.omega.*Ip-(pr.h+pr.q2)*Iss;%(1-h).*Iss;
Qmdot=pr.q1*(Ia*pr.z+(1-pr.p2)*pr.z.*Ip+Inm*pr.z+Ism)-pr.g4*Qm;
Qsdot=pr.q2*(pr.p2.*Ip*pr.z+Ins*pr.z+Iss)-pr.h*Qs;
%}