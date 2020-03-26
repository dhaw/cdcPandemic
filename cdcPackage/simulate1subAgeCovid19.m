function f=simulate1subAgeCovid19(beta,p,sigma,omega,q,gamma,h,mu,tvec,Dvec)
solvetype=2;
if solvetype==2
    lt=length(tvec);
    t0=tvec(1);%tvec(1)=start time for while simulation
    y0=[S0;zn;NNbar-S0];
    for i=1:lt
        D=Dvec(:,:,i);
        tend=tvec(i+1);
        soFar=integr8(t0,tend,y0,beta,params,D);
        
        t0=tend;
    end
elseif solvetype==1
    error('Final size calculations not possible')
elseif solvetype==3
    error('Code not written yet')
    %f=stochSim(y0,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,tau,alpha);
end
end

function f=integr8(plotTau)
    [tout,yout]=ode45(@(t,y)integr8covid(t,y,beta,gamma,nbar,NN0,D,phi1,phi2,seedvec);
    if plotTau>0
        plotEpi(tout,yout,plotTau);
    end
end


function f=plotEpi(tout,yout,plotTau)
solvetype=2;
if solvetype==2
    %[tout,yout]=ode45(@(t,y)integr8covid(t,y,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,alpha,seedvec),[t0,tend],y0);
    %Incidence curve in here:
    %
    if tau==plotTau
        figure
        fs=12; lw=2;
        Y=yout(:,nbar+1:2*nbar); %Z=yout(:,2*nbar+1:3*nbar);
        %Y=sum(Y,2);
        %Y1=Y(:,minNind)+Y(:,minNind+n)+Y(:,minNind+2*n)+Y(:,minNind+3*n); Y1=Y1/NN(minNind);
        %Y2=Y(:,maxNind)+Y(:,maxNind+n)+Y(:,maxNind+2*n)+Y(:,maxNind+3*n); Y2=Y2/NN(maxNind);
        Ysum=sum(Y,2);
        Yall=Y(:,1:n)+Y(:,n+1:2*n)+Y(:,2*n+1:3*n)+Y(:,3*n+1:end);
        Yall=abs(Yall);%****cheat ;)
        %
        %Unlogged plots:
        hold on
        %plot(tout,Y1,'--','linewidth',lw,'color',[.165,.31,.431]);%[.165,.31,.431][.447,.553,.647]
        %plot(tout,Y2,'-','linewidth',lw,'color',[.165,.31,.431]);
        plot(tout,Ysum,'k','linewidth',lw);
        plot(tout,Yall);
        %}
        %{
        %Logged plots:
        semilogy(tout,Ysum,'k','linewidth',lw);
        hold on
        semilogy(tout,Yall);
        %}
        xlabel('Time (days)','FontSize',fs);
        ylabel('Prevalence','FontSize',fs);
        set(gca,'FontSize',fs);
        maxY=max(Ysum);
        %axis([0,tend,0,maxY])
        axis ([0,tend,0,maxY])
        %legend('Min','Max','location','NE')
        grid on
        grid minor
        hold off
    end
    %}
    rec0=yout(end,2*nbar+1:end)+yout(end,nbar+1:2*nbar);%******** Truncation - add infectious to rercovered?
    f=rec0';
%}
end
end


function f=integr8covid(t,y,beta,gamma,nbar,NN0,D,phi1,phi2,seedvec)
%phi=phi1-phi2*cos(pi*t/180);%Seasonality****
phi=phi1;
%%
S=y(1:nbar);
E=y(nbar+1:2*nbar);
Ia=y(2*nbar+1:3*nbar);
Ip=y(3*nbar+1:4*nbar);
Iny(4*nbar+1:5*nbar);
Isy(5*nbar+1:6*nbar);
Q=y(6*nbar+1:7*nbar);
Hb=y(7*nbar+1:8*nbar);
Hv=y(8*nbar+1:9*nbar);
Ho=y(9*nbar+1:10*nbar);
DE=y(10*nbar+1:11*nbar);
%R=y(11*nbar+1:end);
%%
seed1=seedvec.*S./NN0;
Sfoi=phi*(beta*S.*(D*(I./NN0))+seed1);
%%
Sdot=-Sfoi;
Edot=Sfoi-sigma*E;
Iadot=(1-p(1))*sigma*E-gamma(1)*Ia;
Ipdot=p(1)*sigma*E-gamma(1)*Ip-omega*Ip;
Indot=(1-p(2))*omega*Ip-(gamma(2)+q(1))*In;
Isdot=p(2)*omega*Ip-(gamma(4)+q(2)+h(1)+mu(1))*Is;
Qdot=q(1)*In+q(2)*Is-(gamma(3)+h(2))*Q;
Hbdot=h(1)*Is+h(2)*Q-(gamma(5)+h(3)+mu(2))*Hb;
Hvdot=h(3)*Hb-(h(4)+mu(3))*Hv;
Ho=h(4)*Hv-gamma(6)*Ho;
DEdot=mu(1)*Is+mu(2)*Hb+mu(3)*Hv;
Rdot=gamma(1)*Ia_gamma(2)*In_gamma(3)*Q+gamma(4)*Is+gamma(5)*Hb+gamma(6)*Ho;
f=[Sdot;Edot;Iadot;Ipdot;Indot;Isdot;Qdot;Hbdot;Hvdot;Ho;DEdot;Rdot];
end


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