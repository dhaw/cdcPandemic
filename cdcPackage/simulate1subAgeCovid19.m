function [DEout,Rout]=simulate1subAgeCovid19(sigma,omega,gamma,hosp,mu,pvec,qvec,beta,tvec,Dvec,n,nbar,NNbar,NN0,phi1,phi2,seedvec,S0,tau)
solvetype=2;
plotTau=1;%Plot this year
zn=zeros(nbar,1);
if solvetype==2
    lt=length(tvec);
    t0=tvec(1);%tvec(1)=start time for while simulation
    y0=[S0;repmat(zn,11,1);NNbar-S0];
    toutAll=[];
    Sout=[];
    DEout=zeros(nbar,lt);
    Rout=DEout;
    for i=1:lt-1
        D=Dvec(:,:,i);
        tend=tvec(i+1);
        [tout,Sclass,DEcum,Rcum,y0]=integr8(sigma,omega,gamma,hosp,mu,pvec,qvec,beta,nbar,NN0,D,phi1,phi2,seedvec,t0,tend,y0);
        toutAll=[toutAll;tout];
        Sout=[Sout;Sclass];
        DEout(:,i)=DEcum;
        Rout(:,i)=Rcum;
        if tau==plotTau
            plotEpi(tout,Sout,n);
        end
        t0=tend;
    end
elseif solvetype==1
    error('Final size calculations not possible')
elseif solvetype==3
    error('Code not written yet')
    %f=stochSim(y0,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,tau,alpha);
end
end

function [tout,Sclass,DEcum,Rcum,y0new]=integr8(sigma,omega,gamma,hosp,mu,pvec,qvec,beta,nbar,NN0,D,phi1,phi2,seedvec,t0,tend,y0)
%ncomps=13;%Number of compartments
    [tout,yout]=ode45(@(t,y)integr8covid(t,y,sigma,omega,gamma,hosp,mu,pvec,qvec,beta,nbar,NN0,D,phi1,phi2,seedvec),[t0,tend],y0);
    Sclass=yout(:,1:nbar);%Incidence
    DEcum=yout(end,11*nbar+1:12*nbar);%Deaths
    Rcum=yout(end,12*nbar+1:end);
    y0new=yout(end,:)';
end


function f=plotEpi(tout,Y,n)
solvetype=2;
tend=360;%For plot only
na=size(Y,2)/n;
if solvetype==2
    Y=-diff(Y,1);
    tdiff=diff(tout,1);
    Y=Y./tdiff;
    tout(1)=[];
    figure
    fs=12; lw=2;
    if na==4
        Yall=[sum(Y(:,1:n),2),sum(Y(:,n+1:2*n),2),sum(Y(:,2*n+1:3*n),2),sum(Y(:,3*n+1:end),2)];
    elseif na==5
        Yall=[sum(Y(:,1:n),2),sum(Y(:,n+1:2*n),2),sum(Y(:,2*n+1:3*n),2),sum(Y(:,3*n+1:4*n),2),sum(Y(:,4*n+1:end),2)];
    else
        error('Number of age groups not recognised for plotting')
    end
    %
    %Unlogged plots:
    hold on
    plot(tout,Yall,'linewidth',lw);
    %}
    %{
    %Logged plots:
    hold on
    semilogy(tout,Yall);
    %}
    xlabel('Time (days)','FontSize',fs);
    ylabel('Incidence','FontSize',fs);
    set(gca,'FontSize',fs);
    maxY=max(max(Yall));
    axis ([0,tend,0,maxY])
    if na==4
        legend('0-4','5-19','20-64','65+','location','NE')
    elseif na==5
        legend('0-4','5-17','18-49','50-64','65+','location','NE')
    end
    grid on
    grid minor
    box on
    hold off
%}
end
end


function f=integr8covid(t,y,sigma,omega,gamma,h,mu,p,q,beta,nbar,NN0,D,phi1,phi2,seedvec)
%phi=phi1-phi2*cos(pi*t/180);%Seasonality****
phi=phi1;
%%
S=y(1:nbar);
E=y(nbar+1:2*nbar);
Ia=y(2*nbar+1:3*nbar);
Ip=y(3*nbar+1:4*nbar);
Inm=y(4*nbar+1:5*nbar);
Ism=y(5*nbar+1:6*nbar);
Ins=y(6*nbar+1:7*nbar);
Iss=y(7*nbar+1:8*nbar);
Qm=y(8*nbar+1:9*nbar);
Qs=y(9*nbar+1:10*nbar);
H=y(10*nbar+1:11*nbar);
%H2=y(11*nbar+1:12*nbar);
I=Ia+Ip+Inm+Ism+Ins+Iss;%All infectious
%DE=y(10*nbar+1:11*nbar);
%R=y(11*nbar+1:end);
%%
%if t<30
    seed1=seedvec.*S./NN0;
%else
    %seed1=0;
%end
Sfoi=phi*(beta*S.*(D*(I./NN0))+seed1);
%%
Sdot=-Sfoi;
Edot=Sfoi-sigma*E;
Iadot=(1-p(1))*sigma*E-gamma(1)*Ia;
Ipdot=p(1)*sigma*E-omega*Ip;
Inmdot=(1-p(2))*(1-p(3))*omega*Ip-gamma(2)*Inm;
Ismdot=(1-p(2))*p(3)*omega*Ip-(gamma(2)+q(1))*Ism;
Insdot=p(2)*(1-p(4))*omega*Ip-gamma(2)*Ins;%*(1-h).*Ins;
Issdot=p(2)*p(4)*omega*Ip-(gamma(2)+q(2))*Iss;%(1-h).*Iss;
Qmdot=q(1)*Ism-gamma(4)*Qm;
Qsdot=q(2)*Iss-gamma(5)*Qs;
Hdot=gamma(2)*h.*(Ins+Iss)+gamma(5)*h.*Qs-gamma(3)*H-mu.*H;
%H2dot=h(:,2).*Hb-(mu(3))*Hv-h(:,3).*Hv;
%Ho=h(:,3).*Hv-gamma(6)*Ho;
DEdot=mu.*H;
%Rdot=gamma(1)*Ia+gamma(2)*(Inm+Ism+Ins+Iss)+gamma(3)*H+gamma(4)*Qm+gamma(5)*Qs;
Rdot=gamma(1)*Ia+gamma(2)*(Inm+Ism)+gamma(2)*(1-h).*Ins+gamma(2)*(1-h).*Iss+gamma(3)*H+gamma(4)*Qm+gamma(5)*(1-h).*Qs;
f=[Sdot;Edot;Iadot;Ipdot;Inmdot;Insdot;Ismdot;Issdot;Qmdot;Qsdot;Hdot;DEdot;Rdot];
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