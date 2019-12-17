function [c,z]=cdcDecisionSeverityW5(xsto)%ydata,ydataNX)

x=(1:.05:2);%xi
lx=length(x);
chainLength=10000;%Match cdcMCMC
burn=1000;
inc=100;
c=cell(lx,1);
z=zeros(lx,floor((chainLength-burn)/inc),2);
%ydatai=ydata/x(i);
%ydataNXi=ydayaNX/x(i);
%[xsto,~,~,~]=cdcMCMC(ydataNX);%m, x22
for i=1:lx
    [cdf,Nvec]=cdcCDFW5(xsto,x(i));
    c{i}=[cdf,Nvec];
    zhold=zeros(floor((chainLength-burn)/inc),2);
    qmean=1./[393.0256 999.3054 702.6867 406.0680 177.1958]*x(i);
    for j=1:floor((chainLength-burn)/inc)%burn+1:inc:length(xsto)
        jj=inc*(j-1)+1;
        [~,~,z1]=cdcPandemicSimulationW5(xsto(jj,:),4:7,0,0,0,243);
        z1=sum(z1.*qmean);
        [~,~,z2]=cdcPandemicSimulationW5(xsto(jj,:),4:7,0,0,0,1000);
        z2=sum(z2.*qmean);
        zhold(j,1)=z1;
        zhold(j,2)=z2;
    end
    
    z(i,:,:)=zhold;
end
%}

%{
x=(.8:.05:1.2);%xi
lx=length(x);
burn=1000;
inc=100;
chainLength=100000;%Match cdcMCMC
c=cell(lx,1);
z=zeros(lx,floor((chainLength-burn)/inc),2);
for i=1:lx
    %ydatai=ydata/x(i);
    %ydataNXi=ydayaNX/x(i);
    theta3=cdcLSQfit(ydata/x(i));%Not R0
    [xsto,~,~,~]=cdcMCMC(ydataNX,theta3([1,2]),x(i));%m, x22
    [cdf,Nvec]=cdcCDF(xsto,theta3,x(i));
    c{i}=[cdf,Nvec];
    zhold=zeros(floor((chainLength-burn)/inc),2);
    qmean=1./[143.4400,364.7100,148.2000,64.6700]/x(i);
    for j=1:floor((chainLength-burn)/inc)%burn+1:inc:length(xsto)
        jj=inc*(j-1)+1;
        [~,~,z1]=cdcPandemicSimulation2(xsto(jj,:),4:7,0,0,0,243,theta3(1),theta3(2));
        z1=sum(z1.*qmean);
        [~,~,z2]=cdcPandemicSimulation2(xsto(jj,:),4:7,0,0,0,1000,theta3(1),theta3(2));
        z2=sum(z2.*qmean);
        zhold(j,1)=z1;
        zhold(j,2)=z2;
    end
    z(i,:,:)=zhold;
end
%}