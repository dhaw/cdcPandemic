function [hosp,death]=processHospitalData(H,n,na)
%H=matrix: named rows and age category columns
%4 age groups:
n1=ones(n,1);

if na==4
    %Fluscape 4 groups:
    Hout=[H(:,1),mean(H(:,2:4),2),mean(H(:,5:13),2),mean(H(:,14:end),2)];
elseif na==3
    %Sonoma 3 groups:
    Hout=[mean(H(:,1:4),2),mean(H(:,5:13),2),mean(H(:,14:end),2)];%Proxy - 14/15
end
hosp=kron(Hout(11,:)',n1);
%h2=kron(Hout(13,:)',n1);
%f=[h1,h2];
Dout=Hout(13,:).*Hout(14,:)+(1-Hout(13,:)).*Hout(15,:);
death=kron(Dout',n1);