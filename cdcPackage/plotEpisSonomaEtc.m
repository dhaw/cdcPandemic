function f=plotEpisSonomaEtc(f1,g1,f2,g2,f3,g3,f4,g4,f5,g5,f6,g6,f7,g7,f8,g8,f9,g9,f10,g10,f11,g11)%,data)

td1=197.4349;%233.7842;
tdreal=76;%79;
tsip=tdreal+700;%43;%38;
titstr='R_0=2.5';
thresh1=916;
thresh2=1416;

figure
fs=10; lw=2;
cmap=parula(11);%lines(7);
col1=cmap(1,:);
red=75;
%maxy=max([g1;g2;g3]);
maxy=100;
hold on
plot([1,f1(end)],[thresh1,thresh1],'--','linewidth',1,'color',[.5,.5,.5])
plot([1,f1(end)],[thresh2,thresh2],'--','linewidth',1,'color',[.5,.5,.5])

tvec=[tdreal,tsip,tsip,tdreal];
inBetween=[0,0,maxy,maxy];
fill(tvec,inBetween,[.5,.5,.5],'facealpha',.2);
plot([tsip,tsip],[0,maxy],'-','linewidth',1,'color',[.5,.5,.5])
plot([tdreal,tdreal],[0,maxy],'-','linewidth',1,'color',[.5,.5,.5])
plot([tsip,tdreal],[0,0],'-','linewidth',1,'color',[.5,.5,.5])
plot([tsip,tdreal],[maxy,maxy],'-','linewidth',1,'color',[.5,.5,.5])
%
h11=plot(f11-td1+tdreal,g11,'-','linewidth',lw,'color',cmap(11,:));
h10=plot(f10-td1+tdreal,g10,'-','linewidth',lw,'color',cmap(10,:));
h9=plot(f9-td1+tdreal,g9,'-','linewidth',lw,'color',cmap(9,:));
h8=plot(f8-td1+tdreal,g8,'-','linewidth',lw,'color',cmap(8,:));
h7=plot(f7-td1+tdreal,g7,'-','linewidth',lw,'color',cmap(7,:));
h6=plot(f6-td1+tdreal,g6,'-','linewidth',lw,'color',cmap(6,:));
h5=plot(f5-td1+tdreal,g5,'-','linewidth',lw,'color',cmap(5,:));
h4=plot(f4-td1+tdreal,g4,'-','linewidth',lw,'color',cmap(4,:));
h3=plot(f3-td1+tdreal,g3,'-','linewidth',lw,'color',cmap(3,:));
h2=plot(f2-td1+tdreal,g2,'-','linewidth',lw,'color',cmap(2,:));
h1=plot(f1-td1+tdreal,g1,'-','linewidth',lw,'color',cmap(1,:));

%h4=plot(61:61+length(data)-1,data,'k-','linewidth',lw);
%}
%{
h4=plot(f3-td2+tdreal,g3,'-','linewidth',lw,'color',col2);
h5=plot(f4-td2+tdreal,g4,'--','linewidth',lw,'color',col2);
h6=plot(f6-td2+tdreal,g6,':','linewidth',lw,'color',col2);
%}
set(gca,'fontsize',fs)
axis([0,300,0,maxy])%([60,tsip,0,50])%([0,720,0,maxy])
xlabel('Days (from 1st Jan)')
ylabel('Hospitalisations')%(strcat('Hosp. (',num2str(red),'% reduction)'))
title(titstr)
legend([h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11],'75%','70%','65%','60%','55%','50%','45%','40%','35%','30%','25%','location','NW')%,'Unmitigated','P_1','P_2'
%legend([h1,h2,h3],'75%','50%','25%','location','NW')%h4, 'Data',
grid on
grid minor
box on
