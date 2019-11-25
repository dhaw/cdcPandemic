function f=plot2dlhood(yvec)

p1=2.5:.01:3;
p2=.5:.01:1;
lp1=length(p1); lp2=length(p2);
X=zeros(lp1,lp2);
for i=1:lp1
    for j=1:lp2
        X(i,j)=lhoods_david([p1(i),p2(j)],yvec);
    end
end
fs=12;
figure
imagesc(p1,p2,X)
set(gca,'fontsize',fs)