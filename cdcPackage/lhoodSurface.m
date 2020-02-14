function f=lhoodSurface(ydataNX)
R0=(1.3:.01:1.4);
lr=length(R0);
gamma=(.3:.01:.4);
lg=length(gamma);
z=zeros(lr,lg);
parfor i=1:lr
    ri=R0(i);
    zi=zeros(lg,1);
    gi=gamma;
    for j=1:lg
        zi(j)=-cdcLhoodsW5([ri,gi(j)],ydataNX);
    end
    z(i,:)=zi;
end
f=z;
figure
surf(R0,gamma,z');