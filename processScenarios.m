function f=processScenarios(F,NNbar,tI,tH,tD,tV)
if length(F)==8
    scenarioType='RL';
elseif length(F)==24
    scenarioType='HP';
else
    error('How many scenarios are there?')
end
numScen=length(F);
for i=1:numScen
    Fi=F{i};%ndays*na*nruns
    [ndays,na,nruns]=size(Fi);
    na=na/4;%****
    nweeks=floor(ndays/7);
    NNrep=repmat([sum(NNbar),NNbar],nweeks,1);
    Iout=zeros(nweeks,na+1,nruns);
    %nweeks redefined soon
    
    Hout=Iout;
    Dout=Iout;
    Vout=Iout;
    FI=Fi(:,1:na,:);
    FH=Fi(:,na+1:2*na,:);
    FD=Fi(:,2*na+1:3*na,:);
    
    %Deaths to non-cumulative:
    %FD=[zeros(1,na,nruns);diff(FD,1)];
    
    FV=Fi(:,end-na+1:end,:);%Should've separated after rearranging. Never mind. 
    for j=1:nruns%Change
        %{
        FIj=FI(:,:,j);
        FIj=[sum(FIj,2),FIj];
        FIj=squeeze(sum(reshape(FIj',na+1,7,[]),2))';
        Iout(:,:,j)=FIj./NNrep*100;
        bins=[0:.1:19.9,100];
        binsCum=[0:.2:39.8,100];
        tHold=tI;
        scentype='SymIllness';
        %}
        %{
        FHj=FH(:,:,j);
        FHj=[sum(FHj,2),FHj];
        FHj=squeeze(sum(reshape(FHj',na+1,7,[]),2))';
        Iout(:,:,j)=FHj./NNrep*1000000;
        bins=[0:.25:149.75,1000];
        binsCum=[0:1:599,1000];
        tHold=tH;
        scentype='Hosp';
        %}
        %
        FDj=FD(:,:,j);
        FDj=[sum(FDj,2),FDj];
        FDj=squeeze(sum(reshape(FDj',na+1,7,[]),2))';
        Iout(:,:,j)=FDj./NNrep*100000;
        bins=[0:.025:14.975,100];
        binsCum=[0:.1:59.9,100];
        tHold=tD;
        scentype='Death';
        %}
        %{
        FVj=FV(:,:,j);
        FVj=[sum(FVj,2),FVj];
        FVj=squeeze(sum(reshape(FVj',na+1,7,[]),2))';
        Iout(:,:,j)=FVj./NNrep*100;
        bins=[0:.01:3.99,100];
        binsCum=[0:.02:7.98,100];
        tHold=tV;
        scentype='AntiviralTX';
        %}
    end
    nweeks=size(Iout,1);
    nruns=size(Iout,3);
    %%
    %Incidence:
    %Non-cumulative:
    nbins=length(bins)-1;
    Ibinned=zeros(nbins,nweeks+2,na+1);%Cum, peak, W1-W52
    for w=1:nweeks
        for j=1:na+1
           hij=histogram(Iout(w,j,:),bins);
           Ibinned(:,w+2,j)=hij.Values;
        end
    end
    %Cumulative:*****
    bins2=binsCum;
    for a=1:na+1
        ha=cumsum(Iout(:,a,:),1);
        ha=histogram(ha,bins2);
        Ibinned(:,1,a)=ha.Values;
    end
    %Peaks:
    peaks=zeros(nruns,na);
    pweek=zeros(nweeks,na+1);
    for j=1:nruns
        for a=1:na+1
            [pval,pw]=max(Iout(:,a,j));
            peaks(j,a)=pval;
            pweek(pw,a)=pweek(pw,a)+1;
        end
    end
    psum=sum(pweek,1);
    pweek=pweek./repmat(psum,nweeks,1);
    for a=1:na+1
        ha=histogram(peaks(:,a),bins);
        Ibinned(:,2,a)=ha.Values;
    end
    Isum=sum(Ibinned,1);
    Isum=repmat(Isum,nbins,1,1);
    Ibinned=Ibinned./Isum;
    Ibinned(end+1:end+9,:,:)=0;
    %Summary stats:
    means=mean(Iout,3);
    prctiles=prctile(Iout,[50,2.5,5,25,75,95,97.5],3);
    for a=1:na+1
        Ibinned(end-8,3:end,a)=means(:,a);
        Ibinned(end-7:end-1,3:end,a)=permute(prctiles(:,a,:),[3,1,2]);
        Ibinned(end,3:end,a)=pweek(:,a)/sum(pweek(:,a));
        Ibinned(end,1:2,a)=[NaN,NaN];
        %mean/prctiles for cum/peak:
        Isum=sum(Iout(:,a,:),1);
        Ibinned(end-8,1,a)=nanmean(Isum,3);
        Ibinned(end-7:end-1,1,a)=prctile(Isum,[50,2.5,5,25,75,95,97.5],3);
        Ibinned(end-8,2,a)=nanmean(peaks(:,a));%pweek to peaks
        Ibinned(end-7:end-1,2,a)=prctile(peaks(:,a),[50,2.5,5,25,75,95,97.5]);
    end
    Ibinned=reshape(permute(Ibinned,[1,3,2]),[(na+1)*(nbins+9),nweeks+2,1]);
    Ibinned=array2table(Ibinned);
    %
    tHold(:,4)=Ibinned(:,1);
    tHold(:,6:6+nweeks)=Ibinned(:,2:end);
    
    %tHold(:,6)=Ibinned(:,2);
    %tHold(:,7:7+33)=array2table(zeros((na+1)*(nbins+9),34));%Change for 2009
    %tHold(:,7+34:end)=Ibinned(:,3:20);
    tHold(:,7+nweeks:end)=array2table(zeros((na+1)*(nbins+9),58-nweeks-6));
    
    if i<10
        sceno=strcat('0',num2str(i));
    else
        sceno=num2str(i);
    end
    filename=strcat(scenarioType,sceno,scentype,'IMP.csv');
    writetable(tHold,filename)
    %%
    %{
    %Hosp:
    bins=[0:.1:59.9,100];
    bins2=2*bins;
    
    %Death:
    bins=[0:.01:5.99,100];
    bins2=2*bins;
    
    %Anti-v:
    bins=[0:.1:39.9,100];
    bins2=2*bins;
    %}
end