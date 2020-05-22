function vaxparams=prepVax(data)
%Columns 1 to 6 as table
tab=data(data.Season=='0910',:);
tab=unstack(tab,{'VC','VE'},'k');
%Manually extract from october:
fromInd=6;%Include one month of zero vax
toInd=13;
vRate=table2array(tab(fromInd:toInd,4:7));
vEff=table2array(tab(fromInd:toInd,8:11));

%Columns are age groups
vRate=vRate';
vEff=vEff';
vaxparams=zeros(size(vEff,1),size(vEff,2),2);
vaxparams(:,:,1)=vRate;
vaxparams(:,:,2)=vEff;
vaxparams=cat(1,vaxparams(1:3,:,:),vaxparams(3,:,:),vaxparams(4,:,:));%Approximate 5 age-groups

%October - 273
%November - 304
%December - 334
%January - 365
%February - 396
%March - 424
%April - 455
%numdays=[31,30,31,31,28,31,30]';
%vRate=repmat(kron