function f=heMakeDs(NN,Z)
baseRate=1;
ln=length(NN);
NNsum=sum(NN);
NNrel=NN./NNsum;%Column
D=repmat(NNrel',ln,1)*baseRate;
D(1:ln-1,1:ln-1)=D(1:ln-1,1:ln-1)+Z;
f=D;
