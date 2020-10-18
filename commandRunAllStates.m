[Z2ca,realZ2ca]=subWave2Z1(NNbarCA',xca,ycax,xsto,vaxparams);
[Z2co,realZ2co]=subWave2Z1(NNbarCO',xco,ycox,xsto,vaxparams);
[Z2ct,realZ2ct]=subWave2Z1(NNbarCT',xct,yctx,xsto,vaxparams);
[Z2ga,realZ2ga]=subWave2Z1(NNbarGA',xga,ygax,xsto,vaxparams);
[Z2md,realZ2md]=subWave2Z1(NNbarMD',xmd,ymdx,xsto,vaxparams);
[Z2mn,realZ2mn]=subWave2Z1(NNbarMN',xmn,ymnx,xsto,vaxparams);
[Z2nm,realZ2nm]=subWave2Z1(NNbarNMnew',xnmNew,ynmNewx,xsto,vaxparams);
[Z2ny,realZ2ny]=subWave2Z1(NNbarNY',xny,ynyx,xsto,vaxparams);
[Z2or,realZ2or]=subWave2Z1(NNbarOR',xor,yorx,xsto,vaxparams);
[Z2tn,realZ2tn]=subWave2Z1(NNbarTN',xtn,ytnx,xsto,vaxparams);

realZ2all=[realZ2ca,realZ2co,realZ2ct,realZ2ga,realZ2md,realZ2mn,realZ2nm,realZ2ny,realZ2or,realZ2tn];cellZ2={Z2ca,Z2co,Z2ct,Z2ga,Z2md,Z2mn,Z2nm,Z2ny,Z2or,Z2tn};
plotz2(cellZ2,realZ2all);

%[Z2mn,realZ2mn]=subWave2Z1(NNbarMN',xmn,ymnx,xsto);

[Z2ca,realZ2ca]=subWave2Z1(NNbarCA',xca,yca,ycax,X1us30,vaxparams);
[Z2ct,realZ2ct]=subWave2Z1(NNbarCT',xct,yct,yctx,X1us30,vaxparams);
[Z2co,realZ2co]=subWave2Z1(NNbarCO',xco,yco,ycox,X1us30,vaxparams);
[Z2ga,realZ2ga]=subWave2Z1(NNbarGA',xga,yga,ygax,X1us30,vaxparams);
[Z2md,realZ2md]=subWave2Z1(NNbarMD',xmd,ymd,ymdx,X1us30,vaxparams);
[Z2ny,realZ2ny]=subWave2Z1(NNbarNY',xny,yny,ynyx,X1us30,vaxparams);
[Z2or,realZ2or]=subWave2Z1(NNbarOR',xor,yor,yorx,X1us30,vaxparams);
[Z2tn,realZ2tn]=subWave2Z1(NNbarTN',xtn,ytn,ytnx,X1us30,vaxparams);

[Z2ca,realZ2ca]=subWave2Z1(NNbarCA',xca,yca,ycax,x0,vaxparams);
[Z2ct,realZ2ct]=subWave2Z1(NNbarCT',xct,yct,yctx,x0,vaxparams);
[Z2co,realZ2co]=subWave2Z1(NNbarCO',xco,yco,ycox,x0,vaxparams);
[Z2ga,realZ2ga]=subWave2Z1(NNbarGA',xga,yga,ygax,x0,vaxparams);
[Z2md,realZ2md]=subWave2Z1(NNbarMD',xmd,ymd,ymdx,x0,vaxparams);
[Z2mn,realZ2mn]=subWave2Z1(NNbarMN',xmn,ymn,ymnx,x0,vaxparams);
[Z2nm,realZ2nm]=subWave2Z1(NNbarNMnew',xnmNew,ynmNew,ynmNewx,x0,vaxparams);
[Z2ny,realZ2ny]=subWave2Z1(NNbarNY',xny,yny,ynyx,x0,vaxparams);
[Z2or,realZ2or]=subWave2Z1(NNbarOR',xor,yor,yorx,x0,vaxparams);
[Z2tn,realZ2tn]=subWave2Z1(NNbarTN',xtn,ytn,ytnx,x0,vaxparams);