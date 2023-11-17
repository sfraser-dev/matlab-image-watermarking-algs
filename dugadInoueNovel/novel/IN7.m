function IN7();

L_BERStore=[];
L_NCStore=[];
H_BERStore=[];
H_NCStore=[];
D_BERStore=[];
D_NCStore=[];
V_BERStore=[];
V_NCStore=[];
Comb_BERStore=[];
Comb_NCStore=[];
totalWmLenEmbeddedStore=[];
wmLenRecStore=[];
PSNRStore=[];
wPSNRStore=[];
TPEStore=[];
noOfIts=1;	% number of iterations

for itcount=1:noOfIts,
	% Debugging 
	global PRINT_MINI_INFO
	PRINT_MINI_INFO=0; 	% 0=don't print; 1=do print

	% Which subbands to use?
	DO_L=0;			% 0=don't use; 1=use
	DO_H=1; 		% 0=don't use; 1=use
	DO_D=0;			% 0=don't use; 1=use
	DO_V=1;			% 0=don't use; 1=use

	% Get input image and its size
	I=double(imread('lena_256.tif'));
	[rows,cols]=size(I);

	% Determine dyadlength of input image (for FWT2_PO) and
	% calculate the amount of wavelet levels taken for a 
	% requested coarselevel
	[len,po2]=dyadlength(I(:,1));
	cl=5;
	%qmf=MakeONFilter('Daubechies',20);	% 4 6 8 10 12 14 16 18 20
	qmf=MakeONFilter('Daubechies',8);
	%qmf=MakeONFilter('Symmlet',10);	% 4 5 6 7 8 9 10
	%qmf=MakeONFilter('Symmlet',4);	
	%qmf=MakeONFilter('Coiflet',1);		% 1 2 3 4 5
	%qmf=MakeONFilter('Battle',1);		% 1 3 5
	%qmf=MakeONFilter('Beylkin');		% No par
	%qmf=MakeONFilter('Haar');		% No par
	%qmf=MakeONFilter('Vaidyanathan');	% No par
	
	numberWavLev=po2-cl;

	% Determine the size of the LL band (DC/low freq component)
	DCsize=cols;
	for i=1:numberWavLev,
	        DCsize=DCsize/2;
	end; clear i;

	% Take the forward wavelet transform
	fwt=FWT2_PO(I,cl,qmf);   

	% Get the coarsest subands (excluding LL)
	if DO_L==1, L=fwt(1:DCsize,1:DCsize); end			% dc part
	if DO_H==1, H=fwt(1:DCsize,DCsize+1:DCsize*2); end		% horiz
	if DO_D==1, D=fwt(DCsize+1:DCsize*2,DCsize+1:DCsize*2); end	% diag
	if DO_V==1, V=fwt(DCsize+1:DCsize*2,1:DCsize); end		% vert
	
	% Set the two quantization thresholds
	% Rough estimates for using LL subband to embed watermark.
	%T1_L=1800; T2_L=2200; 	% Baboon, Daub20, cl4
	%T1_L=2000; T2_L=2400;	% Lena, Symm10, cl4 
	%T1_L=700; T2_L=1400;	% Lena, Symm10, cl4 
	% Estimates for using detail subbands to embed watermark
	T1=120; T2=200; X1=20; X2=10;
	%T1=115; T2=200; X1=20; X2=10;
	%T1=110; T2=210; X1=10; X2=5;
	%T1=110; T2=210; X1=35; X2=25;
	%T1=100; T2=200; X1=25; X2=20;

	% Histograms
	noOfBins=10;
	if DO_L==1, getHist(L,noOfBins,'L'); end
	if DO_H==1, getHist(H,noOfBins,'H'); end
	if DO_D==1, getHist(D,noOfBins,'D'); end
	if DO_V==1, getHist(V,noOfBins,'V'); end

	% Set the seed values for the random no. generators
	% seedVal1 used in generation of watermark
	% seedVal2 used in detection of watermark
	seedVal1=round( sum(100*clock) * rand(1,1) );
	%seedVal1=100;
	seedVal2=seedVal1;
	fprintf('seedVal1=%d, seedVal2=%d\n',seedVal1,seedVal2);

	rand('seed',seedVal1); 
	binMat=round( rand(DCsize,DCsize) );  % 1s & 0s
	
	printMini(binMat,'binMat');
	if DO_L==1, printMini(L,'L'); end
	if DO_H==1, printMini(H,'H'); end
	if DO_D==1, printMini(D,'D'); end
	if DO_V==1, printMini(V,'V'); end

	% Insert the same watermark into each subband and store
	% insertion locations in the mask
	if DO_L==1, [L_,L_MASK,L_Wo]=INSERT(L,binMat,T1_L,T2_L,X1_L); end
	if DO_H==1, [H_,H_MASK,H_Wo]=INSERT(H,binMat,T1,T2,X1); end
	if DO_D==1, [D_,D_MASK,D_Wo]=INSERT(D,binMat,T1,T2,X1); end
	if DO_V==1, [V_,V_MASK,V_Wo]=INSERT(V,binMat,T1,T2,X1); end
	totalWmLenEmbedded=0;
	if DO_L==1, totalWmLenEmbedded=totalWmLenEmbedded + sum(L_MASK(:)); end
	if DO_H==1, totalWmLenEmbedded=totalWmLenEmbedded + sum(H_MASK(:)); end
	if DO_D==1, totalWmLenEmbedded=totalWmLenEmbedded + sum(D_MASK(:)); end
	if DO_V==1, totalWmLenEmbedded=totalWmLenEmbedded + sum(V_MASK(:)); end
	fprintf('totalWmLenEmbedded=%d\n',totalWmLenEmbedded);
	totalWmLenEmbeddedStore=[totalWmLenEmbeddedStore totalWmLenEmbedded];
	
	if DO_L==1, printMini(L_,'L_'); end
	if DO_H==1, printMini(H_,'H_'); end
	if DO_D==1, printMini(D_,'D_'); end
	if DO_V==1, printMini(V_,'V_'); end

	% Replace the original subbands with the watermarked subbands
	if DO_L==1, fwt(1:DCsize,1:DCsize) = L_; end
	if DO_H==1, fwt(1:DCsize,DCsize+1:DCsize*2) = H_; end
	if DO_D==1, fwt(DCsize+1:DCsize*2,DCsize+1:DCsize*2) = D_; end
	if DO_V==1, fwt(DCsize+1:DCsize*2,1:DCsize) = V_; end
	
	% Inverse wavelet transform
	WI = IWT2_PO(fwt,cl,qmf);

	% Tidy up 
	% Also removing position files (masks) here! Also removing 
	% the original watermark files as these will be reselected
	% from the image sized watermark in the recovery process.
	if DO_L==1, clear L L_ L_MASK L_Wo; end
	if DO_H==1, clear H H_ H_MASK H_Wo; end
	if DO_D==1, clear D D_ D_MASK D_Wo; end
	if DO_V==1, clear V V_ V_MASK V_Wo; end
	clear fwt binMat


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Attack the watermarked image
	%WIA=jpegComp(WI,5);
	WIA=jpegComp(WI,10);
	%WIA=jpegComp(WI,15);
	%WIA = MED2_M(WI,3);
	%WIA = MED2_M(WI,5);
	%WIA = AV2_M(WI,3);
	%WIA = AV2_M(WI,5);
	%WIA=double(imnoise(uint8(WI),'gaussian',0,0.006));
	%WIA=gausAtt(WI,0,375);	
	%WIA=double(imnoise(uint8(WI),'salt & pepper',0.015));
	%WIA=crop(WI,60,190,60,190);
	%WIA=subSamp(WI,0.5);
	%WIA=subSamp(WI,0.6);
	%WIA=subSamp(WI,0.7);
	%WIA = WI;
	%WIA = I;
	%WIA=stirmarkAtt(WI);		% requires StirMark download
	%WIA=stirmarkAttNoGeo(WI);	% requires StirMark download
	%WIA=stirmarkAtt1p0(WI);	% requires StirMark download
	%WIA=stirmarkAtt1p0NoGeo(WI);	% requires StirMark download

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



	% Recover the watermark
	fwt=FWT2_PO(WIA,cl,qmf);

	% Get the coarsest subands (excluding LL)
	if DO_L==1, L=fwt(1:DCsize,1:DCsize); end		       % dc part
	if DO_H==1, H=fwt(1:DCsize,DCsize+1:DCsize*2); end             % horiz
	if DO_D==1, D=fwt(DCsize+1:DCsize*2,DCsize+1:DCsize*2); end    % diag
	if DO_V==1, V=fwt(DCsize+1:DCsize*2,1:DCsize); end             % vert

	% Recover the watermarks from each subband
	if DO_L==1, [L_Wo,L_Wr] = RECOVER(L,T1_L,T2_L,seedVal2,X2_L); end
	if DO_H==1, [H_Wo,H_Wr] = RECOVER(H,T1,T2,seedVal2,X2); end
	if DO_D==1, [D_Wo,D_Wr] = RECOVER(D,T1,T2,seedVal2,X2); end
	if DO_V==1, [V_Wo,V_Wr] = RECOVER(V,T1,T2,seedVal2,X2); end


	% Statistics
	if DO_L==1, [L_MATCHES,L_OUTOF,L_BER,L_NC]=STATS(L_Wo,L_Wr,'L'); end
	if DO_H==1, [H_MATCHES,H_OUTOF,H_BER,H_NC]=STATS(H_Wo,H_Wr,'H'); end
	if DO_D==1, [D_MATCHES,D_OUTOF,D_BER,D_NC]=STATS(D_Wo,D_Wr,'D'); end
	if DO_V==1, [V_MATCHES,V_OUTOF,V_BER,V_NC]=STATS(V_Wo,V_Wr,'V'); end

	if DO_L==1, L_BERStore=[L_BERStore L_BER]; L_NCStore=[L_NCStore L_NC]; end
	if DO_H==1, H_BERStore=[H_BERStore H_BER]; H_NCStore=[H_NCStore H_NC]; end
	if DO_D==1, D_BERStore=[D_BERStore D_BER]; D_NCStore=[D_NCStore D_NC]; end
	if DO_V==1, V_BERStore=[V_BERStore V_BER]; V_NCStore=[V_NCStore V_NC]; end

	% Concatenating the WM's for an overall NC value.
	combWo=[];
	combWr=[];
	if DO_L==1, combWo=[combWo L_Wo]; combWr=[combWr L_Wr]; end
	if DO_H==1, combWo=[combWo H_Wo]; combWr=[combWr H_Wr]; end
	if DO_D==1, combWo=[combWo D_Wo]; combWr=[combWr D_Wr]; end
	if DO_V==1, combWo=[combWo V_Wo]; combWr=[combWr V_Wr]; end
	[Comb_MATCHES,Comb_OUTOF,Comb_BER,Comb_NC]=STATS(combWo,combWr,'Combined');
	Comb_BERStore=[Comb_BERStore Comb_BER];
	Comb_NCStore=[Comb_NCStore Comb_NC];
	wmLenRec=length(combWr);
	wmLenRecStore=[wmLenRecStore wmLenRec];

	fprintf('\nPSNR for original and watermarked images=%.2f dB\n',getPSNR(I,WI));
	TPE=WatsonMetric(I,WI,16);
	%TPE=0.0;
	fprintf('WatsonMetric=%f\n',TPE);
	[PSNR,wPSNR]=psnrMetric(I,WI,4);
	fprintf('PSNR=%.2f, wPSNR=%.2f\n',PSNR,wPSNR);

	PSNRStore=[PSNRStore PSNR];
	wPSNRStore=[wPSNRStore wPSNR];
	TPEStore=[TPEStore TPE];
end; clear itcount; 

% Show the different images
fprintf('\n\ndugadRun1\n');
qmfDug=MakeONFilter('Daubechies',8); % 8-tap daubechies filter (p=8 (vanishing moments); len=2p);
DUG_WI=dugadRun1(I,seedVal1,qmfDug);
GrayImage([I WI DUG_WI WIA]); title('Original, Novel Watermarked, Dugad Watermarked and Attacked');

if DO_L==1, avL_BER=sum(L_BERStore)/noOfIts; avL_NC=sum(L_NCStore)/noOfIts; end
if DO_H==1, avH_BER=sum(H_BERStore)/noOfIts; avH_NC=sum(H_NCStore)/noOfIts; end
if DO_D==1, avD_BER=sum(D_BERStore)/noOfIts; avD_NC=sum(D_NCStore)/noOfIts; end
if DO_V==1, avV_BER=sum(V_BERStore)/noOfIts; avV_NC=sum(V_NCStore)/noOfIts; end
avComb_BER=sum(Comb_BERStore)/noOfIts;
avComb_NC=sum(Comb_NCStore)/noOfIts;
avTotalWmLenEmbedded=sum(totalWmLenEmbeddedStore)/noOfIts;
avWmLenRec=sum(wmLenRecStore)/noOfIts;
avPSNR=sum(PSNRStore)/noOfIts;
avWPSNR=sum(wPSNRStore)/noOfIts;
avTPE=sum(TPEStore)/noOfIts;
MIN_NC = min(Comb_NCStore);
MIN_LEN = min(wmLenRecStore);
MAX_NC = max(Comb_NCStore);
MAX_LEN = max(wmLenRecStore);
PFP=[];
for i=1:length(wmLenRecStore),
	PFP=[PFP falsePosCalc(Comb_NCStore(i),wmLenRecStore(i))];
end; clear i;
PFP_MIN=min(PFP);
PFP_MAX=max(PFP);
PFP_AVERAGE_CASE=falsePosCalc(avComb_NC,round(avWmLenRec));
PFP_WORST_CASE=falsePosCalc(MIN_NC,MIN_LEN);
%PFP_MAX
%PFP_AVERAGE_CASE
%PFP_WORST_CASE

fprintf('\n\n\nnoOfIts=%d\n',noOfIts);
if DO_L==1, fprintf('avL_BER=%f, avL_NC=%f\n',avL_BER,avL_NC); end
fprintf('T1=%d, T2=%d, X1=%d, X2=%d\n',T1,T2,X1,X2);
if DO_H==1, fprintf('avH_BER=%f, avH_NC=%f\n',avH_BER,avH_NC); end
if DO_D==1, fprintf('avD_BER=%f, avD_NC=%f\n',avD_BER,avD_NC); end
if DO_V==1, fprintf('avV_BER=%f, avV_NC=%f\n',avV_BER,avV_NC); end
fprintf('avComb_BER=%f, avComb_NC=%f\n',avComb_BER,avComb_NC);
fprintf('avTotalWmLenEmbedded=%d\n',avTotalWmLenEmbedded);
fprintf('avWmLenRec=%d\n',avWmLenRec);
fprintf('avPSNR=%f\n',avPSNR);
fprintf('avWPSNR=%f\n',avWPSNR);
fprintf('avTPE=%f\n',avTPE);
fprintf('MIN_NC combined = %f\n',MIN_NC);
fprintf('MIN_LEN = %d\n',MIN_LEN);
fprintf('MAX_NC combined = %f\n',MAX_NC);
fprintf('MAX_LEN = %d\n',MAX_LEN);
fprintf('PFP_MIN = %.20f\n',PFP_MIN);
fprintf('PFP_MAX = %.20f\n',PFP_MAX);
fprintf('PFP_AVERAGE_CASE= %.20f\n',PFP_AVERAGE_CASE);
fprintf('PFP_WORST_CASE= %.20f\n',PFP_WORST_CASE);

fp=fopen('DATA_IN7.dat','a');
fprintf(fp,'noOfIts=%d\n',noOfIts);
if DO_L==1, fprintf(fp,'avL_BER=%f, avL_NC=%f\n',avL_BER,avL_NC); fprintf('T1_L=%d, T2_L=%d\n',T1_L,T2_L); end
fprintf(fp,'T1=%d, T2=%d, X1=%d, X2=%d\n',T1,T2,X1,X2);
if DO_H==1, fprintf(fp,'avH_BER=%f, avH_NC=%f\n',avH_BER,avH_NC); end
if DO_D==1, fprintf(fp,'avD_BER=%f, avD_NC=%f\n',avD_BER,avD_NC); end
if DO_V==1, fprintf(fp,'avV_BER=%f, avV_NC=%f\n',avV_BER,avV_NC); end
fprintf(fp,'avComb_BER=%f, avComb_NC=%f\n',avComb_BER,avComb_NC);
fprintf(fp,'avTotalWmLenEmbedded=%d\n',avTotalWmLenEmbedded);
fprintf(fp,'avWmLenRec=%d\n',avWmLenRec);
fprintf(fp,'avPSNR=%f\n',avPSNR);
fprintf(fp,'avWPSNR=%f\n',avWPSNR);
fprintf(fp,'avTPE=%f\n',avTPE);
fprintf(fp,'MIN_NC combined = %f\n',MIN_NC);
fprintf(fp,'MIN_LEN = %d\n',MIN_LEN);
fprintf(fp,'MAX_NC combined = %f\n',MAX_NC);
fprintf(fp,'MAX_LEN = %d\n',MAX_LEN);
fprintf(fp,'PFP_MIN = %.20f\n',PFP_MIN);
fprintf(fp,'PFP_MAX = %.20f\n',PFP_MAX);
fprintf(fp,'PFP_AVERAGE_CASE= %.20f\n',PFP_AVERAGE_CASE);
fprintf(fp,'PFP_WORST_CASE= %.20f\n\n\n',PFP_WORST_CASE);
fclose(fp);


%%%%%
function printMini(M,str);
global PRINT_MINI_INFO
if PRINT_MINI_INFO == 1
	fprintf(str);
	M(1:3,1:3)
end

%%%%%
function [SCALE_WM,MASK,Wo]=INSERT(SCALE,binMat,T1,T2,X);
MASK = ((abs(SCALE)) > T1) & ((abs(SCALE)) < T2);
T1=T1+X;
T2=T2-X;
Ck = SCALE .* MASK;
T2Pos= Ck > 0 & binMat==1; 	T2Pos=T2Pos*T2;
T1Pos= Ck > 0 & binMat==0;	T1Pos=T1Pos*T1;
negT2Pos= Ck < 0 & binMat==1;	negT2Pos=negT2Pos*(-1)*T2;
negT1Pos= Ck < 0 & binMat==0;	negT1Pos=negT1Pos*(-1)*T1;
SCALE_WM = SCALE;
SCALE_WM = SCALE_WM .* (MASK==0);
SCALE_WM = SCALE_WM + T2Pos + T1Pos + negT2Pos + negT1Pos;
LOCS = T2Pos + T1Pos + negT2Pos + negT1Pos;
[Rv,Cv]=find(LOCS);
Wo=[];
for i=1:length(Rv),
        Wo=[Wo LOCS(Rv(i),Cv(i))];
end; clear i;
Wo=abs(Wo)==T2;



%%%%%
function [Wo,Wr] = RECOVER(SCALE,T1,T2,seedVal2,X2);
%LOCS = SCALE.*SCALE_MASK;
LOCS = ((abs(SCALE)) >= T1+X2) & ((abs(SCALE)) <= T2-X2); 
LOCS = LOCS .* SCALE;
% Regenerate watermark
rand('seed',seedVal2); binmat=round(rand(size(SCALE)));
[Rv,Cv]=find(LOCS);
Wo=[];
Wr=[];
for i=1:length(Rv),
	Wo=[Wo binmat(Rv(i),Cv(i))];
	Wr=[Wr LOCS(Rv(i),Cv(i))];
end; clear i;
CENTRE=(T1+T2)/2;
Wr= abs(Wr) >= CENTRE;



%%%%%
function [MATCHES,OUTOF,BER,NC]=STATS(O,R,str);
matches=sum(O==R);
non_matches = sum(O~=R);
nc = (matches - non_matches) /length(O);
fprintf('Subband %s MATCHES   = %d/%d\n',str,matches,length(O));
fprintf('Subband %s BER       = %f\n',str,100*(1-(matches/length(O))));
fprintf('Subband %s NC        = %f\n',str,nc);
fprintf('Subband %s Vector NC = %f\n',str,nc);
MATCHES=matches;
OUTOF=length(O);
BER=100*(1-(matches/length(O)));
NC=nc;

%%%%%
function NC=get_VectorNC(v,v_);
v=ShapeAsRow(v);
v_=v_(:);
%NC=(v*v_) / (sqrt(sum(v.^2))*sqrt(sum(v_.^2)));
NC=(v*v_) /  sqrt( sum(v.^2) * sum(v_.^2) );



%%%%%
function JC=jpegComp(I,Q);
imwrite(uint8(I), 'I.jpg','Quality',Q);
JC=double(imread('I.jpg')); 
!rm I.jpg

%%%%%
function out=gausAtt(I,M,V);
out = I + sqrt(V)*randn(size(I)) + M;

%%%%%
function WIA=crop(WI,v1,v2,v3,v4);
WIA=zeros(size(WI));
WIA(v1:v2,v3:v4)=WI(v1:v2,v3:v4);


%%%%%
function AT=stirmarkAtt(I);
pnmwrite(uint8(I),'stir.ppm');
!pnmnoraw stir.ppm > stirAscii.ppm
!StirMark stirAscii.ppm stirAscii.jpg
AT=imread('stirAscii.jpg');
AT=rgb2gray(AT);
AT=double(AT);
!rm stir.ppm stirAscii.ppm stirAscii.jpg

%%%%%
function AT=stirmarkAttNoGeo(I);
pnmwrite(uint8(I),'stir.ppm');
!pnmnoraw stir.ppm > stirAscii.ppm
!StirMark -i0 -o0 stirAscii.ppm stirAscii.jpg
AT=imread('stirAscii.jpg');
AT=rgb2gray(AT);
AT=double(AT);
!rm stir.ppm stirAscii.ppm stirAscii.jpg

%%%%%
function WIA = stirmarkAtt1p0(WI);
imwrite(uint8(WI),'WI.tif');
!tifftopnm WI.tif > WI.pgm
!stirmark1.0 WI.pgm stirred.pgm
!pnmtotiff stirred.pgm > stirred.tif
WIA = double(imread('stirred.tif'));
!rm WI.tif WI.pgm stirred.pgm stirred.tif

%%%%%
function WIA = stirmarkAtt1p0NoGeo(WI);
imwrite(uint8(WI),'WI.tif');
!tifftopnm WI.tif > WI.pgm
!stirmark1.0 -i0 -o0 WI.pgm stirred.pgm
!pnmtotiff stirred.pgm > stirred.tif
WIA = double(imread('stirred.tif'));
!rm WI.tif WI.pgm stirred.pgm stirred.tif

%%%%%
function WIATT=subSamp(WI,subFact);
[WIR,WIC]=size(WI);
subSampIm=imresize(WI,subFact);
WIATT=imresize(subSampIm,[WIR,WIC]);

%%%%%
function psnr=getPSNR(O,A);
difference=(O-A).^2;
difference=difference(:);
mse=sum(difference) / length(difference);
psnr = 10 * log10( (255^2) / mse);

%%%%%
function [Hist,Vals]=myHist(I,noOfBins);
I=abs(I);
mn=min(I(:));
mx=max(I(:));
%step=floor(mx/mn);
step=mx/noOfBins;
Hist=[];
Vals=[];
mask = I < mn;
tot=sum(mask(:));
Hist = [Hist tot];
Vals=[Vals 0];
for i=mn:step:mx
	mask = I >= i & I < (i+step);
	tot = sum(mask(:));
	Hist=[Hist tot];
	Vals=[Vals i];
end; clear i;

%%%%%
function getHist(Scale,noOfBins,ScaleStr);
[Hist,Vals]=myHist(Scale,noOfBins);
%figure, bar(Vals,Hist);
fprintf('Histogram for scale %s:\n',ScaleStr);
fprintf('Hist\n');
for i=1:length(Hist),
	fprintf('%6d ',Hist(i));
end; clear i;
fprintf('\n');
fprintf('Vals\n');
for i=1:length(Vals),
	fprintf('%6d ',round(Vals(i)));
end; clear i;
fprintf('\n');
