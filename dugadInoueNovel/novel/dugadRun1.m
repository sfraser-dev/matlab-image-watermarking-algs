function WI=dugad(IMAGEIN,SEEDIN,qmf);

CorrStore=[];
SStrore=[];
NCStore=[];
wmInLenStore=[];
wmOutLenStore=[];
PSNRStore=[];
wPSNRStore=[];
TPEStore=[];
noOfIts=1;

noOfFailures=0;
for itcount=1:noOfIts,
	T1=40;
	T2=50;

	seedVal1=SEEDIN;
	seedVal2=seedVal1;
	fprintf('seed=%d\n',seedVal1);
	% seedVal2 is the seedValue used to reconstruct the image sized watermark.
	% hopefully, only the key holder will know this value. Only if
	% seedVal1=seedVal2 should the watermark be detected.

	alpha=0.2;

	% FWT2_PO computes (image_dyadlength - cl) wavelet levels,which
	% in the case of 256x256 image, is (8-cl).
	cl=5;
	%qmf=MakeONFilter('Daubechies',8);	% 4 6 8 10 12 14 16 18 20
	
	I=double(IMAGEIN);

	[theLen,po2]=dyadlength(I(1,:));
	noWavLevs=po2-cl;	% how many wavelet levels to compute.

	% setting the random number generator
	randn('seed',seedVal1);
	WM=randn(size(I));

	% inserting the watermark
	[WI,wmInsertLen]=insertWM(I,WM,qmf,alpha,T1,cl,noWavLevs);

	TPE=WatsonMetric(I,WI,16);
	%TPE=0;
	fprintf('WatsonMetric=%f\n',TPE);    
	[PSNR,wPSNR]=psnrMetric(I,WI,4); % 3rd argument can be 2 or 4 (give same results)
	fprintf('PSNR=%.2f, wPSNR=%.2f\n',PSNR,wPSNR);
end; clear itcount;


%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%

%%%%
function [WI,wmInsertLen]=insertWM(I,WM,qmf,alpha,T1,cl,noWavLevs);
[rows,cols]=size(I);
DC_r=rows;
DC_c=cols;
for i=1:noWavLevs,
	DC_r=DC_r/2;
	DC_c=DC_c/2;
end; clear i
FWT=FWT2_PO(I,cl,qmf);
DC=FWT(1:DC_r,1:DC_c);
FWT(1:DC_r,1:DC_c)=0;

MASK=(abs(FWT) > T1);
%figure, GrayImage(MASK*255); title('The Mask');
WM=WM.*MASK;
WM=WM*alpha;
WM=WM.*abs(FWT);
%figure, GrayImage(WM); title('The Watermark');

FWT=FWT+WM;
FWT(1:DC_r,1:DC_c)=DC;
WI=IWT2_PO(FWT,cl,qmf);
wmInsertLen=sum(MASK(:));


%%%%%
function MASK=recoverMASK(WI,qmf,T2,cl,noWavLevs);
[rows,cols]=size(WI);
DC_r=rows;
DC_c=cols;
for i=1:noWavLevs,
        DC_r=DC_r/2;
        DC_c=DC_c/2;
end; clear i
FWT=FWT2_PO(WI,cl,qmf);
FWT(1:DC_r,1:DC_c)=0;
MASK=(abs(FWT) > T2);


%%%%%
function [FWTcoefs,wmRec]=recoverWM(WI,MASKrec,qmf,seedVal,cl,noWavLevs);
[rows,cols]=size(WI);
FWT=FWT2_PO(WI,cl,qmf);

randn('seed',seedVal);
WMregen=randn(size(WI));
[rVec,cVec]=find(MASKrec~=0);
FWTcoefs=[];
wmRec=[];
for i=1:length(rVec),
        FWTcoefs=[FWTcoefs FWT(rVec(i),cVec(i))];
	wmRec=[wmRec WMregen(rVec(i),cVec(i))];
end; clear i;



%%%%%
function Corr=print_CORR(v,v_);
v=ShapeAsRow(v);
v_=v_(:);
Corr=(v*v_) / (length(v));
fprintf('CORR=%f\n',Corr);

%%%%%
function NC=print_NC(v,v_);
v=ShapeAsRow(v);
v_=v_(:);
% OLD !!!  NC=(v*v_) / (sqrt(sum(v.^2))*sqrt(sum(v_.^2)));
NC=(v*v_) /  sqrt( sum(v.^2) * sum(v_.^2) );
fprintf('Vector NC=%f\n',NC);

v=ShapeAsRow(v);
v_=ShapeAsRow(v_);
meerNC=sum(v.*v_) / sqrt(sum(v.*v)*sum(v_.*v_));
fprintf('meerNC=%f\n',meerNC);


%%%%%
function S=detThresh(VI,alpha,M);
tot=sum(abs(VI(:)));
S= (alpha*tot) / (2*M);

%%%%%
function J=jpegComp(I,qf);
imwrite(uint8(I),'tmp.jpg','Quality',qf);
J=double(imread('tmp.jpg'));
!rm tmp.jpg

%%%%%
function psnr=getPSNR(O,A);
difference=(O-A).^2;
difference=difference(:);
mse=sum(difference) / length(difference);
psnr = 10 * log10( (255^2) / mse);

%%%%%
function WIA=crop(WI,v1,v2,v3,v4);
WIA=zeros(size(WI));
WIA(v1:v2,v3:v4)=WI(v1:v2,v3:v4);

%%%%%
function WIATT=subSamp(WI,subFact);
[WIR,WIC]=size(WI);
subSampIm=imresize(WI,subFact);
WIATT=imresize(subSampIm,[WIR,WIC]);


