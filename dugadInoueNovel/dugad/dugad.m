function dugad();

T1=40;
T2=60;

seedVal1=393;
seedVal2=393;
% seedVal2 is the seedValue used to reconstruct the image sized watermark.
% hopefully, only the key holder will know this value. Only if
% seedVal1=seedVal2 should the watermark be detected.

alpha=0.2;

% FWT2_PO computes (image_dyadlength - cl) wavelet levels,which
% in the case of 256x256 image, is (8-cl).
cl=5;
qmf=MakeONFilter('Daubechies',8);

I=double(imread('lena_256.tif'));
[theLen,po2]=dyadlength(I(1,:));
noWavLevs=po2-cl;	% how many wavelet levels to compute.

% setting the random number generator
randn('seed',seedVal1);
WM=randn(size(I));

% inserting the watermark
WI=insertWM(I,WM,qmf,alpha,T1,cl,noWavLevs);

% Attacking the watermark
WI_Att=jpegComp(WI,10);
%WI_Att=double(imnoise(uint8(WI),'gaussian',0,0.0075));
%WI_Att=zeros(size(WI));WI_Att(50:205,50:205)=WI(50:205,50:205); % Crop attack
%WI_Att=WI; % This is the case of no attack to the watermarked image.

% recovering the positions/coefficients in the wavelet transform whos
% magnitude is greater than T2.
MASKrec=recoverMASK(WI_Att,qmf,T2,cl,noWavLevs);

% Recovering the watermark from possibly corrupted image.
[FWTcoefs,wmRec]=recoverWM(WI_Att,MASKrec,qmf,seedVal2,cl);

% Correlation, threshold determination and watermark length stats.
Corr=print_CORR(FWTcoefs,wmRec);
print_NC(FWTcoefs,wmRec);
wmLen=length(FWTcoefs);
fprintf('length of FWTcoefs=%d\n',wmLen);
S=detThresh(FWTcoefs,alpha,length(FWTcoefs));
fprintf('S=%f\n',S);
GrayImage([I WI WI_Att]); 
title('Original image : Watermarked Image : Attacked Watermarked Image');
fprintf('PSNR=%.2f\n',getPSNR(I,WI));

% Storing the results to file.
fp=fopen('dugad.dat','a');
fprintf(fp,'alpha=%.2f T1=%.2f T2=%.2f\ncorr=%.2f\nS=%.2f\n',alpha,T1,T2,Corr,S);
fprintf(fp,'wmLength=%d\n\n',length(FWTcoefs));
fclose(fp);



%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%

%%%%
function WI=insertWM(I,WM,qmf,alpha,T1,cl,noWavLevs);
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
function print_NC(v,v_);
v=ShapeAsRow(v);
v_=v_(:);
%NC=(v*v_) / (sqrt(sum(v.^2))*sqrt(sum(v_.^2)));
NC=(v*v_) /  sqrt( sum(v.^2) * sum(v_.^2) );
fprintf('Vector NC=%f\n',NC);


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



