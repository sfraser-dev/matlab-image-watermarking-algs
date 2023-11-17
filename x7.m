function xxxxxx();

ALPHA=0.2;

% BEWARE of changing "slideWinSize" 
% (in case differnet pixels in window have same values)
slideWinSize=3; 	

% Reading input image and taking its WT.
% Grabbing the DC componet from it.
cl=4;
Iin=MyReadImage('Lenna');
[I_noRows,I_noCols]=size(Iin);
[len,po2]=dyadlength(Iin(:,1));
qmf=MakeONFilter('Daubechies',8);
numberWavLev=po2-cl;
DCsize=I_noCols;
for i=1:numberWavLev,
        DCsize=DCsize/2;
end; clear i;
fwt=FWT2_PO(Iin,cl,qmf);
DC=fwt(1:DCsize,1:DCsize);

% Cropping the input image (I is now the DC component)
I=DC;
[nrows,ncols]=size(I);
HEAD=floor((ncols/slideWinSize))*slideWinSize;
Icrop=I(:,1:HEAD);
[nrowsCrop,ncolsCrop]=size(Icrop);

% creating binary watermark to file all locations
spaceForWM=(HEAD/slideWinSize)*nrows;
wm=round(rand(1,spaceForWM));


% Quantizing the cropped image.
IcropQ=[];
rowstore=[];
i=1;
for r=1:nrowsCrop,
	for c=1:slideWinSize:ncolsCrop
		vec = Icrop(r,c:c+(slideWinSize-1));
		vecQ = quantizeNEW(vec,ALPHA,wm(i));
		rowstore=[rowstore vecQ];
		i=i+1;
	end; clear c;
	IcropQ = [IcropQ ; rowstore];
	rowstore=[];
end; clear r i;

% Replacing original image with its cropped part
% (in the appropriate place)
dcQ=I;
dcQ(:,1:HEAD)=IcropQ;

% putting the DC component back with the other wavelet coefficients
fwt(1:DCsize,1:DCsize)=dcQ;

% getting the watermarked image
WI=IWT2_PO(fwt,cl,qmf);


%%%%%%%%
% Attack
%%%%%%%%
%WIA = jpegComp(WI,5);
%WIA = stirmarkAttNoGeo(WI);
%WIA = stirmarkAtt1p0(WI);  
%WIA = stirmarkAtt(WI);
%WIA = crop(WI,100,200);
%WIA = MED2_M(WI,9);
%WIA = double(imnoise(uint8(WI),'Gaussian',0,0.02));
WIA = WI;
%WIA = Iin;

WatermarkImg=WI;



%% Recovery
fwtRec=FWT2_PO(WIA,cl,qmf);
DCrec=fwtRec(1:DCsize,1:DCsize);
WI=DCrec;	% WI is now the recovered DC component from possibly watermarked image.

% Taking the crop of possibly watermarked image.
[nrowsWI,ncolsWI]=size(WI);
HEAD=floor((ncolsWI/slideWinSize))*slideWinSize;
WIcrop=WI(:,1:HEAD);
[nrowsWICrop,ncolsWICrop]=size(WIcrop);

% recovering the watermark from the cropped image section.
recBitStream=[];
for r=1:nrowsWICrop,
        for c=1:slideWinSize:ncolsWICrop
		vec = WIcrop(r,c:c+(slideWinSize-1));
		recBit=dequantizeNEW(vec,ALPHA);
		recBitStream=[recBitStream recBit];
	end; clear c;
end; clear r;

[wm(:) recBitStream(:)]
STATS(wm,recBitStream);
fprintf('psnr=%.2f\n',getPSNR(Iin,WatermarkImg));
%fprintf('tpe=%f\n',WatsonMetric(Iin,WatermarkImg,8));
GrayImage([Iin WatermarkImg WIA]); title('Iin, WatermarkImg and WIA');
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT=quantizeNEW(x,ALPHA,Wo);
 
sorted=sort(x);
[R,C]=find( sorted((length(sorted)+1)/2) == x);
F2=x(C);
F2Pos=C;
[F1,F1Pos]=min(x);
[F3,F3Pos]=max(x);
 
WM_BIT=Wo;
 
S = ALPHA * (abs(F3) - abs(F1)) / 2;
if WM_BIT == 1,
        L=F1+S;
else
        L=F1;
end
while ((L+(2*S)) < F2),
        L=L+(2*S);
end
if (F2-L) < (L+(2*S)-F2),
        F2new = L;
else
        F2new = L+(2*S);
end
 
OUT=x;
OUT(F2Pos)=F2new;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT=dequantizeNEW(x,ALPHA);
IN=x; clear x;
sorted=sort(IN);
[R,C]=find( sorted((length(sorted)+1)/2) == IN);
F2=IN(C);
F2Pos=C;
[F1,F1Pos]=min(IN);
[F3,F3Pos]=max(IN);
 
S= ALPHA * (abs(F3) - abs(F1)) / 2;
L=F1;
X=0;
while L < F2,
        L=L+S;
        X=X+1;
end
if abs(L-S-F2) < abs(L-F2),
        DQ=mod(X+1,2);
else
        DQ=mod(X,2);
end
% Dummy value.
OUT=DQ;


%%%%%
function JC=jpegComp(I,Q);
imwrite(uint8(I), 'I.jpg','Quality',Q);
JC=double(imread('I.jpg'));
!rm I.jpg
 
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
function WIA = stirmarkAtt1p0(WI);
imwrite(uint8(WI),'WI.tif');
!tifftopnm WI.tif > WI.pgm
!stirmark1.0 WI.pgm stirred.pgm
!pnmtotiff stirred.pgm > stirred.tif
WIA = double(imread('stirred.tif'));
!rm WI.tif WI.pgm stirred.pgm stirred.tif
 
%%%%%
function AT=stirmarkAttNoGeo(I);
pnmwrite(uint8(I),'stir.ppm');
!pnmnoraw stir.ppm > stirAscii.ppm
!StirMark -i0 -o0 -d0 stirAscii.ppm stirAscii.jpg
AT=imread('stirAscii.jpg');
AT=rgb2gray(AT);
AT=double(AT);
!rm stir.ppm stirAscii.ppm stirAscii.jpg
 
%%%%%
function WIA=crop(WI,v1,v2);
WIA=zeros(size(WI));
WIA(v1:v2,v1:v2)=WI(v1:v2,v1:v2);
 
%%%%%
function psnr=getPSNR(O,A);
difference=(O-A).^2;
difference=difference(:);
mse=sum(difference) / length(difference);
psnr = 10 * log10( (255^2) / mse);

%%%%%
function NCORR = STATS(Wo,Wr);
MATCHES=sum(Wo==Wr);
fprintf('MATCHES=%d/%d\n',MATCHES,length(Wr));
BER=1 - (MATCHES/length(Wr));
fprintf('BER=%f\n',BER);
NCORR=sum(Wo==Wr);
NCORR=(NCORR-(sum(Wo~=Wr))) / length(Wr);
fprintf('NC=%f\n',NCORR);
