%% Inoue Method B.
%%  - 3 level WT.
%%  - Same WM into all 3 detail bands at 3rd wavelet level.
%%  - Position file used in recovery (along with T1 + T2).
%%  - The greater the difference between T1 and T2, the
%%    greater the robustness but the poorer the quality
%%    of the marked image.


function IN2();

% Debugging 
global PRINT_MINI_INFO
PRINT_MINI_INFO=0; 	% 0=don't print; 1=do print

% Which subbands to use?
DO_H=1; 		% 0=don't use; 1=use
DO_D=0;			% 0=don't use; 1=use
DO_V=0;			% 0=don't use; 1=use

% Get input image and its size
I=double(imread('lena_256.tif'));
[rows,cols]=size(I);

% Determine dyadlength of input image (for FWT2_PO) and
% calculate the amount of wavelet levels taken for a 
% requested coarselevel
[len,po2]=dyadlength(I(:,1));
cl=5;
qmf=MakeONFilter('Daubechies',16);
numberWavLev=po2-cl;

% Determine the size of the LL band (DC/low freq component)
DCsize=cols;
for i=1:numberWavLev,
        DCsize=DCsize/2;
end; clear i;

% Take the forward wavelet transform
fwt=FWT2_PO(I,cl,qmf);   

% Get the coarsest subands (excluding LL)
if DO_H==1, H=fwt(1:DCsize,DCsize+1:DCsize*2); end		% horiz
if DO_D==1, D=fwt(DCsize+1:DCsize*2,DCsize+1:DCsize*2); end	% diag
if DO_V==1, V=fwt(DCsize+1:DCsize*2,1:DCsize); end		% vert

% Set the two quantization thresholds
T1=50; T2=120;

% Set the seed values for the random no. generators
seedVal=2888;

rand('seed',seedVal); 
binMat=round( rand(DCsize,DCsize) );  % 1s & 0s

printMini(binMat,'binMat');
if DO_H==1, printMini(H,'H'); end
if DO_D==1, printMini(D,'D'); end
if DO_V==1, printMini(V,'V'); end

% Insert the same watermark into each subband and store
% insertion locations in the mask
if DO_H==1, [H_,H_MASK,H_Wo]=INSERT(H,binMat,T1,T2); end
if DO_D==1, [D_,D_MASK,D_Wo]=INSERT(D,binMat,T1,T2); end
if DO_V==1, [V_,V_MASK,V_Wo]=INSERT(V,binMat,T1,T2); end

if DO_H==1, printMini(H_,'H_'); end
if DO_D==1, printMini(D_,'D_'); end
if DO_V==1, printMini(V_,'V_'); end

% Replace the original subbands with the watermarked subbands
if DO_H==1, fwt(1:DCsize,DCsize+1:DCsize*2) = H_; end
if DO_D==1, fwt(DCsize+1:DCsize*2,DCsize+1:DCsize*2) = D_; end
if DO_V==1, fwt(DCsize+1:DCsize*2,1:DCsize) = V_; end

% Inverse wavelet transform
WI = IWT2_PO(fwt,cl,qmf);

% Tidy up
if DO_H==1, clear H H_; end
if DO_D==1, clear D D_; end
if DO_V==1, clear V V_; end
clear fwt binMat



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attack the watermarked image
%WIA = MED2_M(WI,5);
WIA=jpegComp(WI,5);
%WIA=stirmarkAtt(WI);
%WIA=double(imnoise(uint8(WI),'Gaussian',0,0.01));
%WIA=crop(WI,100,200);
%WIA = WI;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Recover the watermark
fwt=FWT2_PO(WIA,cl,qmf);

% Get the coarsest subands (excluding LL)
if DO_H==1, H=fwt(1:DCsize,DCsize+1:DCsize*2); end             % horiz
if DO_D==1, D=fwt(DCsize+1:DCsize*2,DCsize+1:DCsize*2); end    % diag
if DO_V==1, V=fwt(DCsize+1:DCsize*2,1:DCsize); end             % vert

% Recover the watermarks from each subband
if DO_H==1, H_Wr = RECOVER(H,H_MASK,T1,T2); end
if DO_D==1, D_Wr = RECOVER(D,D_MASK,T1,T2); end
if DO_V==1, V_Wr = RECOVER(V,V_MASK,T1,T2); end

% Show the different images
GrayImage([I WI WIA]); title('Original, Watermarked and Attacked');

% Statistics
if DO_H==1, STATS(H_Wo,H_Wr,'H'); end
if DO_D==1, STATS(D_Wo,D_Wr,'D'); end
if DO_V==1, STATS(V_Wo,V_Wr,'V'); end
fprintf('\nPSNR for original and watermarked images=%.2f dB\n',getPSNR(I,WI));



%%%%%
function printMini(M,str);
global PRINT_MINI_INFO
if PRINT_MINI_INFO == 1
	fprintf(str);
	M(1:3,1:3)
end

%%%%%
function [SCALE_WM,MASK,Wo]=INSERT(SCALE,binMat,T1,T2);
MASK = (abs(SCALE) > T1) & (abs(SCALE) < T2);
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
function Wr = RECOVER(SCALE,SCALE_MASK,T1,T2);
LOCS = SCALE.*SCALE_MASK;
[Rv,Cv]=find(LOCS);
Wr=[];
for i=1:length(Rv),
	Wr=[Wr LOCS(Rv(i),Cv(i))];
end; clear i;
CENTRE=(T1+T2)/2;
Wr= abs(Wr) >= CENTRE;



%%%%%
function STATS(O,R,str);
matches=sum(O==R);
non_matches = sum(O~=R);
nc = (matches - non_matches) /length(O);
fprintf('Subband %s MATCHES = %d/%d\n',str,matches,length(O));
fprintf('Subband %s BER     = %f\n',str,100*(1-(matches/length(O))));
fprintf('Subband %s NC      = %f\n',str,nc);


%%%%%
function JC=jpegComp(I,Q);
imwrite(uint8(I), 'I.jpg','Quality',Q);
JC=double(imread('I.jpg')); 
!rm I.jpg

%%%%%
function WIA=crop(WI,v1,v2);
WIA=zeros(size(WI));
WIA(v1:v2,v1:v2)=WI(v1:v2,v1:v2);

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
function psnr=getPSNR(O,A);
difference=(O-A).^2;
difference=difference(:);
mse=sum(difference) / length(difference);
psnr = 10 * log10( (255^2) / mse);

