%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Victor Iquise,  Maribel Madueno                   %
%                         Computer Vision Group                        %
%                               16/03/01                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computes perceptual error. The perceptual error is computed using the
% Watson model which takes into account three factors: contrast 
% sensitivity, luminance masking and contrast masking.
% For a description of the Watson model see section 2.3 of the paper 
% "A comparison of image quality models and metrics based on human 
% visual sensitivity".

% Copyright (c) 2001,University of Geneva
% 
%  Permission to use, copy, modify, and distribute this software and its
%  documentation for any non-commercial purpose and without fee is hereby
%  granted (GPL), provided that the above copyright notice appear in all
%  copies and that both that copyright notice and this permission notice
% appear in supporting documentation. This software is provided "as is" 
% without express or implied warranty. The authors shall not be held
% liable in any event for incidental or consequential damages in
% connection with, or arising out of, the furnishing, performance, or
% use of this program.
% 
% If you use the Checkmark software package for your research, please cite:
%
% Shelby Pereira, Sviatoslav Voloshynovskiy, Maribel Madueño, Stéphane Marchand-Maillet
% and Thierry Pun, Second generation benchmarking and application oriented evaluation,
% In Information Hiding Workshop, Pittsburgh, PA, USA, April 2001.
%
%  http://cui.unige.ch/~vision/Publications/watermarking_publications.html
%
%  
%
% See the also the "Copyright" file provided in this package for
% copyright information about code used in the Checkmark package.
%
function [TPE,NLPE1,NLPE2,LPE]=WatsonMetric(I1,MI,BlockSize,varargin)

%[TPE,NLPE1,NLPE2,MAXLPE]=WatsonMetric(I,MI,BlockSize)
%[TPE,NLPE1,NLPE2,MAXLPE]=WatsonMetric(I,MI,BlockSize,LPETh1,LPETh2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BlockSize: Block size for calculating the local perceptual error. It 
%            should be multiple of 8.
% I: Original Image. 
% MI: Modified Image.
% LPETh1: First local perceptual error threshold.
% LPETh2: Second local perceptual error threshold.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TPE: Global perceptual error.
% NLPE1: Number of blocks of size BlockSize greaters than the first 
%        local perceptual error threshold.
% NLPE2: Number of blocks of size BlockSize greaters than the second 
%        local perceptual error threshold.
% MAXLPE: Maximum local perceptual error.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local perceptual error threshold.
if nargin==5
    LPETh1=varargin{1};
    LPETh2=varargin{2};    
else
    LPETh1=0.1;
    LPETh2=0.12;
end
I1=getlum255(I1);MI=getlum255(MI);
% Generates the base thresholds.
tij=VisibilityTh;  

% Compute contrast masking thresholds.
CIJK=blkproc(I1,[8 8],'dct2');
MIJK=blkproc(CIJK,[8 8],'ContrastMasking',tij);

% Compute quantization errors.
UIJK=blkproc(MI,[8 8],'dct2');
EIJK=abs(CIJK-UIJK);

% Compute perceptual error i.e. dijk.        
DIJK=EIJK./MIJK;   

% Compute local perceptual errors.
LPE=blkproc(DIJK,[BlockSize BlockSize],'mean2');

%figure(3);imagesc(LPE);

% Compute NLPE1, NLPE2, MAXLPE.
NLPE1=sum(LPE(:)>LPETh1);
NLPE2=sum(LPE(:)>LPETh2);
MAXLPE=max(LPE(:));

% Compute global perceptual error.
TPE=mean2(DIJK);

