%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Maribel Madueno                            %
%                         Computer Vision Group                        %
%                               19/03/01                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculates PSNR in dB for two images of the same size and same format
% (0:255).

% Copyright (c) 2001, University of Geneva
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
function [PSNR,wPSNR]=psnrMetric(a,b,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a:
% b:
% type:
% type=1 -> psnr
% type=2 -> nsG
% type=3 -> sgG
% type=4 -> psnr,nsG
% type=5 -> psnr,sgG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wPSNR=[];
if (a==b)
    out=inf;
else
    if size(a,3)==3
        a=[a(:,:,1) a(:,:,2) a(:,:,3)]; 
        b=[b(:,:,1) b(:,:,2) b(:,:,3)];
    end    
    if type>=1 & type<=5
        if type==2 | type==4
            statistics='nsG';
        else
            statistics='sgG';
        end
        NVF=nvf(a,statistics,150);
        NVF=NVF/max(NVF(:));
        c=NVF.*(a-b).^2;
    end
    c=(a-b).^2;
    PSNR=10*log10(255^2*prod(size(a))/sum(c(:)));
    c=NVF.*c;
    wPSNR=10*log10(255^2*prod(size(a))/sum(c(:)));        
end   

