%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Victor Iquise and  Maribel Madueno           %
%                         Computer Vision Group                        %
%                               16/03/01                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Copyright (c) 2001, Shelby Pereira, University of Geneva
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
function [mijk]=ContrastMasking(cijk,tij)

% Contrast masking degree i.e. wij.
wij=0.7;

% Compute the mean luminance in DCT domain.
coo=128*8;
    
% Compute luminance masking threshods i.e. tijk.
cook=cijk(1,1);

if cook==0, cook=8; end 
% the factor 2.2 was added to increase the luminance masking
tijk=3.2*tij*power(cook/coo,0.65);
if cook<600   
  tijk=tijk*2.2;  % very dark areas can be marked a bit more
end

% Compute contrast masking thresholds i.e. mijk.
cijk_tijk=abs(cijk).^wij.*tijk.^(1-wij);
cijk_tijk(1,1)=tijk(1,1);

mijk=tijk;
index=find(tijk<cijk_tijk);
mijk(index)=cijk_tijk(index);
