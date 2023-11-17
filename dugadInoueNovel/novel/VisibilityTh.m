%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generates the base thresholds.

% Copyright (c) 2001,Victor Iquise, Maribel Madueno, Shelby Pereira, University of Geneva
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
function [tij]=VisibilityTh()

% Viewing distance in centimeters.
viewing_distance=72; 

% Matlab horizontal and vertical image size in centimeters.
horizontal_screen_image_size=9.4;  
vertical_screen_image_size=8.8; 

% In 9.4 centimeters there is 344 pixels according to xv.
% In 8.8 centimeters there is 342 pixels according to xv.
vertical_pixel_size=vertical_screen_image_size/344;
horizontal_pixel_size=horizontal_screen_image_size/342;

% Compute pixel size in radians.
Wx=vertical_pixel_size/viewing_distance;
Wy=horizontal_pixel_size/viewing_distance;

% Convert pixel size in radians to size in degrees.
Wx=Wx*180/pi;
Wy=Wy*180/pi;

% DCT block size.
block_size=8;

% Compute spatial frequencies i.e. fij.
fio=0:block_size-1;
fio=fio/(2*block_size);
foj=fio;
fio=fio/Wx;
foj=foj/Wy;
[Foj,Fio]=meshgrid(fio,foj);
fij=sqrt(Fio.^2+Foj.^2);

% Compute angular parameters i.e. tetaij.
tetaij=2*Fio.*Foj;
tetaij(2:block_size,2:block_size)=...
  tetaij(2:block_size,2:block_size)./fij(2:block_size,2:block_size).^2;
tetaij=asin(tetaij);

% Determine obliqueness effect values i.e. rij.
r=0.7;
rij=r + (1 - r)*cos(tetaij).^2;

% Compute the mean luminance i.e. L in candels per square meter
% according to an exponential law.
%L=128; % mean luminance in graylevels  

L=240; % mean luminance in graylevels  


L=power(10, (7/255*L-4)); % mean luminance in candels per square meter

% Compute luminance functions i.e. Tmin, fmin and K.
LT=13.45; So=94.7; at=0.649;
Lk=300; ko=3.125; ak=0.0706;
Lf=300; fo=6.78; af=0.182;

% Minimal threshold function i.e. Tmin.
if L > LT
  Tmin=L/So;   
else
  Tmin=L/So*power((LT/L), 1-at);   
end

% Parabola stepness function i.e. K.
if L > Lk
  parabola_stepness=ko;
else
  parabola_stepness=ko*power(L/Lk, ak);   
end

% minimal frequency function i.e. fmin
if L > Lf
  fmin=fo;
else
  fmin=fo*power(L/Lf, af);   
end

% Compute luminance thresholds i.e. tij.
% fij(1,1) is 0 but to don't have log(0) we make fij(1,1)=1.
fij(1,1)=1;
tij=log10(Tmin./rij)+parabola_stepness*power(log10(fij/fmin), 2);
tij=power(10, tij);

% Approximate luminance threshold for DC coefficient
tij(1,1)=min([ tij(1,2) tij(2,1) ]);

% According to our experiments we can multiply the thresholds until a
% factor of 3.5 so that they are not yet visibly after luminance
% masking corrections.
%tij=3.5*tij;

tij=1*tij;

