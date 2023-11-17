% function lum=getlum255(imin)
% outputs the luminance of a color image
% normalized between 0 and 255, assumes the input is an RGB image
% if the image is not color, the original image is returned.

function lum=getlum255(imin)


if iscolor(imin)
  tmp1=rgb2ntsc(imin);
  lum=double(255*tmp1(:,:,1));
else
  lum=imin;
end
