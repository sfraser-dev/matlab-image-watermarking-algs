% Modified 7/02/2001

% Copyright (c) 2001,S. Voloshynovskiy, University of Geneva
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
function [NVF,varargout]=nvf(x,statistic,varargin)

% NVF=nvf_ng(x,statistic)
% NVF=nvf_ng(x,statistic,nhood)
% NVF=nvf_ng(x,statistic,Th0)
% NVF=nvf_ng(x,statistic,nhood,Th0)
% [NVF,LocalVar]=nvf_ng(x,statistic)
% [NVF,LocalVar]=nvf_ng(x,statistic,nhood)
% [NVF,LocalVar]=nvf_ng(x,statistic,Th0)
% [NVF,LocalVar]=nvf_ng(x,statistic,nhood,Th0)

% Calculation of the texture masking fuction based on Noise Visibility
% Function developed according non-stationary Gaussian model.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x:
% nhood: Local window size for the calculation of the local mean and the
%        local variance.
% Th0: Th0 is contrast sensitivity (Th|range:50-1000) Th=150;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NVF: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=cell(1,2);

if nargin==2
    % Default value for Th0.
    Th0=150;
    % Default value for nhood.
    nhood=[3 3];
elseif nargin==3
    if length(varargin{1})>1
        nhood=varargin{1};
        Th0=150;
    else
        Th0=varargin{1};
        nhood=[3 3];
    end
elseif nargin==4
    nhood=varargin{1};
    Th0=varargin{2};
else
    errordgl('The number of parameters is not correct',...
             'Input arguments');
end

switch statistic
case 'nsG'

    % Local mean and local variance calculation.
    LocalMean=filter2(ones(nhood),x)/prod(nhood);
    LocalVar=filter2(ones(nhood),x.^2)/prod(nhood)-LocalMean.^2;

    % Contrast sensitivity.
    Th=Th0/max(LocalVar(:));

    % NVF calculation.
    NVF=1./(1+Th*LocalVar);

    % Normalization of NVF to be in the range of (0;1).
    %NVF=NVF/max(NVF(:));    
    
case 'sgG'
    
end

if nargout>1
    varargout{1}=LocalMean;
    if nargout==3
        varargout{2}=LocalVar;
    end    
end
        