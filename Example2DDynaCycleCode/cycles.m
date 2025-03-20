function c = cycles(m)
% CYCLES Shades of colors to represent velocity during seismic cycles.
%   CYCLES(M), is an M-by-3 matrix that defines a colormap.
%
%   cycles, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(cycles)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

if nargin < 1, m = size(get(gcf,'colormap'),1); end

cpt=[
-12.    0       0       0       
-11.0   0       0       0
-9.3    106     135     196
-8.8    135     164     224
-3.0    247     236     044
-1.0    239     064     035
+0.0    128     021     023
+1.0    050     021     023
];

x=-12+(0:(m-1))'/(m-1)*13.0;

r=interp1(cpt(:,1),cpt(:,2),x)/255;
g=interp1(cpt(:,1),cpt(:,3),x)/255;
b=interp1(cpt(:,1),cpt(:,4),x)/255;

c = [r g b]; 

