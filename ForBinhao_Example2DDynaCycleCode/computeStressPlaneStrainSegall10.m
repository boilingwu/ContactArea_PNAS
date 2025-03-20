function [s22,s23,s33]=computeStressPlaneStrainSegall10(x2,x3,q2,q3,slip,opening,width,dip,mu,nu)
% COMPUTESTRESSPLANESTRAINSEGALL10
% compute the stress field due to a two-dimensional edge dislocation
% with shear and opening components.
%
% INPUT:
%   x2      horizontal component of observation points
%   x3      depth component of observation points (must be positive)
%   y2      horizontal component of updip end
%   y3      depth component of updip end (must be positive)
%   slip    shear component of dislocation
%   opening normal component of dislocation
%   width   down-dip dimension of the dislocation
%   dip     dip angle measured positive down from x2 axis (degrees)
%   nu      Poisson's ratio
%
% OUTPUT:
%   s22     horizontal stress
%   s23     shear stress
%   s33     vertical stress
%
% Authors: P. Cervelli, P. Segall, (2010)

assert(min(x3(:))>=0,'x3 must be positive.')
assert(q3>=0,'x3 must be positive.')

% The vertical axis is reckoned NEGATIVE down.
x= x2;
z=-x3;
y1= q2;
y2=-q3;

% make sure vertical coordinates are negative
z=-abs(z);	Z=-abs(y2);

% define Constants
bx=slip*sind(dip);
bz=slip*cosd(dip);

% components of slip vector
b1=slip*cosd(dip)-opening*sind(dip);
b2=slip*sind(dip)+opening*cosd(dip);

a1=[bx -bx]*mu/(2*pi*(1-nu));
a2=[bz -bz]*mu/(2*pi*(1-nu));

X=[y1 y1+width*cosd(dip)];
Z=[Z Z-width*sind(dip)];

Sxx=zeros(size(x));
Sxz=zeros(size(x));
Szz=zeros(size(x));

%Compute stresses
for i=1:2
    zmZ=z-Z(i);
    zmZ2=zmZ.*zmZ;
    zpZ=z+Z(i);
    zpZ2=zpZ.*zpZ;
    xmX=x-X(i);
    xmX2=xmX.*xmX;
    r12=zmZ2+xmX2;
    r22=zpZ2+xmX2;
    r122=r12.*r12;
    r222=r22.*r22;
    r223=r222.*r22;
    b1=xmX.*(3*zmZ2+xmX2)./(r122);
    b2=xmX.*(3*zpZ2+xmX2)./(r222);
    b3=4.0*Z(i)*z.*xmX.*(3.0*zpZ2-xmX2)./(r223);
    b4=zmZ.*(zmZ2-xmX2)./(r122);
    b5=zpZ.*(zpZ2-xmX2)./(r222);
    a3=(3.0*z+Z(i)).*zpZ2.*zpZ-6.0*z.*zpZ.*xmX2-xmX2.*xmX2;
    b6=2.0*Z(i).*a3./(r223);
    
    c1=xmX.*(zmZ2-xmX2)./(r122);
    c2=xmX.*(zpZ2-xmX2)./(r222);
    a4=(2.0*Z(i)-z).*zpZ2+(3.0*z+2.0*Z(i)).*xmX2;
    c3=4.0*Z(i).*xmX.*a4./(r223);
    c4=zmZ.*(zmZ2+3.0*xmX2)./(r122);
    c5=zpZ.*(zpZ2+3.0*xmX2)./(r222);
    a5=zmZ.*zpZ2.*zpZ-6.0*z.*zpZ.*xmX2+xmX2.*xmX2;
    c6=2.0*Z(i)*a5./(r223);
    
    Sxx=Sxx+a1(i)*(c1-c2+c3)+a2(i)*(-c4+c5+c6);
    Sxz=Sxz+a1(i)*(b4-b5-c6)+a2(i)*(-c1+c2-b3);
    Szz=Szz+a1(i)*(-b1+b2+b3)+a2(i)*(-b4+b5-b6);
end

s22=-Sxx;
s23=+Sxz;
s33=-Szz;

end