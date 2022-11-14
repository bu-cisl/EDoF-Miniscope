%Code designed to calculate the instintaneous Seidal coefficient aberration
%as a function of image height. Note, Seidal coefficients assume a radial
%geometry in the overarching lens system and this code only analyzes the
%five principal seidal coefficients. Higher order terms may be significant.

function [W] = Seidal_Coeff(u0,v0,xx_norm,yy_norm,Wd,W040,W131,W222,W220,W311, vis_flag)
%{
    u0,v0: Normalized height in IMAGE PLANE
    xx_norm,yy_norm: Normalized pupil grid from meshgrid
    Wd: Defocus term
    W040: Spherical Aberration
    W131: Coma
    W222: Astigmatism
    W220: Field Curvature
    W311: Distortion
%}
    beta = atan2(v0,u0); %Rotation angle on image plane
    h_im = sqrt(u0.^2+v0.^2); %Image height in polar
    
    %Alignment unit vector of Pupil to Image Point 
    %(this will calculate beyond desired values)and 
    %computes through a conversion into polar coordinates
    Xr = xx_norm*cos(beta) + yy_norm*sin(beta);
    Yr = -xx_norm*sin(beta)+yy_norm*cos(beta);
    
    %Rotational Pupil Coordinates
    rho = sqrt(Xr.^2+Yr.^2);
    
    %Calculate Defocus
    W = Wd*rho.^2 + ...
        W040*rho.^4 + ...
        W131*h_im*rho.^2.*Xr + ...
        W222*h_im^2*Xr.^2 +...
        W220*h_im^2*rho.^2 +...
        W311*h_im^3*Xr;
    
    %Visualize the wavefront distortion
    if vis_flag
       ap = ones(size(xx_norm)); 
       ap(rho > 1) = 0;
       figure
       h = surfc(xx_norm,yy_norm,W);
       xlabel('Normalized Pupil x')
       ylabel('Normalized Pupil y')
       title(['Seidal Aberrations at ' num2str(u0) ' ,' num2str(v0)])
       h.set('LineStyle','none')
    end
end