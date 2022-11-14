%The Goal of this function is to model propagation through an ideal medium. The
%Parameters are:
%{
    Volume = 3D stack of 2D imaging planes seperated by axial distance dz
    Depth = 1D depth information for provided stack
    Pupil = designed mask
    Kernel = defocus kernel (Paraxial, angular spectrum, aberrated, etc...
    iF = shorthand handle for inverse fourier transform
    F = same for Fourier transform
    ab = pupil plane aberrations. Since I play with these pretty readily, I
    explicityl separated them from the defocus kernel.
%}
function [image, vol] = propagate(volume, depth, pupil, kernel, iF, F, ab)
    acrr= @(x) iF(conj(F(x)).*F(x));
    vol = zeros(size(volume));
    for j = 1:length(depth)
        vol(:,:,j) = real(iF(F(volume(:,:,j)).*acrr(pupil.*kernel(depth(j)).*ab)));
        %vol(:,:,j) =
        %abs(iF(F((volume(:,:,j))).*pupil.*kernel(depth(j)).*ab)).^2;
    end
    image = sum(vol,3);

end