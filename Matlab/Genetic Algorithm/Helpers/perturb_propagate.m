%The Goal of this function is to model misalignment and scattering on the
%effect of a pupil mask while imaging through a preset volume. The
%Parameters are:
%{
    Volume = 3D stack of 2D imaging planes seperated by axial distance dz
    Depth = 1D depth information for provided stack
    Pupil = designed mask
    Kernel = defocus kernel (Paraxial, aberrated, with scattering...
    iF = shorthand handle for inverse fourier transform
    F = same for Fourier transform
    Aperture = pupil aperture for system
    axmis = axial misalignment for generated mask (meters)
    latmis = laterial misalignment for generated mask (pixels)
    showFlag = shows modified pupil
%}
function [image, vol, pupil] = perturb_propagate(volume, depth, pupil, kernel, iF, F, Aperture, axmis, latmis,showFlag)
    pupil = circshift(pupil,[1,latmis]).*Aperture; %Misalign
    acrr= @(x) iF(conj(F(x)).*F(x));
    if showFlag
       figure
       imagesc(real(pupil))
       title('misaligned pupil','FontSize',14)
    end
    vol = zeros(size(volume));
    for j = 1:length(depth)
        vol(:,:,j) = real(iF(F((volume(:,:,j))).*acrr(pupil.*kernel(depth(j)))));
        %vol(:,:,j) = abs(iF(F(volume(:,:,j)).*pupil.*kernel(depth(j)))).^2;
    end
    image = sum(vol,3);

end