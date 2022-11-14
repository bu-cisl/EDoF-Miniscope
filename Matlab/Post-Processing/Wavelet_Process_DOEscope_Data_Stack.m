%Process collected data of a certain frame using wavelet transform & Hard
%Thresholding
addpath('Export_Fig')
%save_dir = [pwd '\'];
save_dir = 'D:\DOEscope Data Collection\Figure 5 Thin Brain Slice\Alberto Redo GRIN\';

%SaveDir = 'D:/DOEscope Data Collection/DOE Acute Slices/'; %Save directory (leave blank to save in folder)
%addpath(SaveDir) %Incase our target file is NOT in our current directory

%Define System Parameters from Zemax - Note some parameters might be
%slightly different from optimization files. Hardware choices were refined
%after the mask was optimized. Previous different pixel sizes and NA were
%imagined.
lambda = 509e-9; %Center wavelength of OBJECT
n2 = 1.33; %Index of sample
NA0 = 0.35; %NA of GRIN Obj lens
mag = 9.41; %System magnification from optimized ZEMAX
pix = 2.9e-6; %Sensor array pixel size - my camera here is Sony IMX290
dx = pix/mag; %Image Space discretization

%Useful handles
F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));
acrr= @(x) iF(conj(F(x)).*F(x)); %autocorrelation using fourier transform (much faster)


%Target File Settings
% dir = 'D:/DOEscope Data Collection/DOE Acute Slices/';
% filename = 'Big_Test';
% datatype = '.tif';
% save = 'Wavelet_'; %Adds to beg orginal name
%dir = 'D:\DOEscope Data Collection\Pycromanager Captures\Alberto Redo DOE\Full resolution\';
%dir = 'D:\DOEscope Data Collection\Pycromanager Captures\Alberto Redo GRIN\Full resolution\';
dir = 'D:\DOEscope Data Collection\Figure 5 Thin Brain Slice\Alberto Redo GRIN\Full resolution\';
filename = 'Acquisition_NDTiffStack';
datatype = '.tif';
save = 'Wavelet_'; %Adds to beg orginal name

info = imfinfo([dir filename datatype]); %Get image data
w = info(1).Width; %Pixels for first image in stack (if stack)
h = info(1).Height;

%Generate grids
[xx,yy] = meshgrid([-floor(w/2):(floor(w/2))-1]*dx, [-floor(h/2):(floor(h/2))-1]*dx); %Image plane representation of image plane -> res/2 for inco
[uu,vv] = meshgrid([-floor(w/2):(floor(w/2))-1]*1/(w*dx), [-floor(h/2):(floor(h/2))-1]*1/(h*dx)); %du defined by 1/FoV, Fourier plane
NAx = uu*lambda; %converting to NA space (alottable angles) -> unitless! easy of scaling and design
NAy = vv*lambda;
NA = sqrt(NAx.^2+NAy.^2); %NA @ any given point

%Define Wavelet: Second Derivative of Gaussian is a peaked wavelet function
%that is summable to 0. Ergo it will extract features of a certain size and
%reject dc or slowly varing background
r = sqrt(xx.^2+yy.^2);
rho = 8*dx;
filt = 2/(sqrt(3)*pi^(1/4)).*(1-1/2*(r/(rho)).^2).*exp(-(r.^2./(2.*(rho).^2))); %Second derivative of gaussian
filt = filt./max(filt(:)); %Normalize to 1

%Apply wavelet filtering!
% raw = imread([dir filename datatype], 1); %Raw frame for normalization
% res = real(iF(F(double(raw)/double(max(raw(:)))).*F(filt))); %Normalize image to 1 and apply filter

%Normalize to original range
%res = round(res./max(res(:))*double(max(raw(:)))); %Convert processed result to same range as raw frame for scaled analysis

%Now lets process the stack and save it!
z = length(info); %Number slices
num_vids = 1;
dz = ceil(z/num_vids);

f_filt = F(filt); %Normally we divide by the totla sum of the filter to ensure the filtering process does not affect the max val, but here the filter is summable to 0 so.....
%f_filt = f_filt./max(abs(f_filt(:)));

%We run the risk of maxing out our dynamic range through convolving with
%the kernel. To keep the wavelet results tangible, I'll normalize the total
%fluorescence to it exactly matches the dynamic range. This won't adjust
%calcium fluorescent statistics (all linear transforms) but make the stacks easier to deal with.
tic
k = 1;
for j = 1:num_vids
    for zi = 1:dz
        fprintf(['Processing Frame ' num2str(k) '\n'])
        raw = double(imread([dir filename datatype], k)); %Raw frame & normalize by max dynamic range value
        res = real(iF(F(raw).*f_filt))/25; %Convolve - I know the maximum expected value the wavelet transform will amplify target signals, so I can rescale using that

        if info(1).BitDepth == 8
           imwrite(uint8(res),[save_dir save filename num2str(j) datatype],'WriteMode','append','Compression','none'); %Write to a new file
        elseif info(1).BitDepth == 16
           imwrite(uint16(res),[save_dir save filename num2str(j) datatype],'WriteMode','append','Compression','none'); %Write to a new file
        end
        k = k+1;
        %pause(0.01)
    end
end
toc