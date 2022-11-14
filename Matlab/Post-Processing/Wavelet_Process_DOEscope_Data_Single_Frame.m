%Process collected data of a certain frame using wavelet transform & Hard
%Thresholding
addpath('Export_Fig')

SaveFlag = 0; %If 1, will save result and wavelet, if 0 it will just show the filtered result for inspection
SaveDir = 'C:\Users\joe19\Desktop\'; %Save directory (leave blank to save in folder)
addpath(SaveDir) %In case our target file is NOT in our current directory (saves us the need of changing the matlab directory, just add full path)

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
% dir = 'D:\DOEscope Data Collection\Alberto Redo\';
% filename = 'DOE Redo binned2';
% datatype = '.tiff';
dir = 'D:\DOEscope Data Collection\Figure 3 Fluorescent Fiber\Used Data\';
filename = 'Redo_DOE_Fiber_5';
datatype = '.tiff';
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
rho = 16*dx;
filt = 2/(sqrt(3)*pi^(1/4)).*(1-1/2*(r/(rho)).^2).*exp(-(r.^2./(2.*(rho).^2))); %Second derivative of gaussian
filt = filt./max(filt(:)); %Normalize to 1\

%Hard thresholding param
%dthres = 0.075; %Frac of maximum boundary

%Apply wavelet filtering!
raw = imread([dir filename datatype], 1); %Raw frame for normalization
res = real(iF(F(double(raw)/double(max(raw(:)))).*F(filt))); %Normalize image to 1 and apply filter
%Perform Min-Max image normalization since wavelet transform produces pos
%and neg values
res = (res - min(res(:)))/(max(res(:))-min(res(:)));
%res(res < 0) = 0;
%res = sqrt(res.^2);
%Normalize to original range
res = round(res./max(res(:))*double(max(raw(:)))); %Convert processed result to same range as raw frame for scaled analysis

%Apply hard thresholding
% res(res < 0) = 0; %Wavelet filter can cause more varying bg to assume negative values (not a signal)
% res(res <= dthres*max(res(:))) = 0; %Hard Thresh

figure
imagesc(res)
title(['Filtered w/ rho = ' num2str(rho/dx) 'pix'],'fontsize',16)
%colormap gray
%Custom Colormap
% vec = [      100;      90;       0]; %Corresponding points - this is in REVERSE 
% hex = ['#f4771e';'#41332d';'#63cd4a']; %Hex for low,mid,high (orange,black,green) - interactive colorbar to adjust visually
% dec = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255; %Convert to decimal
% N = 256; %Number of steps on colormap
% map = interp1(vec,dec,linspace(100,0,N),'pchip');
%Custom Colormap
map = colorcet('fire');
colormap(map)
colormap('hot')
%colormap(map)

colorbar
truesize
set(gca,'xtick',[])
set(gca,'xticklabel',[])

if SaveFlag == 1
    pause %Hit enter in the command window to resume - take this time to play with interactive colorbar!
    %After resuming, export fig
    export_fig([SaveDir 'manual_' filename datatype],'-native')
    
    %Also save greyscale image for reference
    if info.BitDepth == 8
       imwrite(uint8(res),[SaveDir save filename datatype],'WriteMode','overwrite','Compression','none'); %Write to a new file
    elseif info.BitDepth == 16
       imwrite(uint16(res),[SaveDir save filename datatype],'WriteMode','overwrite','Compression','none'); %Write to a new file
    end

close()
end

%Now Raw Frame
figure
imagesc(raw)
title('Raw Frame','fontsize',16)
%Custom Colormap
map = colorcet('fire');
colormap(map)
colorbar
axis('image')
truesize
set(gca,'xtick',[])
set(gca,'xticklabel',[])

%After resuming, export fig
if SaveFlag == 1
    pause %Hit enter in the command window to resume - take this time to play with interactive colorbar!
    export_fig([SaveDir 'manual_raw_' filename datatype],'-native')
    close()
end

figure
imagesc(res)
set(gca,'visible','off')
colormap gray
axis equal
export_fig([SaveDir 'fibbers_' num2str(round(rho/dx)) '_' filename datatype],'-native')

% figure
% imagesc(filt)
% set(gca,'visible','off')
% axis equal
% colormap hot
% 
% figure
% imagesc(real(F(filt)))
% set(gca,'visible','off')
% axis equal
% colormap hot
%%
h = surf(filt(1025-50:1025+50,1225-50:1225+50));  
set(h,'linestyle','none'); 
colormap jet
axis('square')
xlabel('Pixels','fontsize',12)
ylabel('Pixels','fontsize',12)
zlabel('Amplitude','fontsize',12)