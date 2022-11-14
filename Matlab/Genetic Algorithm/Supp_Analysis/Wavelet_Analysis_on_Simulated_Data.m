%Lets just use Params
%Load and generate grids
load Objects/Params.Mat
addpath('Objects','Helpers')
load Neur
lambda = Params.lambda;
n2 = Params.n2;
targ_size = Params.targ_size; %Target particle size
targ_depth = Params.targ_depth; %Target optimization depth
show_depth = Params.show_depth; %Additional obj simulation depth. We penalize light here
maxDepth = Params.maxDepth; %Max depth in final images
depth_step = Params.depth_step; %dz
depth = 0:depth_step:show_depth; %Depth vector for obj simulation
neurDepth = 0:depth_step:maxDepth;
N = Params.N; % Fixed - number of pixels on array
saveDir = Params.name;

F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));
acrr= @(x) iF(conj(F(x)).*F(x)); %autocorrelation using fourier transform (much faster)
maxnorm = @(x) x./max(abs(x(:)));
minmaxnorm = @(x) (x - min(x(:)))./(max(x(:))-min(x(:)));

%Physical parameters of system considering magnification and sensor size
%(assume not pixel limited)
NA0 = Params.NA0; %NA of objective lens - sets resolution
mag = Params.mag; %Leave as 1 else need to rescale all generated masks (though can easily improve code to generalize)
pix = Params.pix; %Pixel Pitch
dx = pix/mag; %Effective pixel size - magnification in system akin to having unit mag with smaller pixels
%innerCircle = Params.innerCircle; %Radius of inner Ring used for alignment

%Load the appropriate grids
ap = Params.ap; %This doesn't take up any more space, simply adds a shorter pointer to this quantity
xx = Params.xx;
yy = Params.yy;
uu = Params.uu;
vv = Params.vv;
NA = Params.NA;
eva = Params.eva;

%Define initial phase set
k = 2*pi*n2/lambda;
defo = @(d) exp(1i*k*d.*sqrt(1-(NA/n2).^2).*eva).*eva; %Angular spectrum kernel in a new media

%Define on-axis aberrations
Wd = 0;
W040 = Params.W040;
W131 = 0;
W222 = 0;
W220 = 0;
W311 = 0;
f_Num = Params.fnum;
f0 = 1/(2*(lambda*f_Num)); %This quantity is conserved due to the relationship between cutoff frequency and NA
W = Seidal_Coeff(0,0,uu./f0,vv./f0,Wd,W040,W131,W222,W220,W311,0); %On axis aberrations
ab = exp(-1i*k*W*lambda).*ap;

load mask

%% Simuate neurological
Neural.N = N/2; %Make EVEN
Neural.dx = dx; %Real space 'pixel' size
Neural.depth = neurDepth*3; %Long Range
%Neural.targetsize = targ_size;
Neural.targetsizes = [4,6]*1e-6;
Neural.ratios = [1,0]; %By setting a ratio to 0, that size is the spacing between neurons
Neural.targetsig = 1;
Neural.M = 3; %Number of targets
Neural.backgroundsize = 1e-6;
Neural.backgroundsig = Neural.targetsig*0.75;
Neural.ratio = 25; %background to foreground per slice
Neural.rand = 0; %Number of random targets added to scene (keeps background consistent)
Neural.noise = 'None';
Neural.var1 = 0.25; %Changes dymanics of selected noise (see function)
Neural.var2 = 0.1;

%Generate stack representing neural volume
Neur = generate_neural_volNoOverlap(Neural);
%Neur = generate_neural_vol2D(Neural);
Neur = padarray(Neur,[N/4,N/4],'both');
% figure('units','normalized','outerposition',[0 0 1 1])
% imagesc(sum(Neur,3))
% title(['Total depth = ' num2str(max(Neural.depth))])

[immask,~] = propagate_pseudoscattering(Neur,neurDepth,mask,defo,iF,F,ab ,Params.ls);

% Gaussian Filter
r = sqrt(xx.^2+yy.^2);
gaus = exp(-(r.^2./(2.*(N/4*dx).^2))); %Second derivative of gaussian
gaus = gaus./max(gaus(:)); %Normalize to 1

immask = immask.*gaus;
%% Make the wavelet 
rho = 8*dx;
filt = 2/(sqrt(3)*pi^(1/4)).*(1-1/2*(r/(rho)).^2).*exp(-(r.^2./(2.*(rho).^2))); %Second derivative of gaussian
filt = filt./max(filt(:)); %Normalize to 1

%Make a gaussian 'neuron'
slice = zeros(N,N);
slice(r.^2 <= (Params.targ_size/2).^2) = 1;
gaus = exp(-(r.^2./(2.*(Params.targ_size/4).^2))); %Second derivative of gaussian
gaus = gaus./max(gaus(:)); %Normalize to 1
slice = slice .* gaus;
%% Lets analyze filter
figure;
imagesc((-N/4:N/4-1),(-N/4:N/4-1),filt(N/2-N/4:N/2+N/4-1,N/2-N/4:N/2+N/4-1))
title('Optimized Filter: Real Space','fontsize',18)
xlabel('Side (pix)','fontsize',14)
ylabel('Side (pix)','fontsize',14)
axis('square')
colormap hot
colorbar

figure;
imagesc((-N/4:N/4-1),(-N/4:N/4-1),maxnorm(real(F(filt))))
title('Optimized Filter: Fourier','fontsize',18)
xlabel('Freq. (1/N*pix)','fontsize',14)
ylabel('Freq. (1/N*pix)','fontsize',14)
axis('square')
colormap hot
colorbar

figure;
imagesc((-N/4:N/4-1),(-N/4:N/4-1),sum(Neur(N/2-N/4:N/2+N/4-1,N/2-N/4:N/2+N/4-1,:),3))
title('Simulated Volume Sources','fontsize',18)
xlabel('Side (pix)','fontsize',14)
ylabel('Side (pix)','fontsize',14)
axis('square')
colormap hot
caxis([0 1])

figure;
imagesc((-N/4:N/4-1),(-N/4:N/4-1),minmaxnorm(immask(N/2-N/4:N/2+N/4-1,N/2-N/4:N/2+N/4-1)))
title('Raw Simulated Volume','fontsize',18)
xlabel('Side (pix)','fontsize',14)
ylabel('Side (pix)','fontsize',14)
axis('square')
colormap hot
caxis([0 1])
colorbar

tmp = real(iF(F(immask).*F(filt)));
figure;
imagesc((-N/4:N/4-1),(-N/4:N/4-1),minmaxnorm(tmp(N/2-N/4:N/2+N/4-1,N/2-N/4:N/2+N/4-1)))
title('Filtered Image','fontsize',18)
xlabel('Side (pix)','fontsize',14)
ylabel('Side (pix)','fontsize',14)
axis('square')
colormap hot
caxis([0 1])
colorbar

tmp = real(iF(F(immask).*F(1-filt)));
figure;
imagesc((-N/4:N/4-1),(-N/4:N/4-1),minmaxnorm(tmp(N/2-N/4:N/2+N/4-1,N/2-N/4:N/2+N/4-1)))
title('Rejected Signal','fontsize',18)
xlabel('Side (pix)','fontsize',14)
ylabel('Side (pix)','fontsize',14)
axis('square')
colormap hot
caxis([0 1])
colorbar

%% Now lets see the support vs filter!
tmp = maxnorm(real(log10(real(F(immask)))));
laterplot = tmp(N/2-N/4:N/2+N/4-1,N/2);
figure
imagesc((-N/4:N/4-1),(-N/4:N/4-1),(tmp(N/2-N/4:N/2+N/4-1,N/2-N/4:N/2+N/4-1)))
title('Log10 Image Spectrum','fontsize',18)
xlabel('Freq. (1/N*pix)','fontsize',14)
ylabel('Freq. (1/N*pix)','fontsize',14)
axis('square')
colormap hot
caxis([0.2 1])

tmp = maxnorm(abs(F(filt)));
h_ax = gca;
h_ax.Colormap = hot;
h_ax_c = axes('position', get(h_ax, 'position'), 'Color', 'none');
contour(h_ax_c, tmp(N/2-N/4:N/2+N/4-1,N/2-N/4:N/2+N/4-1),[0,0.2,0.6,0.9],'linewidth',1.25)
h_ax_c.Color = 'none';
% h_ax_c.Colormap = autumn;
h_ax_c.XTick = [];
h_ax_c.YTick = [];
set(h_ax_c,'box','off')
axis square

%% Now we can do a cross-section
figure
hold on
plot((-N/4:N/4-1),laterplot,'color','r','linewidth',1.25)
tmp2 =log10(abs(F(slice)));
plot((-N/4:N/4-1),tmp2(N/2-N/4:N/2+N/4-1,N/2)./max(tmp2(:)),'linestyle','--','color','k','linewidth',1.25)
plot((-N/4:N/4-1),tmp(N/2-N/4:N/2+N/4-1,N/2),'linewidth',1.25)
title('Filter: Fourier Analysis','fontsize',18)
xlabel('Freq. (1/N*pix)','fontsize',14)
ylabel('Norm. Amplitude','fontsize',14)
legend('Log10 Image Spectrum','Log10 Neuron Support','Designed Filter','location','southeast')