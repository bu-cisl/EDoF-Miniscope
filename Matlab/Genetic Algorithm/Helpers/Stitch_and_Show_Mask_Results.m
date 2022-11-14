function mask = Stitch_and_Show_Mask_Results(x,weights,apoType,erodeRate,W040,saveDir)
    %Normalize weights in case inputs off
    weights = weights./sum(weights);
    %Load and generate grids
    load Objects/Params.Mat
    lambda = Params.lambda;
    n2 = Params.n2;
    targ_size = Params.targ_size; %Target particle size
    targ_depth = Params.targ_depth; %Target optimization depth
    show_depth = Params.show_depth; %Additional obj simulation depth. We penalize light here
    maxDepth = Params.maxDepth; %Max depth in final images
    depth_step = Params.depth_step; %dz
    depth = 0:depth_step:show_depth; %Depth vector for obj simulation
    neurDepth = 0:depth_step:maxDepth;
    N_pixels = Params.N; % Fixed - number of pixels on array
    
    F = @(x) fftshift(fft2(ifftshift(x)));
    iF = @(x) fftshift(ifft2(ifftshift(x)));
    acrr= @(x) iF(conj(F(x)).*F(x)); %autocorrelation using fourier transform (much faster)

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
    gam = @(g) exp(1i*pi*lambda/n2*g^4*(uu.^2+vv.^2).^2); %Spherical defocus 
    axi = @(q) exp(1i*pi*q*(sqrt((uu.^2+vv.^2)))); %Rotationally symmetric
    
    %Define on-axis aberrations
    Wd = 0;
    W040 = W040;
    W131 = 0;
    W222 = 0;
    W220 = 0;
    W311 = 0;
    f_Num = Params.fnum;
    f0 = 1/(2*(lambda*f_Num)); %This quantity is conserved due to the relationship between cutoff frequency and NA
    W = Seidal_Coeff(0,0,uu./f0,vv./f0,Wd,W040,W131,W222,W220,W311,0); %On axis aberrations
    ab = exp(-1i*k*W*lambda).*ap;
    
    %Develop mask
    mask = zeros(N_pixels);
    numSections = length(weights);
    
    if apoType == 1 %Concetric circles
        cutoffs = cumsum(sqrt(pi*(NA0*weights).^2/pi));
        if numSections == 1
            mask = (axi(x(1)).*defo(x(2)).*gam(x(3)));
        else
            for i = 1:numSections
                phase = (axi(x(1+(3*(i-1)))).*defo(x(2+(3*(i-1)))).*gam(x(3+(3*(i-1)))));
                if i == 1 %First section lower bounded by 0
                    mask(NA>0 & NA < cutoffs(i)) = phase(NA>0 & NA < cutoffs(i));
                elseif i < numSections
                    mask(NA> cutoffs(i-1) & NA < cutoffs(i)) = phase(NA> cutoffs(i-1) & NA < cutoffs(i));
                elseif i == numSections %Last section bounded by maxNAa
                    mask(NA > cutoffs(i) & NA < NA0) = phase(NA > cutoffs(i) & NA < NA0);
                end
            end
        end
        mask = binarize(mask,0).*ap; %Binarize final result

        figure('units','normalized','outerposition',[0 0 1 1])
        imagesc(mask)
        title('Generated Mask')
        axis('square')
        saveas(gcf,[saveDir '/mask.jpg'])
        close()
        
    elseif apoType == 2
        cutoffs = cumsum(2*pi*weights);
        %Define Angular Grid
        ang = atan(-yy./xx);
        ang(xx < 0 & yy < 0) = ang(xx < 0 & yy < 0)+pi; %Normalize to 2pi rotation
        ang(xx < 0 & yy >= 0) = ang(xx < 0 & yy >= 0)+pi;
        ang(xx >= 0 & yy > 0) = ang(xx >= 0 & yy > 0)+2*pi;
        mask = zeros(N_pixels);
        
        for i = 1:numSections
            phase = (axi(x(1+(3*(i-1)))).*defo(x(2+(3*(i-1)))).*gam(x(3+(3*(i-1)))));
            if i == 1 %First section lower bounded by 0
                mask(ang>0 & ang < cutoffs(i)) = phase(ang>0 & ang < cutoffs(i));
            elseif i < numSections
                mask(ang >= cutoffs(i-1) & ang < cutoffs(i)) = phase(ang >= cutoffs(i-1) & ang < cutoffs(i));
            elseif i == numSections %Last section bounded by maxNAa
                mask(ang >= cutoffs(i-1) & ang < 2*pi) = phase(ang >= cutoffs(i-1) & ang < 2*pi);
            end
        end
        mask = binarize(mask,0).*ap; %Binarize final result
        
        figure('units','normalized','outerposition',[0 0 1 1])
        imagesc(mask)
        title('Generated Mask')
        axis('square')
        saveas(gcf,[saveDir '/mask.jpg'])
        close()
        
    elseif apoType == 3
        phase = zeros(N_pixels);
        for i = 1:numSections
            phase = phase + weights(i)*((axi(x(1+(3*(i-1)))).*defo(x(2+(3*(i-1)))).*gam(x(3+(3*(i-1))))));
        end
        
        mask = binarize(phase,0).*ap; %Binarize final result
    end
    
%     %Erode and Dilate mask if desired (define minimum feature size)
%     if erodeRate == 0
%         mask = mask;
%     else
%         se = strel('disk',erodeRate);
%         dil = imdilate(mask,se);
%         mask = imerode(dil,se);
% %         disk = zeros(N_pixels);
% %         disk(xx.^2+yy.^2 <= erodeRate*dx.^2) = 1;
% %         mask = conv2(disk,mask,'same');
%         mask = binarize(mask./max(abs(mask(:))),0).*ap;
%     end
%     
    %mask(pupilGrid <= innerCircle) = -1; %Apply optical marker
    
    %Erode Mask
    mask = filter_mask_rings(mask,Params.erodeRate); %Use a custom filtering operation since erode/dilate is no bueno in polar coordinates
    
    
    %Find Axial PSF
    slice = zeros(N_pixels,N_pixels);
    slice(xx.^2+yy.^2 <= targ_size.^2) = 1; %10um sized object is target
    volume = repmat(slice,1,1,length(depth)); %on-axis object throughout depth
    
    %Load in neural volume. If the size doesn't match the pregenerated set,
    %remake a new one with the same params
    load neur.mat;
    if size(Neur,1) ~= size(mask,1)
        Neural.N = Params.N/2; %Make EVEN
        Neural.dx = dx; %Real space 'pixel' size
        Neural.depth = neurDepth; %Long Range
        %Neural.targetsize = targ_size;
        Neural.targetsizes = [4,6]*1e-6;
        Neural.ratios = [1,0]; %By setting a ratio to 0, that size is the spacing between neurons
        Neural.targetsig = 1;
        Neural.M = 2; %Number of targets
        Neural.backgroundsize = 1e-6;
        Neural.backgroundsig = Neural.targetsig*0.6;
        Neural.ratio = 0; %background to foreground per slice
        Neural.rand = 0; %Number of random targets added to scene (keeps background consistent)
        Neural.noise = 'None';
        Neural.var1 = 0.002; %Changes dymanics of selected noise (see function)
        Neural.var2 = 0.0004;

        %Generate stack representing neural volume
        Neur = generate_neural_volNoOverlap(Neural);
        %Neur = generate_neural_vol2D(Neural);
        Neur = padarray(Neur,[Params.N/4,Params.N/4],'both');
    end
    %Plot result with or without scattering
    if Params.ScatterFlag == 0
        [im,~] = propagate(Neur,depth,ap,defo,iF,F,ab); %Propagate through volume
        [~,vmask] = propagate(volume,depth,mask,defo,iF,F,ab);
        [immask,~] = propagate(Neur,depth,mask,defo,iF,F,ab);
    else
        [im,~] = propagate_pseudoscattering(Neur,neurDepth,ap,defo,iF,F,ab ,Params.ls); %Now with scattering! 
        [immask,~] = propagate_pseudoscattering(Neur,neurDepth,mask,defo,iF,F,ab ,Params.ls); %Now with scattering! 
        [~,vmask] = propagate_pseudoscattering(volume, depth, mask, defo, iF, F, ab ,Params.ls); %Now with scattering! 
    end
    
    %A = corrcoef(reshape(vmask,N_pixels*N_pixels,length(depth)));
    
    leng = (-N_pixels/2:N_pixels/2-1)*dx;

    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,3,1)
    imagesc(leng*10^6,leng*10^6,sum(Neur,3))
    title(['Vol Max Depth = ' num2str(max(neurDepth))],'FontSize',18)
    xlabel('Length, \mu m')
    ylabel('Length, \mu m')
    axis([-N_pixels/4 N_pixels/4 -N_pixels/4 N_pixels/4]*dx*10^6)
    axis('square')
    subplot(1,3,2)
    imagesc(leng*10^6,leng*10^6,im)
    title('Bare Ap Image','FontSize',18)
    xlabel('Length, \mu m')
    ylabel('Length, \mu m')
    axis([-N_pixels/4 N_pixels/4 -N_pixels/4 N_pixels/4]*dx*10^6)
    axis('square')
    subplot(1,3,3)
    imagesc(leng*10^6,leng*10^6,immask)
    title('Optimized Mask Image','FontSize',18)
    xlabel('Length, \mu m')
    ylabel('Length, \mu m')
    axis([-N_pixels/4 N_pixels/4 -N_pixels/4 N_pixels/4]*dx*10^6)
    axis('square')
    saveas(gcf,[saveDir '/res1.jpg'])
    close()
    
    %Compare to section
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    imagesc(real(mask))
    title(['Optimal Mask max depth ' num2str(targ_depth)],'FontSize',18)
    axis('square')
    %axis([[N_pixels/2+1-N_pixels/4 N_pixels/2+1+N_pixels/4 N_pixels/2+1-N_pixels/4 N_pixels/2+1+N_pixels/4]])
%     subplot(1,3,2)
%     imagesc(depth,depth,A)
%     title('Correlation between 10um objects','FontSize',18)
%     axis('square')
%     xlabel('depth, m')
%     ylabel('depth, m')
%     colorbar
    subplot(1,2,2)
    imagesc(depth*10^6,leng*10^6,squeeze(vmask(N_pixels/2+1,:,:)))
    title('Axial iPSF (Target Sized Particle)','FontSize',18)
    axis('square')
    xlabel('length, \mu m')
    ylabel('depth, \mu m')
    %colorbar
    saveas(gcf,[saveDir '/res2.jpg'])
    close()
    
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    imagesc(real(mask))
    title(['Vol Max Depth = ' num2str(targ_depth)],'FontSize',18)
%     xlabel('Length, m')
%     ylabel('Length, m')
    %axis([N_pixels/2+1-N_pixels/4 N_pixels/2+1+N_pixels/4 N_pixels/2+1-N_pixels/4 N_pixels/2+1+N_pixels/4])
    axis('square')
    subplot(1,2,2)
    imagesc(depth*10^6,leng*10^6,squeeze(vmask(N_pixels/2+1,:,:)))
    title('Axial iPSF (Target Sized Particle)','FontSize',18)
    axis('square')
    xlabel('length, \mu m')
    ylabel('depth, \mu m')
    colorbar
    saveas(gcf,[saveDir '/res3.jpg'])
    close()
    
    %Lets Predict Strehl Ratio and Plot it
    max_lens = max(max(abs(iF(F(slice).*acrr(ap))))); %Lens focal spot - only signal contained within boundary of orig obj (since we are not comparing iPSF)
    %max_lens(slice == 0) = 0;
    %max_lens = sum(max_lens,'all');
    
    %Now lets do this!
    tmp = max(max(vmask));
    %tmp(volume == 0) = 0; %Filter by obj size
    %tmp = sum(sum(tmp,1),2); %Sum across lateral dimensions
    strehl = tmp./max_lens; %BINARY strehl (approx)
    figure('units','normalized','outerposition',[0 0 1 1])
    hold on
    plot(depth*10^6,squeeze(strehl(1,1,:)),'linewidth',2)
    line([depth(1)*10^6 depth(end)*10^6], [max(squeeze(strehl(1,1,:)))/2 max(squeeze(strehl(1,1,:)))/2],'linewidth',2,'color','r','linestyle','--')
    xlabel('depth \mu m','fontsize',14)
    ylabel('Strehl','fontsize',14)
    title('Strehl Ratio of Obj per Depth','fontsize',16)
    legend('Mask DoF','Cutoff','fontsize',16)
    saveas(gcf,[saveDir '/Strehl.jpg'])
    close()
    
    vmask = vmask./max(vmask(:));
    vmask(vmask >= 0.5) = 1;
    vmask(vmask < 0.5) = 0;
    figure('units','normalized','outerposition',[0 0 1 1])
    imagesc(depth*10^6,leng*10^6,squeeze(vmask(N_pixels/2+1,:,:)))
    title('Recoverable Depth','FontSize',18)
    axis('square')
    xlabel('length, \mu m')
    ylabel('depth, \mu m')
    colorbar
    saveas(gcf,[saveDir '/recoverable_depth.jpg'])
    close()
    
    %Save the mask for easier loading in the future - since these arrays
    %are small, we have little to worry about! (unlike ND2P)
    save mask.mat mask
    
    %Save PSFs at some key depths - need to go back and rectify so its
    %evenly spaced over designed depth range
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,6,1)
    imagesc(real(mask))
    title('GA Mask')
    axis('square')
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    subplot(1,6,2)
    imagesc(log10(abs(iF(mask.*defo(0e-6))).^2))
    title('Log iPSF @ 0um')
    axis('square')
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    colormap hot
    subplot(1,6,3)
    imagesc(log10(abs(iF(mask.*defo(38e-6))).^2))
    title('@ 38um')
    axis('square')
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    colormap hot
    subplot(1,6,4)
    imagesc(log10(abs(iF(mask.*defo(78e-6))).^2))
    title('@ 78um')
    axis('square')
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    colormap hot
    subplot(1,6,5)
    imagesc(log10(abs(iF(mask.*defo(118e-6))).^2))
    title('@ 118um')
    axis('square')
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    colormap hot
    subplot(1,6,6)
    imagesc(log10(abs(iF(mask.*defo(158e-6))).^2))
    title('@ 160um')
    axis('square')
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    colormap hot
    saveas(gcf,[saveDir '/PSFs.jpg'])
    close()
    
end