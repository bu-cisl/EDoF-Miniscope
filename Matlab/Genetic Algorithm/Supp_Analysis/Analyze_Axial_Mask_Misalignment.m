function Analyze_Axial_Mask_Misalignment(saveStep,maxAx,axStep)
    %{
        The Input parameters are:
            1.) maxAx: Maximum Axial shift we will test in MICRONS
            2.) axStep: Axial stepsize in MICRONS (we will round to closest
            factor)
    %}

%Where as the last script we can use a relaxed sampling condition, to
%misalign the phase mask as a REAL SPACE object to defocus it, we need to
%match our sampling for that object. I.E we need a 0.55 NA across a 1.8mm
%FoV. Since our fourier space allows for complex components for intensity
%spatial frequencies, we can't just ignore this phenomenon

%% Physical grids and parameters
    %Setup the appropriate grids and load appropriate quantities from our
    %parameter file
    addpath('Helpers','Objects')
    load Objects/Params.Mat
    %Params.ScatterFlag = 0;
     %Load and generate grids

    lambda = Params.lambda;
    n2 = Params.n2;
    targ_size = Params.targ_size;
    targ_depth = Params.targ_depth;
    maxDepth = Params.maxDepth;
    depth_step = Params.depth_step;
    depth = 0:depth_step:targ_depth;
    %neurDepth = 0:depth_step:maxDepth; %Fix the max depth, 91 to match loaded Neural volume dimensions
    N = Params.N; % Fixed - number of pixels on array
    saveDir = Params.name;
    
    %Load our handles
    F = @(x) fftshift(fft2(ifftshift(x)));
    iF = @(x) fftshift(ifft2(ifftshift(x)));
    acrr= @(x) iF(conj(F(x)).*F(x)); %autocorrelation using fourier transform (much faster)
    maxi = @(x) x./max(x(:)); %Normalize an N-dim array
    
    %Physical parameters of system considering magnification and sensor size
    %(assume not pixel limited)
    NA0 = Params.NA0; %NA of objective lens - sets resolution
    mag = Params.mag; %Leave as 1 else need to rescale all generated masks (though can easily improve code to generalize)
    pix = Params.pix; %Pixel Pitch
    dx = pix/mag; %Effective pixel size
    
    %However we cannot use the same grids in this simulation!! We need to
    %match our imaging requirements across the FoV of our mask. This way we
    %match our grid dimensions as we go across the different steps.
    
    %First lets generate the REAL SPACE array for our phase mask
    d_phase = 1.8e-3; %1.8e-3 equates to 0.55 NA for our GRIN lens since we have an INCOHERENT problem, we need to simulate TWICE this amount since we will eventually correlate this grid to the FOURIER grid for our PSF
    N = 2*round(d_phase/dx); %Make it closest even number for ease & correct dft discretization
    
    %Make the appropriate NA grids - we don't need the spatial grids - will
    %upsample our pregenerated mask array to the correct values
    [NAx,NAy] = meshgrid([-N/2:(N/2)-1]*lambda/(N*dx)); %Directly calculate NA grids
    
    %Define Aperture
    NA = sqrt(NAx.^2+NAy.^2);
    ap = logical(NA < NA0); %Logical arrays take 1/4 the space - going to be sure to use for these big arrays (tips from a later project!)

    %Define evanescence
    eva = logical((NA/Params.lambda).^2 < (Params.n2/Params.lambda).^2);

    %Male our defocus handles
    k = 2*pi*n2/lambda;
    defo = @(d) exp(1i*k*d.*sqrt(1-(NA/n2).^2).*eva).*eva; %Angular spectrum kernel in a new media
    
    %Load the mask & upsample it so it makes physical sense (i.e. real
    %space discretization matches our NA grid)
    [xx,yy] = meshgrid(linspace(min(Params.xx(:)),max(Params.xx(:)),N)); %Need for interp
    load mask.mat
    mask = binarize(interp2(Params.xx,Params.yy,mask,xx,yy),0); %Since interpolation does, well, interpolation need to rebinarize
    
    %% Phase Abberations and 3D Volume Matrix
    %Now lets prepare some useful matricies
    %First we need our on-axis aberrations - GRIN has prominent spherical
    %and minimal others
    Wd = 0;
    W040 = Params.W040; %Match the spherical aberration we previously optimized over
    W131 = 0;
    W222 = 0;
    W220 = 0;
    W311 = 0;
    f_Num = Params.fnum;
    f0 = 1/(2*(lambda*f_Num));
    W = Seidal_Coeff(0,0,NAx./(lambda*f0),NAy./(lambda*f0),Wd,W040,W131,W222,W220,W311,0); %On axis aberrations
    ab = exp(-1i*k*W*lambda).*ap;
    clear NAx NAy %Not useful & big
    
    %Next we need our 3D Volume
    slice = zeros(N,N);
    slice(xx.^2+yy.^2 <= targ_size.^2) = 1; %Make an extended volume of our target size
    volume = repmat(slice,1,1,length(depth)); %on-axis object throughout depth
    clear xx yy
    
   %% Axial Shifts
    %Now lets do the Axial shift
    axShifts = (0:axStep:maxAx)*10^-6; %Determine the vector of axial shifts
    axErr = zeros(1,length(axShifts));
    numSliceXZ = zeros(1,length(axShifts)); %Since our cost function isn't super intuitive, we can also count the number of slices that have a discerable signal as we move. This is effectively our immediate depth of field (numslices*dz)
    numSliceYZ = zeros(1,length(axShifts));
    subdir = 'Axial Shifts';
    mkdir([saveDir '/' subdir])
    for i = 1:length(axShifts)
        tmp = iF(F(mask).*defo(-1*axShifts(i))).*ap; %Here, the mask is the OBJECT and we backpropagate it to the fourier plane to then act as a FILTER (why the grids need to match)
        
        %Now we can actually propagate
        if Params.ScatterFlag == 0
            [~,vol] = propagate(volume,depth,tmp,defo,iF,F,ab); %Propagate
        else
            [~,vol] = propagate_pseudoscattering(volume,depth,tmp,defo,iF,F,ab,Params.ls); %Propagate
        end
        
        %Calculate the cost
        vol = vol./max(vol(:));
        vol(vol > 0.5) = 1; %Set a threshold where are values above are equally penalized (else small breaches in desired conditions are ignored)
        vol(vol <= 0.5) = 0;
        cost = -1*sum(vol.*volume,'all');
        axErr(i) = cost;
        
        %Now determine the number of recoverable slices (or DoF)
        xz = squeeze(vol(N/2+1,:,:)); %Take cross section as 2D slice in middle of 3D PSF (I know its not really the
        yz = squeeze(vol(:,N/2+1,:));
        numSliceXZ(i) = max(sum(xz,2)); %Here we sum across depth to get the total number of 'active slices' (since we already applied a boolean 0 or 1 per depth). Then we take the maximum as the best case scenario.
        numSliceYZ(i) = max(sum(yz,2)); %Here we sum across depth to get the total number of 'active slices' (since we already applied a boolean 0 or 1 per depth). Then we take the maximum as the best case scenario.
        
        %Save the images at each iter
        figure('visible','off');
        %figure
        subplot(1,2,1)
        imagesc(real(tmp))
        axis('square')
        title(['Axial Shift: ' num2str(round(axShifts(i)*10^6)) '\mu m'],'FontSize',16)
        colorbar
        caxis([-1 1])
        subplot(1,2,2)
        imagesc(depth*10^6,leng,xz)
        title(['xz x-sec:' num2str(round(axShifts(i)*10^6)) '\mu m'],'FontSize',16)
        axis('square')
        xlabel('z, \mu m')
        ylabel('x, \mu m')
        colorbar
        filebase = ['xzax'  num2str(i)];
        name = [saveDir '/' subdir '/' filebase '.svg'];
        saveas(gcf,name)
        close()
                
        figure('visible','off');
        subplot(1,2,1)
        imagesc(real(tmp))
        axis('square')
        title(['Axial Shift: ' num2str(round(axShifts(i)*10^6)) '\mu m'],'FontSize',16)
        colorbar
        caxis([-1 1])
        subplot(1,2,2)
        imagesc(depth*10^6,leng*10^6,yz)
        title(['yz x-sec:' num2str(round(axShifts(i)*10^6)) '\mu m'],'FontSize',16)
        axis('square')
        xlabel('z, \mu m','FontSize',14)
        ylabel('x, \mu m','FontSize',14)
        colorbar
        
        filebase = ['yzax'  num2str(i)];
        name = [saveDir '/' subdir '/' filebase '.svg'];
        saveas(gcf,name)
        close()  
    end
    
    figure('visible','off');
    %figure
    t = tiledlayout(1,1);
    ax1 = axes(t);
    %First set of axis - DoF vs. Lateral Displacement
    plot(ax1,axShifts*10^6,numSliceYZ*depth_step*10^6,'LineStyle','-','Color','r','linewidth',2)
    hold on
    plot(ax1,axShifts*10^6,numSliceXZ*depth_step*10^6,'LineStyle','--','Color','b','linewidth',2)
    %title(ax1,'DoF vs. Lat. Shift','fontsize',18, 'Position',[0 0.5])
    xlabel(ax1,'Axial Shift (\mu m)','fontsize',16)
    ylabel(ax1,'DoF (\mu m)','fontsize',16)
    %axis('square')
    axis([0 round(max(axShifts(:))*10^6) 0 max(numSliceXZ(:)*depth_step*10^6)])
    legend(ax1,'yz DoF','xz DoF','fontsize',12)
    
    %Second Axis - Percent DoF vs. NA Displacement
    ax2 = axes(t);
    %plot(ax2,latShifts./max(latShifts(:))*0.55*(max(latShifts(:))*dp)/1.8e-3,numSliceYZ./max(numSliceYZ(:))*100,'LineStyle','-','Color','r','linewidth',1.5) %Convert to % per NA (scaled to max aperture size)
    %plot(ax2,latShifts./max(latShifts(:))*0.55*(max(latShifts(:))*dp)/1.8e-3,numSliceXZ./max(numSliceXZ(:))*100,'LineStyle','--','Color','b','linewidth',1.5)
    xlabel(ax2,'Axial Shift (Wavelengths)','fontsize',16)
    ylabel(ax2,'% EDoF Retained','fontsize',16)
    axis([0 round(max(axShifts(:)./lambda)) min(numSliceYZ(:))/max(numSliceYZ(:))*100 100])
    %axis('square')
    ax2.XAxisLocation = 'top';
    ax2.YAxisLocation = 'right';
    ax2.Color = 'none';
    ax1.Box = 'off';
    ax2.Box = 'off';
    hold off
    
    filebase = 'RecSlicesPerAxShift';
    name = [saveDir '/' subdir '/' filebase '.svg'];
    saveas(gcf,name)
    close()

end