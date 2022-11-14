% The purpose of this script is to analyze the tolerance of the designed
% phase mask. We will analyze laterial misalignment (along one dimension
% since the mask is circularly symmetric), axial misalignment, chromatic aberration and etch
% error.

function Analyze_Lateral_Mask_Misalignment(saveStep,maxLat,latStep)
    %{
        The Input parameters are:
            1.) maxLat: Maximum lateral shift we will test in PIXELS
            2.) latStep: Lateral stepsize in PIXELs
            3.) maxAx: Maximum Axial shift we will test in MICRONS
            4.) axStep: Axial stepsize in MICRONS (we will round to closest
            factor)
    %}

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
    dp = Params.dp;
    
    %Load the appropriate grids
    ap = Params.ap; %This doesn't take up any more space, simply adds a shorter pointer to this quantity
    xx = Params.xx;
    yy = Params.yy;
    uu = Params.uu;
    vv = Params.vv;
    NA = Params.NA;
    eva = Params.eva;
    
    k = 2*pi*n2/lambda;
    defo = @(d) exp(1i*k*d.*sqrt(1-(NA/n2).^2).*eva).*eva; %Angular spectrum kernel in a new media

%% Phase Error and 3D Volume Matrix
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
    W = Seidal_Coeff(0,0,uu./f0,vv./f0,Wd,W040,W131,W222,W220,W311,0); %On axis aberrations
    ab = exp(-1i*k*W*lambda).*ap;
    
    %Next we need our 3D Volume
    slice = zeros(N,N);
    slice(xx.^2+yy.^2 <= targ_size.^2) = 1; %Make an extended volume of our target size
    volume = repmat(slice,1,1,length(depth)); %on-axis object throughout depth
    
    %Now lets load the mask and test it!
    load mask.mat
    figure
    imagesc(mask)
    mask(NA > NA0) = 1; %Basically assign all values outside the pupil to 'bare aperture' so when we circshift, it exposure bare aperture
    
    %Step 1: Test Lateral tolerence (Previous experience indicates this
    %will be the most sensitive parameter)
    %First, make sure we don't overlap our pupil with our maximum shift
    if maxLat >= N
        maxLat = N-1;
    end

%% Perform Lateral Shifts. Save Results
    leng = (-N/2:N/2-1)*dx;    
    latShifts = 0:latStep:maxLat; %Determine the vector of laterial shifts
    latErr = zeros(1,length(latShifts));
    numSliceXZ = zeros(1,length(latShifts)); %Since our cost function isn't super intuitive, we can also count the number of slices that have a discerable signal as we move. This is effectively our immediate depth of field (numslices*dz)
    numSliceYZ = zeros(1,length(latShifts));
    
    %Lets do laterial shifts first
    j = 0; %Used for saving stacks
    subdir = 'Lateral Shifts';
    mkdir([saveDir '/' subdir])
    for i = 1:length(latShifts)
        tmp = circshift(mask,latShifts(i),2).*ap; %Shift the mask and apodize
        if Params.ScatterFlag == 0
            [~,vol] = propagate(volume,depth,tmp,defo,iF,F,ab); %Propagate
        else
            [~,vol] = propagate_pseudoscattering(volume,depth,tmp,defo,iF,F,ab,Params.ls); %Propagate
        end
        
        %Calculate the cost
        vol = vol./max(vol(:));
        vol(vol > 0.5) = 1; %Set a threshold where are values above are equally penalized (else small breaches in desired conditions are ignored)
        vol(vol <= 0.5) = 0;
        cost = -1*sum(vol.*volume,'all'); %Here we do the effective DoF as defined by the 50% cutoff in the Strel ratio
        latErr(i) = cost;
        
        %Now determine the number of recoverable slices (or DoF)
        xz = squeeze(sum(vol,1)); %Take projection as 2D slice across one dim of the vol
        yz = squeeze(sum(vol,2));
        xz(xz > 0) = 1; %Filter scaling by projection.
        yz(yz > 0) = 1;
        numSliceXZ(i) = max(sum(xz,2)); %Here we sum across depth to get the total number of 'active slices' (since we already applied a boolean 0 or 1 per depth). Then we take the maximum as the best case scenario.
        numSliceYZ(i) = max(sum(yz,2)); %Here we sum across depth to get the total number of 'active slices' (since we already applied a boolean 0 or 1 per depth). Then we take the maximum as the best case scenario.
        
        %Save the images at each iter
        figure('visible','off');
        subplot(1,2,1)
        imagesc(real(tmp))
        axis('square')
        title(['Lat Shift: ' num2str(round(latShifts(i)*dp*10^6)) '\mu m'],'FontSize',16)
        colorbar
        caxis([-1 1])
        subplot(1,2,2)
        imagesc(depth*10^6,leng*10^6,xz)
        title(['xz x-sec:' num2str(round(latShifts(i)*dp*10^6)) '\mu m'],'FontSize',16)
        axis('square')
        xlabel('z, \mu m')
        ylabel('x, \mu m')
        colorbar
        
        filebase = ['xzlat'  num2str(i)];
        name = [saveDir '/' subdir '/' filebase '.svg'];
        saveas(gcf,name)
        close()
        
        %Save full 3D PSF to view in FIJI
        vol = vol./max(vol(:));
        if mod(i,saveStep) == 0
            %Make a folder to hold each stack
            dirBase = ['LatShift ' num2str(j*saveStep*dp) 'mm'];
            name = [saveDir '/' subdir '/' dirBase];
            mkdir(name)
            
            %Iterate through stack and save as tiff
            for k = 1:size(vol,3)
                imwrite(squeeze(vol(:,:,k)),[name '/' 'LatShift' num2str(j*saveStep) 'mm.tif'],'WriteMode','append','Compression','none');
            end
            
            %Iterate for next folder
            j = j+1;
            
        end        
                
        figure('visible','off');
        subplot(1,2,1)
        imagesc(real(tmp))
        axis('square')
        title(['Lat Shift: ' num2str(round(latShifts(i)*dp*10^6)) '\mu m'],'FontSize',16)
        colorbar
        caxis([-1 1])
        subplot(1,2,2)
        imagesc(depth*10^6,leng*10^6,yz)
        title(['yz x-sec:' num2str(round(latShifts(i)*dp*10^6)) '\mu m'],'FontSize',16)
        axis('square')
        xlabel('z, \mu m','FontSize',14)
        ylabel('x, \mu m','FontSize',14)
        colorbar
        
        filebase = ['yzlat'  num2str(i)];
        name = [saveDir '/' subdir '/' filebase '.svg'];
        saveas(gcf,name)
        close() 
    end
    
    %Save results
%     %figure('visible','off');
%     figure
%     plot(latShifts*dp*10^6,latErr)
%     title('Cost vs. Lat. Shift')
%     xlabel('Lat Shift \mu m')
%     ylabel('Cost')
%     filebase = 'CostPerLatShift';
%     name = [saveDir '/' subdir '/' filebase '.svg'];
%     saveas(gcf,name)
%     close()
    
    figure('visible','off');
    %figure
    t = tiledlayout(1,1);
    ax1 = axes(t);
    %First set of axis - DoF vs. Lateral Displacement
    plot(ax1,latShifts*dp*10^6,numSliceYZ*depth_step*10^6,'LineStyle','-','Color','r','linewidth',2)
    hold on
    plot(ax1,latShifts*dp*10^6,numSliceXZ*depth_step*10^6,'LineStyle','--','Color','b','linewidth',2)
    %title(ax1,'DoF vs. Lat. Shift','fontsize',18, 'Position',[0 0.5])
    xlabel(ax1,'Lat Shift (\mu m)','fontsize',16)
    ylabel(ax1,'DoF (\mu m)','fontsize',16)
    %axis('square')
    legend(ax1,'yz DoF','xz DoF','fontsize',12)
    
    %Second Axis - Percent DoF vs. NA Displacement
    ax2 = axes(t);
    %plot(ax2,latShifts./max(latShifts(:))*0.55*(max(latShifts(:))*dp)/1.8e-3,numSliceYZ./max(numSliceYZ(:))*100,'LineStyle','-','Color','r','linewidth',1.5) %Convert to % per NA (scaled to max aperture size)
    %plot(ax2,latShifts./max(latShifts(:))*0.55*(max(latShifts(:))*dp)/1.8e-3,numSliceXZ./max(numSliceXZ(:))*100,'LineStyle','--','Color','b','linewidth',1.5)
    xlabel(ax2,'Lat Shift (NA)','fontsize',16)
    ylabel(ax2,'% EDoF Retained','fontsize',16)
    axis([0 0.55*(max(latShifts(:))*dp)/1.8e-3 min(numSliceYZ(:))/max(numSliceYZ(:))*100-0.1 100])
    %axis('square')
    ax2.XAxisLocation = 'top';
    ax2.YAxisLocation = 'right';
    ax2.Color = 'none';
    ax1.Box = 'off';
    ax2.Box = 'off';
    hold off
    
    filebase = 'RecSlicesPerLatShift';
    name = [saveDir '/' subdir '/' filebase '.svg'];
    saveas(gcf,name)
    close()
 
end