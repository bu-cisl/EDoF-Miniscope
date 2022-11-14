function cost = Generate_Cost_Multiap_Selected_PointSources(x,fval,plotFlag,generateFlag)
%% Generate Simulation Parameters
    %addpath('Helpers','Objects');
    %Normalize weights in case inputs off
 
    %Load and generate grids
    
    load Objects/Params.Mat
    lambda = Params.lambda;
    n2 = Params.n2;
    targ_size = Params.targ_size; %Target particle size
    targ_depth = Params.targ_depth; %Target optimization depth
    %show_depth = Params.show_depth; %Additional obj simulation depth. We penalize light here
    maxDepth = Params.maxDepth; %Max depth in final images
    depth_step = Params.depth_step; %dz
    depth = 0:depth_step:targ_depth; %Depth vector for obj simulation
    neurDepth = 0:depth_step:maxDepth;
    N_pixels = Params.N;% Fixed - number of pixels on array
    saveDir = Params.name;
    
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
    W040 = Params.W040/3;
    W131 = 0;
    W222 = 0;
    W220 = 0;
    W311 = 0;
    f_Num = Params.fnum;
    f0 = 1/(2*(lambda*f_Num)); %This quantity is conserved due to the relationship between cutoff frequency and NA
    W = Seidal_Coeff(0,0,uu./f0,vv./f0,Wd,W040,W131,W222,W220,W311,0); %On axis aberrations
    ab = exp(-1i*k*W*lambda).*ap;
    
    %Develop mask
    numSections = ceil(length((x))/3); %We can use the number of entries in x to determine how many foci we optimized over
    
    %Mask is summation of all submasks!
    phase = zeros(N_pixels);
    for i = 1:numSections
        phase = phase + (axi(x(1+(3*(i-1)))).*defo(x(2+(3*(i-1)))).*gam(x(3+(3*(i-1)))));
    end

    mask = binarize(-1*phase,0).*ap; %Binarize final result

    %Erode Mask
    mask = filter_mask_rings(mask,Params.erodeRate); %Use a custom filtering operation since erode/dilate is no bueno in polar coordinates
    
    %% Finish Sim
    slice = zeros(N_pixels,N_pixels);
    slice(xx.^2+yy.^2 < targ_size.^2) = 1; %Set our target object
    volume = repmat(slice,1,1,length(depth)); %on-axis object throughout depth
    
    %Discern 3D Profile over desired range
    if Params.ScatterFlag == 0
        [~,vol] = propagate(volume,depth,mask,defo,iF,F,ab); %Propagate through volume
    else
        [~,vol] = propagate_pseudoscattering(volume, depth, mask, defo, iF, F, ab ,Params.ls); %Now with scattering! 
    end
    
    vol = vol./max(vol(:));
    vol(vol > 0.5) = 1; %Set a threshold where are values above are equally penalized (else small breaches in desired conditions are ignored)
    vol(vol <= 0.5) = 0;
    cost = -1*sum(vol(:,:,1:floor(targ_depth/depth_step)+1),'all') + 5*sum(vol(:,:,floor(targ_depth/depth_step)+2:end),'all');
    %cost = -1*sum(vol,'all');
     
    %% Test Mask In Neural Volume vs Bare Aperture
    if plotFlag
        if generateFlag
            Neural.N = N_pixels/2; %Make EVEN
            Neural.dx = dx; %Real space 'pixel' size
            Neural.depth = neurDepth; %Long Range
            %Neural.targetsize = targ_size;
            Neural.targetsizes = [4,6]*1e-6;
            Neural.ratios = [1,0]; %By setting a ratio to 0, that size is the spacing between neurons
            Neural.targetsig = 1;
            Neural.M = 3; %Number of targets
            Neural.backgroundsize = 1e-6;
            Neural.backgroundsig = Neural.targetsig*0.6;
            Neural.ratio = 0; %background to foreground per slice
            Neural.rand = 0; %Number of random targets added to scene (keeps background consistent)
            Neural.noise = 'None';
            Neural.var1 = 0.2; %Changes dymanics of selected noise (see function)
            Neural.var2 = 0.004;

            %Generate stack representing neural volume
            Neur = generate_neural_volNoOverlap(Neural);
            %Neur = generate_neural_vol2D(Neural);
            Neur = padarray(Neur,[N_pixels/4,N_pixels/4],'both');
            figure('units','normalized','outerposition',[0 0 1 1])
            imagesc(sum(Neur,3))
            title(['Total depth = ' num2str(max(Neural.depth))])
        else
            load neur.mat; % I presaved a neural volume for universal comparison
        end

        %% Propagate
        %Since we don't calculate this per iteration, need it here in this
        %routine
        %Discern 3D Profile over desired range
        if Params.ScatterFlag == 0
            [~,vol] = propagate(volume,depth,mask,defo,iF,F,ab); %Propagate through volume
        else
            [~,vol] = propagate_pseudoscattering(volume, depth, mask, defo, iF, F, ab ,Params.ls); %Now with scattering! 
        end
        
        if Params.ScatterFlag == 0
            [im,~] = propagate(Neur,neurDepth,ap,defo,iF,F,ab);
            [immask,~] = propagate(Neur,neurDepth,mask,defo,iF,F,ab);
        else
            [im,~] = propagate_pseudoscattering(Neur,neurDepth,ap,defo,iF,F,ab ,Params.ls); %Now with scattering! 
            [immask,~] = propagate_pseudoscattering(Neur,neurDepth,mask,defo,iF,F,ab ,Params.ls); %Now with scattering! 
        end
        %% Plot Results
        leng = (-N_pixels/2:N_pixels/2-1)*dx;

        figure('units','normalized','outerposition',[0 0 1 1]) %Make image full screen
        subplot(1,3,1)
        imagesc(leng*10^6,leng*10^6,sum(Neur,3))
        title(['Vol Max Depth = ' num2str(max(neurDepth))],'FontSize',18)
        xlabel('Length, m')
        ylabel('Length, m')
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
        imagesc(leng,leng,immask)
        title('Optimized Mask Image','FontSize',18)
        xlabel('Length, m')
        ylabel('Length, m')
        axis([-N_pixels/4 N_pixels/4 -N_pixels/4 N_pixels/4]*dx)
        axis('square')
        
        %Save the figure('units','normalized','outerposition',[0 0 1 1]) - first check if image exists. We need to do this
        %since this code runs for each focus we optimize over. So we if
        %only declare one focus, great! Else we want to see the effect of
        %each independent focus.
        saveme = false; %This is not a cry for help embedded in my code. I'm just literal with my variable names.
        filebase = 'ImageCompareFoci';
        while saveme == false
            name = [saveDir '/' filebase num2str(i) '.jpg'];
            if ~isfile(name)
               saveas(gcf,name)
               saveme = true;
            else
                i = i+1;
            end
        end
        close()
        
        %Compare to cross section
        figure('units','normalized','outerposition',[0 0 1 1])
        subplot(1,2,1)
        imagesc(real(mask))
        title(['Optimal Mask max depth ' num2str(max(depth))],'FontSize',18)
        axis('square')
%         subplot(1,3,2)
%         imagesc(depth,depth,A)
%         title('Correlation between 10um objects','FontSize',18)
%         axis('square')
%         xlabel('depth, m')
%         ylabel('depth, m')
%         colorbar
        subplot(1,2,2)
        imagesc(depth*10^6,leng*10^6,squeeze(vol(N_pixels/2+1,:,:))./max(vol(:)))
        title('Axial iPSF','FontSize',18)
        axis('square')
        xlabel('length, \mu m')
        ylabel('depth, \mu m')
        colorbar
        saveme = false; %This is not a cry for help embedded in my code. I'm just literal with my variable names.
        filebase = 'AxialInfoFoci';
        while saveme == false
            name = [saveDir '/' filebase num2str(i) '.jpg'];
            if ~isfile(name)
               saveas(gcf,name)
               saveme = true;
            else
                i = i+1;
            end
        end
        close()

        %Plot FWHM as a function of depth - useful since radially symmetric
        %design
        cross = zeros(1,length(depth));
        for j = 1:length(depth)
            temp = vol(:,N_pixels/2+1,j);
            temp = temp./max(temp(:));
            cross(1,j) = length((find(temp >= 0.5)))*dx; %FWHM
        end
        
        figure('units','normalized','outerposition',[0 0 1 1])
        plot(depth,cross)
        title('Maxthresh FWHM/Depth')
        xlabel('Depth')
        ylabel('FWHM (m)')
        saveme = false; %This is not a cry for help embedded in my code. I'm just literal with my variable names.
        filebase = 'RecoverableDepthFoci';
        while saveme == false
            name = [saveDir '/' filebase num2str(i) '.jpg'];
            if ~isfile(name)
               saveas(gcf,name)
               saveme = true;
            else
                i = i+1;
            end
        end
        close()
        
        
        %Now lets look at the reference no ap case
        if Params.ScatterFlag == 0
            [~,vol] = propagate(volume,depth,ap,defo,iF,F,ab); %Propagate through volume
        else
            [~,vol] = propagate_pseudoscattering(volume, depth, ap, defo, iF, F, ab ,Params.ls); %Now with scattering! 
        end
        
                %Compare to cross section
        figure('units','normalized','outerposition',[0 0 1 1])
        subplot(1,2,1)
        imagesc(real(ap))
        title(['Aperture max depth ' num2str(max(depth))],'FontSize',18)
        axis('square')
%         subplot(1,3,2)
%         imagesc(depth,depth,A)
%         title('Correlation between 10um objects','FontSize',18)
%         axis('square')
%         xlabel('depth, m')
%         ylabel('depth, m')
%         colorbar
        subplot(1,2,2)
        imagesc(depth*10^6,leng*10^6,squeeze(vol(N_pixels/2+1,:,:))./max(vol(:)))
        title('Axial iPSF','FontSize',18)
        axis('square')
        xlabel('length, \mu m')
        ylabel('depth, \mu m')
        colorbar
        saveme = false; %This is not a cry for help embedded in my code. I'm just literal with my variable names.
        filebase = 'BareApReference';
        while saveme == false
            name = [saveDir '/' filebase num2str(i) '.jpg'];
            if ~isfile(name)
               saveas(gcf,name)
               saveme = true;
            else
                i = i+1;
            end
        end
        close()
        cost = cost
    end
end