%% Make Parameters for simulation
%Add appropriate subdirectories:
%{
    Helpers: functions that assist in the execution of the genetic
    algorithm
    Objects: Generated data to assist in the running of the genetic
    algorithm
    Supp_Analysis: Other scripts that assist in the simulation of factors
    shown in the paper supplement (misalignment, wavelet analysis, etc)
%}
addpath('Helpers','Objects','Supp_Analysis')

%%%%%%%%%%%%%%%%%%%%%%%%%%% User Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.lambda = 509e-9; %Center wavelength
Params.targ_size = 5e-6; %Target particle size
Params.targ_depth = 150e-6; %Optimization depth 
Params.show_depth = 120e-6; %Visualization depth for plotting cross sections. Also determines extra room cost is determined within
Params.maxDepth = 150e-6; %Total simulation depth
Params.n2 = 1; %Index of refraction for brain (mean)
Params.depth_step = 2e-6; %Depth resolution
Params.NA0 = 0.55; %NA of lens
Params.f = 0.0017; %Focal length
Params.fnum = 0.96; %F number of lens 
Params.mag = 8.71; %System magnification
Params.pix = 2.2e-6; %Sensor array pixel size - I increased this by a factor of two so the magnified pixel size is ~lambda/4 in accordance with the scattering BPM
Params.dx = Params.pix/Params.mag; %Image Space discretization
Params.W040 = 29.6; %Seidel aberration: spherical (in wavelength deviations)
Params.erodeRate = 1; %Remove rings less than this size (in pixels)

%Make a base folder name to save all the results
filebase = 'Optimize_DOE';

%Turn Scatteing on or off
Params.ScatterFlag = 1;

%Setup scattering properties
Params.ls = 100e-6; %scattering length in meters!

%From the aforementioned parameters we may deduce the number of pixels
%required to properly capture the defocused PSF at the maximum depth -This
%is a rough geometric argument with a slight scaling factor to account for
%the wave nature of light
Params.N = 2*floor(1.5*(Params.maxDepth*tan(asin(Params.NA0/Params.n2)))/(Params.pix/Params.mag)); %Ensure this number is closest even number (Lei likes an even array)

%Set a minimum value for optimization - do this so it fits our
%pre-generated neural volume
if Params.N < 1000
    Params.N = 1000;
end

%% Now Run simulation and save important results!

%Make a new directory to save files
%filebase = 'Optimization in Different Media Water ls 250um';
if Params.ScatterFlag == 1
    filebase = [filebase ' Scattering'];
end
makeFlag = false;
i = 1;

%Check if its exists and if not, make it 
while makeFlag == false
    name = [filebase num2str(i)];
    if ~exist(name, 'dir')
       mkdir(name)
       makeFlag = true;
    else
        i = i+1;
    end
end
Params.name = name;

%Lets make our real in fourier space grids here so we don't need to make
%them each iteration. We could preallocate other quantities but the
%generation time versus the save/load time is trivial.
%Generate grids
[xx,yy] = meshgrid([-floor(Params.N/2):(floor(Params.N/2))-1]*Params.dx); %Image plane representation of image plane -> res/2 for inco
[uu,vv] = meshgrid([-floor(Params.N/2):(floor(Params.N/2))-1]*1/(Params.N*Params.dx)); %du defined by 1/FoV, Fourier plane
R = sqrt(xx.^2+yy.^2); %Radial spatial grid
%Rmax = maxi(R(Params.NA <= Params.NAs)); %Cutoff from Fourier space recipriocals
NAx = uu*Params.lambda; %converting to NA space (alottable angles) -> unitless! easy of scaling and design
NAy = vv*Params.lambda;

%Define Aperture
NA = sqrt(NAx.^2+NAy.^2);
ap = zeros(Params.N);
ap(NA <= Params.NA0) = 1;

%Convert NA space into the physical grid for the pupil plane (grin lens
%is 1.8mm in diameter, this is to enforce the innerCircle hyperparam)
pupilGrid = 1.8e-3*NA/Params.NA0; %We know NA0 is maximum size of GRIN, which is 1.8mm


%Setup Evanescence cutoff for free space propagation - this limits
%kernel from blowing up (i.e becoming physically invalid) and works in
%conjunction with the aperture in terms of the detected signal
eva = logical((NA/Params.lambda).^2 < (Params.n2/Params.lambda).^2);

Params.NA = NA;
Params.xx = xx;
Params.yy = yy;
Params.uu = uu;
Params.vv = vv;
Params.eva = eva;
Params.ap = ap;
Params.dp = pupilGrid(1) - pupilGrid(2); %Discretization on pupil plane

save Objects/Params.Mat Params

%% Run the GA
%{
    For more information on the matlab genetic algorithm: https://www.mathworks.com/help/gads/genetic-algorithm.html
%}

%User input
%Where we want light
PosRange = [0  100e-6];

%Weighting of the foci
maskvars = zeros(1,3*length(PosRange)/2);
NAmax = Params.NA0;

tic
%modify phase functions to make these values relatively similar in
%scale/sensitivity.
LBphase = [1*Params.lambda/Params.n2,-100e-6,-9e-5];
UBphase = [75*Params.lambda/Params.n2,-10e-6,0];
for i = 1:length(PosRange)/2
    %For each foci, lets define our positive and negative points
    %Write_Axial_Points_One_Foci(PosRange(2*i-1:2*i));
    ObjectiveFunction = @Optimize_Mask_Multiap_Selected_PointSources; %Uses info generated by Generate_Postive_and_Negative_Axial_Points
    nvars = 3; %Phase arguments - since only optimize one at a time, few nvars so rapid optimization

    %Setup options for GA
    options = optimoptions('ga','CrossoverFcn',{@crossoverscattered}); %Children weighted mean of parents
    options.FunctionTolerance = 0; %Stops if improvement less than this amount
    options.MaxStallTime = inf;
    options.MaxTime = inf; %In seconds, so 3600 = 1 hr.
    options.MaxGenerations = 3;
    options.PopulationSize = 1*nvars; %Decrease population to create more generations - default allows only 9 gen/hr
    options.EliteCount = ceil(0.2*options.PopulationSize); %Percent population carried onto the next generation
    options.CrossoverFraction = 0.1; %Need more mutant children to push out of local minima
    options.MutationFcn = {@mutationuniform,0.1}; %Default is gaussian with shrinkage (bad as less variation as generations progress).
    options.PlotFcn = @gaplotbestf;
    [optvars,cost] = dsGA_First_Step(1, 1, NAmax, 0, ObjectiveFunction, LBphase, UBphase, options);
    cost = Generate_Cost_Multiap_Selected_PointSources(optvars,cost,1,1); %Will plot optimal mask
    
    %Save result
    maskvars(3*(i-1)+1:3*i) = optvars;
    
    %Save the GA result
    fh = findobj( 'Type', 'Figure', 'Name', 'Genetic Algorithm'); %Find the image with the cost function (automatically generated by alg)
    saveas(fh,[name '/GA_Cost_Foci' num2str(i) '.jpg'])
end
toc

weights = ones(1,length(PosRange)/2)/(length(PosRange)/2); %Weighted summation of parts
apoType = 3; %Phase combination type
erodeRate = Params.erodeRate; %Fixed in cost function
W040 = Params.W040; %Spherical ab to test on - usually just optimization parameter

Stitch_and_Show_Mask_Results(maskvars, weights, apoType, erodeRate, W040, Params.name);
save([name '/optvars.mat'], 'optvars')
a = 'done'

%% Perform Perturbation Analysis on Result
saveStep = 1;
maxLat = 450; %Shift in pixels
latStep = 5; %Shift in pixels
maxAx = 5000; %Shift in microns
axStep = 100; %Shift in microns

%Run perturbation analysis
Analyze_Lateral_Mask_Misalignment(saveStep,maxLat,latStep)
%Analyze_Axial_Mask_Misalignment(saveStep,maxAx,axStep) %This is much more computationally intensive

% %% Test wavelet transform
% rho = 8; %Waist in pixels
% filt = 2/(sqrt(3)*pi^(1/4)).*(1-1/2*(r/(rho)).^2).*exp(-(r.^2./(2.*(rho).^2))); %Second derivative of gaussian
% filt = filt./max(filt(:)); %Normalize to 1
