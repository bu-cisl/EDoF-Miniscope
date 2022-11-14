%Simplifies Neural Volume generation by eliminating overlap all together.

%M is total particles per slice assigned with particles of size x_n based
%on the user set ratios

%Note if the user inputs M = 0, this will generate 3D stack for generated
%just raw background with no 'neurons'

function stack = generate_neural_volNoOverlap(System)
    tarrads = ceil(System.targetsizes/(2*System.dx)); %Radius of objs on grid in pixels
    bgrad = ceil(System.backgroundsize/(2*System.dx));
    stack = zeros(System.N,System.N,length(System.depth)); %Preallocate
    [xgrid,ygrid] = meshgrid([-System.N/2:System.N/2-1]); %Unitless pixel grids
    maxSize = max(tarrads); %Largest Neuron we want to generate. We will use this to discretize grid further
    div = floor(System.N/(2*maxSize)); % Number of divisions in 1D
    if System.M*length(System.depth) > div^2 %We are asking for more particles than this supports
        error('Number of particles not supported. Increase pixel # or decrease max particle size.');
    else
        total_particles = randperm(div^2,System.M*length(System.depth)); %Randomly select desired number
    end
   
    for i = 1:length(System.depth)
        
        %Generate non overlapping target particles
        temp_tar = zeros(System.N,System.N);
        if System.M ~= 0
            for j = 1:System.M
                %Select square on divided grid
                num = randperm(length(total_particles),1); %Choose entry
                [xsec,ysec] = ind2sub([div,div],total_particles(num)); %Find coordinates on new grid
                total_particles(total_particles  == total_particles(num)) = []; %Remove from array
                %Choose object size
                x = rand; %Random selector
                select = cumsum(System.ratios); %Determine 'cutoff' points for random assignment between zero and one
                ind = find(select >= x); %Find all values greater than rand variable
                indsel = ind(1); %First largest is closest value
                randlim = 1*(maxSize - tarrads(indsel)); %How much random movement is allows for Neuron
                temp_tar(sqrt((xgrid-(xsec-div/2-1/2)*2*maxSize+randi([-randlim randlim])).^2 + (ygrid-(ysec-div/2-1/2)*2*maxSize+randi([-randlim randlim])).^2) < tarrads(indsel)) = System.targetsig; %Generate points
            end
        end
        %Generate background particles
        if System.M == 0
            temp_bg = zeros(System.N,System.N);
            centroidsXbg = randperm(System.N-2*bgrad,floor(1*System.ratio)); %Rand select non-overlapping centroids
            centroidsYbg = randperm(System.N-2*bgrad,floor(1*System.ratio));
            centroidsXbg = centroidsXbg - round(System.N/2)+bgrad; %Normalize to grid between -N/2:N/2-1
            centroidsYbg = centroidsYbg - round(System.N/2)+bgrad;
            for j = 1:floor(System.ratio)
              temp_bg(sqrt((xgrid-centroidsXbg(j)).^2 + (ygrid-centroidsYbg(j)).^2) < bgrad) = System.backgroundsig; %Generate points
            end

        else
            temp_bg = zeros(System.N,System.N);
            centroidsXbg = randperm(System.N-2*bgrad,floor(System.M*System.ratio)); %Rand select non-overlapping centroids
            centroidsYbg = randperm(System.N-2*bgrad,floor(System.M*System.ratio));
            centroidsXbg = centroidsXbg - round(System.N/2)+bgrad; %Normalize to grid between -N/2:N/2-1
            centroidsYbg = centroidsYbg - round(System.N/2)+bgrad;
            for j = 1:floor(System.M*System.ratio)
              temp_bg(sqrt((xgrid-centroidsXbg(j)).^2 + (ygrid-centroidsYbg(j)).^2) < bgrad) = System.backgroundsig; %Generate points
            end
        end
       stack(:,:,i) = max(temp_bg,temp_tar);
    end
    
    %Add Noise
    if strcmp(System.noise,'Gaus')
        stack = imnoise(stack,'gaussian',System.var1,System.var2);
    elseif strcmp(System.noise,'Pois')
        stack = imnoise(stack,'poisson');
    elseif strcmp(System.noise,'SP')
        stack = imnoise(stack,'salt & pepper',System.var1);
    end
end