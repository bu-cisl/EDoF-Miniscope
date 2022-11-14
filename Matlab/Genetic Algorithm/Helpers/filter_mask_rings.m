%The goal of this code is to remove any rings below a certain size. Also
%Binarization sometimes form inperfect rings (they vary off their radial
%position) so we take care of that too.
function filtmask = filter_mask_rings(mask,erodeRate)
    %First develop a pixel grid for later use
    N = size(mask,1);
    [px,py] = meshgrid([-floor(N/2):(floor(N/2))-1]);
    pr = floor(sqrt(px.^2+py.^2)); %Radial pixel grid
    
    %Get a cross section of the mask and determine edge
    %locations/thicknesses
    cross = mask(N/2+1,:);
    edges = abs(cross - circshift(cross,1))/2; %Shift the mask by one pixel and subtract. This artificially increases the measured ring size by 1. I tried matlab edge(cross) command but it didn't pick up on 1pix edges
    edges = edges(1:N/2+1); %Since radial symmetry, only need have the edges
    edgeloc = abs(find(edges(1:end-1) > 0)-N/2-1); %Find location of pixels and normalize to grid. The last entry in edges is an artificial edge outside the aperture created by the circshift
    edgeloc = [edgeloc 0];
    edgesize = (edgeloc(1:end-1) - edgeloc(2:end)); %Find size of each fringe
    
    %Method 1: Now id rings that are too small and remove them
    edgelist = find(edgesize <= erodeRate); 
    for i = 1:length(edgelist)-1
        mask(find(pr >= edgeloc(edgelist(i)) & pr <= edgeloc(edgelist(i))+edgesize(i)+1)) = 1;
    end
    
%     %Method 2: Instead Rebuild Mask ring by ring - inflate ring size to
%     %fill need
%     edgelist = find(edgesize > erodeRate);
%     blank = zeros(N);
%     for i = 1:length(edgelist)
%         blank(find(pr <= edgeloc(edgelist(i)) & pr >= edgeloc(edgelist(i))-edgesize(edgelist(i))-erodeRate)) = cross(-edgeloc(i)+N/2+1);
%     end
    
    %filtmask = blank;
    filtmask = mask;
    
    