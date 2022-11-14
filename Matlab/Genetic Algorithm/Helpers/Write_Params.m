%The goal of this function is to write a universal struct that holds the
%desired parameters for all GA simulations. Those parameters are:

%{
    lambda = center wavelength 
    n2 = submersion media
    N = number of pixels in each dimension
    targ_size = target particle size
    max_depth = target depth range for optimization
    depth_step = stepsize
    NA0 = max NA
    f = focal length of system
    fnum = system F_number
    mag = system magnification
    pix = pixel size
%}

function Params = Write_Params(lambda,n2,N,targ_size,targ_depth,show_depth,maxDepth,depth_step,NA0,f,fnum,mag,pix,W040,erodeRate,innerCircle,saveDir)
    %Generate Struct
    Params.lambda = lambda;
    Params.n2 = n2;
    Params.N = N;
    Params.targ_size = targ_size;
    Params.show_depth = show_depth;
    Params.targ_depth = targ_depth;
    Params.maxDepth = maxDepth;
    Params.depth_step = depth_step;
    Params.NA0 = NA0;
    Params.f = f;
    Params.fnum = fnum;
    Params.mag = mag;
    Params.pix = pix;
    Params.W040 = W040;
    Params.erodeRate = erodeRate;
    Params.innerCircle = innerCircle;
    Params.saveDir = saveDir;
end