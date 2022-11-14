function [initialPop,popCost] = dsGA_First_Step(maxIter, numSections, NAmax, offset, ObjectiveFunction, LBphase, UBphase, options)
   %Generate lower and upper bound
    nvars = 3*numSections+numSections-1; %Define population by number of sections
    dNA = (NAmax)/(numSections-1); %One less cutoff point than total sections
    LB = [];
    UB = [];

    for i = 1:numSections
        if i < numSections
            LB = [LB LBphase 0+dNA*(i-1)+offset];
            UB = [UB UBphase 0+dNA*(i)];
        else
            LB = [LB LBphase];
            UB = [UB UBphase];
        end
    end 

    %Setup Genetic Algorithm and run!
    initialPop = zeros(maxIter,length(LB));
    initialVal = zeros(1,maxIter);
    for i = 1:maxIter
        rng('shuffle') %Shuffle seed to ensure random start point
        [ipop,ival] = ga(ObjectiveFunction,nvars, [], [], [], [], LB, UB, [], options);
        initialPop(i,:) = ipop; %Save pop
        initialVal(i) = ival; %Save cost
    end
    popCost = initialVal;
end