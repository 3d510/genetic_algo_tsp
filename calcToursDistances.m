function totalDist = calcToursDistances(pop, popSize, dmat, n)
% this function calculates the tour length of all members of the population
% pop. 
       for p = 1:popSize
            d = dmat(pop(p,n),pop(p,1)); % Closed Path
            for k = 2:n
                d = d + dmat(pop(p,k-1),pop(p,k));
            end
            totalDist(p) = d;
        end
end

