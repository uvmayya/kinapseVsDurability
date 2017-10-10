function [datacell,stabInd] = includeLconfZoneStabInd(inputDatacell)
%INCLUDECONF Summary of this function goes here
%   Detailed explanation goes here
% built from the function confinement

Sm=10; % maximum segment length (could try 10 or 12, but 10 or 15 didn't make a huge difference)
diffCoeff=2; % a value chosen to give the best chance of picking real positionally confined regions. It also avoids getting abusrdly high L values
% 6.64 mean value from all tracks from all datasets pooled together (without ccl21 & lo aCD3) for that cell-type; separate for bilayers
% 8.09 for hCD8 uniform coated; 5.48 for hCD4 uniform coated; 7.09 for mCD8 uniform coated; 1.87 for hCD8 bilayer
% mean value of D is used rather than individual D value for the track.
% Individual D value will shrink or expand the permissible R for the same probability psi. 
% However this will make it difficult to have uniformity in deciding on positional confinement.

datacell=remove_shortpaths(inputDatacell,40); % anything below 40 time-steps will be removed
stabInd=zeros(length(datacell),1); % column vector

%cp=6;
for cp=1:length(datacell)
    t_len = length(datacell{cp}(:,3));
    x_p  = datacell{cp}(:,3); % x coordinates are the third column
    y_p  = datacell{cp}(:,4); % y coordinates are the fourth column
    
    allL = cell(1,t_len); % to store L from each segment for each track-point
    for t=1:t_len, allL{t}=[]; end % initiation of each cell-entry
    avL = zeros(t_len,1); % column vector is best to include into datacell
    
    for i = 4:Sm % vary the segment length (Sm needs to be optimized); Ricky had it start at 4, but paper says 3-steps & (k+j) makes 3 steps
        for j = 1:t_len-i % for every starting point; t_len-i ensures that k+j (segment) doesn't go beyond t_len in the loop below
            R = 0;
            for k = 1:i % as per the current segment length; Ricky had it start at 2
                d = (x_p(k+j)-x_p(j)).^2+(y_p(k+j)-y_p(j)).^2; % d is distance-squared here
                if (d > R)
                    R = d;
                end
            end
            lPsi = 0.2048-(2.5117*diffCoeff*i/R); % as per definition it should be R^2; however here R is same as distance-squared
            if (exp(lPsi) <= .1)
                L = -lPsi - 1;
            elseif (exp(lPsi) > .1)
                L = 0;
            end
            for k = 0:i, allL{k+j}(end+1)=L; end % every point on the segment gets the value L in this iteration. 
        end
    end
    for t = 1:t_len, avL(t) = mean(allL{t}); end
    datacell{cp}(:,15)=avL; % save average L %42
    
    confZone = zeros(t_len,1); % column vector is best to include into datacell
    t=1;
    % scan through the track 
    while(t<=t_len)
        % identify putative confZone
        if avL(t)>3 % equivalent to probability of 0.007 for L of 4 and 0.02 for L of 3; 3 or 4 didn't make a difference
            startPhase=t;
            while(avL(t)>3 && t<t_len)% to go till the end of the putative confZone
                t=t+1; 
            end
            if t<t_len, endPhase=t-1;
            else endPhase=t; end
            
            %to confirm and home in on the true positionally confZone 
            meanTop10dist=35; % initialize to high distance so that the while loop below runs atleast once
            RsqByT=(meanTop10dist^2)/(endPhase-startPhase+1); % this needs to be a constant or low to keep the same log-psi
            while(RsqByT>0.66 && endPhase-startPhase>10) % at least for 5 minutes (10 steps) & happens to have larger confinement zone on time-relative terms
                vectorLength=(endPhase-startPhase+1)*(endPhase-startPhase)/2; % based on combination formula of factorials
                dist=zeros(vectorLength,1); % initialize 
                distInd=1; 
                for j=endPhase:-1:startPhase+1
                    for i=j-1:-1:startPhase
                        dx = x_p(i) - x_p(j);
                        dy = y_p(i) - y_p(j);
                        dist(distInd) = norm([dx,dy]);% Eucledian distance
                        distInd=distInd+1;                        
                    end
                end
                dist=sort(dist, 'descend');
                meanTop10dist=mean(dist(1:10)); % gives diameter of the confined zone
                RsqByT=(meanTop10dist^2)/(endPhase-startPhase+1); % this needs to be a constant or low to keep the same log-psi
                if RsqByT<=0.66 % turned out to be very cir
                    confZone(startPhase:endPhase)=1;
                    break;  % you have identified the zone with positional stability
                else % clip-off points from either or both ends to see if positional confinement improves in the next iteration
                    if avL(startPhase)>1.5*avL(endPhase)
                        endPhase=endPhase-1; 
                        % in some instances, at the ends of tracks or as the clipping process is ongoing, one end may have higher L than the other. 
                        % in such instances clip only the lower end 
                    elseif avL(startPhase)<1.5*avL(endPhase)
                        startPhase=startPhase+1;
                    else % clip both ends   
                        startPhase=startPhase+1;
                        endPhase=endPhase-1; 
                    end    
                end    
            end
        end
        t=t+1;   
    end
    datacell{cp}(:,16)=confZone; % save average L % 43
    datacell{cp}(:,17)=nnz(confZone)/length(confZone); % positional stability index % 44
    stabInd(cp)=nnz(confZone)/length(confZone);
end        




end

