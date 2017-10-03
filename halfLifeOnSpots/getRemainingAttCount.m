function remAttCount=getRemainingAttCount(datacell, startFrame)
% to measure number of remaining attached cells from the frame specified in startFrame
% this will not count any new cells that attach after the startFrame. 

datacell=remove_shortpaths(datacell,10);
noOfTracks = size(datacell, 2);  % could have used length(datacell)
% cell-array of track data is called 'datacell' by Willie
% the no. of columns in the cell-array corresponds to noOfTracks

endFrame=0;
for cp=1:noOfTracks
    if datacell{cp}(1,2)>endFrame
        endFrame=datacell{cp}(1,2);
    end
end

remAttCount=zeros(endFrame-startFrame+1,1);
for fr=startFrame:endFrame
    for cp=1:noOfTracks
        if datacell{cp}(1,1)<=startFrame
            if datacell{cp}(1,2)>=fr
                remAttCount(fr-startFrame+1)=remAttCount(fr-startFrame+1)+1;
            end    
        end
    end
end

end