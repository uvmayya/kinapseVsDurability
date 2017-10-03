function datacell2=selectArrestedPortionsOnSpots(datacell,flurThresh,minLength,proxThresh)
% selects portions of tracks that were attached on spots and didn't have another attached cell on the same spot. 
% adopted from the function file filterOutArrestedPortions in 041414/anals/072415_motilityAnals

% removes arrested/attached portions of tracks based on either one or a set of three conditions.
% three conditions (that are combinations of 2 conditions) are used mainly because some arrested cells inexplicably don't have IRM

%---------------------------------
datacell1={}; % stores arrested portions of tracks that are on the spots
for cp=1:length(datacell)
    cond=zeros(size(datacell{cp}(:,1))); % col vector
    for pos=1:length(datacell{cp}(:,1))
        cond1=(datacell{cp}(pos,5)>5 && datacell{cp}(pos,6)>flurThresh); % irm and flur 
        cond2=(datacell{cp}(pos,5)>5 && datacell{cp}(pos,10)<1.5); % irm and speed4
        cond3=(datacell{cp}(pos,6)>flurThresh && datacell{cp}(pos,10)<1.5); % flur and speed4
        cond(pos)=(cond1+cond2+cond3)>=1; % cond is 1 if arrested on the spot       
    end
    startPos=[]; endPos=[]; 
    if cond(1)==0, startPos(end+1)=2; endPos(end+1)=2;
    elseif cond(1)==1, startPos(end+1)=1; endPos(end+1)=1; end
    for pos=2:length(cond)
        if cond(pos)==1 && cond(pos-1)==1
            endPos(end)=pos; % extend the existing arrested sub-track length
        elseif cond(pos)==1 && cond(pos-1)==0
            startPos(end+1)=pos; % start a new arrested sub-track
            endPos(end+1)=pos;
        end
        % don't need to bother about 0 followed by 0 in cond
    end
    
    for n=1:length(startPos)
        if endPos(n)>startPos(n)+minLength % to ensure that arrested portion is more than minimum frames in length
           datacell1{end+1}=datacell{cp}(startPos(n):endPos(n),:);
           datacell1{end}(:,1)= datacell{cp}(1,1)+startPos(n)-1; % update start & end frames
           datacell1{end}(:,2)= datacell{cp}(1,1)+endPos(n)-1;
        end    
    end
end

%-----------------------------------------

% build the proximity condition cell-array first
proxCond=cell(1,length(datacell1)); % cell-array of empty matrices
for cp=1:length(datacell1)
    cond=ones(size(datacell1{cp}(:,1))); % col vector
    proxCond{cp}=cond; % stores 
end    
for cp1=1:length(datacell1) % this portion is copied from findNeighbor2 function under tiam_attachmentRate package
    for pos1=1:length(datacell1{cp1}(:,1));
        loc1=datacell1{cp1}(pos1,3:4); % x,y location of the cell1 in the start frame
        fr=datacell1{cp1}(1,1)+pos1-1;  
        for cp2=cp1+1:length(datacell1) % all the subsequent track entries
            if fr>=datacell1{cp2}(1,1) && fr<=datacell1{cp2}(1,2) && cp1~=cp2 % frame of cell1 within the track of cell2
                pos2=fr-datacell1{cp2}(1,1)+1;
                loc2=datacell1{cp2}(pos2,3:4); % x,y location of cell2 in the relevant frame
                dist=pdist([loc1;loc2],'euclidean');
                if dist<=proxThresh && dist~=0  % 0 dist will mean the same cell 
                    proxCond{cp1}(pos1)=0; % there is a close neighboring attached cell
                    proxCond{cp2}(pos2)=0;
                end
            end    
        end
    end
end 

datacell2={}; % stores arrested portions of tracks that are on spots and don't have another arrested neighbor
for cp=1:length(datacell1)
    cond=proxCond{cp};
    startPos=[]; endPos=[]; 
    if cond(1)==0, startPos(end+1)=2; endPos(end+1)=2;
    elseif cond(1)==1, startPos(end+1)=1; endPos(end+1)=1; end
    for pos=2:length(cond)
        if cond(pos)==1 && cond(pos-1)==1
            endPos(end)=pos; % extend the existing arrested sub-track length
        elseif cond(pos)==1 && cond(pos-1)==0
            startPos(end+1)=pos; % start a new arrested sub-track
            endPos(end+1)=pos;
        end
        % don't need to bother about 0 followed by 0 in cond
    end
    
    for n=1:length(startPos)
        if endPos(n)>startPos(n)+minLength % to ensure that arrested portion is more than minimum frames in length
           datacell2{end+1}=datacell1{cp}(startPos(n):endPos(n),:);
           datacell2{end}(:,1)= datacell1{cp}(1,1)+startPos(n)-1; % update start & end frames
           datacell2{end}(:,2)= datacell1{cp}(1,1)+endPos(n)-1;
        end    
    end
end


end


