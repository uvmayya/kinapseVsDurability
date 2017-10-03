function [samplingEff,meanOverlap] = getSampEffOverlapForTrack(track, videocell, spotImg, stFrame, halfCropSize, detectParams, umPerPixel)
% adopted from the function with the same name in 041414/anals/072415_motilityAnals
% modified from function 'getSamplingMapFromTracks' in the tiam_attachmentRate function


overlap=zeros(size(track(:,1:3))); % initialize as an array with 3 columns
meanOverlap=[0,0,0];
%for cp = 1 : size(datacell, 2)
    % initialize the map for every track
    map=zeros(size(videocell{1}));
    totalArea=0; % this stores the total area the cell could potentially sample
    trackLength=size(track, 1);
    cellArea=zeros(trackLength,1);
    
	for pos = 1 : trackLength % consider every frame
    %for pos = 1 : 5: trackLength  % consider every 5th frame  
        frame = stFrame + pos -1;        
        center_x = round(track(pos, 3)/umPerPixel); % rounded to nearest integer
        center_y = round(track(pos, 4)/umPerPixel);
        
        halfcropXsize = halfCropSize; %size influences how watershed works (because the size determines the gradient)
        halfcropYsize = halfCropSize;
        if center_x<=halfcropXsize
            halfcropXsize=center_x-1;
        elseif center_x+halfcropXsize>=size(videocell{frame},2)
            halfcropXsize=size(videocell{frame},2)-center_x; % column value corresponds to x-coordinate of the pixel
        end % to account for centroids that are too close the edges of the frames
        if center_y<=halfcropYsize
            halfcropYsize=center_y-1;
        elseif center_y+halfcropYsize>=size(videocell{frame},1)
            halfcropYsize=size(videocell{frame},1)-center_y; % row value corresponds to y-coordinate of the pixel
        end % to account for centroids that are too close the edges of the frames

        if(length(size(videocell{1,frame})) == 3)
            I=rgb2gray(videocell{frame});
        else I=videocell{frame};
        end 

        I = imcrop(I, [uint16(center_x-halfcropXsize) uint16(center_y-halfcropYsize) halfcropXsize*2 halfcropYsize*2]);
        % the first pixel in the cropped image corresponds to center_y-halfcropsize and center_x-halfcropsize in the original image
        % actual size of the cropped image is twice the halfcropsize plus 1
        % thus centroid pixel(halfcropYsize+1,halfcropXsize+1) is the absolute center pixel and there are equal number of pixels on all sides
        % unit16 converts to unsigned 16-bit integer type
        
        [Ie,thresh] = edge(I,'canny'); %[temp,thresh] = edge(Ifilt,'canny');
        %Ie = edge(I, 'canny', 0.7*thresh, 0.8); %0.7*thresh for exp1
        Ie = edge(I, 'canny', 1.0*thresh, 0.8); % 1.3*thresh for exp5
        %Ie = edge(I, 'canny', 1.0*thresh, 2.2); % this didn't distinguish nve and memory sampling
        %Ie=edge(I, 'canny', detectParams(2), 0.8); % using the input edge value didn't work
        Icht = im2uint8(Ie); % needed for CircularHough_Grd function
            
        % structural elements
        se90 = strel('line', 2, 90);
        se0 = strel('line', 2, 0);
        disk3 = strel('disk', 4); % 4/3 for exp1 and 3 for exp5; 4/3 worked well overall and didn't cause cells to merge after dilation
        disk1 = strel('disk', 1);
        disk2 = strel('disk', 2);

        % take canny image and dilate, fill, erode, and dilate
        Id = imdilate(Ie, [se90, se0]);
        If = imfill(Id, 'holes');
        Ir = imerode(If, disk3);
        img = imdilate(Ir, disk2);
        centerBlob = imerode(img, disk3);
            
        % circular hough transform
        %[accum, circen, cirrad] = CircularHough_Grd(test, [rad_min rad_max], gradientThresh, searchRadius);  % function call format
        accum=CircularHough_Grd(Icht, [round(detectParams(3)/detectParams(1)) round(detectParams(4)/detectParams(1))], detectParams(5), detectParams(6)); % [7 15]
        %accum1=imhmax(accum,300); % 1000 appeared to worked well; 
        %ultimately I didn't need to suppress local maxima as there was no difference with and without suppression 
        imgDist=imimposemin(-accum,centerBlob); % forces the minimum to be on the centerBlob
        
        imgDist(~img)=-inf; % '~' inverts the image; sets zero values to negative infinity which ensures that the background doesn't get segmented by watershed
        imgLabel=watershed(imgDist)>0; 
        %watershed assigns 0 values to all the boundary pixels and positive integers (labels) to regions
        % but by having '>0' only non-zero values are stored in imgLabel.
        % This allows me to do the .* operation (because of compatible data-type) later on so that I can use bwlabel function. This will eliminate the need to handle background as a component.
        bwLabel=img.*imgLabel;
        imgLabel=bwlabel(bwLabel);
        
        stats=regionprops(imgLabel,'Area', 'BoundingBox', 'Centroid', 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength', 'Perimeter');
        % the component/label with the largest area appears to be the 1st component, which is also the background;
        % but at times the background is split into two components/labels

        dist=zeros(size([stats.Area]));
        for i=1:length(dist)
            dist(i)=pdist([stats(i).Centroid;halfcropXsize+1,halfcropYsize+1],'euclidean');
            % distance between the centroid of the object and center pixel of the cropped box
            % centroid of the background component also tends to be very close to center pixel of the box
        end
        %save('label','imgLabel', 'stats', 'dist', 'I', 'img');
        
        if length(dist) >= 2 % if there are two or more components, pick the foreground component that is closest to the center pixel of the cropped box
            [temp,index]=sort([stats.Area], 'descend'); %sorts in the descending order of area
            % index holds the indices of the elements that were sorted 
            for i=1:length(dist) % to pick based on proximity to the center of the cropped image
                if dist(index(i))<10 % closer to the center of the cropped image   
                    imgLabel(imgLabel~=index(i))=0; % set the rest to zero to make it a binary image
                    imgLabel(imgLabel~=0)=1; % set the object pixels to 1 to enable faithfull counting
                    break
                end
            end   
        end
        
        cellArea(pos)=nnz(imgLabel); % counts the number of pixels within the DIC cell boundary
        totalArea=totalArea+cellArea(pos);
        
        %update map counts for that section; imgLabel has 1 whereever the cell associated pixel is with the rest being zero
        map(uint16(center_y-halfcropYsize):uint16(center_y+halfcropYsize), uint16(center_x-halfcropXsize):uint16(center_x+halfcropXsize)) = map(uint16(center_y-halfcropYsize):uint16(center_y+halfcropYsize), uint16(center_x-halfcropXsize):uint16(center_x+halfcropXsize)) + imgLabel;
        
        % for overlap calculations
        imgLabel=logical(imgLabel);
        spotImgCrop=imcrop(spotImg,[uint16(center_x-halfcropXsize) uint16(center_y-halfcropYsize) halfcropXsize*2 halfcropYsize*2]);
        spotFilt = medfilt2(spotImgCrop); 
        threshLevel = graythresh(spotFilt); %doing min-max normalization might help with higher EM values, but I didn't implement it
        spotBW = im2bw(spotFilt,threshLevel);
        cellArea=nnz(imgLabel); % counts the number of pixels within the DIC cell boundary, i.e. cell area
        overlap(pos,1)=nnz(imgLabel.*spotBW)/cellArea; % .* is same as AND operation on the two masks
        overlap(pos,2)=nnz(imgLabel.*spotBW)*umPerPixel^2/78.54; % overlap area relative to 10um spot area
        overlap(pos,3)=cellArea*umPerPixel^2/78.54; % cell area given as a fraction of 10um spot area
    end
    sampledArea=nnz(map); %the non zero pixel values represent all the pixels visited by the cell
    samplingEff=sampledArea/totalArea;
    %avgArea=totalArea/trackLength;
    meanOverlap(1)=mean(overlap(:,1)); % fractional overlap
    meanOverlap(2)=mean(overlap(:,2)); % overlap relative to spot area
    meanOverlap(3)=mean(overlap(:,3)); % cell area relative to spot area
   
    % display information
	%fprintf('Sampling efficiency calculated for cell-track: %d/%d \n',cp,size(datacell,2));
    
%end    

end

