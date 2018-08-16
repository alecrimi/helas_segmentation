function L = segment_HeLas (Input, min_area, max_area, size_se1, size_se2, threshold )

% Marker-Controlled Watershed Segmentation
% The markers are given by the regionalmax, specify the minimum size of the
% area

if nargin < 5
   size_se1 = 1; 
   size_se2 = 3; 
end
% Assure the image is unsigned int
Image = uint8(Input);

if nargin < 6
   threshold = graythresh(Image); 
end

% Opening-closing by reconstruction 
se1 = strel('disk', size_se2);
Ie = imerode(Image, se1);
Iobr = imreconstruct(Ie, Image);
se2 = strel('disk', size_se1);
Iobrd = imdilate(Iobr, se2);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);

% Compute Regional maxima of opening-closing by reconstruction
fgm = imregionalmax(Iobrcbr);
se3 = strel(ones(size_se2, size_se2));
fgm2 = imclose(fgm, se3);
fgm3 = imerode(fgm2, se3);
fgm4 = bwareaopen(fgm3, size_se2 * 4);

% Dilate again to fill holes in the segmentation
se4 = strel('disk',size_se1);
fgm4 = imdilate(fgm4,se4);

%{
% Compute the Watershed Transform of the Segmentation Function.
bw = imbinarize(Iobrcbr, threshold);
imagesc(fgm4);
D = bwdist(fgm4);
DL = watershed(D);
bgm = DL == 0;
figure; imagesc(bgm);
gradmag2 = imimposemin(gradmag, bgm | fgm4);
%}

% Keep only region within the area range
L = bwlabel(fgm4,4);
n_connected_areas = max(max(L));
%figure; imagesc(L);
for kk = 1 : n_connected_areas
    sum_count = sum(sum( L == kk )); 
    
    if ( (sum_count < min_area) || (sum_count > max_area)  )
        [r c v] = find(L == kk);
        L(r,c) = 0;
    end
end
%figure; imagesc(L);
