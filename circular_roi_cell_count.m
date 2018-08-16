function res_count = count_cells(im)

rows = 8;
columns = 12;

res_count = zeros(rows,columns);

imag = imread(im);
imshow(imag);

for kk = 0 : columns -1 
 for jj = 0 : rows -1

    % Select ROI0
    %h =  imellipse;
    shiftx = 1130 * kk;
    shifty = 1150 * jj;
    h = imellipse(gca,[(780 + shiftx) (680+shifty) 650 650]);
    % Mask data
    BW = createMask(h);
    %imag = rgb2gray(imag);
    sel = imag .* uint16(BW);

    %{
    [r,c] = find(sel);
    start_y = min(r);
    end_y = max(r);
    start_x = min(c);
    end_x = max(c);

    ROI = imag(start_y:end_y,start_x:end_x);
    ROI = medfilt2(ROI);
    %}
    ROI = medfilt2(sel);
    ROI = imadjust(ROI, [0.5 1],[] );
    %ROI = imadjust(ROI );
    level = graythresh(ROI) ;

    BW = imbinarize(ROI,0.7); %0.30
    J = imcomplement(BW);
    [L,counts] = bwlabel(J);
    %figure; imagesc(L);

    res_count(jj+1,kk+1) = counts - 1; %subtract 1 because of the background 
 end
end
