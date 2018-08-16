% Test script for segmenting HeLas with multple selection of files
clear all;
close all; 
%%%%%%%%% Setup variables %%%
min_area = 40; %Min area per cell
max_area = 100; %Max area per cell
legnth = 1000; %Length ROI
width = 700; %Length ROI
resolution = 4 ; %Resolution (um/pixel)

%Load file
[FileNames,PathName] = uigetfile('*.tif','Select the stack you want to process','MultiSelect','on');

% Check if the selected files are more than 1
if (iscell(FileNames))
    n_files = length(FileNames);
else
    n_files = 1 ;
end

init_points = zeros(n_files,2);
%sum_int = zeros(n_files,1);

for ff = 1 : n_files
    disp('#######################################')
    disp('Loading images');
    if (n_files > 1)
        fname = cell2mat(FileNames(ff));
    else
        fname = FileNames;
    end
    
    info = imfinfo(fname); 
    num_images =  numel(info);
    stack = zeros(info(1).Width, info(1).Height, num_images);
    for k = 1 :  num_images
        stack(:,:,k) = imread(fname, k);
    end
    disp('Image loaded')
    original_overlay = max(stack, [], 3);
  
    disp('#######################################')
    disp('Select the NordEast corner of the Region of Interest')
    figure('Name',fname,'NumberTitle','off'); imagesc(original_overlay); colormap(hsv); grid on;
    imcontrast;
    %rect = round(getrect);
    [x, y] = ginput(1);
    init_points(ff,1) = x;
    init_points(ff,2) = y;
end

for kk = 1 : n_files
    disp('#######################################')
    disp('Loading images');
    if (n_files > 1)
        fname = cell2mat(FileNames(kk));
    else
        fname = FileNames;
    end
    
    info = imfinfo(fname); 
    num_images =  numel(info);
    stack = zeros(info(1).Width, info(1).Height, num_images);
    for k = 1 :  num_images
        stack(:,:,k) = imread(fname, k);
    end
    disp('Image loaded')

    rect(1) = round( init_points(kk,1));
    rect(2) = round( init_points(kk,2));
    rect(3) = legnth; 
    rect(4) = width;
    Im_test_cropped = stack(rect(2):rect(2)+rect(4),rect(1)-rect(3):rect(1),:);
    segmented_stack = zeros(size(Im_test_cropped));
    fprintf('The selected region is the square starting from the point (NordEast) %d,%d and long %d and %d \n', rect(2), rect(1),  rect(4), rect(3) );

    center_list = []; %Incremental, this has to be improved
    mean_list = [];

    % Stack of images
    for nn = 1 : num_images 
        segmented_stack(:,:,nn) = segment_HeLas(Im_test_cropped(:,:,nn),min_area,max_area);
    %end

    % Compute intensity within cells (connected elements
    %overall_mask = sum(segmented_stack,3);
    %overall_imm = sum(Im_test_cropped,3);
    %L = bwlabel(overall_mask,4);
    %n_detected_cells = max(max(L));

        L = bwlabel(segmented_stack(:,:,nn),4);
        n_detected_cells = max(max(L));

        for jj = 1 : n_detected_cells
            [r c v] = find(L == jj);
            mean_list(end+1) = mean(mean(Im_test_cropped(r,c,nn))); %mean(mean(overall_imm(r,c)));%Im_test_cropped(round(mean(r)),round(mean(c)),nn);%
            center_list(end+1) = mean(c);
        end
    end 
    % Flip x-axis
    c_list = (rect(2)+center_list)*resolution;
    c_list_flipped = max(c_list) - c_list;

    figure('Name',fname,'NumberTitle','off');
    plot( c_list_flipped, mean_list,'b*','MarkerSize',3);
    xlabel('Depth in \mum') % x-axis label
    ylabel('Mean intensity value inside segmented cells') % y-axis label
    box on
    axis([ 0 1600*resolution  0 1500])

    %overlay_stack_seg  = max( segmented_stack,3);
    overlay_stack_seg = max(segmented_stack, [], 3);
    overlay_stack_mask = overlay_stack_seg > 1 ;
    %print(strcat('plot_',fname),'-dpng');
    saveas(gcf,strcat(PathName,'analyzed/p_',fname,'.png'));

    % Mask
    %figure; imagesc(overlay_stack_mask);
    %original_overlay = max(stack,3); %sum(Im_test_cropped,3);
    original_overlay = max(stack, [], 3);
    global_overlay = original_overlay;%zeros(size(original_overlay));
    global_overlay(rect(2):rect(2)+rect(4),rect(1)-rect(3):rect(1),:) = global_overlay(rect(2):rect(2)+rect(4),rect(1)-rect(3):rect(1),:) + 10000*overlay_stack_mask;

    RGB = cat(3, original_overlay, original_overlay, original_overlay);
    RGB(:,:,1) = global_overlay; 
    figure('Name',fname,'NumberTitle','off'); imshow(uint8(RGB));
    hold on; quiver(1800,2000,250,0,'ShowArrowHead','off','color',[1 1 1],'linewidth',10) %Scalebar
    saveas(gcf,strcat(PathName,'analyzed/s_',fname,'.png'));

    % Save points
    data = [ c_list_flipped' mean_list'];
    % First column center, second column intensity
    csvwrite([PathName 'analyzed/d_' fname '.csv'],data)

    % Histogram
    figure('Name',fname,'NumberTitle','off');
    hist(c_list_flipped);
    xlabel('Depth in \mum') % x-axis label
    ylabel('Number of cells') % y-axis label
    box on
    axis([ 0 1600*resolution  0 1500])
    saveas(gcf,strcat(PathName,'analyzed/h_',fname,'.png'));
 
    %Cumulative intensities
    sorted_intensities = sortrows(data,1);
    interval = 1 : max(data(:,1))/10 : max(data(:,1));
    n_bins = length(interval);
    bins = zeros(n_bins,1);
    interval(end + 1) = max(data(:,1));
    for hh = 1 : n_bins     
        pos = find(sorted_intensities(:,1)>interval(hh) & sorted_intensities(:,1)<interval(hh+1) );
        bins(hh) = sum(sorted_intensities(pos,2));
    end
    figure; bar(interval(1:end-1)/2,bins);
    xlim([-200 max(data(:,1))])
    xlabel('Depth in \mum') % x-axis label
    ylabel('Cumulative sum of intensities from cells') % y-axis label
    saveas(gcf,strcat(PathName,'analyzed/ci_',fname,'.png'));
end
