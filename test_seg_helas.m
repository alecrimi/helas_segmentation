% Test script for segmenting HeLas
clear all;
close all; 

%%%%%%%%% Setup variables %%%
min_area = 40; %Min area per cell
max_area = 100; %Max area per cell
legnth = 1000; %Length ROI
width = 700; %Length ROI
resolution = 4 ; %Resolution (um/pixel)

%Load file
fname = 'USExpII_HeLa_580kHz_40%e_30min_LTT-l_8%AA_CK8-18-a649_647ex_Quadruple_lp10_rightL_1-6x_15umStep(417).tif';
disp('#######################################')
disp('Loading images');
info = imfinfo(fname); 
num_images =   150;%numel(info);
stack = zeros(info(1).Width, info(1).Height, num_images);
for k = 1 :  num_images
    stack(:,:,k) = imread(fname, k);
end
disp('Image loaded')

% Select region of Interest
%disp('Select Region of Interest')
disp('Select the NordEast corner of the Region of Interest')
figure; imagesc(stack(:,:,1)); colormap(hsv); grid on;
imcontrast;
%rect = round(getrect);
[x, y] = ginput(1);
rect(1) = round(x);
rect(2) = round(y);
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
print(strcat('plot_',fname),'-dpng');
saveas(gcf,strcat('../analyzed/plot_',fname,'.png'));

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
saveas(gcf,strcat('../analyzed/segmented_',fname,'.png'));

% Save points
data = [ c_list_flipped' mean_list'];
% First column center, second column intensity
csvwrite(['../analyzed/datapoints_' fname '.csv'],data)

% Histogram
figure('Name',fname,'NumberTitle','off');
hist(c_list_flipped);
xlabel('Depth in \mum') % x-axis label
ylabel('Number of cells') % y-axis label
box on
axis([ 0 1600*resolution  0 1500])
saveas(gcf,strcat('../analyzed/histogram_',fname,'.png'));


%%%%%%%%%%%%%%%%%%%%%OLD STUFF, KEEP IT TEMPORARILY%%%%%%%%%%%%%%%%%%%%%%
%{
% Compute intensity within cells (connected elements
%overall_mask = sum(segmented_stack,3);
%overall_imm = sum(Im_test_cropped,3);
%L = bwlabel(overall_mask,4);
%n_detected_cells = max(max(L));
overlay_stack = sum( Im_test_cropped .* segmented_stack,3 );
average_int = mean(overlay_stack); 

figure('Name',fname,'NumberTitle','off');
hold on
     for kk = 1  : length(average_int)
         if ( average_int(kk) ~= 0  )   
             plot(kk,average_int(kk),'b*','MarkerSize',3);
         end
     end
hold off
xlabel('Position in tube from bottom to top') % x-axis label
ylabel('Mean intensity value inside segmented cells') % y-axis label
box on
%}

%{
figure('Name',fname,'NumberTitle','off'); imagesc(original_overlay );
print(strcat('maxProj_',fname),'-dpng');

% Make a truecolor all-green image.
green = cat(3, zeros(size(original_overlay)), ones(size(original_overlay)), zeros(size(original_overlay)));
hold on
h = imagesc(green);
hold off
set(h, 'AlphaData', global_overlay)
% figure('Name',fname,'NumberTitle','off'); imagesc(original_overlay );
print(strcat('maxProjoverlay_',fname),'-dpng');
%}
