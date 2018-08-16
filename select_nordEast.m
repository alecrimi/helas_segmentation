% Select region of Interest
function [ x ,y ] = select_NordEast(fname)

disp('#######################################')
first_pic = imread(fname, 1);
disp('Select the NordEast corner of the Region of Interest')
 figure('Name',fname,'NumberTitle','off');; imagesc(first_pic); colormap(hsv); grid on;
imcontrast;
%rect = round(getrect);
[x, y] = ginput(1);

end
