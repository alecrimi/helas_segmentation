%Load file
[FileNames,PathName] = uigetfile('*.csv','Select the stack you want to process','MultiSelect','on');

% Check if the selected files are more than 1
if (iscell(FileNames))
    n_files = length(FileNames);
else
    n_files = 1 ;
end

%barr_data = zeros(n_files,10);
data =[];

for ff = 1 : n_files
    % Select initial point for ROI    
    if (n_files > 1)
        data_exp = csvread([PathName cell2mat(FileNames(ff))]);
    else
        data_exp = csvread(FileNames);
    end
    data = [ data ; data_exp];
end

figure;

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
    saveas(gcf,strcat(PathName,'analyzed/tot_ci'.png'));
