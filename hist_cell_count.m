%function hist_cell_count
close all
%Load file
[FileNames,PathName] = uigetfile('*.csv','Select the stack you want to process','MultiSelect','on');

% Check if the selected files are more than 1
if (iscell(FileNames))
    n_files = length(FileNames);
else
    n_files = 1 ;
end
n_bins = 10;
barr_data = zeros(n_files,n_bins);
bins = zeros(n_files,n_bins);
for ff = 1 : n_files
    %%%%%% Data for cell count
    % Select initial point for ROI    
    if (n_files > 1)
        data_exp = csvread([PathName cell2mat(FileNames(ff))]);
    else
        data_exp = csvread(FileNames);
    end
    [count, centers] = hist(data_exp(:,1));
    barr_data(ff,:) = count;
    %%%%%%% Data for intensity count
    sorted_intensities = sortrows(data_exp,1);
    interval = 1 : max(data_exp(:,1))/10 : max(data_exp(:,1));
    interval(end + 1) = max(data_exp(:,1));
    for hh = 1 : n_bins     
            pos = find(sorted_intensities(:,1)>interval(hh) & sorted_intensities(:,1)<interval(hh+1) );
            bins(ff,hh) = sum(sorted_intensities(pos,2));
    end

end

figure;
A = [100: 165: 1650];
boxplot(barr_data,A);
hold on
[n1,~]=size(barr_data);
plot(ones(n1,1)*[1:length(barr_data)],barr_data,'r*')
hold off
xlabel('Interval Depth in \mum') % x-axis label
ylabel('Number of cells') % y-axis label
saveas(gcf,strcat(PathName,'analyzed/cc_',fname,'.png'));


figure; boxplot(bins, A);
hold on
plot(ones(n1,1)*[1:n_bins],bins,'r*')
hold off
xlabel('Interval Depth in \mum') % x-axis label
ylabel('Cumulative sum of intensities from cells') % y-axis label
saveas(gcf,strcat(PathName,'analyzed/ci_',fname,'.png'));
