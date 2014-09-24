%% Get response tc for each trial and avg.
% Show pearson corr image for the field of view, measure on avg movie

datapath = '/Users/sasha/Documents/Wilson lab/QuasAr imaging/2014-09-12 (Flies)/'
fname = 'Sq_camera';

clear data;
close all;

% This specifies the series of runs to average
search_dirs = '*_Fly4_oAOTF400_p5s_AL';
search_dirs_savename = strrep(search_dirs, '*', '_wild_');
dirs = dir([datapath '/' search_dirs]);
nRun = length(dirs);

%dir_prefix = {'155958_Fly2_oAOTF400','155842_Fly2_oAOTF400','155929_Fly2_oAOTF400','160023_Fly2_oAOTF400','155631_Fly2_oAOTF400', '160817_Fly2_oAOTF400_0p125s','160838_Fly2_oAOTF400_0p125s','160935_Fly2_oAOTF400_0p125s','161000_Fly2_oAOTF400_0p125s','161023_Fly2_oAOTF400_0p125s' }; % 0.25 stim

%dir_prefix = {'155448_Fly2_oAOTF400', '155038_Fly2_oAOTF400', '161430_Fly2_rAOTF250_0p5s', '161453_Fly2_rAOTF250_0p5s', '161620_Fly2_rAOTF250_0p5s'}; % 0.50 stim

% dir_prefix = {'160148', '160217', '160334', '160244', '160123'}; % 1.00 stim
%search_dirs_savename = [ dir_prefix{1} '_Fly2_oAOTF400' ];
%nRun = length(dir_prefix);

STIM_PERIOD = 0.5;

Xsize = 320;
%Ysize = 240;
Ysize = 320;

dx = 2;
dy = dx;

data = zeros(Xsize/dx,Ysize/dy,4000);

t = [0:3999]./500.0;
cnt = 0;
for i=1:nRun
    path2 = [datapath '/' dirs(i).name];
    %path2 = [datapath '/' dir_prefix{i} '_Fly2_oAOTF400'];
    %path2 = [datapath '/' dir_prefix{i} ];
    
    filepath = [path2 '/' fname '.bin']

    [mov,nframe] = readBinMov( filepath, Xsize, Ysize );

    cur_data = squeeze(mean(mean(reshape(double(mov), [dx, Xsize/dx, dy, Ysize/dy, nframe]),3),1));
    
    %cur_data = double(mov);
    
    if size(cur_data,3) ~= 4000
        continue
    end
    
    disp(['Now serving #: ' num2str(i)]);
    
    if( i == 1 )

        % Get an ROI
        f = figure;
        subplot(1,3,1)
        refimg = squeeze(mean(cur_data,3));
        imshow(refimg, [], 'InitialMagnification', 'fit')
        hold on;

        [ysize, xsize] = size(refimg);
        npts = 1;
        colorindex = 0;
        order = get(gca,'ColorOrder');
        nroi = 1;
        intens = [];
        [x, y] = meshgrid(1:xsize, 1:ysize);
         
        %  
        [xv, yv] = (getline(gca, 'closed'));
        inpoly = inpolygon(x,y,xv,yv);
        
        %draw the bounding polygons and label them
        currcolor = order(1,:);
        plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
        text(mean(xv),mean(yv),num2str(1),'Color',currcolor,'FontSize',12);
    end
        
    % Get photo bleach correct and calculate df/f movie
    clear data_pb_rem;
    intens = squeeze(mean(mean(cur_data)));
    [intensN, pbleach] = rem_pbleach(intens, 400);
    
    [ysize, xsize, nframes] = size(cur_data);
    cur_data_pb_rem = cur_data./repmat(reshape(pbleach,[1,1,nframes]),[ysize,xsize,1]);   % photobleaching correction without bkg subtraction
    
    baseline = repmat(squeeze(mean(cur_data_pb_rem,3)), [1 1 size(cur_data_pb_rem,3)]);
    cur_data_norm = (cur_data_pb_rem-baseline)./baseline;
    
    itrace = squeeze(sum(sum(cur_data_norm.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
    
    subplot(1,3,2:3) % plot the trace
    hold on;
    plot(t,itrace,'Color', rgb('Bisque'), 'LineWidth', 1);
    
    data = data + cur_data_norm;
    cnt = cnt + 1;
end

data = data ./ cnt;

% Plot average intensity
avg_tc = squeeze(sum(sum(data.*repmat(inpoly, [1, 1, size(data,3)]))))/sum(inpoly(:));

subplot(1,3,2:3) % plot the trace
hold on;
plot(t,avg_tc,'Color', rgb('Brown'), 'LineWidth', 3);
xlabel('Time (s)', 'FontSize', 14);
ylabel('dF/F', 'FontSize', 14);
title(['Number of trials: ' num2str(cnt) '         Stim: ' num2str(STIM_PERIOD) ' s'], 'FontSize', 14);
ylim([-0.01 0.01]);
xlim([2.0 6.0]);
set(gca, 'FontSize', 14);

saveas(f, [datapath '/' search_dirs_savename '_tc_all_stim_' num2str(STIM_PERIOD) '.eps']);
saveas(f, [datapath '/' search_dirs_savename '_tc_all_stim_' num2str(STIM_PERIOD) '.png']);
saveas(f, [datapath '/' search_dirs_savename '_tc_all_stim_' num2str(STIM_PERIOD) '.fig']);

%%
BEGIN_TC = 1800;
%END_TC = 3600;
END_TC = 2600;

%fill([PRE_STIM PRE_STIM PRE_STIM+STIM_LEN PRE_STIM+STIM_LEN], [-0.05 0.05 0.05 -0.05], rgb('Bisque'));

%baseline = squeeze(mean(intens));
%intens_norm = (intens-baseline)./baseline;

f2 = figure;
rho = corr(avg_tc(BEGIN_TC:END_TC), reshape(data(:,:,BEGIN_TC:END_TC), [size(data,1)*size(data,2) size(data(:,:,BEGIN_TC:END_TC),3) ])' );
corr_img = reshape(rho', [size(data,1),  size(data,2)]);
imagesc( corr_img );
axis image;
colorbar;

saveas(f2, [datapath '/' search_dirs_savename '_corr_stim_' num2str(STIM_PERIOD) '.eps']);
saveas(f2, [datapath '/' search_dirs_savename '_corr_stim_' num2str(STIM_PERIOD) '.png']);
saveas(f2, [datapath '/' search_dirs_savename '_corr_stim_' num2str(STIM_PERIOD) '.fig']);

%% play movie of response

baseline = repmat(squeeze(mean(data_pb_rem,3)), [1 1 size(data_pb_rem,3)]);
data_pb_rem_norm = (data_pb_rem-baseline)./baseline;

my_video = VideoWriter(['01_response.avi']);
my_video.FrameRate = 50;
my_video.Quality = 100;
open(my_video);

f=figure;
for i=BEGIN_TC:END_TC
    imagesc(data_pb_rem_norm(:,:,i));
    colormap gray; colorbar; axis image;
    %if( i == BEGIN_TC ) c1 = caxis(); else caxis(c1); end
    caxis([-0.01 0.01]);
    title(['Frame: ' num2str(i)]);
    %pause(0.1);
    F = getframe(f);
    writeVideo(my_video,F);
end
close(my_video);
close(f);
