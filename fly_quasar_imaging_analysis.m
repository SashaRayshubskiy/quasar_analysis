%% clicky on data
dfile = [ '03_640nm_0OD_100ms_10001frm_0x_135z_QuasAr1_20xobj_AL' ];

clear data;
data = ReadTifFast([dfile '.tif']);
%data = data(:,:,2:end);
BEGIN_FRAME = 1000;
data = data(:,:,BEGIN_FRAME:4000);

baseline = repmat(squeeze(mean(data,3)), [1 1 size(data,3)]);
data_norm = (data-baseline)./baseline;

% Photobleach correct the data.
clear data_pb_rem;
intens = squeeze(mean(mean(data)));
[intensN, pbleach] = rem_pbleach(intens, 1000);
[ysize, xsize, nframes] = size(data);
data_pb_rem = data./repmat(reshape(pbleach,[1,1,nframes]),[ysize,xsize,1]);   % photobleaching correction without bkg subtraction

clicky_flies1(data);

%% Plot GCaMP, red laser, GCaMP data: modified from clicky

dfiles = { [ '01_488nm_.7V_0OD_100ms_101frm_0x_135z_GFP_20xobj_AL' ], ...
           [ '02_488nm_.7V_0OD_100ms_101frm_0x_135z_GFP_20xobj_AL_post 12s red light' ], ...
           [ '04_488nm_.7V_0OD_100ms_101frm_0x_135z_GFP_20xobj_AL_post 40s red light' ], ...
           [ '05_488nm_.7V_0OD_100ms_101frm_0x_135z_GFP_20xobj_AL_post 120s red light' ], ...
           [ '06_488nm_.7V_0OD_100ms_101frm_0x_135z_GFP_20xobj_AL_post 300s red light' ], ...
           [ '07_488nm_.7V_0OD_100ms_101frm_0x_135z_GFP_20xobj_AL_post 600s red light' ], ...
           };

colorindex = 0;

order = get(gca,'ColorOrder');

for i=1:size(dfiles,2)
    
    clear tmp tmp1;
    tmp = ReadTifFast([dfiles{i} '.tif']);

    tmp0 = tmp(:,:,2:end);
    
    baseline = repmat(squeeze(mean(tmp0,3)), [1 1 size(tmp0,3)]);
    tmp1 = (tmp0-baseline)./baseline;
    
    if i == 1
        % show the image
        refimg = squeeze(mean(tmp0(:,:,:),3));
        
        figure;
        subplot(1,3,1);    
        imshow(refimg, [], 'InitialMagnification', 'fit');        
        hold on;
        
        [ysize, xsize] = size(refimg(:,:,1));
        [x, y] = meshgrid(1:xsize, 1:ysize);
        
        [xv, yv] = (getline(gca, 'closed'));
        inpoly = inpolygon(x,y,xv,yv); 
        plot(xv, yv, 'b', 'Linewidth', 1);
    end

    % Get time course
    currcolor = order(1+mod(colorindex,size(order,1)),:);
    
    nframes = size( tmp1, 3 );
    
    itrace = squeeze(sum(sum(tmp1.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
    
    subplot(1,3,2:3) % plot the trace
    
    hold on;
    
    plot(itrace,'Color',currcolor);
    
    colorindex = colorindex + 1;
end

legend(['Before 640nm'],['After 12 sec'],['After 40 sec'],['After 120 sec'],['After 300 sec'],['After 600 sec']);

%% Get the frequency of the 'response'
BEGIN_TC = 2500; 
END_TC   = 2700; 
FR = 500;

[roi_pts, intens] = clicky(data);
figure;
hold on;
[P,m] = pwelch( intens(BEGIN_TC:END_TC), 128, 120, 128, FR );
hold on;
plot( m, 10*log10(P) );

%% Create a Fourier image on the onset and offset of the response.

fourier_image = zeros( size(data,1), size(data,2) );

[P,m] = pwelch( squeeze(data(50, 50, BEGIN_TC:END_TC)), 128, 120, 128, FR );

m_range = find( m > 58.0 & m < 62.0 );

for x = 1:size(data,1)
    for y = 1:size(data,2)

        [P,m] = pwelch( squeeze(data(x, y, BEGIN_TC:END_TC)), 128, 120, 128, FR );
        Power = 10*log10(P);

        fourier_image(x,y) = mean(Power(m_range));

    end
    disp(['Processing x: ' num2str(x)]);
end

f = figure;
imagesc(fourier_image);

%% Calculate the Pearson correlation
BEGIN_TC = 2500 - BEGIN_FRAME;
END_TC = 3500 - BEGIN_FRAME;

[roi_pts, intens] = clicky(data_pb_rem);

%baseline = squeeze(mean(intens));
%intens_norm = (intens-baseline)./baseline;

f = figure;
rho = corr(intens(BEGIN_TC:END_TC), reshape(data_pb_rem(:,:,BEGIN_TC:END_TC), [size(data,1)*size(data,2) size(data(:,:,BEGIN_TC:END_TC),3) ])' );
corr_img = reshape(rho', [size(data,1),  size(data,2)]);
imagesc( corr_img );
axis image;
colorbar;

%%
clicky_flies1(data_pb_rem);

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
