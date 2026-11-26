clear
close all

patternDims = [14 14];
family = 'DICT_5X5_250';
% imageSize = patternDims*6e2;
checkerSize = 100;
markerSize = 75;

imageFinalSize=[6702 6702];

real_size=[56 56]; % size of slide in mm
target_size=[28 28]; % required size of grid in mm

imageSize=imageFinalSize.*target_size./real_size; % required size
imageSizeUse=patternDims.*round(imageSize./patternDims); % scale to closest multiple of patternDims

pad_rows = 0.5*(imageFinalSize(1) - imageSizeUse(1));
pad_cols = 0.5*(imageFinalSize(2) - imageSizeUse(2));

% pad_rows = ceil(0.05*(imageFinalSize(1) - imageSizeUse(1)));
% pad_cols = ceil(0.05*(imageFinalSize(2) - imageSizeUse(2)));

I = generateCharucoBoard(imageSizeUse,patternDims,family,checkerSize,markerSize);

I2 = padarray(I, [pad_rows, pad_cols], 255, 'both');

% I3 = insertText(I2,[0.5*imageFinalSize(2) abs(0.25*(imageFinalSize(1)-imageSizeUse(1)))],['ChArUco ' family ' ' num2str(patternDims(1)) 'x' num2str(patternDims(2))],...
%     'font','Sans Serif Collection','fontsize',200,'TextBoxColor','white','AnchorPoint','Center');

BW=logical(I2(:,:,1));
% BW=1-logical(I3(:,:,1));
imagesc2(BW)

imwrite(BW,['ChArUco' family '_' num2str(patternDims(1)) 'x' num2str(patternDims(2)) '.tif'],'tif')

