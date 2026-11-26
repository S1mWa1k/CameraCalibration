clear
close all

%% Physical Parameters
nRows          = 36;     % Number of dot rows (vertical)
nCols          = 36;     % Number of dot columns (horizontal)
pitch_mm       = 0.75;   % Center-to-center dot spacing in mm (= 2 * dotDia_mm)
margin_mm      = 2.0;    % White border around outer dot edges in mm
dpi            = 2400;   % Output resolution (dots per inch)

ringOuterRatio = 1.2;
ringInnerRatio = 0.7;

% Specify [row, col] positions to convert into rings (1-based; top-left is [1, 1])
ringPositions = [
    17, 17;
    17, 18;
    17, 20;
    20, 17;
    20, 20
];

%% Physical-to-Pixel Conversion
pxPerMM   = dpi / 25.4;                 % Exact conversion (1 inch = 25.4 mm)
dotDia_mm = pitch_mm / 2;               % Edge-to-edge gap equals dot diameter

% Keep geometry in floating-point pixels to avoid cumulative spacing drift
pitch    = pitch_mm  * pxPerMM;
dotDia   = dotDia_mm * pxPerMM;

dotRad   = dotDia / 2;
outerRad = (ringOuterRatio * dotDia) / 2;
innerRad = (ringInnerRatio * dotDia) / 2;

% Overall physical and pixel dimensions
totalHeight_mm = 2 * margin_mm + (nRows - 1) * pitch_mm + dotDia_mm;
totalWidth_mm  = 2 * margin_mm + (nCols - 1) * pitch_mm + dotDia_mm;
imgHeight      = round(totalHeight_mm * pxPerMM);
imgWidth       = round(totalWidth_mm  * pxPerMM);

% Center the grid symmetrically inside the rounded pixel canvas
gridHeightPx = (nRows - 1) * pitch;
gridWidthPx  = (nCols - 1) * pitch;
startY       = (imgHeight - gridHeightPx) / 2 + 0.5;
startX       = (imgWidth  - gridWidthPx)  / 2 + 0.5;
rowCenters   = startY + (0:nRows-1) * pitch;
colCenters   = startX + (0:nCols-1) * pitch;

fprintf('Dot diameter: %.2f mm (%.2f px) | Pitch: %.2f mm (%.2f px)\n', ...
    dotDia_mm, dotDia, pitch_mm, pitch);
fprintf('Total image:  %.2f x %.2f mm (W x H) | %d x %d px at %d DPI\n', ...
    totalWidth_mm, totalHeight_mm, imgWidth, imgHeight, dpi);

%% Generate Image
img     = true(imgHeight, imgWidth);    % Logical array (true = white, false = black)
halfWin = ceil(max(dotRad, outerRad)) + 1;

for r = 1:nRows
    cy = rowCenters(r);
    for c = 1:nCols
        cx = colCenters(c);

        xIdx = max(1, floor(cx - halfWin)) : min(imgWidth,  ceil(cx + halfWin));
        yIdx = max(1, floor(cy - halfWin)) : min(imgHeight, ceil(cy + halfWin));
        [X, Y] = meshgrid(xIdx, yIdx);
        dist = hypot(X - cx, Y - cy);

        isRing = ~isempty(ringPositions) && ismember([r, c], ringPositions, 'rows');

        patch = img(yIdx, xIdx);
        if isRing
            patch(dist <= outerRad) = false;
            patch(dist <= innerRad) = true;
        else
            patch(dist <= dotRad) = false;
        end
        img(yIdx, xIdx) = patch;
    end
end

%% Display and Save with Physical DPI Metadata
figure;
imagesc(img, [0 1]);
colormap(gray);
axis image off;
title(sprintf('%dx%d Grid (R x C) | %.2f mm pitch | %d DPI', nRows, nCols, pitch_mm, dpi));
return
% PNG stores physical resolution in pixels per meter
pxPerMeter = round(dpi / 0.0254);
filename   = sprintf('dot_grid_%dx%d_%ddpi.png', nRows, nCols, dpi);
imwrite(img, filename, ...
    'ResolutionUnit', 'meter', ...
    'XResolution', pxPerMeter, ...
    'YResolution', pxPerMeter);

% Or save as TIFF, which stores DPI directly:
% imwrite(img, sprintf('dot_grid_%dx%d_%ddpi.tif', nRows, nCols, dpi), 'Resolution', dpi);