function [BW,maskedImage] = segmentImage(X)
%segmentImage Segment image using auto-generated code from Image Segmenter app
%  [BW,MASKEDIMAGE] = segmentImage(X) segments image X using auto-generated
%  code from the Image Segmenter app. The final segmentation is returned in
%  BW, and a masked image is returned in MASKEDIMAGE.

% Auto-generated by imageSegmenter app on 09-Oct-2023
%----------------------------------------------------


% Normalize input data to range in [0,1].
Xmin = min(X(:));
Xmax = max(X(:));
if isequal(Xmax,Xmin)
    X = 0*X;
else
    X = (X - Xmin) ./ (Xmax - Xmin);
end

% Graph cut
foregroundInd = [111374 111393 111403 111413 111419 114710 114795 116298 116305 116326 116339 116360 116368 121291 121506 121508 121509 121512 121513 127943 132955 136265 139577 139603 142723 142726 142728 142729 145310 145312 145313 145314 145316 145318 145320 145321 146066 146068 147735 147736 147738 148607 148614 148622 148627 148629 148638 148641 ];
backgroundInd = [59008 62350 67368 70703 70705 70707 70709 70712 70714 76909 76911 77399 79068 85757 90733 90739 90746 90749 90762 90768 90774 91979 98682 102430 114141 123798 132172 135523 145569 147603 150951 157647 162672 170679 177742 180723 184445 189093 196168 202485 207893 214595 217551 229677 230938 235956 239302 242648 249793 252625 254831 255975 256008 256012 256018 263215 266574 273288 ];
L = superpixels(X,1688);
BW = lazysnapping(X,L,foregroundInd,backgroundInd);

% Active contour
iterations = 50;
BW = activecontour(X, BW, iterations, 'Chan-Vese');

% Clear borders
BW = imclearborder(BW);

% Create masked image.
maskedImage = X;
maskedImage(~BW) = 0;
end

