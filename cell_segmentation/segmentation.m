
lenas_orig = importdata('/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/lenas.txt');

%% Takki_leFP_003_A1
% Read in histology image with imread function %
imgRGB = imread(['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{1},'.png']);

% convert image datatype to double %
imgRGB = im2double(imgRGB);

% These numbers are manually adjusted for each image until desired contrast is obtained %
imgRGB_smooth_adj = imadjust(imgRGB, [.2 .3 0; .6 .7 1],[]);

%Applying Kmeans color based segmentation on the above resulting image
%following the technique mentioned in https://www.mathworks.com/help/images/color-based-segmentation-using-k-means-clustering.html%
he = imgRGB_smooth_adj;
lab_he = rgb2lab(he); %convert image from rgb color space to lab color space%
ab = lab_he(:,:,2:3);
ab = im2single(ab);
nColors = 8; %user defined number, number of colors visually seen by the user in a image%
pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',3); %apply Kmeans

%resulting clusters from Kmeans
mask1 = pixel_labels==1; %mask* is a binary segmented image per color
cluster1 = he .* double(uint8(mask1)); %cluster* is a colored segmented image per color
mask2 = pixel_labels==2;
cluster2 = he .* double(uint8(mask2));
mask3 = pixel_labels==3;
cluster3 = he .* double(uint8(mask3));
mask4 = pixel_labels==4;
cluster4 = he .* double(uint8(mask4));
mask5 = pixel_labels==5;
cluster5 = he .* double(uint8(mask5));
mask6 = pixel_labels==6;
cluster6 = he .* double(uint8(mask6));
mask7 = pixel_labels==7;
cluster7 = he .* double(uint8(mask7));
mask8 = pixel_labels==8;
cluster8 = he .* double(uint8(mask8));

%one of these five cluster/binary images is the segmented nuclei image
%Usually nuclei in histology images have a distinct dark color and the
%nuclei are clearly segmented from the background

%imwrite(cluster1,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{1},'imgRGB_smooth_adjus_cluster1.png'])
%imwrite(cluster2,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{1},'imgRGB_smooth_adjus_cluster2.png'])
%imwrite(cluster3,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{1},'imgRGB_smooth_adjus_cluster3.png'])
%imwrite(cluster4,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{1},'imgRGB_smooth_adjus_cluster4.png'])
%imwrite(cluster5,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{1},'imgRGB_smooth_adjus_cluster5.png'])
%imwrite(cluster6,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{1},'imgRGB_smooth_adjus_cluster6.png'])
%imwrite(cluster7,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{1},'imgRGB_smooth_adjus_cluster7.png'])
%imwrite(cluster8,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{1},'imgRGB_smooth_adjus_cluster8.png'])
%save(['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/',lenas_orig{1},'/',lenas_orig{1},'.mat'],'mask6');
%% Takki_leFP_003_B1
% Read in histology image with imread function %
imgRGB = imread(['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{2},'.png']);

% convert image datatype to double %
imgRGB = im2double(imgRGB);

% These numbers are manually adjusted for each image until desired contrast is obtained %
imgRGB_smooth_adj = imadjust(imgRGB, [.2 .3 0; .6 .7 1],[]);

%Applying Kmeans color based segmentation on the above resulting image
%following the technique mentioned in https://www.mathworks.com/help/images/color-based-segmentation-using-k-means-clustering.html%
he = imgRGB_smooth_adj;
lab_he = rgb2lab(he); %convert image from rgb color space to lab color space%
ab = lab_he(:,:,2:3);
ab = im2single(ab);
nColors = 7; %user defined number, number of colors visually seen by the user in a image%
pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',3); %apply Kmeans

%resulting clusters from Kmeans
mask1 = pixel_labels==1; %mask* is a binary segmented image per color
cluster1 = he .* double(uint8(mask1)); %cluster* is a colored segmented image per color
mask2 = pixel_labels==2;
cluster2 = he .* double(uint8(mask2));
mask3 = pixel_labels==3;
cluster3 = he .* double(uint8(mask3));
mask4 = pixel_labels==4;
cluster4 = he .* double(uint8(mask4));
mask5 = pixel_labels==5;
cluster5 = he .* double(uint8(mask5));
mask6 = pixel_labels==6;
cluster6 = he .* double(uint8(mask6));
mask7 = pixel_labels==7;
cluster7 = he .* double(uint8(mask7));

%one of these five cluster/binary images is the segmented nuclei image
%Usually nuclei in histology images have a distinct dark color and the
%nuclei are clearly segmented from the background

%imwrite(cluster1,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{2},'imgRGB_smooth_adjus_cluster1.png'])
%imwrite(cluster2,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{2},'imgRGB_smooth_adjus_cluster2.png'])
%imwrite(cluster3,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{2},'imgRGB_smooth_adjus_cluster3.png'])
%imwrite(cluster4,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{2},'imgRGB_smooth_adjus_cluster4.png'])
%imwrite(cluster5,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{2},'imgRGB_smooth_adjus_cluster5.png'])
%imwrite(cluster6,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{2},'imgRGB_smooth_adjus_cluster6.png'])
%imwrite(cluster7,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{2},'imgRGB_smooth_adjus_cluster7.png'])
%save(['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/',lenas_orig{2},'/',lenas_orig{2},'.mat'],'mask5');
%% Takki_leFP_004_C1
% Read in histology image with imread function %
imgRGB = imread(['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{3},'.png']);

% convert image datatype to double %
imgRGB = im2double(imgRGB);

% These numbers are manually adjusted for each image until desired contrast is obtained %
imgRGB_smooth_adj = imadjust(imgRGB, [.2 .3 0; .6 .7 1],[]);

%Applying Kmeans color based segmentation on the above resulting image
%following the technique mentioned in https://www.mathworks.com/help/images/color-based-segmentation-using-k-means-clustering.html%
he = imgRGB_smooth_adj;
lab_he = rgb2lab(he); %convert image from rgb color space to lab color space%
ab = lab_he(:,:,2:3);
ab = im2single(ab);
nColors = 7; %user defined number, number of colors visually seen by the user in a image%
pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',3); %apply Kmeans

%resulting clusters from Kmeans
mask1 = pixel_labels==1; %mask* is a binary segmented image per color
cluster1 = he .* double(uint8(mask1)); %cluster* is a colored segmented image per color
mask2 = pixel_labels==2;
cluster2 = he .* double(uint8(mask2));
mask3 = pixel_labels==3;
cluster3 = he .* double(uint8(mask3));
mask4 = pixel_labels==4;
cluster4 = he .* double(uint8(mask4));
mask5 = pixel_labels==5;
cluster5 = he .* double(uint8(mask5));
mask6 = pixel_labels==6;
cluster6 = he .* double(uint8(mask6));
mask7 = pixel_labels==7;
cluster7 = he .* double(uint8(mask7));

%one of these five cluster/binary images is the segmented nuclei image
%Usually nuclei in histology images have a distinct dark color and the
%nuclei are clearly segmented from the background

imwrite(cluster1,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{3},'imgRGB_smooth_adjus_cluster1.png'])
imwrite(cluster2,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{3},'imgRGB_smooth_adjus_cluster2.png'])
imwrite(cluster3,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{3},'imgRGB_smooth_adjus_cluster3.png'])
imwrite(cluster4,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{3},'imgRGB_smooth_adjus_cluster4.png'])
imwrite(cluster5,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{3},'imgRGB_smooth_adjus_cluster5.png'])
imwrite(cluster6,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{3},'imgRGB_smooth_adjus_cluster6.png'])
imwrite(cluster7,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{3},'imgRGB_smooth_adjus_cluster7.png'])
save(['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/',lenas_orig{3},'/',lenas_orig{3},'.mat'],'mask3');
%% Takki_leFP_004_D1
% Read in histology image with imread function %
imgRGB = imread(['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{4},'.png']);

% convert image datatype to double %
imgRGB = im2double(imgRGB);

% These numbers are manually adjusted for each image until desired contrast is obtained %
imgRGB_smooth_adj = imadjust(imgRGB, [.2 .3 0; .6 .7 1],[]);

%Applying Kmeans color based segmentation on the above resulting image
%following the technique mentioned in https://www.mathworks.com/help/images/color-based-segmentation-using-k-means-clustering.html%
he = imgRGB_smooth_adj;
lab_he = rgb2lab(he); %convert image from rgb color space to lab color space%
ab = lab_he(:,:,2:3);
ab = im2single(ab);
nColors = 7; %user defined number, number of colors visually seen by the user in a image%
pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',3); %apply Kmeans

%resulting clusters from Kmeans
mask1 = pixel_labels==1; %mask* is a binary segmented image per color
cluster1 = he .* double(uint8(mask1)); %cluster* is a colored segmented image per color
mask2 = pixel_labels==2;
cluster2 = he .* double(uint8(mask2));
mask3 = pixel_labels==3;
cluster3 = he .* double(uint8(mask3));
mask4 = pixel_labels==4;
cluster4 = he .* double(uint8(mask4));
mask5 = pixel_labels==5;
cluster5 = he .* double(uint8(mask5));
mask6 = pixel_labels==6;
cluster6 = he .* double(uint8(mask6));
mask7 = pixel_labels==7;
cluster7 = he .* double(uint8(mask7));

%one of these five cluster/binary images is the segmented nuclei image
%Usually nuclei in histology images have a distinct dark color and the
%nuclei are clearly segmented from the background

imwrite(cluster1,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{4},'imgRGB_smooth_adjus_cluster1.png'])
imwrite(cluster2,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{4},'imgRGB_smooth_adjus_cluster2.png'])
imwrite(cluster3,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{4},'imgRGB_smooth_adjus_cluster3.png'])
imwrite(cluster4,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{4},'imgRGB_smooth_adjus_cluster4.png'])
imwrite(cluster5,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{4},'imgRGB_smooth_adjus_cluster5.png'])
imwrite(cluster6,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{4},'imgRGB_smooth_adjus_cluster6.png'])
imwrite(cluster7,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/perSpotCellCount/',lenas_orig{4},'imgRGB_smooth_adjus_cluster7.png'])
save(['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/',lenas_orig{4},'/',lenas_orig{4},'.mat'],'mask7');
