imPath = '/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/H&EbinaryImages/';
wellPath = '/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/';
sgeID = importdata('/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/lenas.txt');

%% Takki_leFP_003_A1

d = dir(fullfile(imPath, '*.mat'));
fname = sgeID{1};
BW = load(fullfile(imPath, fname));
jsonname = fullfile(wellPath, sgeID{1}, '/scalefactors_json.json');
w = jsondecode(fileread(jsonname));
R = ceil(w.spot_diameter_fullres/2);
tbl = readtable(fullfile(wellPath, sgeID{1}, '/tissue_positions_list.txt'));
nSpots = size(tbl, 1);
count = zeros(nSpots, 1);
mask = zeros(size(BW.mask6));
crow = table2array(tbl(:, 5));
ccol = table2array(tbl(:, 6));

for i = 1:nSpots
    mask(crow(i), ccol(i)) = 1;
end

mask = bwdist(mask) <= R;
mask = bwlabel(mask);
BW = bwlabel(BW.mask6);
parfor i = 1:nSpots
    idx = mask(crow(i), ccol(i));
    %tmpBW = BW;
    %tmpBW(mask~=idx) = 0;
    %[~, c] = bwlabel(tmpBW);
    spot = BW(mask==idx & BW>0);
    c = length(unique(spot));
    count(i) = c;
end

tbl = [tbl array2table(count)];
tbl.Properties.VariableNames = {'barcode','tissue','row','col','imagerow','imagecol','count'};

mkdir(fullfile('/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/results/Histology/', fname));
writetable(tbl, fullfile('/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/results/Histology/', fname, '/tissue_spot_counts.csv'), 'Delimiter', ',');

%% Takki_leFP_003_B1

d = dir(fullfile(imPath, '*.mat'));
fname = sgeID{2};
BW = load(fullfile(imPath, fname));
jsonname = fullfile(wellPath, sgeID{2}, '/scalefactors_json.json');
w = jsondecode(fileread(jsonname));
R = ceil(w.spot_diameter_fullres/2);
tbl = readtable(fullfile(wellPath, sgeID{2}, '/tissue_positions_list.txt'));
nSpots = size(tbl, 1);
count = zeros(nSpots, 1);
mask = zeros(size(BW.mask5));
crow = table2array(tbl(:, 5));
ccol = table2array(tbl(:, 6));

for i = 1:nSpots
    mask(crow(i), ccol(i)) = 1;
end

mask = bwdist(mask) <= R;
mask = bwlabel(mask);
BW = bwlabel(BW.mask5);
parfor i = 1:nSpots
    idx = mask(crow(i), ccol(i));
    %tmpBW = BW;
    %tmpBW(mask~=idx) = 0;
    %[~, c] = bwlabel(tmpBW);
    spot = BW(mask==idx & BW>0);
    c = length(unique(spot));
    count(i) = c;
end

tbl = [tbl array2table(count)];
tbl.Properties.VariableNames = {'barcode','tissue','row','col','imagerow','imagecol','count'};

mkdir(fullfile('/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/results/Histology/', fname));
writetable(tbl, fullfile('/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/results/Histology/', fname, '/tissue_spot_counts.csv'), 'Delimiter', ',');


%% Takki_leFP_004_C1


d = dir(fullfile(imPath, '*.mat'));
fname = sgeID{3};
BW = load(fullfile(imPath, fname));
jsonname = fullfile(wellPath, sgeID{3}, '/scalefactors_json.json');
w = jsondecode(fileread(jsonname));
R = ceil(w.spot_diameter_fullres/2);
tbl = readtable(fullfile(wellPath, sgeID{3}, '/tissue_positions_list.txt'));
nSpots = size(tbl, 1);
count = zeros(nSpots, 1);
mask = zeros(size(BW.mask3));
crow = table2array(tbl(:, 5));
ccol = table2array(tbl(:, 6));

for i = 1:nSpots
    mask(crow(i), ccol(i)) = 1;
end

mask = bwdist(mask) <= R;
mask = bwlabel(mask);
BW = bwlabel(BW.mask3);
parfor i = 1:nSpots
    idx = mask(crow(i), ccol(i));
    %tmpBW = BW;
    %tmpBW(mask~=idx) = 0;
    %[~, c] = bwlabel(tmpBW);
    spot = BW(mask==idx & BW>0);
    c = length(unique(spot));
    count(i) = c;
end

tbl = [tbl array2table(count)];
tbl.Properties.VariableNames = {'barcode','tissue','row','col','imagerow','imagecol','count'};

mkdir(fullfile('/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/results/Histology/', fname));
writetable(tbl, fullfile('/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/results/Histology/', fname, '/tissue_spot_counts.csv'), 'Delimiter', ',');


%% Takki_leFP_004_D1


d = dir(fullfile(imPath, '*.mat'));
fname = sgeID{4};
BW = load(fullfile(imPath, fname));
jsonname = fullfile(wellPath, sgeID{4}, '/scalefactors_json.json');
w = jsondecode(fileread(jsonname));
R = ceil(w.spot_diameter_fullres/2);
tbl = readtable(fullfile(wellPath, sgeID{4}, '/tissue_positions_list.txt'));
nSpots = size(tbl, 1);
count = zeros(nSpots, 1);
mask = zeros(size(BW.mask7));
crow = table2array(tbl(:, 5));
ccol = table2array(tbl(:, 6));

for i = 1:nSpots
    mask(crow(i), ccol(i)) = 1;
end

mask = bwdist(mask) <= R;
mask = bwlabel(mask);
BW = bwlabel(BW.mask7);
parfor i = 1:nSpots
    idx = mask(crow(i), ccol(i));
    %tmpBW = BW;
    %tmpBW(mask~=idx) = 0;
    %[~, c] = bwlabel(tmpBW);
    spot = BW(mask==idx & BW>0);
    c = length(unique(spot));
    count(i) = c;
end

tbl = [tbl array2table(count)];
tbl.Properties.VariableNames = {'barcode','tissue','row','col','imagerow','imagecol','count'};

mkdir(fullfile('/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/results/Histology/', fname));
writetable(tbl, fullfile('/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/results/Histology/', fname, '/tissue_spot_counts.csv'), 'Delimiter', ',');

