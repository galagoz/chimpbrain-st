%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we must run it on Matlab2020 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lenas_orig=importdata('/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/lenas.txt');

for x=1:length(lenas_orig)
    imgRGB = imread(['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/folded_regions/',lenas_orig{x},'.png']);
    imgRGB = im2double(imgRGB);
    tbl = readtable(['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/',lenas_orig{x},'/tissue_positions_list.csv']);
    w = jsondecode(fileread(['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/',lenas_orig{x},'/scalefactors_json.json']));
    mask=zeros(size(imgRGB,1),size(imgRGB,2));
    % R = ceil(w.spot_diameter_fullres*w.tissue_lowres_scalef/2);
    R=1; % we must use 1, because it's low-resolution
    crow = round(table2array(tbl(:, 5))*w.tissue_lowres_scalef);
    ccol = round(table2array(tbl(:, 6))*w.tissue_lowres_scalef);
    for k = 1:length(table2array(tbl(:, 5)))
        mask(crow(k), ccol(k)) = 1;
    end
    mask = bwdist(mask) <= R;
    mask = bwlabel(mask);
    % extract the color from the map with labelled regions
    mask_supervised=zeros(size(imgRGB,1),size(imgRGB,2),size(imgRGB,3));
    folded_label=zeros(length(table2array(tbl(:,5))),1);
    for i=1:length(table2array(tbl(:,5)))
        ind=imgRGB(crow(i),ccol(i),:);
        ind3(i,:)=ind(:);
        ind2=mask(crow(i),ccol(i));
        [hang lie]=find(mask==ind2);
        ind3(ind3<0.4)=0;
        %     ind4=isempty(find(ind3(i,:)==0));
        %     if ind4==0
        %         ind3(i,:)=0;
        %     end
        for j=1:length(hang)
            mask_supervised(hang(j),lie(j),1)=ind3(i,1);
            mask_supervised(hang(j),lie(j),2)=ind3(i,2);
            mask_supervised(hang(j),lie(j),3)=ind3(i,3);
        end
        if sum(ind3(i,:))==0
            folded_label(i,1)=1;
        end
    end
    imwrite(mask_supervised,['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/folded_regions/images_alignment/',lenas_orig{x},'.png']);
    
    % save a list with barcode and its color, so that we can identify the regions by colors
    data=importdata(['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/folded_regions/data/',lenas_orig{x},'.csv']);
    barcode_filtered=data.textdata(2:end,1);
    [c ia ib]=intersect(barcode_filtered, table2array(tbl(:,1)),'stable');
    folded_label_filtered=folded_label(ib);
    fid=fopen(['/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/folded_regions/data/',lenas_orig{x},'_fold.csv'],'w');
    fprintf(fid,'%s\n','barcode,imagecol,imagerow,fold');
    for kk=1:length(barcode_filtered)
        fprintf(fid,'%s,',barcode_filtered{kk});
        fprintf(fid,'%s,',num2str(data.data(kk,1)));
        fprintf(fid,'%s,',num2str(data.data(kk,2)));
        fprintf(fid,'%s\n',num2str(folded_label_filtered(kk)));
    end
    fclose(fid);
    clearvars -except lenas lenas_orig;
end


