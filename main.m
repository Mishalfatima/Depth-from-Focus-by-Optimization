% clear all
% 
 data = '/home/mishal/Desktop/CV/PA1/PA1_dataset1_balls/'; 
 cd(data)
% 
% % Step 1: Image alignment
 files = dir('*.jpg');
  n0_of_files = length(files)-1;% Since filenames start at zero
 
 ref = imread(strcat(data,int2str(n0_of_files),'.jpg'));  % The image with the nearest focus is used as first ref
 transform = 'affine';
 
 mkdir('warp');
 for i = n0_of_files:-1:0
     fname = strcat(data,int2str(i),'.jpg');
     im = imread(fname);
     
     %Extracting SURF descriptors
     [dr,lr] = iat_surf(ref);
     [d , l] = iat_surf(im);
     [~, ~, refInd, imInd] = iat_match_features_mex(dr,d,.7);
     ptsA = lr(refInd,1:2);
     ptsB = l(imInd,1:2);
     
     %Using Ransac to determine transformation
     [inliers, ransacWarp] = iat_ransac(iat_homogeneous_coords(ptsA'),iat_homogeneous_coords(ptsB'),...
     transform,'tol',.05, 'maxInvalidCount', 10);
 
     %Warping image
     [wimage, support] = iat_inverse_warping(im, ransacWarp, transform, 1:size(ref,2),1:size(ref,1));
     im_name = strcat('warp/',int2str(i),'-warped.jpg');
     imwrite(uint8(wimage),strcat(im_name));
     ref = uint8(wimage);
 end

%step 2: Finding Focus Measure

data = '/home/mishal/Desktop/CV/PA1/PA1_dataset1_balls/'; 
dataset = 'PA1_dataset1_balls/';
cd(data)

files = dir('*.jpg');
no_of_files = length(files);

Focus_Measure = [];
p="/home/mishal/Desktop/CV/PA1/PA1_dataset1_balls/warp/"

path=strcat(p,"0-warped.jpg");


img=imread(path);
[ROWS,COLS,~]=size(img);

%Calling OTF function
H=OTF(ROWS,COLS);

for k= 1:no_of_files
    
    %Reading warped images
    path=strcat(p,int2str(k-1),"-warped.jpg");
    img=imread(path);
    [ROWS,COLS,~] = size(img);
    
    H=OTF(ROWS,COLS);
    
    %Finding spectrum of each intensity image
    spectrum=fft2(img);
    FI=fftshift(spectrum);
    max1 = max(FI);
    max2 = max(max1);
    scale = 1.0/max2;
    spectrum = FI.*scale;
    
    %Finding Focus Measure by Multiplying with OTF
    FHI = spectrum.*H';
    HI = ifft2(FHI);
    max1 = max(HI);
    max2 = max(max1);
    scale = 1.0/max2;
    FM = HI.*scale;
   
    %Concatenating Focus Measures together
    Focus_Measure = cat(4,Focus_Measure,FM);
    
    
end

%Finding Mf and Mp
[Rows,Cols,ch] = size(img);
Mp = zeros(Rows,Cols);        % Matrix of best focus values
Mf = zeros(Rows,Cols);        % Matrix of frames per pixel having the best focus

for i = 1:Rows
    for j = 1:Cols
        a = zeros(1,2); id = 0;
        [a(1,1),a(1,2)] = max(Focus_Measure(i,j,1,:));
        [a(2,1),a(2,2)] = max(Focus_Measure(i,j,2,:));
        [a(3,1),a(3,2)] = max(Focus_Measure(i,j,3,:));
        [Mp(i,j),id] = max(a(:,1));
        Mf(i,j) = a(id,2);
    end
end

figure;
RGB = label2rgb(Mf);
mkdir('output');
imwrite(RGB,'output/depthmap.jpg');
imwrite(uint8(Mp),'output/edges.jpg');
save('output/Mf.mat','Mf');

raw=[]
%Concatenating all original warped images together
for i = 0 : no_of_files
    path=strcat(p,int2str(k-1),"-warped.jpg");
    A = imread(path);
    raw = cat(4, raw, A);
    
end

%All-in-focus image  using raw images and Mf(This step will be repeated for Graphcuts)
[r,c] = size(Mf);
ch = 3;
stitched = zeros(r,c,ch);

for i = 1:r
    for j = 1:c
        f = Mf(i,j);
        for chan = 1:ch
            stitched(i,j,chan) = raw(i,j,chan,f);
        end
    end
end

figure;

imwrite(uint8(stitched),'output/out-01.jpg');

figure;

imwrite(uint8(Mp*256),'output/edges-01.jpg');
save('output/Mf.mat','Mf');

%Step 3: Graph cuts

[rows, cols, ~]=size(Mf)

Data_Cost = zeros(rows, cols, no_of_files+1);    
Smoothing_Cost = zeros(no_of_files+1,no_of_files+1);

%First calculate Data_Cost
for i = 1:no_of_files+1
    Data_Cost(:,:,i) = min(2,abs(Mf-i));      % Thresholding data cost, t=2

end

%Then calculate Smoothing_Cost
for i = 1:no_of_files+1
    for j = 1:no_of_files+1
        if(i~=j) 
            Smoothing_Cost(i,j) = min(abs(i-j),10); % Thresholding smoothing cost, t=12
            Smoothing_Cost(j,i) = Smoothing_Cost(i,j); 
        end
    end
end

[gch] = GraphCut('open', Data_Cost, Smoothing_Cost);

[gch L] = GraphCut('expand',gch);
gch = GraphCut('close', gch);

L_n = L+1;
Graphcut = label2rgb(L_n);

% cd ..
% cd(dataset)
save('output/L.mat','L_n');
imwrite(Graphcut,'output/graph_cut_depth.jpg');
%  
% %% Step 4: All-in-focus image by using result of GraphCuts
[r,c] = size(L_n);
ch = 3;
out_gc = zeros(r,c,ch);

for i = 1:r
    for j = 1:c
        f = L_n(i,j);
        for chan = 1:ch
            out_gc(i,j,chan) = raw(i,j,chan,f);
        end
    end
end

figure;

imwrite(uint8(out_gc),'output/out-GraphCut.jpg');

% %% Step 5: Depth Refinement
Z = double(out_gc);
mask_size = 3;
mask_ctr = ceil(mask_size/2);
mask_start = mask_ctr - 1;
W = ones(mask_size);
W(mask_ctr, mask_ctr) = 4;
sum_W = sum(W(:));
W = W./sum_W; 

image_padded = padarray(Z,[mask_start, mask_start],'symmetric','both');
[row, col, ~] = size(Z);
row_pad = row + 2 * mask_start;
col_pad = col + 2 * mask_start;

for x = mask_ctr:row_pad - mask_start - 1
    for y = mask_ctr:col_pad - mask_start - 1
        %% To make a 3x3 weighted mask into a 1x9 mask
        patch_selected = image_padded(x - mask_start: x + mask_start, y - mask_start: y + mask_start);
        a1 = W.*patch_selected;
        med = median(a1(:))*sum_W; % the5th value is the weighted median
        Z(x,y) = med;
    end
end

imwrite(uint8(Z), 'output/out_median.jpg');



