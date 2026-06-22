 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Written By Azadeh Jafari											%%
%%2017									                    	    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
commandwindow

            %% *****read MRI data*****
%%We should add NIFTI toolbox to matlab that was impossible for my version
%%of matlab so I add them here

nii = load_nii('NAME.nii');
% view_nii(nii);%see Image in 1st:coronal 2th:sagittal 3rd:transverse VIEW
map=nii.img;

prompt = 'please Enter slice number : ';
slice_number = input(prompt);

if isempty(slice_number)
    slice_number = 'ENTER SLICE NUMBER';
end

map = imrotate(map,-90);%rotate image clockwise to see better
conn = 4;

figure
imshow(map(:,:,slice_number),[]);%choose slice number 160 from MRICron transverse plane(X:saggital Y:frontal Z:transverse)
title('Original Quantitative MRI','Color', 'm');
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);%show full screen

            %% ******Preprocessing*****
%% SKULL STRIPPING ALGORITHM FOR BRAIN MRI by ROY & MAJI 2015

%% Step1: skull Stripping the brain map with morphological theoeriem
skull_Stripped_Image= skull_stripping_preprocess(map,slice_number);
%skull_Stripped_Image=imcrop(map(:,:,slice_number),[]);
% figure;
% imshow(skull_Stripped_Image,[]);
% skull_Stripped_Image=map(:,:,slice_number);
picture =  double(skull_Stripped_Image);

%% step2 :we can apply median filter (optional)
% medfilt_skull_Stripped_Image = medfilt2(skull_Stripped_Image);
% figure;
% imshow(medfilt_skull_Stripped_Image,[]);
% title('skull-striped filtered MRI','Color', 'r');

%% *****labeling with clustering Kmeans and segmentation result with Kmeans*****

%%based on principal of graph cut @(IMAGE SEG:A SURVEY OF GRAPH CUT METHODS by fALIU YI 2012) 
%%We should dedicate initial labeling to pixels based on k-means and then use classification algorithm
[ind, cntr] = kmeans( skull_Stripped_Image(:),2 );

%%Distance: distances from each point to every centroids
%%ind:cluster indices of each observation(size of idx=r*c)
%%Cntr: clusters centroid locations
%%sumd: within-cluster sums of point-to-centroid distances for each cluster

cntr = round(cntr'); %Round to the nearest decimal or integer
f = reshape(ind, size(skull_Stripped_Image));%specify indices to each pixel of images
f(~isnan(f)) = cntr(1);%all the pixels that are not Nan must be equal to center(1) and...
%this checks robustness of the algorithm to the initial labeling, by providing a constant initial labeling
cntr = sort(cntr);
kmean_ind = ind; 
kmean_ind(isnan(kmean_ind)) = 0;
% figure, imshow( label2rgb( reshape(kmean_ind,size(skull_Stripped_Image)),[0 1 1; 1 0 1; 1 1 0],'k') ), title('segmentation result by kmeans method')


%% *****create custom prabability distribiution function fitted to a training sample data****
%%based on AN EXPERIMENTAL COMPARISON OF MIN CUT/MAX FLOW FOR ENERGY MINIMIZATION by YURI BOYKOV 2004 page2/34

%%Explanation:data penalties D(·) indicate individual label-preferences of pixels based on ...
%%observed intensities and pre-specifed likelihood function
%%call function handle by @,,

%%For regional term of the cost function, we should have a PD based on
%%seeds (Guassian Kernel), so that we can calculate the weight accordingly
PD1 = @(input) .4*normpdf(input,1687,266.5) + .6*normpdf(input,1212.6,206);
PD2 = @(input) 1; % logarithm of 1 is 0. so when it is summed with the logarithm of probabilities obtained from seeds, it will have no effect
prior_likelihood = {PD1, PD2};

%% *****calculate camera noise*****
%%based on GRAPH CUT AND EFFICIENT N-D IAMGE SEG by YURI BOYKOV 2006 PAGE
%%10 and IMAGE SEG:A SURVEY OF GRAPH CUT METHOD by FALIU YI 2012

%%camera_noise is the intensity diffrence between two adjacent nodes (n=8) that is expected to be normal...
%%if they belong to the same region, we calculate diffrences between 2 adjacent node for all pixels...
%%and then mean of them can use as camera noise for all the Image

%%Finglish: baraye be dast avardane camera noise ekhtelafe intensitye har pixel ba
%%hamsayegie 8taei'ash ra hesab kardim va ruye kole tasvir miangin
%%gereftim.

mask = [-1, -1, -1; -1, 8, -1; -1, -1, -1]/8;% second derivative mask in pactical
filtered_image = conv2(picture, mask, 'same');%same=central part of the convolution
camera_noise = round( mean( abs(filtered_image(~isnan(filtered_image))) )/2 );


%% *****get init seeds for object and background and implement alpha expansion move for label changing****

[segind, Labels, final_alpha_gcost, seeds] = initseeds(picture, cntr, f, conn,...
                                                            camera_noise, prior_likelihood);
% init_seeds(:,:,slice_number)
%%init_seed(:,:,slice_number) use when matrix of seeds is Entered before in
%%previous run of code and you want to compare with changes in constant or type of qmri or
%%cmri
