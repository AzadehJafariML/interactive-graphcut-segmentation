function skull_stripped=skull_stripping_preprocess(I,slice_number)
%%skull stripping by Morphologic OPening (A SIMPLE SKULL STRIPPING ALGORITHM FOR BRAIN MRI by ROY & MAJI 2015)

I=I(:,:,slice_number);
morphologyopen = imopen(I, strel('disk',55));%open the image I with structure SE
%%SE=Create a disk-shaped structuring element with a radius of 55 pixels.

skull_stripped = I;
mask = morphologyopen>0;
skull_stripped(~mask) = NaN;%instead of assigning intensity of these pixels to be zero, by assigning them as NaN, we delete them
figure, imshow(skull_stripped, []),title('skull_Stripped MRI','Color', 'r');
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1])

%%
%%an other way to skull strip
% function skull_stripped=skull_stripping_preprocess(I,slice_number)
% I=I(:,:,slice_number);
% figure
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
% mask = livewire(I);
% livewired_slice = NaN(size(I)); % using NaN instead of zeros, gets rid of using fmold altogether
% livewired_slice(mask) = I(mask);
% figure, imshow(livewired_slice, []), title('livewired slice')
% 
% skull_stripped = livewired_slice;
% boundary = bwboundaries(mask);
% boundary = boundary{1};
% bw = false(size(I));
% for i = 1:length(boundary)
%     bw( boundary(i,1),boundary(i,2) ) = true;
% end
% bw1 = imdilate( bw, strel('disk',1) );
% dummy = skull_stripped(bw1);
% dummy( dummy<800 ) = NaN;
% skull_stripped(bw1) = dummy;
% mask2 = ~isnan(skull_stripped);
% bw2 = bwareaopen(mask2, 2, 4);
% pruned_livewired_slice(~bw2) = NaN;
% figure, imshow(skull_stripped, []),title('skull_Stripped Quantitative MRI','Color', 'r')
% clear i dummy