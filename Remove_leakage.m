clear all
close all

 original =imread('160.tif');
 se = strel('square',7);
 afterOpening = imerode(original,se);
 CC = bwconncomp(afterOpening);
 number  = CC.NumObjects;
 count=1;
 for i=1:inf
     
   if number ~=2
       se1 = strel('square',3);
       afterOpening = imerode(afterOpening,se1);
       CC = bwconncomp(afterOpening);
       number  = CC.NumObjects;
       count=count+1;
   else
     break
%      
   end
 end
  for j=1:count
      se2 = strel('square',3);
     afterOpening = imdilate(afterOpening,se2);
  end
Z = imsubtract(original,afterOpening);

se1 = strel('square',7);
Z = imerode(Z,se1);
se1 = strel('square',13);
Z = imdilate(Z,se1);
Z = imsubtract(original,Z);

figure, imshow(Z);
figure, imshow(original);
figure, imshow(afterOpening,[]);