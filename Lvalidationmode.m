 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Written By Azade Jafari											%%				  									    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validation Criterion

% clc
% close all
% clear all

%read segmented image by my code
dummy1_cMRI4= load_untouch_nii('LefPutamenHighresBin90.nii');%load the segmented region by my code
mysegment_cMRI4=dummy1_cMRI4.img;
mysegment_cMRI4=flip(mysegment_cMRI4,1);
mysegment_cMRI4=imrotate(mysegment_cMRI4,-90);
i=103;%i=zaribe slice haye joda shode+1
figure
imshow(mysegment_cMRI4(:,:,i),[]);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);%show full screen size
title('Segmented Binary region mask by me of cMRI 4','Color', 'm');



%read Gold standard image
dummy2= load_untouch_nii('LP_GOLDSTANDARD0to197.nii');%gold standard binary image which is obtained by MIPAV
GOLDstandard=dummy2.img;
GOLDstandard=flip(GOLDstandard,1);
GOLDstandard=imrotate(GOLDstandard,-90);
slice_number1=i+49;%havasat bashe ke matlab az 1 shoru mishe va to bayad slice number asli ra+1 koni
figure
imshow(GOLDstandard(:,:,slice_number1),[]);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);%show full screen size
title('Gold standard Binary region mask','Color', 'm');

%%
%DICE
%DICE SIMILARITY COEFFICIENT form1 based on BABALOLA

%slice by slice
k=1;
j=137;
Dice_coefficient_cMRI4=zeros(1,29);
for i=88:116
       common = nnz(mysegment_cMRI4(:,:,i) & GOLDstandard(:,:,j));%eshterak 
       joint = nnz(mysegment_cMRI4(:,:,i)|GOLDstandard(:,:,j));%ejtema
a = sum(common(:));
b = sum(joint(:));
Dice_coefficient_cMRI4(1,k) = (2*a/(a+b))*100; 
k=k+1;           
j=j+1;
end
mean_DSC_cMRI4=mean(Dice_coefficient_cMRI4);
% %CSF of whole slice
% common = nnz(mysegment(:,:,[80,108]) & GOLDstandard(:,:,[137,165]));%eshterak 
% joint = nnz(mysegment(:,:,[80,108])|GOLDstandard(:,:,[137,165]));%ejtema
% a = sum(common(:));
% b = sum(joint(:));
% Dice_coefficient_whole = (2*a/(a+b))*100;

%% show result of DSC as box plot
figure
boxplot(Dice_coefficient_cMRI4,'colors',[1 0 1])
xlabel('DICE Similarity Coefficient(DSC)of Left Putamen of cMRI 4')
ylabel('Percentage')
title('EVALUATION')

figure
plot(Dice_coefficient_cMRI4,'-gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('DICE Similarity Coefficient(DSC)of Left Putamen of cMRI 4')
ylabel('Percentage')
title('EVALUATION')

% %DICE SIMILARITY COEFFICIENTform2
% common = (mysegment(:,:,[80,108]) & GOLDstandard(:,:,[137,165]));
%  a = nnz(common(:));%nnz=number of nonzero elements
%  b = nnz(mysegment(:,:,[80,108]));
%  c = nnz(GOLDstandard(:,:,[137,165]));
%  Dice = 2*a/(b+c);


%%
%AVD

area1=zeros(1,29);
k=1;
for i=88:116
% area(1,k) = bwarea(mysegment(:,:,i));
% area(1,k)=regionprops(mysegment(:,:,i),'Area');
area1(1,k)=nnz(mysegment_cMRI4(:,:,i));
k=k+1;
end
numberofnnz_mysegment=sum(area1,2);
% stats = regionprops3(mysegment(:,:,[80,108]),'Volume');

area2=zeros(1,29);
k=1;
for i=137:165
% area(1,k) = bwarea(mysegment(:,:,i));
% area(1,k)=regionprops(mysegment(:,:,i),'Area');
area2(1,k)=nnz(GOLDstandard(:,:,i));
k=k+1;
end
numberofnnz_GOLDSTANDARD=sum(area2,2);
k=1;
AVD_slices_cMRI4=zeros(1,29);
for i=1:29
AVD_slices_cMRI4(1,k)=(abs(area1(1,i)-area2(1,i))/area2(1,i))*100;
k=k+1;
end
mean_AVD_cMRI4=mean(AVD_slices_cMRI4);
% AVD_whole=((abs(numberofnnz_GOLDSTANDARD-numberofnnz_mysegment))/numberofnnz_GOLDSTANDARD)*100;

%% show result of AVD as box plot
figure
boxplot(AVD_slices_cMRI4,'colors',[0 0 1])
xlabel('Absolute Volumetric Difference(AVD)of Left Putamen of cMRI 4')
ylabel('Percentage')
title('EVALUATION')

figure
plot(AVD_slices_cMRI4,'-cs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Absolute Volumetric Difference(AVD)of Left Putamen of cMRI 4')
ylabel('Percentage')
title('EVALUATION')
%% Accuracy,Sensivity,Precision,Specificity

TP=0;
FN=0;
FP=0;
TN=0;
[r,c,z]=size(GOLDstandard);
val_parameter_cMRI4=zeros(29,4);%24 slice darim va 4 parameter TP FN FP TN ke mikhahim baraye har slice jodagane mohasebe shavad
k=1;
for y=88:116
   for i=1:r
      for j=1:c
          if (GOLDstandard(i,j,y+49)~=0 && mysegment_cMRI4(i,j,y)~=0)
              TP=TP+1;
          elseif (GOLDstandard(i,j,y+49)~=0 && mysegment_cMRI4(i,j,y)==0)
              FN=FN+1;
          elseif (GOLDstandard(i,j,y+49)==0 && mysegment_cMRI4(i,j,y)~=0)
              FP=FP+1;
          elseif(GOLDstandard(i,j,y+49)==0 && mysegment_cMRI4(i,j,y)==0)
              TN=TN+1;
          end
      end
    end
%       val_parameter(k,1)=TP;
%       val_parameter(k,2)=FN;
%       val_parameter(k,3)=FP;
%       val_parameter(k,4)=TN;
      val_parameter_cMRI4(k,1)=(TP/(TP+FN))*100;%Sensivity
      val_parameter_cMRI4(k,2)=(TN/(TN+FP))*100;%Specificity
      val_parameter_cMRI4(k,3)=((TP+TN)/(TP+TN+FP+FN))*100;%Accuracy
      val_parameter_cMRI4(k,4)=(TP/(TP+FP))*100;% Precision
      k=k+1;
      TP=0;
      FN=0;
      FP=0;
      TN=0;
      
end

mean_Sensivity_cMRI4=mean(val_parameter_cMRI4(:,1));
mean_Specificity_cMRI4=mean(val_parameter_cMRI4(:,2));
mean_Accuracy_cMRI4=mean(val_parameter_cMRI4(:,3));
mean_Precision_cMRI4=mean(val_parameter_cMRI4(:,4));

%% show result of TP,TN,FP,FN as box plot
figure
boxplot(val_parameter_cMRI4(:,1),'colors',[1 0 1])
xlabel('Sensivity of Left Putamen of cMRI 4')
ylabel('Percentage')
title('EVALUATION')

figure
plot(val_parameter_cMRI4(:,1),'-rs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Sensivity of Left Putamen of cMRI 4')
ylabel('Percentage')
title('EVALUATION')


figure
boxplot(val_parameter_cMRI4(:,2),'colors',[0 0 1])
xlabel('Specificity of Left Putamen of cMRI 4')
ylabel('Percentage')
title('EVALUATION')

figure
plot(val_parameter_cMRI4(:,2),'-ms',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Specificity of Left Putamen of cMRI 4')
ylabel('Percentage')
title('EVALUATION')


figure
boxplot(val_parameter_cMRI4(:,3),'colors',[1 0 1])
xlabel('Accuracy of Left Putamen of cMRI 4')
ylabel('Percentage')
title('EVALUATION')

figure
plot(val_parameter_cMRI4(:,3),'-ys',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Accuracy of Left Putamen of cMRI 4')
ylabel('Percentage')
title('EVALUATION')


figure
boxplot(val_parameter_cMRI4(:,3),'colors',[0 0 1])
xlabel('Precision of Left Putamen of cMRI 4')
ylabel('Percentage')
title('EVALUATION')

figure
plot(val_parameter_cMRI4(:,3),'-ks',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Precision of Left Putamen of cMRI 4')
ylabel('Percentage')
title('EVALUATION')
%%
%% Hausdorff distance
% [ hd ] = computeHD( GOLDstandard,mysegment, 'Euclidean');
k=1;
dH_cMRI4=zeros(1,29);
for i=88:116
% dH_cMRI4(1,k) = hausdorff( mysegment_cMRI4(:,:,i),GOLDstandard(:,:,i+49));
dH_cMRI4(1,k) = hausdorffDist(mysegment_cMRI4(:,:,i),GOLDstandard(:,:,i+49))
k=k+1;
end

mean_Hausdorff_dist_cMRI4=mean(dH_cMRI4);

figure
boxplot(dH_cMRI4,'colors',[1 0 1])
xlabel('Hausdorff Distance (HD)of Left Putamen of cMRI 4')
ylabel('value')
title('EVALUATION')

figure
plot(dH_cMRI4,'-ks',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Hausdorff Distance (HD)of Left Putamen of cMRI 4')
ylabel('value')
title('EVALUATION')

%% Statistical Analysis
[h1,p1,ci1,stats1] =ttest(Dice_coefficient_cMRI4);%if h=1 it means there is a meaningfull diffrence between values of DSC
[h2,p2,ci2,stats2] =ttest(AVD_slices_cMRI4);%if h=1 it means there is a meaningfull diffrence between values of DSC
[h3,p3,ci3,stats3] =ttest(dH_cMRI4);%if h=1 it means there is a meaningfull diffrence between values of DSC
[h4,p4,ci4,stats4] =ttest(val_parameter_cMRI4(:,1));%if h=1 it means there is a meaningfull diffrence between values of DSC
[h5,p5,ci5,stats5] =ttest(val_parameter_cMRI4(:,2));%if h=1 it means there is a meaningfull diffrence between values of DSC
[h6,p6,ci6,stats6] =ttest(val_parameter_cMRI4(:,3));%if h=1 it means there is a meaningfull diffrence between values of DSC
[h7,p7,ci7,stats7] =ttest(val_parameter_cMRI4(:,4));%if h=1 it means there is a meaningfull diffrence between values of DSC
