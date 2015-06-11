% =========================================================================
% Example code for emotion transfer
% 
% Li He
% EECS, University of Tennessee, Knoxville
% 
% Paper
% Li He, Hairong Qi, Russell Zaretzki, 
% "Image color transfer to evoke different emotions based on color combinations", 
% Signal, Image and Video Processing, Aug 2014
% contact: lihe.now@gmail.com
% =========================================================================

function []=emotion_transfer(imgFile,dest_emotion,refImgFile)

    load pantone.mat;
    clear pcombine_cmyk;

    addpath(genpath('./gmm'))
    addpath(genpath('./optprop'))
    
    if exist('refImgFile','var')
        ref_img = 1;
    else
        ref_img = 0;
    end

    % mode = 'rgb' / 'lab' / 'hsv' / 'lch' /'ycc'
    mode = 'lab';
    im_name=findname(imgFile);
    img=imread(imgFile);
    if (ref_img)
        ref_name=findname(refImgFile);
        rimg=imread(refImgFile);
    end

    disp(['====== Classify ',imgFile,' in ',mode,' space ======']);
    [centers emotion imgc w model]=class_img(img,mode,'em','group');
    if (ref_img)
        disp(['====== Classify ',refImgFile,' in ',mode,' space ======']);
        [rcenters remotion rimgc rw rmodel]=class_img(rimg,mode,'em','group');
        dest_emotion=remotion;
    end

    patch_size=50;
    figure,subplot(2,3,1);
    imshow(img);    
    title('test image');
    subplot(2,3,2);
    [m n d]=size(img);
    cimg=reshape(imgc,m*n,1);
    cimg(find(cimg==1))=0;
    cimg(find(cimg==2))=0.5;
    cimg(find(cimg==3))=1;
    imshow(reshape(cimg,m,n));
    title(['Segmentation result:main->black']);
    subplot(2,3,3);        
    imshow(combination_patch(centers,patch_size));
    title(['Three main colors : ',cell2mat(names{emotion(1)})]);

    % transfer color to another emotion
    if (~ref_img)
        if numel(dest_emotion)==2
            [timg,coeff,score,tcenters]=transfer_img(img,imgc,centers,dest_emotion,w,mode);
        else
            disp('Select best image from 24 images');
            for i=1:24
                dest_emotion(2)=i;
                [timg,coeff,score(i,:),allcenters{i}]=transfer_img(img,imgc,centers,dest_emotion,w,mode);
                oimg(i,:,:,:)=timg;
            end
            fscore=score*[0.7;0.3];
            [~,I]=sort(fscore,'ascend');
            % select the one with highest score
            dest_emotion(2)=I(1);
            tcenters=allcenters{I(1)};

            timg = reshape(oimg(dest_emotion(2),:,:,:),m,n,d);
        end
    
        imwrite(timg,['./result/',im_name,'_',num2str(dest_emotion(1)),...
            '_',num2str(dest_emotion(2)),'.jpg'],'jpg','Quality',100);
    else
        [timg,coeff,score,tcenters]=transfer_img(img,imgc,centers,dest_emotion,w,mode);
        imwrite(timg,['./result/',im_name,'_',ref_name,'.jpg'],'jpg','Quality',100);
    end


    subplot(2,3,4);
    pimg=combination_patch(pcombine_rgb(dest_emotion(1),dest_emotion(2),:,:),patch_size);
    pimg=[pimg combination_patch(tcenters*255,patch_size)];
    imshow(pimg);
    title([{'transfered color combinations'};{'Pantone,Target'}]);
    
    subplot(2,3,5);
    imshow(timg);
    title(['transfered image: ' cell2mat(names{dest_emotion(1)})]);
    subplot(2,3,6);
    timg_single=transfer_img_single(img,dest_emotion,mode);
    imshow(timg_single);
    title('transfered image: single color');

end