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

function [timg score]=transfer_img_single(img,dest_emotion,mode)

    load pantone.mat;
    clear pcombine_cmyk;

    [num_emotion block_size num_color d]=size(pcombine_rgb);
    
    switch(mode)
        case 'rgb'
            pcombine = pcombine_rgb;
            img = double(img);
        case 'lab'
            pcombine = pcombine_lab;
            img = rgb2lab(double(img)/255);
        case 'hsv'
            pcombine = reshape(pcombine_rgb/255,num_emotion*block_size*num_color,d);
            pcombine = rgb2hsv(pcombine);
            pcombine = reshape(pcombine,num_emotion,block_size,num_color,d);
            img = rgb2hsv(im2double(img));
        case 'lch'
            pcombine = rgb2lab(pcombine_rgb/255);
            pcombine = lab2lch(pcombine);
            img = rgb2lab(im2double(img));
            img = lab2lch(img);
        case 'ycc'
            pcombine = rgb2ycc(pcombine_rgb/255);
            img = rgb2ycc(im2double(img));
    end
    [m n dim]=size(img);
    
    % select one color combinations from group. 
    disp(['Target mood: [',cell2mat(names{dest_emotion(1)}),']']);
   
    
    % transfer color
    dest_color=pcombine(dest_emotion(1),dest_emotion(2),1,:);
    dest_color=reshape(dest_color,1,3);
    img = reshape(img,m*n,d);
    center = mean(img);
    center(1)=0;dest_color(1)=0;
    timg=img-repmat(center,m*n,1)+repmat(dest_color,m*n,1);
    % Luminace score
    score=sum(abs(timg(:,1)-img(:,1)))/(m*n);
    
    timg = reshape(timg,m,n,d);
    
    switch(mode)
        case 'rgb'
            timg=timg/255;
        case 'lab'
            timg=lab2rgb(timg);
        case 'hsv'
            timg=hsv2rgb(timg);
        case 'lch'
            timg=lch2lab(timg);
            timg=lab2rgb(timg);
        case 'ycc'
            timg=ycc2rgb(timg);
    end
end

