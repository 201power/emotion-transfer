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

function [centers emotion imgc w model]=class_img(img,mode,fmethod,cmethod)
    
    load pantone.mat;
    clear pcombine_cmyk;
    
    plot3d=0;
    model = 0;
    
    [num_emotion block_size num_color d]=size(pcombine_rgb);

    switch(mode)
        case 'rgb'
            pcombine = pcombine_rgb;
        case 'lab'
            pcombine = pcombine_lab;
            if ((max(max(img)))<=1)
                img = rgb2lab(img);
            else
                img = rgb2lab(double(img)/255);
            end
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
   
    % method for find centers
    switch(fmethod)
        case 'kmean'
            % fixed weights for kmean
            w=[1;0.3;0.1];
            % kmeans
            % us_centers: unsorted centers
            [us_centers nr cid]= kmeans(reshape(double(img),m*n,3),3);
            imgc = cid;
            if (plot3d)
                scatterd(reshape(im2double(img),m*n,3)',cid);
            end
            % put the biggest cluster in first vector in centers
            [tmp,I]=sort(nr,'descend');
            centers=us_centers;
        case 'em'
            x=reshape(im2double(img),m*n,3)';
            [imgc model]=emgm(x,3);
                
            if (plot3d)
                scatterd(x,imgc);
            end
            imgc=imgc';
            centers=model.mu';
            [w,I]=sort(model.weight,'descend');
            w=w';
        case 'kem'
            x=reshape(im2double(img),m*n,3)';
            % kmean before em
            [us_centers nr cid]= kmeans(reshape(double(img),m*n,3),3);
            [imgc model]=emgm(x,us_centers');
            if (plot3d)
                scatterd(x,imgc);
            end
            imgc=imgc';
            centers=model.mu';
            [w,I]=sort(model.weight,'descend');
            w=w';
        case 'ems'
            % em with spatial smooth 
            % kmean to initialize
            [us_centers nr cid]= kmeans(reshape(double(img),m*n,3),3);
            x=reshape(im2double(img),m*n,3)';
            [imgc model]=emgms(x,us_centers',m,n);
            if (plot3d)
                scatterd(x,imgc);
            end
            imgc=imgc';
            centers=model.mu';
            [w,I]=sort(model.weight,'descend');
            w=w';   
    end
    centers=centers(I,:);
    
    imgct=imgc;
    imgct(find(imgc==1))=find(I==1);
    imgct(find(imgc==2))=find(I==2);
    imgct(find(imgc==3))=find(I==3);
    imgc=reshape(imgct,m,n);
    
    % classify image
    kk(1,:,:)=centers;
    kk=repmat(kk,block_size,1);
    kkk(1,:,:,:)=kk;
    kkk=repmat(kkk,num_emotion,1);
    clear kk;

    dist=sum((pcombine-kkk).^2,4);
    score = reshape(dist,num_emotion*block_size,d)*w;
    switch(cmethod)
        case 'single'
            min_score = min(score);
            ind=find(score==min_score);
            ind=ind(1); % if have multiple mins

            emotion(2) = fix((ind-1)/num_emotion)+1;
            emotion(1) = rem(ind-1,num_emotion)+1;
            disp(['Min score equal to: ',num2str(min_score)]);
        case 'group'            
            score = reshape(score,num_emotion,block_size);
            gscore = sum(score,2);
            gmin=min(gscore);
            emotion(1)=find(gmin==gscore);
            smin=min(score(emotion(1),:));
            emotion(2)=find(score(emotion(1),:)==smin);
    end
    disp(['Image is classified as: [',cell2mat(names{emotion(1)}),'] in ',mode]);
    
    switch(mode)
        case 'rgb'
        case 'lab'
            centers=lab2rgb(centers)*255;
        case 'hsv'
            centers=hsv2rgb(centers)*255;
        case 'lch' %problem
            centers=lch2lab(centers);
            centers=lab2rgb(centers)*255;
        case 'ycc'
            centers=ycc2rgb(centers)*255;
    end
end
