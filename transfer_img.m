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

function [timg coeff score tcenters]=transfer_img(img,imgc,centers,dest_emotion,w,mode)

    load pantone.mat;
    clear pcombine_cmyk;

    [num_emotion block_size num_color d]=size(pcombine_rgb);
    draw_process=0;
    
    switch(mode)
        case 'rgb'
            pcombine = pcombine_rgb;
            img = double(img);
        case 'lab'
            pcombine = pcombine_lab;
            img = rgb2lab(double(img)/255);
            centers = rgb2lab(centers/255);
            color_lb=[0;-100;-100];
            color_hb=[100;100;100];
        case 'hsv'
            pcombine = reshape(pcombine_rgb/255,num_emotion*block_size*num_color,d);
            pcombine = rgb2hsv(pcombine);
            pcombine = reshape(pcombine,num_emotion,block_size,num_color,d);
            img = rgb2hsv(im2double(img));
            centers = rgb2hsv(centers);
        case 'lch'
            pcombine = rgb2lab(pcombine_rgb/255);
            pcombine = lab2lch(pcombine);
            img = rgb2lab(im2double(img));
            img = lab2lch(img);
            centers = rgb2lab(centers/255);
            centers = lab2lch(centers);
        case 'ycc'
            pcombine = rgb2ycc(pcombine_rgb/255);
            img = rgb2ycc(im2double(img));
            centers = rgb2ycc(centers/255);
            color_lb=[0;0;0];
            color_hb=[1;1;1];
    end
    [m n dim]=size(img);
    
    % select one color combinations from group. 
    disp(['Target mood: [',cell2mat(names{dest_emotion(1)}),']']);
    centers=reshape(centers,1,num_color,d);
    
    % transfer color
    dest_color=pcombine(dest_emotion(1),dest_emotion(2),:,:);
    centers=reshape(centers,num_color,d);
    dest_color=reshape(dest_color,num_color,d);      
    img = reshape(img,m*n,d);
    imgc = reshape(imgc,m*n,1);
    c1=img(find(imgc==1),:);
    c2=img(find(imgc==2),:);
    c3=img(find(imgc==3),:);
    
    % each dimension separately
    c=centers';
    t=dest_color';
     % find cluster min / max
     if isempty(c3)
         c3=[color_lb';color_hb'];
     end
    cmin = [min(c1)';min(c2)';min(c3)'];
    cmax = [max(c1)';max(c2)';max(c3)'];
    lb=repmat(color_lb,3,1)-cmin;
    hb=repmat(color_hb,3,1)-cmax;
    % reshape to vector
    c=c(:);
    t=t(:);
    w=sqrt(repmat(w,3,1));
    f=@(x)wdist(x,c,t,w);
    % find new centers
    [delta fval]=fmincon(f,t-c,[],[],[],[],lb,hb,[],optimset('Algorithm','interior-point')); % final two values:lower bound, higher bound
    delta=reshape(delta,3,3)';           
    
    coeff=delta;
    
    %% color transfer
    % each pixel in same cluster move same distance
    % don't move the luminance channel
    delta(:,1)=0;
    tcenters=centers+delta;
    %delta=reshape(delta,3,3)';
    timg=img;
    for i=1:m*n
       h = imgc(i);
       timg(i,:)=img(i,:)+delta(h,:);
    end

    
    timg = reshape(timg,m,n,d);
    img = reshape(img,m,n,d);
    

    %% Gradient preservation
    fimg=timg;
    % optimization 
    fimg= reshape(fimg,m,n,d);
    lambda=20;
    % foward
    Dx=[-1 0 1;-2 0 2;-1 0 1];
    % calculate output image for each dimension
    disp('Gradient preservation in three channels...');
    for i=1:3
       % enable padding 
       oimg=iter_solve(img(:,:,i),fimg(:,:,i),Dx,lambda);
       timg(:,:,i)=oimg;
    end
    
    score = calc_score(img,timg,centers,dest_color,delta);
    
    switch(mode)
        case 'rgb'
            timg=timg/255;
        case 'lab'
            timg=lab2rgb(timg);
            tcenters=lab2rgb(tcenters);
        case 'hsv'
            timg=hsv2rgb(timg);
            tcenters=hsv2rgb(tcenters);
        case 'lch'
            timg=lch2lab(timg);
            timg=lab2rgb(timg);
        case 'ycc'
            timg=ycc2rgb(timg);
            tcenters=ycc2rgb(tcenters);
    end
end



% x: initial value: target centers-input centers
% c: centers of input image
% t: target centers: from pantone color combinations
% w: sqrt(weights)
function f=wdist(x,c,t,w)
        f=norm(w.*(c+x-t));
end

function [img]=func(img,Dx,lambda)
    Dy=Dx';
    Dxt=flipdim(Dx,2);
    Dyt=flipdim(Dy,1);
    [m n]=size(img);
    
    img=imfilter(imfilter(img,Dx,'symmetric'),Dxt,'symmetric')+imfilter(imfilter(img,Dy,'symmetric'),Dyt,'symmetric');
    img=lambda*img;
end

% solve equation (3) in gradient preserve paper
function [oimg err]=iter_solve(simg,fimg,Dx,lambda)
    % initialize, output = intermidiate    
    img_new = fimg;
    alpha = 0.01;
    alpha_reduce_rate=0.9;
    maxiter = 15000;
    
    precision=1e-6;
    err_tol=200;
    % f+lambda(Dxt*Dx+Dyt*Dy)*s ->constant 
    timg=fimg+func(simg,Dx,lambda);
    iter = 0;
    total_iter=0;
    while (1)
        img_old = img_new;
        f_prime=img_new+func(img_new,Dx,lambda)-timg;
        img_new = img_old - alpha*f_prime;
        iter = iter + 1;
        err(iter)=sum(sum(abs(img_new-img_old)));
        %disp(['iteration ',num2str(iter), ' err=',num2str(err(iter),3)]);
        if (err(iter)<=err_tol)
            break;
        end
        % check if the error grow, if then reduce alpha
        if (iter>10 && mod(iter,10)==0)
            if (checkerr(err,iter))
                alpha=alpha*alpha_reduce_rate;
                %disp(['  !error growth, reduce alpha by ',num2str((1-alpha_reduce_rate)*100),'%...alpha=',num2str(alpha),', Reset gradient preservation...']);
                fprintf('*');
                total_iter=total_iter+iter;
                iter=1; % reset algorithm. 
                img_new=fimg;
                if (alpha < precision)
                    disp('!too small alpha,exit gradient preservation...');
                    oimg=fimg;
                    return;
                end
            end
            if (mod(iter,40)==0)
                fprintf('.');
            end
            if (mod(iter,1000)==0)
                if (iter>maxiter)
                    disp(['  Reach max iteration ',num2str(maxiter),', exit....']);
                    break;
                end
                disp(['  iteration ',num2str(iter),', error=',num2str(err(iter))]);
            end
        end            
    end
    disp(['  Gradient preservation complete...error=',num2str(err(iter)),',total_iter=',num2str(total_iter+iter)]);
    oimg = img_new;    
end

% check if error growth continuoues in last ten iterations. 
function result=checkerr(err,iter)
    result=1;
    for i=iter-10:iter-1
        if (err(i)>err(i+1))
            result=0;
            return;
        end
    end
end

% score for output image compare to input image 
function score=calc_score(simg,oimg,centers,target,delta)
    [m n]=size(simg);
    %weights of lscore and cscore
    w=[0.5 0.5];
    % Luminace score
    lscore=sum(sum(abs(oimg(:,:,1)-simg(:,:,1))))/(m*n);
    % center movement score
    cscore=sum(sum(abs(centers+delta-target)))/9;
    score=[lscore cscore];
end

function n=dir_vector(points)
    % calculate the normal vector (direction) of plane formed by three points
    % Ref: http://www.jtaylor1142001.net/calcjat/Solutions/VPlanes/VP3Pts.htm
    % points: row - point index, col - 3 axis value
    A=points(1,:);B=points(2,:);C=points(3,:);
    AB=B-A;
    AC=C-A;
    % cross product
    n=cross(AB,AC);
end