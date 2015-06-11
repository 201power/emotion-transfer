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

function [img]=combination_patch(colors,imgsize)
    shrink=0.8; % how much each shrink
    % test dimension of input colors
    if numel(size(colors))>2
        colors=reshape(colors,3,3);
    end
    psize=[imgsize;fix(imgsize*shrink);fix(imgsize*shrink*shrink)];
    img=ones(sum(psize),imgsize,3)/1.25;
    rindex=1;
    for i=1:3
        for j=1:psize(i)
            s=(imgsize-psize(i))/2+1;
            for k=s:s+psize(i)-1
                img(rindex,k,:)=colors(i,:)/255;
            end
            rindex=rindex+1;
        end
    end
    %imwrite(img,'./result/color-comb.jpg','jpg','Quality',100);
end