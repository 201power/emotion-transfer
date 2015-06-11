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

% Table
%    1       2         3      4         5             6          7
%'serene' 'earthy' 'mellow' 'Muted' 'Capricious' 'Spiritual' 'Romantic'
%    8        9         10        11        12        13         14
%'Sensual';'Powerful';'Elegant';'Robust';'Delicate';'Playful';'Energetic';
%    15           16       17        18         19      20         21
%'Traditional';'Classic';'Festive';'Fanciful';'Cool';'Warm';'Luscious&sweet'
%    22           23        24          25         26         27
%'Spicy.Tangy';'Unique 1';'Unique 2';'Unique 3';'Unique 4';'Naturals';

%% emotion transfer without reference image
% only output one image
% [4 7] 4 - muted, 7 - 7th color combination in muted
% output image the \result folder
dest_emotion=[4 7];
emotion_transfer('./photo/p1.jpg',dest_emotion);

% select best output from 24 output images 
% (could be slow!)
% output image the \result folder
dest_emotion=4;
emotion_transfer('./photo/p1.jpg',dest_emotion);

%% emotion transfer with with reference
% output image the \result folder
% p9: input, p10: reference
emotion_transfer('./photo/p9.jpg',[],'./photo/p10.jpg');
