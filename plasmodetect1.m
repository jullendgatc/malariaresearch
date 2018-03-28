% Program Preprocessing Data untuk Paper Deteksi Parasit Plasmodium
% Falciparumd Dari Citra Mikroskopis Sel Darah Merah 
% Created by Jullend Gatc 
% Faculty of Computer Science, Universitas Indonesia.
% 8 Desember 2012 - 2 januari 2013
%==========================================================================
% Program ini menjalankan Prepocessing data pada citra mikroskopis 
% Preprocessing data yaitu : RGBtoGray >> Filtering menggunakan median-
% filter(3x3) >> Histogram.
%==========================================================================
% Cara menjalankan program ini adalah sbb :
% 1. pada command window ketiklah >> medianfilter('nama_file_citra.jpg')
% 2. nama_file_citra disesuaikan dengan nama file citra yang ada pada
% current folder

%==========================================================================
% PREPROCESSING
%==========================================================================
function A = plasmodetect1(x)

bacadata = imread(x);                      % membaca citra
G = rgb2gray(bacadata);                    % mengubah citra rgb menjadi citra grayscale
M = medfilt2(G,[5 5]);                     % Median Filter

%==========================================================================
% Double Threshold
%==========================================================================
figure, imhist (M);
[row,col]=size(M);
T=mean2(M);
err=T;
while err ~=0
    R1=M(M>=T);
    R2=M(M<T);
    miu1=mean2(R1);
    miu2=mean2(R2);
   
    temp=T;
    T=(miu1+miu2)/2;
    err=T-temp;
end
% =========================================================================
for jj=1:row
    for kk=1:col
        if M(jj,kk)>=100
            M(jj,kk)=0;
        else
            M(jj,kk)=255;
        end
    end
end
% =========================================================================
% DETEKSI PARASIT DARI CITRA
% =========================================================================
%    [~, threshold] = edge(M, 'sobel');
%     fudgeFactor = 3.;
%     BWs = edge(M,'sobel', threshold * fudgeFactor);
%     figure, imshow(BWs), title('EDGE DETECTION BINARY MASK');

%     se90 = strel('line', 5, 90);
%     se0 = strel('line', 5, 0);

%     BWsdil = imdilate(M, [se90 se0]);
%     figure, imshow(BWsdil), title('DILATED GRADIENT MASK');

     BWdfill = imfill(M, 'holes');
     figure, imshow(BWdfill);
     title('CITRA BINER DENGAN FILLED HOLES');

%     BWnobord = imclearborder(BWdfill, 8);
%     figure, imshow(BWnobord), title('HAPUS CITRA BORDER');

     seD = strel('diamond',1);
    
     BWfinal = imerode(BWdfill,seD);
     figure, imshow(BWfinal), title('SEGMENTASI CITRA');

     BWoutline = bwperim(BWfinal);
     
    se = strel('disk',10);              % dilasi edge (memperjelas edge)
    BWoutline2 = imdilate(BWoutline,se);

     Segout = bacadata;
     Segout(BWoutline2) = 255;
     
     
     figure, imshow(Segout), title('IDENTIFIKASI PARASIT');   
%==========================================================================
% TAMPILKAN DAN TULIS CITRA
%==========================================================================
% Citra Asli
figure, imshow(bacadata), title('CITRA ASLI');
figure, imshow(bacadata), title('Grayscale');
figure,imshow(BWfinal), title('Hasil Thresholding');
figure, imshow(Segout), title('IDENTIFIKASI PARASIT');

imwrite (Segout, 'DETEKSI PARASIT CITRA 1.jpg');
imwrite (BWfinal, 'PARASIT CITRA 1.jpg');

end % Ending Fungsi