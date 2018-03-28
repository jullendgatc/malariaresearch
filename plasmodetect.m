% Program Preprocessing Data untuk Paper Deteksi Parasit Plasmodium
% Falciparumd Dari Citra Mikroskopis Sel Darah Merah 
% Created by Jullend Gatc 
% Faculty of Computer Science, Universitas Indonesia.
% 8 Desember 2012
%==========================================================================
% Program ini menjalankan Prepocessing data pada citra mikroskopis 
% Preprocessing data yaitu : RGBtoGray >> Filtering menggunakan median-
% filter(3x3) >> Histogram.
%==========================================================================
% Cara menjalankan program ini adalah sbb :
% 1. pada command window ketiklah >>medianfilter('nama_file_citra.jpg')
% 2. nama_file_citra disesuaikan dengan nama file citra yang ada pada
% current folder

%==========================================================================
% PREPROCESSING
%==========================================================================
function A = plasmodetect(x)

bacadata = imread(x);                      % membaca citra
G = rgb2gray(bacadata);                    % mengubah citra rgb menjadi citra grayscale
M = medfilt2(G,[3 3]);                     % Median Filter
%==========================================================================
% Menggunakan Metode Otsu
%==========================================================================
I2 = im2uint8(M(:));
N=256;
[count,x]=imhist(I2,N);
figure, bar(x,count);
xlabel('variasi level keabuan');
ylabel('jumlah piksel pada level keabuan');
title('Histrogram');

% menghitung nilai histrogram ternormalisasi
p=(count/sum(count))';
L=length(x);

% menghitung jumlah kumulatif
p1k=cumsum(p);

% menghitung rerata kumulatif kelas
m=cumsum((1:L).*p);

% menghitung rerata intensitas global
mg=sum((1:L).*p);

% menghitung varians antar kelas
varB=(mg*p1k-m).^2./(p1k.*(1-p1k));

% mendapatkan threshold 
val=max(varB);
idx=mean(find(varB==val));
T=(idx-1)/(N-1);

% menghitung separability measure
varG=sum((((1:L)-mg).^2).*p);
sm=varB(T*255)/varG;

% thresholding
Ii=im2bw(M,T);

hasilotsu = 1-Ii; % inverse citra background putih menjadi hitam

%==========================================================================
% PROSES MORPOLOGI (LUBANG CITRA + SEDIKIT BERSIH-BERSIH)
%==========================================================================
fillholes = imfill(hasilotsu,'holes');
%==========================================================================
% Fugure Semua Image
%==========================================================================
% Citra Asli
%figure, imshow(bacadata), title('CITRA ASLI');
%figure, imshow(bacadata), title('Grayscale');
%figure, imshow(hasilotsu), title('Otsu Method');


figure, imshow(fillholes), title('Fill Holes');
figure, imshow(1-fillholes), title('Fill Holes');

%imwrite (fillholes, 'Fill Holes 30.jpg');
%imwrite (hasilotsu, 'Metode otsu 30.jpg');

end