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
function A = preprocessing(x)
F = imread(x);% membaca citra
F = rgb2gray(F);%mengubah citra rgb menjadi citra grayscale
F = im2double(F);%mengubah presisi citra menjadi duakalinya.
[r,c]=size(F);
for i=1:r
    for j=1:c
        ip=i+1;
        im=i-1;
        jm=j-1;
        jp=j+1;
        if(im<1)
            im=i;
        elseif (ip>r)
            ip=i;
        end
        if(jm<1)
            jm=j;
        elseif (jp>c)
            jp=j;
        end
        W=[F(i,j) F(i,jm) F(i,jp) F(ip,j) F(im,j) F(im,jm) F(ip,jm) F(ip,jp) F(im,jp)];
        m = median(W); %matriks nilai median 
        MF(i,j)=m; %matriks filter spasial median
        
       end
end
% =========================================================================
% Konversi image tipe data double ke image tipe data uint
%I = im2uint8(MF);
% =========================================================================

% Mengubah citra ke biner
%biner(I);

% Menampilkan citra asli
figure
imshow(F);
title('CITRA ASLI');

%figure
%imshow(MF);
%title('CITRA MEDIAN FILTER');
%saveas(MF,'plasmo3filter.jpg');

% Menampilkan citra median filter
%figure
%imshow(I);
%title('MEDIAN FILTER');
%A = I;

% Menampilkan citra setelah proses morfological
% thinimage = bwmorph(I,'skel',Inf);
% figure, imshow(thinimage), title('THINNING IMAGE');
% =========================================================================
   % Mengubah Citra menjadi histogram.
   %GambarHisto = imhist(I); 
   %ch = cumsum(GambarHisto); % Hitung nilai komulatif dari histogram.
   %bar(ch,'BarWidth', 1)
   
   %% image equalization dan komulatif histogram.
   %figure; 
   %J = histeq(I);
   %imshow(J), title('CITRA HISTOGRAM');
   %figure; 
   %histogramEq = imhist(J);  % Image histogram.
   %ch = cumsum(histogramEq); % Hitung nilai komulatif dari histogram.
   %bar(ch,'BarWidth',1)      % Draw bar. 
% =========================================================================
% SEGMENTASI CITRA (PAKE EDGE DETECTION)
% =========================================================================
% Setelah diperoleh gambar dari hasil histogram sebelumnya citra akan
% disegmentasi untuk membedakan sitoplasma dan nukleus yang sudah terkena
% parasit plasmodium falciparum.

% Segmentasi Citra Menggunakan Edge Detection
%[~, threshold] = edge(J, 'sobel');
%fudgeFactor = .4;
%BWs = edge(I,'sobel', threshold * fudgeFactor);
%figure, imshow(BWs), title('EDGE DETECTION BINARY MASK');

%se90 = strel('line', 3, 90);
%se0 = strel('line', 3, 0);

%BWsdil = imdilate(BWs, [se90 se0]);
% figure, imshow(BWsdil), title('DILATED GRADIENT MASK');

%BWdfill = imfill(BWsdil, 'holes');
%figure, imshow(BWdfill);
% title('CITRA BINER DENGAN FILLED HOLES');

%BWnobord = imclearborder(BWdfill, 4);
% figure, imshow(BWnobord), title('HAPUS CITRA BORDER');

%seD = strel('diamond',2);
%BWfinal = imerode(BWnobord,seD);
%BWfinal = imerode(BWfinal,seD);
%figure, imshow(BWfinal), title('SEGMENTASI CITRA');

%BWoutline = bwperim(BWfinal);
%Segout = I;
%Segout(BWoutline) = 255;
%figure, imshow(Segout), title('EDGE DETECTION DARI CITRA ASLI');
%figure, imhist(Segout), title('HISTOGRAM SEGMENTASI');

% =========================================================================
% Image = imresize(I,0.5); %Crop citra menjadi setengahnya
% [m, p] = size(Image);
% L = 256;

% axis([0 L-1 0 L-1]);
% hold on;
% xlabel({'r';'Input Intensity Level'}, 'fontweight','b');
% ylabel({'Output Intensity Level' ; 's = T(r)'}, 'fontweight','b');
% title('Please Enter Two Points as (r1, s1) > (r2, s2)', 'fontweight','b');

% n = 1;
% xy(:, n) = [0, 0];
% for i = 1 : 2
%     [a, b] = ginput(1);
%     r(i) = a;
%     s(i) = b;
%     n = n + 1;
%     xy(:, n) = [a, b];
%     plot(a, b, 'r*');
% end

% xy(:, n + 1)= [255, 255];
% t = 1 : n + 1;
% ts = 1 : n + 1;
% xys = spline(t,xy,ts);
% plot(xys(1,:),xys(2,:), 'g-');
% pause(0.2);
% syms x;

% Defining a Function from (0, 0) to (r(1), s(1)):
% m1 = s(1) / r(1);
% f1 = m1 * x;

% Defining a Function from (r(1), s(1)) to (r(2), s(2)):
% m2 = (s(1) - s(2)) / (r(1) - r(2));
% f2 = m2 * x - m2 * r(1) + s(1);

% Defining a Function from (r(2), s(2)) to (L-1, L-1):
% m3 = (s(2) - (L-1)) / (r(2) - (L-1));
% f3 = m3 * x - m3 * r(2) + s(2);

% Computing Intensity Levels of New Image:
% NewImage = zeros(m, p);
% for i = 1 : m
%     for j = 1 : p
%       if 0 <= Image(i, j) <= r(r)  -->  s=T(Image(i, j))= f1(Image(i, j)):
%         if (Image(i, j) >= 0) && (Image(i, j) <= r(1))
%             T = subs(f1, Image(i, j));
            
%       if r(1) <= Image(i, j) <= r(2) --> s=T(Image(i, j))=f2(Image(i, j)):      
%         elseif (Image(i, j) >= r(1)) && (Image(i, j) <= r(2))
%             T = subs(f2, Image(i, j));
            
%       if r(2) <= Image(i, j) <= (L-1) --> s=T(Image(i, j))=f3(Image(i, j)): 
%         else
%             T = subs(f3, Image(i, j));
%         end
%         NewImage(i, j) = T;
%     end
% end

% Plot Original Image & Contrast Stretched Image:
% hold off;
% subplot(1,2,1); imshow(Image)
% title('Citra Asli', 'fontsize', 14, 'fontweight', 'b')
% subplot(1,2,2); imshow(NewImage,[])
% title('Contrast Stretched', 'fontsize' , 14, 'fontweight', 'b')
%==========================================================================
    % Segmentasi Citra Menggunakan Edge Detection
%==========================================================================
% [~, threshold] = edge(NewImage, 'sobel');
%     fudgeFactor = .;
%     BWs = edge(NewImage,'sobel', threshold * fudgeFactor);
%     figure, imshow(BWs), title('EDGE DETECTION BINARY MASK');

%     se90 = strel('line', 3, 90);
%     se0 = strel('line', 3, 0);

%     BWsdil = imdilate(BWs, [se90 se0]);
%     figure, imshow(BWsdil), title('DILATED GRADIENT MASK');

%     BWdfill = imfill(BWsdil, 'holes');
%     figure, imshow(BWdfill);
%     title('CITRA BINER DENGAN FILLED HOLES');

%     BWnobord = imclearborder(BWdfill, 4);
%     figure, imshow(BWnobord), title('HAPUS CITRA BORDER');

%     seD = strel('diamond',1);
%     BWfinal = imerode(BWnobord,seD);
%     BWfinal = imerode(BWfinal,seD);
%     figure, imshow(BWfinal), title('SEGMENTASI CITRA');

%     BWoutline = bwperim(BWfinal);
%     Segout = Image;
%     Segout(BWoutline) = 255;
%     figure, imshow(Segout), title('EDGE DETECTION DARI CITRA ASLI');
%     figure, imhist(Segout), title('HISTOGRAM SEGMENTASI');
%     figure, imhist(F), title('Histogram Citra Asli');


% =========================================================================
% Double THRESHOLDING (
% =========================================================================
figure, imhist (F);
[row,col]=size(F);
T=mean2(F);
err=T;
while err ~=0
    R1=F(F>=T);
    R2=F(F<T);
    miu1=mean2(R1);
    miu2=mean2(R2);
    
    temp=T;
    T=(miu1+miu2)/2;
    err=T-temp;
end
% =========================================================================
for jj=1:row
    for kk=1:col
        if F(jj,kk)>=75
            F(jj,kk)=0;
        else
            F(jj,kk)=255;
        end
    end
end
figure, imshow(F);
title('Image Thresholding','fontsize', 14);

end