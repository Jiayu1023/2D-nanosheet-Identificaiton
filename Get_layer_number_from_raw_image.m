close all;
clear all;
clc;

xx=imread('F:/Calculation Results/WSe2 270nm Example.bmp'); % Raw image
RGB=double(xx)./255;

SubstrateRGB270=[0.630660265400170 0.282595577548276 0.295851511688653]; % color of 270 nm SiO2/Si substrate, under specified condition
SubstrateXYZ270=[0.1827 0.1272 0.0823];% color of 270 nm SiO2/Si substrate, under specified condition

SRGB=[0.6527 0.5847 0.3488]; % color of substrate, raw image

sizew=size(xx);
w1=sizew(1);
w2=sizew(2);
DeltaE=zeros(w1,w2);
Layernumber_MoS2=zeros(w1,w2);

DeltaRGB=[0.11 0.11 0.11]; % extract substrate

SR=0;
SG=0;
SB=0;
count=0;
RGBII=zeros(w1,w2,3);

for j=w1:-1:1
    for jj=w2:-1:1
        if abs(RGB(j,jj,1)-SRGB(1))<=DeltaRGB(1)&&abs(RGB(j,jj,2)-SRGB(2))<=DeltaRGB(2)&&abs(RGB(j,jj,3)-SRGB(3))<=DeltaRGB(3)
            SR=SR+RGB(j,jj,1);
            SG=SG+RGB(j,jj,2);
            SB=SB+RGB(j,jj,3);
            count=count+1;
            RGB(j,jj)=0;
            RGBII(j,jj)=1;
        end
    end
end

SR=SR./count;
SG=SG./count;
SB=SB./count;

Show=ones(w1,w2,3);
Show(:,:,1)=RGB(:,:,1)+RGBII(:,:,1).*SR;
Show(:,:,2)=RGB(:,:,2)+RGBII(:,:,2).*SG;
Show(:,:,3)=RGB(:,:,3)+RGBII(:,:,3).*SB;

SXYZSum=rgb2xyz([SR SG SB]);

ShowXYZ=rgb2xyz(Show);

ShowXYZ(:,:,1)=(SubstrateXYZ270(1)-SXYZSum(1)).*ones(w1,w2)+ShowXYZ(:,:,1);
ShowXYZ(:,:,2)=(SubstrateXYZ270(2)-SXYZSum(2)).*ones(w1,w2)+ShowXYZ(:,:,2);
ShowXYZ(:,:,3)=(SubstrateXYZ270(3)-SXYZSum(3)).*ones(w1,w2)+ShowXYZ(:,:,3);

Lab=rgb2lab(double(uint8(xyz2rgb(ShowXYZ).*255))./255);
SLab=[42.3339 37.1177 16.0223]; % color of 270 nm SiO2/Si substrate, under specified condition

DeltaE=sqrt((Lab(:,:,1)-SLab(1).*ones(w1,w2)).^2+(Lab(:,:,2)-SLab(2).*ones(w1,w2)).^2+(Lab(:,:,3)-SLab(3).*ones(w1,w2)).^2);
Layernumber_WSe2=(2913068890455131.*DeltaE.^6)./1208925819614629174706176 - (8384126003118983.*DeltaE.^5)./37778931862957161709568 + (3023233943789855.*DeltaE.^4)./590295810358705651712 + (8266497270349265.*DeltaE.^3)./295147905179352825856 + 4083586185560143./633825300114114700748351602688; % WSe2 270 Eye-A

figure(1)
pcolor(Layernumber_WSe2);
shading('interp');
colormap(jet);
caxis([0 4]);
colorbar;
view([0 -90]);
