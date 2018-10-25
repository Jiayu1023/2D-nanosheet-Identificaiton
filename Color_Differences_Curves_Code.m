close all;
clc;
clear all;

Load=load('H:\Data\Final_nk_4121_20180725.txt'); % n and k of materials
lambda=Load(:,1);
MoS2_n=Load(:,2);
MoS2_k=Load(:,3);
WS2_n=Load(:,4);
WS2_k=Load(:,5);
MoSe2_n=Load(:,6);
MoSe2_k=Load(:,7);
WSe2_n=Load(:,8);
WSe2_k=Load(:,9);
Air_n=Load(:,10);
Air_k=Load(:,11);
Si_n=Load(:,12);
Si_k=Load(:,13);
SiO2_n=Load(:,14);
SiO2_k=Load(:,15);

SiO2_d=300;
Layer=10;
Layernumber=0:1:Layer;
Layernumber=Layernumber';

WS2_d=zeros(Layer,1);
MoS2_d=zeros(Layer,1);
WSe2_d=zeros(Layer,1);
MoSe2_d=zeros(Layer,1);

WS2_d=1.0+0.7.*(Layernumber-1);
MoS2_d=1.0+0.7.*(Layernumber-1);
WSe2_d=1.0+0.7.*(Layernumber-1);
MoSe2_d=1.0+0.7.*(Layernumber-1);

WS2_d(1)=0;
MoS2_d(1)=0;
WSe2_d(1)=0;
MoSe2_d(1)=0;

Air_d=0.3; %10.1038/ncomms5966  10.1021/nl3026357 air layer between substrate and nanosheet

gamma=2.2;

Air_ni=Air_n+i.*Air_k;
WS2_ni=WS2_n+i.*WS2_k;
MoS2_ni=MoS2_n+i.*MoS2_k;
WSe2_ni=WSe2_n+i.*WSe2_k;
MoSe2_ni=MoSe2_n+i.*MoSe2_k;
SiO2_ni=SiO2_n+i.*SiO2_k;
Si_ni=Si_n+i.*Si_k;

n4=Air_ni;
n3=MoS2_ni;
n2=Air_ni;
n1=SiO2_ni;
n0=Si_ni;

Con=zeros(4121,Layer+1);
for yg=1:1:Layer+1
    r43=(n4-n3)./(n4+n3);
    r32=(n3-n2)./(n3+n2);
    r21=(n2-n1)./(n2+n1);
    r10=(n1-n0)./(n1+n0);
    beta1=2.*pi.*n3.*MoS2_d(yg)./lambda;
    MoS2_r=(r43+r32.*exp(2.*i.*beta1))./(1+r43.*r32.*exp(2.*i.*beta1));
    
    beta2=2.*pi.*n2.*Air_d./lambda;
    Air_r=(MoS2_r+r21.*exp(2.*i.*beta2))./(1+MoS2_r.*r21.*exp(2.*i.*beta2));
    
    beta3=2.*pi.*n1.*SiO2_d./lambda;
    SiO2_r=(Air_r+r10.*exp(2.*i.*beta3))./(1+Air_r.*r10.*exp(2.*i.*beta3));
    
    R_SiO2=(abs(SiO2_r)).^2;
    
    Con(:,yg)=R_SiO2;
end

Con_2=zeros(4121,Layer+1);
for yg_2=1:1:Layer+1
    r21_2=(Air_ni-SiO2_ni)./(Air_ni+SiO2_ni);
    r10_2=(SiO2_ni-Si_ni)./(SiO2_ni+Si_ni);
    beta1_2=2.*pi.*SiO2_ni.*SiO2_d./lambda;
    SiO2_r_2=(r21_2+r10_2.*exp(2.*i.*beta1_2))./(1+r21_2.*r10_2.*exp(2.*i.*beta1_2));
    
    R_2=(abs(SiO2_r_2)).^2;
    
    Con_2(:,yg_2)=R_2;
end

Contrast=zeros(4121,Layer+1);
for yg_3=1:1:Layer+1
    Contrast(:,yg_3)=(Con_2(:,yg_3)-Con(:,yg_3))./Con_2(:,yg_3);
end

figure(1)
surf(Layernumber,lambda,Contrast);
shading('interp'),colormap(jet),colorbar;
view(0,90);
axis([0 Layer 388 800]);
xlabel('Layer number'),ylabel('Wavelength (nm)');
set(gca,'tickdir','out');

D65=load('H:\Data\A_interp_388_800.txt'); % Illuminant data
I_lambda=D65(:,2);

CIE=load('H:\Data\Camera xyz Interp 388_800.txt'); % Color matching functions data

CIEINDX=1;
x_=CIE(3.*(CIEINDX)-2,:)';
y_=CIE(3.*(CIEINDX)-1,:)';
z_=CIE(3.*(CIEINDX),:)';

R_lambda=Con;

R_read=zeros(Layer+1,1);
R_linear=zeros(Layer+1,1);
G_read=zeros(Layer+1,1);
G_linear=zeros(Layer+1,1);
B_read=zeros(Layer+1,1);
B_linear=zeros(Layer+1,1);

Y_k=0;
sum_Y_k=0;
for gail_Y_k=1:1:4121
    Y_k=(I_lambda(gail_Y_k).*y_(gail_Y_k).*0.1);
    sum_Y_k=sum_Y_k+Y_k;
    Y_k=0;
end
k=100./sum_Y_k;

X=zeros(Layer+1,1);
sum_X=zeros(Layer+1,1);
Y=zeros(Layer+1,1);
sum_Y=zeros(Layer+1,1);
Z=zeros(Layer+1,1);
sum_Z=zeros(Layer+1,1);
x=zeros(Layer+1,1);
y=zeros(Layer+1,1);
z=zeros(Layer+1,1);


for jy=1:1:Layer+1
    
    
    for gail_X=1:1:4121
        X(jy)=(I_lambda(gail_X).*R_lambda(gail_X,jy).*x_(gail_X).*0.1);
        sum_X(jy)=sum_X(jy)+X(jy);
        X(jy)=0;
    end
    
    for gail_Y=1:1:4121
        Y(jy)=(I_lambda(gail_Y).*R_lambda(gail_Y,jy).*y_(gail_Y).*0.1);
        sum_Y(jy)=sum_Y(jy)+Y(jy);
        Y(jy)=0;
    end
    
    for gail_Z=1:1:4121
        Z(jy)=(I_lambda(gail_Z).*R_lambda(gail_Z,jy).*z_(gail_Z).*0.1);
        sum_Z(jy)=sum_Z(jy)+Z(jy);
        Z(jy)=0;
    end
    
    sum_X(jy)=k.*sum_X(jy);
    sum_Y(jy)=k.*sum_Y(jy);
    sum_Z(jy)=k.*sum_Z(jy);
    
    x(jy)=sum_X(jy)./(sum_X(jy)+sum_Y(jy)+sum_Z(jy));
    y(jy)=sum_Y(jy)./(sum_X(jy)+sum_Y(jy)+sum_Z(jy));
    z(jy)=sum_Z(jy)./(sum_X(jy)+sum_Y(jy)+sum_Z(jy));
    
    Matrix_XYZ=zeros(3,1);
    Matrix_RGB=zeros(3,1);
    Matrix_XYZ=[sum_X(jy);sum_Y(jy);sum_Z(jy)];
    Matrix_Constant=[3.241 -1.5374 -0.4986;-0.9692 1.876 0.0416;0.0556 -0.204 1.057];
    Matrix_RGB=(Matrix_Constant)*(Matrix_XYZ);
    R_linear(jy)=Matrix_RGB(1)./100;
    G_linear(jy)=Matrix_RGB(2)./100;
    B_linear(jy)=Matrix_RGB(3)./100;
    
    if R_linear(jy)<=0.00304&R_linear(jy)>0
        R_read(jy)=1.*(12.92.*R_linear(jy));
    elseif R_linear(jy)>1
        R_read(jy)=1;
    elseif R_linear(jy)<=0
        R_read(jy)=0;
    else
        R_read(jy)=1.*((1+0.055).*(R_linear(jy).^(1./gamma))-0.055);
    end
    
    if G_linear(jy)<=0.00304&G_linear(jy)>0
        G_read(jy)=1.*(12.92.*G_linear(jy));
    elseif G_linear(jy)>1
        G_read(jy)=1;
    elseif G_linear(jy)<=0
        G_read(jy)=0;
    else
        G_read(jy)=1.*((1+0.055).*(G_linear(jy).^(1./gamma))-0.055);
    end
    
    if B_linear(jy)<=0.00304&B_linear(jy)>0
        B_read(jy)=1.*(12.92.*B_linear(jy));
    elseif B_linear(jy)>1
        B_read(jy)=1;
    elseif B_linear(jy)<=0
        B_read(jy)=0;
    else
        B_read(jy)=1.*((1+0.055).*(B_linear(jy).^(1./gamma))-0.055);
    end
    
    
end

R=R_read;
G=G_read;
B=B_read;

R_read=round(255.*R_read);
G_read=round(255.*G_read);
B_read=round(255.*B_read);

figure(2)
la=0:1:Layer;
la=la';
RGB=zeros(Layer+1,3);
for ljy=1:1:Layer+1
    plot(la(ljy),0,'s','markersize',40,'color',[R(ljy) G(ljy) B(ljy)],'MarkerFaceColor',[R(ljy) G(ljy) B(ljy)]);
    RGB(ljy,1)=R(ljy);
    RGB(ljy,2)=G(ljy);
    RGB(ljy,3)=B(ljy);
    hold on;
end
set(gca,'color',[0 0 0]);
xlabel('Layer number');
axis([0 Layer -0.02 0.02]);
set(gca,'tickdir','out');
box off;
grid off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lab=zeros(Layer+1,3);
for ljy=1:1:Layer+1
    Lab(ljy,:)=rgb2lab([R(ljy) G(ljy) B(ljy)]);
end
L=Lab(:,1);
a=Lab(:,2);
b=Lab(:,3);
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lab_2=zeros(Layer+1,3);
L_2=zeros(Layer+1,1);
a_2=zeros(Layer+1,1);
b_2=zeros(Layer+1,1);

R_lambda=Con_2;

R_read=zeros(Layer+1,1);
R_linear=zeros(Layer+1,1);
G_read=zeros(Layer+1,1);
G_linear=zeros(Layer+1,1);
B_read=zeros(Layer+1,1);
B_linear=zeros(Layer+1,1);

Y_k=0;
sum_Y_k=0;
for gail_Y_k=1:1:4121
    Y_k=(I_lambda(gail_Y_k).*y_(gail_Y_k).*0.1);
    sum_Y_k=sum_Y_k+Y_k;
    Y_k=0;
end
k=100./sum_Y_k;

X=zeros(Layer+1,1);
sum_X=zeros(Layer+1,1);
Y=zeros(Layer+1,1);
sum_Y=zeros(Layer+1,1);
Z=zeros(Layer+1,1);
sum_Z=zeros(Layer+1,1);
x=zeros(Layer+1,1);
y=zeros(Layer+1,1);
z=zeros(Layer+1,1);


for jy=1:1:Layer+1
    
    
    for gail_X=1:1:4121
        X(jy)=(I_lambda(gail_X).*R_lambda(gail_X,jy).*x_(gail_X).*0.1);
        sum_X(jy)=sum_X(jy)+X(jy);
        X(jy)=0;
    end
    
    for gail_Y=1:1:4121
        Y(jy)=(I_lambda(gail_Y).*R_lambda(gail_Y,jy).*y_(gail_Y).*0.1);
        sum_Y(jy)=sum_Y(jy)+Y(jy);
        Y(jy)=0;
    end
    
    for gail_Z=1:1:4121
        Z(jy)=(I_lambda(gail_Z).*R_lambda(gail_Z,jy).*z_(gail_Z).*0.1);
        sum_Z(jy)=sum_Z(jy)+Z(jy);
        Z(jy)=0;
    end
    
    sum_X(jy)=k.*sum_X(jy);
    sum_Y(jy)=k.*sum_Y(jy);
    sum_Z(jy)=k.*sum_Z(jy);
    
    x(jy)=sum_X(jy)./(sum_X(jy)+sum_Y(jy)+sum_Z(jy));
    y(jy)=sum_Y(jy)./(sum_X(jy)+sum_Y(jy)+sum_Z(jy));
    z(jy)=sum_Z(jy)./(sum_X(jy)+sum_Y(jy)+sum_Z(jy));
    
    Matrix_XYZ=zeros(3,1);
    Matrix_RGB=zeros(3,1);
    Matrix_XYZ=[sum_X(jy);sum_Y(jy);sum_Z(jy)];
    Matrix_Constant=[3.241 -1.5374 -0.4986;-0.9692 1.876 0.0416;0.0556 -0.204 1.057];
    Matrix_RGB=(Matrix_Constant)*(Matrix_XYZ);
    R_linear(jy)=Matrix_RGB(1)./100;
    G_linear(jy)=Matrix_RGB(2)./100;
    B_linear(jy)=Matrix_RGB(3)./100;
    
    if R_linear(jy)<=0.00304&R_linear(jy)>0
        R_read(jy)=1.*(12.92.*R_linear(jy));
    elseif R_linear(jy)>1
        R_read(jy)=1;
    elseif R_linear(jy)<=0
        R_read(jy)=0;
    else
        R_read(jy)=1.*((1+0.055).*(R_linear(jy).^(1./gamma))-0.055);
    end
    
    if G_linear(jy)<=0.00304&G_linear(jy)>0
        G_read(jy)=1.*(12.92.*G_linear(jy));
    elseif G_linear(jy)>1
        G_read(jy)=1;
    elseif G_linear(jy)<=0
        G_read(jy)=0;
    else
        G_read(jy)=1.*((1+0.055).*(G_linear(jy).^(1./gamma))-0.055);
    end
    
    if B_linear(jy)<=0.00304&B_linear(jy)>0
        B_read(jy)=1.*(12.92.*B_linear(jy));
    elseif B_linear(jy)>1
        B_read(jy)=1;
    elseif B_linear(jy)<=0
        B_read(jy)=0;
    else
        B_read(jy)=1.*((1+0.055).*(B_linear(jy).^(1./gamma))-0.055);
    end
    
end

R=R_read;
G=G_read;
B=B_read;

R_read=round(255.*R_read);
G_read=round(255.*G_read);
B_read=round(255.*B_read);

for ljy=1:1:Layer+1
    Lab_2(ljy,:)=rgb2lab([R(ljy) G(ljy) B(ljy)]);
end

L_2=Lab_2(:,1);
a_2=Lab_2(:,2);
b_2=Lab_2(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_E=zeros(Layer+1,1);
for lalala=1:1:Layer+1
    delta_E(lalala)=sqrt((L(lalala)-L_2(lalala)).^2+(a(lalala)-a_2(lalala)).^2+(b(lalala)-b_2(lalala)).^2);
end

figure(3)
plot(Layernumber,delta_E,'-s','markerfacecolor', [0 0 0],'markersize',7,'color','k','linewidth',1.5);
xlabel('Layer number'),ylabel('Color difference');

delta_E=delta_E';
TransLab=(rgb2lab(RGB));

