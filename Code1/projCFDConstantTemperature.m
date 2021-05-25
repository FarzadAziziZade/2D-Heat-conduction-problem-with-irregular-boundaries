clc;
clear all;
close all;
format long
digits(64);
disp('Please choose method that you want to use Jacobi=1  PSOR/Point Gauss Seidel=2 LSOR/Line Gauss Seidel=3.')
choose=input('Enter the number of method = ');
disp('Note that if you use "LSOR/Line Gauss Seidel" or "PSOR/Point Gauss Seidel" you need to set "w" as input.')
X_number=input('number of x points, X_number=');
b2=input('Elongation in x cordinate, X Beta (1 < X Beta < inf ) =');
Y_number=input('number of y points, Y_number=');
b=input('Elongation in y cordinate, Y Beta (1 < Y Beta < inf ) =');
conv=input('Convergence limit =');

for j=1:Y_number
    g_y(j)=(j-1)/(Y_number-1);
end
for i=1:X_number
    g_x(i)=(i-1)/(X_number-1);
end
for i=1:X_number
    x_down(i)=32.*((1+b2)*(((b2+1)/(b2-1))^((g_x(i)-0.5)/0.5))+(1-b2))/((2*(1+(((b2+1)/(b2-1))^((g_x(i)-0.5)/0.5)))));
    x_up(i)=(16*tan((pi/180)*50)+32-64*cot((pi/180)*70)).*((1+b2)*(((b2+1)/(b2-1))^((g_x(i)-0.5)/0.5))+(1-b2))/((2*(1+(((b2+1)/(b2-1))^((g_x(i)-0.5)/0.5)))))-16*tan((pi/180)*50);
end
for i=1:X_number
h(i)=sqrt((x_up(i)-x_down(i))^2+((((16-64)/(-16*tan((pi/180)*50)-32+64*cot((pi/180)*70)))*(x_up(i)+16*tan((pi/180)*50))+16)-sqrt(x_down(i)))^2);
end
for i=1:X_number
    for j=1:Y_number
    y(i,j)=h(i).*((1+b)*(((b+1)/(b-1))^((g_y(j)-0.5)/0.5))+(1-b))/((2*(1+(((b+1)/(b-1))^((g_y(j)-0.5)/0.5)))));
    end
end
for i=1:X_number
% teta(i)=180+acotd((-x_down(i)+x_up(i))/(((((16-64)/(-16*tan((pi/180)*50)-32+64*cot((pi/180)*70)))*(x_up(i)+16*tan((pi/180)*50))+16)-sqrt(x_down(i)))));
CCC(i)=(-x_down(i)+x_up(i))/h(i);
SSS(i)=(((((16-64)/(-16*tan((pi/180)*50)-32+64*cot((pi/180)*70)))*(x_up(i)+16*tan((pi/180)*50))+16)-sqrt(x_down(i))))/h(i);
end
for i=1:X_number
    for j=1:Y_number
%         x_point(i,j)=x_down(i)+y(i,j)*cosd(teta(i));
%         y_point(i,j)=sqrt(x_down(i))+y(i,j)*sind(teta(i));
        x_point(i,j)=x_down(i)+y(i,j)*CCC(i);
        y_point(i,j)=sqrt(x_down(i))+y(i,j)*SSS(i);
    end
end

% boundary functions
x_down2=linspace(0,32,200);
y_down=sqrt(x_down2);
x_up2=linspace(-16*tan((pi/180)*50),32-64*cot((pi/180)*70),2);
y_up=((16-64)/(-16*tan((pi/180)*50)-32+64*cot((pi/180)*70)))*(x_up2+16*tan((pi/180)*50))+16;
x_left=linspace(-16*tan((pi/180)*50),0,2);
y_left=((0-16)/(0+16*tan((pi/180)*50)))*(x_left);
x_right=linspace(32-64*cot((pi/180)*70),32,2);
y_right=((sqrt(32)-64)/(64*cot((pi/180)*70)))*(x_right-32)+sqrt(32);
for i=1:length(x_down2)
    xp(i,1)=x_down2(i);
    yp(i,1)=y_down(i);
    zp(i,1)=0;
end
for i=1+length(x_down2):length(x_down2)+length(x_up2)
    xp(i,1)=x_up2(i-length(x_down2));
    yp(i,1)=y_up(i-length(x_down2));
    zp(i,1)=0;
end
for i=1+length(x_down2)+length(x_up2):length(x_down2)+length(x_up2)+length(x_left)
    xp(i,1)=x_left(i-length(x_down2)-length(x_up2));
    yp(i,1)=y_left(i-length(x_down2)-length(x_up2));
    zp(i,1)=0;
end
for i=1+length(x_down2)+length(x_up2)+length(x_left):length(x_down2)+length(x_up2)+length(x_left)+length(x_right)
    xp(i,1)=x_right(i-length(x_down2)-length(x_up2)-length(x_left));
    yp(i,1)=y_right(i-length(x_down2)-length(x_up2)-length(x_left));
    zp(i,1)=0;
end
%temp of points
temp=zeros(X_number,Y_number);
temp(:,1)=900;
temp(:,Y_number)=300;
temp(1,:)=1200;
temp(X_number,:)=600;
temp(1,1)=(temp(1,2)+temp(2,1))/2;
temp(1,Y_number)=(temp(1,Y_number-1)+temp(2,Y_number))/2;
temp(X_number,1)=(temp(X_number-1,1)+temp(X_number,2))/2;
temp(X_number,Y_number)=(temp(X_number,Y_number-1)+temp(X_number-1,Y_number))/2;
tempk=zeros(X_number,Y_number);
tempk(:,1)=900;
tempk(:,Y_number)=300;
tempk(1,:)=1200;
tempk(X_number,:)=600;
tempk(1,1)=temp(1,1);
tempk(1,Y_number)=temp(1,Y_number);
tempk(X_number,1)=temp(X_number,1);
tempk(X_number,Y_number)=temp(X_number,Y_number);

%initializing

for i=1+1:X_number-1 
    for j=1+1:Y_number-1
       temp(i,j)=200;
       tempk(i,j)=200;
    end 
end
iteration=0;
nn=1;
error=0;
err=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choose==1
%Jacobi

while nn 
    for i=1+1:X_number-1 
        for j=1+1:Y_number-1
delx22=x_point(i,j)-x_point(i-1,j);
dely22=y_point(i,j)-y_point(i-1,j);
delx2=sqrt(delx22^2+dely22^2);

delx11=x_point(i+1,j)-x_point(i,j);
dely11=y_point(i+1,j)-y_point(i,j);
delx1=sqrt(delx11^2+dely11^2);

dely222=y_point(i,j)-y_point(i,j-1);
delx222=x_point(i,j)-x_point(i,j-1);
dely2=sqrt(delx222^2+dely222^2);

dely111=y_point(i,j+1)-y_point(i,j);
delx111=x_point(i,j+1)-x_point(i,j);
dely1=sqrt(delx111^2+dely111^2);
ratio=(delx1/dely1)*((delx1+delx2)/(dely1+dely2));
alfaa=delx1/delx2;
betaa=dely1/dely2;
gama=(1+betaa)*ratio+(1+alfaa);
            tempk(i,j)=(1/gama)*(temp(i+1,j)+alfaa*temp(i-1,j)+ratio*(temp(i,j+1)+betaa*temp(i,j-1)));
            error=error+abs((tempk(i,j)-temp(i,j)));
        end
    end
    for i=1+1:X_number-1 
        for j=1+1:Y_number-1
            temp(i,j)=tempk(i,j);
        end
    end
    iteration=iteration+1;
    err(iteration)=error;
    if error<conv
    nn=0;
    end
error=0;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choose==2
%PSOR and Point Gauss Seidel

w=input('for the soul to converge (0<w<2) , for Under relaxation (0<w<1) , for Point Gauss Seidel Method (w=1).   w =');
while nn 
    for i=1+1:X_number-1 
        for j=1+1:Y_number-1
delx22=x_point(i,j)-x_point(i-1,j);
dely22=y_point(i,j)-y_point(i-1,j);
delx2=sqrt(delx22^2+dely22^2);

delx11=x_point(i+1,j)-x_point(i,j);
dely11=y_point(i+1,j)-y_point(i,j);
delx1=sqrt(delx11^2+dely11^2);

dely222=y_point(i,j)-y_point(i,j-1);
delx222=x_point(i,j)-x_point(i,j-1);
dely2=sqrt(delx222^2+dely222^2);

dely111=y_point(i,j+1)-y_point(i,j);
delx111=x_point(i,j+1)-x_point(i,j);
dely1=sqrt(delx111^2+dely111^2);
ratio=(delx1/dely1)*((delx1+delx2)/(dely1+dely2));
alfaa=delx1/delx2;
betaa=dely1/dely2;
gama=(1+betaa)*ratio+(1+alfaa);
            tempk(i,j)=temp(i,j)+(w/gama)*(temp(i+1,j)+alfaa*tempk(i-1,j)+ratio*(temp(i,j+1)+betaa*tempk(i,j-1))-gama*temp(i,j));
            error=error+abs((tempk(i,j)-temp(i,j)));
        end
    end
    for i=1+1:X_number-1 
        for j=1+1:Y_number-1
            temp(i,j)=tempk(i,j);
        end
    end
    iteration=iteration+1;
    err(iteration)=error;
    if error<conv
    nn=0;
    end
error=0;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choose==3
%LSOR and Line Gauss Seidel

w=input('for the soul to converge (0<w<1) , for Under relaxation (0<w<1) , for Line Gauss Seidel Method (w=1).   w =');

ratio=zeros(X_number,Y_number);
ratio2=zeros(X_number,Y_number);
alfaa=zeros(X_number,Y_number);
betaa=zeros(X_number,Y_number);
gama=zeros(X_number,Y_number);
delx2=zeros(X_number,Y_number);
delx1=zeros(X_number,Y_number);
dely2=zeros(X_number,Y_number);
dely1=zeros(X_number,Y_number);
for j=1+1:Y_number-1      
    for i=1+1:X_number-1 
delx22=x_point(i,j)-x_point(i-1,j);
dely22=y_point(i,j)-y_point(i-1,j);
delx2(i,j)=sqrt(delx22^2+dely22^2);

delx11=x_point(i+1,j)-x_point(i,j);
dely11=y_point(i+1,j)-y_point(i,j);
delx1(i,j)=sqrt(delx11^2+dely11^2);

dely222=y_point(i,j)-y_point(i,j-1);
delx222=x_point(i,j)-x_point(i,j-1);
dely2(i,j)=sqrt(delx222^2+dely222^2);

dely111=y_point(i,j+1)-y_point(i,j);
delx111=x_point(i,j+1)-x_point(i,j);
dely1(i,j)=sqrt(delx111^2+dely111^2);
alfaa(i,j)=delx1(i,j)/delx2(i,j);
betaa(i,j)=dely1(i,j)/dely2(i,j);
ratioo=(delx1(i,j)/dely1(i,j))*((delx1(i,j)+delx2(i,j))/(dely1(i,j)+dely2(i,j)));
gama(i,j)=(1+betaa(i,j))*ratioo+(1+alfaa(i,j));
        ratio(i,j)=ratioo;            
        ratio2(i,j)=-gama(i,j);
    end    
end

while nn 
        for j=1+1:Y_number-1
            for i=1:X_number-2          
            BB(i,1)=(1-w)*ratio2(i+1,j)*temp(i+1,j)-w*ratio(i+1,j)*(temp(i+1,j+1)+betaa(i+1,j)*tempk(i+1,j-1));
            end
            BB(1,1)=BB(1,1)-w*alfaa(2,j)*tempk(1,j);
            BB(X_number-2,1)=BB(X_number-2,1)-w*tempk(X_number,j);            
            NN=X_number-2;
            for i=1:NN-1    
                AA(i,i+1)=w;
                AA(i,i)=ratio2(i+1,j);    
                AA(i+1,i)=w*alfaa(i+2,j);
            end            
            AA(NN,NN)=ratio2(NN+1,j);
RR(1)=0;
PP=zeros(1,NN);
QQ=zeros(1,NN-1);
RR=zeros(1,NN);
YY=zeros(1,NN-1);
for i=1:NN

      PP(i)=AA(i,i);

end
for i=1:NN-1
    QQ(i)=AA(i,i+1);

end
for i=1:NN-1
    RR(i+1)=AA(i+1,i);

end
YY(1)=QQ(1)/PP(1);

for i=2:NN-1
    YY(i)=QQ(i)/(PP(i)-RR(i)*YY(i-1));
end
WW(1)=BB(1)/PP(1);
for i=2:NN
    WW(i)=(BB(i)-RR(i)*WW(i-1))/(PP(i)-RR(i)*YY(i-1));
end

XX(NN+1)=WW(NN);
for i=NN-1:-1:1
    XX(i+1)=WW(i)-YY(i)*XX(i+2);
end
            for i=1:X_number-2
                tempk(i+1,j)=XX(i+1);
                error=error+abs((tempk(i+1,j)-temp(i+1,j)));
            end
        end        
    for i=1+1:X_number-1 
        for j=1+1:Y_number-1
            temp(i,j)=tempk(i,j);
        end
    end
    iteration=iteration+1;
    err(iteration)=error;
    if error<conv
    nn=0;
    end
error=0;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % corners
temp(1,1)=(temp(1,2)+temp(2,1)+temp(2,2))/3;
temp(1,Y_number)=(temp(1,Y_number-1)+temp(2,Y_number)+temp(2,Y_number-1))/3;
temp(X_number,1)=(temp(X_number-1,1)+temp(X_number,2)+temp(X_number-1,2))/3;
temp(X_number,Y_number)=(temp(X_number,Y_number-1)+temp(X_number-1,Y_number)+temp(X_number-1,Y_number-1))/3;
%tamp at x=16 or y=32
x16=0;
y32=0;
for j=1:Y_number
    for i=1:X_number
        if x_point(i,j)>=16 & x_point(i-1,j)<=16
            x16=x16+1;
            tempx16(x16,1)=((temp(i,j)*(16-x_point(i-1,j))+(temp(i-1,j)*(x_point(i,j)-16))))/(x_point(i,j)-x_point(i-1,j));
            tempx16(x16,2)=((y_point(i,j)*(16-x_point(i-1,j))+(y_point(i-1,j)*(x_point(i,j)-16))))/(x_point(i,j)-x_point(i-1,j));
%             if j~=Y_number & j~=1
%             tempx16(x16,1)=(temp(i,j)+temp(i-1,j)+(temp(i,j+1)+temp(i-1,j+1)+temp(i,j-1)+temp(i,j-1))/2)/4;
%             tempx16(x16,2)=(y_point(i,j)+y_point(i-1,j)+(y_point(i,j+1)+y_point(i-1,j+1)+y_point(i,j-1)+y_point(i,j-1))/2)/4;
%             end
%             if j==Y_number
%             tempx16(x16,1)=(temp(i,j)+temp(i-1,j)+(temp(i,j-1)+temp(i,j-1))/2)/3;
%             tempx16(x16,2)=(y_point(i,j)+y_point(i-1,j)+(y_point(i,j-1)+y_point(i,j-1))/2)/3;                
%             end
%             if j==1
%             tempx16(x16,1)=(temp(i,j)+temp(i-1,j)+(temp(i,j+1)+temp(i-1,j+1))/2)/3;
%             tempx16(x16,2)=(y_point(i,j)+y_point(i-1,j)+(y_point(i,j+1)+y_point(i-1,j+1))/2)/3;                
%             end
        end
        if y_point(i,j)>=32 & y_point(i,j-1)<=32
            y32=y32+1;   
            tempy32(y32,1)=((temp(i,j)*(32-y_point(i,j-1))+(temp(i,j-1)*(y_point(i,j)-32))))/(y_point(i,j)-y_point(i,j-1));
            tempy32(y32,2)=((x_point(i,j)*(32-y_point(i,j-1))+(x_point(i,j-1)*(y_point(i,j)-32))))/(y_point(i,j)-y_point(i,j-1));
%             if i~=X_number & i~=1
%             tempy32(y32,1)=(temp(i,j)+temp(i,j-1)+(temp(i+1,j)+temp(i+1,j)+temp(i-1,j)+temp(i-1,j))/2)/4;
%             tempy32(y32,2)=(x_point(i,j)+x_point(i,j-1)+(x_point(i+1,j)+x_point(i+1,j)+x_point(i-1,j)+x_point(i-1,j))/2)/4;
%             end
%             if i==X_number
%             tempy32(y32,1)=(temp(i,j)+temp(i,j-1)+(temp(i-1,j)+temp(i-1,j))/2)/4;
%             tempy32(y32,2)=(x_point(i,j)+x_point(i,j-1)+(x_point(i-1,j)+x_point(i-1,j))/2)/3;               
%             end
%             if i==1
%             tempy32(y32,1)=(temp(i,j)+temp(i,j-1)+(temp(i+1,j)+temp(i+1,j))/2)/4;
%             tempy32(y32,2)=(x_point(i,j)+x_point(i,j-1)+(x_point(i+1,j)+x_point(i+1,j))/2)/3;                
%             end
        end        
    end 
end
%final output matrix
kk=0;
for j=1:Y_number
    for i=1:X_number
        kk=kk+1;
        output(kk,3)=temp(i,j);
        output(kk,2)=y_point(i,j);
        output(kk,1)=x_point(i,j);
    end 
end
fprintf('ZONE= I= %1.0f,  J= %g \n',X_number,Y_number)
display('Variables= "x"   ,   "y"   ,   "T"')
output
display('Variable= "temp at line x=16"')
tempx16(:,1)
display('Variable="temp at line y=32"')
tempy32(:,1)
for i=1:length(err)
     numberofiter(i,1)=i;
     numberofiter(i,2)=err(i);
end
numberofiter
display('Variables= "iteration"  ,  "error"')
fprintf('Solution coverged after %1.0f iterations \n',length(err))
%%%%plot part
%Mesh plot
subplot(2,3,1)
line(x_point,y_point,'color','green')
hold on
for i=1:X_number
    line(x_point(i,:),y_point(i,:),'color','green')
    hold on
end
title(' Mesh ')
xlabel(' x ')
ylabel(' y ')

%Boundary plot
subplot(2,3,2)
plot(x_down2,y_down,'color','red')
hold on
plot(x_up2,y_up,'color','red')
hold on
plot(x_left,y_left,'color','red')
hold on
plot(x_right,y_right,'color','red')
hold on
title(' Boundaries ')
xlabel(' x ')
ylabel(' y ')
%Temp plot
subplot(2,3,3)
line(x_point,y_point,'color','green')
hold on
for i=1:X_number
    line(x_point(i,:),y_point(i,:),'color','green')
    hold on
end
title(' Temperature contour ')
xlabel(' x ')
ylabel(' y ')
hold on
surf(x_point,y_point,temp)
colormap jet
view(2)
hold on
subplot(2,3,6)
title(' Temperature contour ')
xlabel(' x ')
ylabel(' y ')
hold on
s=surf(x_point,y_point,temp);
colormap jet
colorbar
s.EdgeColor= 'none';
view(2)
hold on
%Temp plot at x=16 and y=32
% x_down2=linspace(0,32,200);
% y_down=sqrt(x_down2);
% x_up2=linspace(-16*tan((pi/180)*50),32-64*cot((pi/180)*70),2);
% y_up=((16-64)/(-16*tan((pi/180)*50)-32+64*cot((pi/180)*70)))*(x_up2+16*tan((pi/180)*50))+16;
% x_left=linspace(-16*tan((pi/180)*50),0,2);
% y_left=((0-16)/(0+16*tan((pi/180)*50)))*(x_left);
% x_right=linspace(32-64*cot((pi/180)*70),32,2);
% y_right=((sqrt(32)-64)/(64*cot((pi/180)*70)))*(x_right-32)+sqrt(32);
subplot(2,3,4)
yx16=linspace(sqrt(16),((sqrt(32)-64)/(64*cot((pi/180)*70)))*(16-32)+sqrt(32),length(tempx16));
xy32=linspace(((32-16)/((16-64)/(-16*tan((pi/180)*50)-32+64*cot((pi/180)*70))))-16*tan((pi/180)*50),((32-sqrt(32))/(((sqrt(32)-64)/(64*cot((pi/180)*70)))))+32,length(tempy32));
for nn=1:length(tempx16)
    xx16(nn)=16;
end
for nn=1:length(tempy32)
    yy32(nn)=32;
end
plot(tempx16(:,2),tempx16(:,1),'*')
hold on
title('Temperature at x=16')
xlabel(' y ')
ylabel(' Temperature ')
subplot(2,3,5)
plot(tempy32(:,2),tempy32(:,1),'*')
hold on
title('Temperature at y=32')
xlabel(' x ')
ylabel(' Temperature ')
subplot(2,3,2)
line(xy32,yy32)
hold on
line(xx16,yx16)
% Export to temperature_field.plt
fid=fopen('temperature_field.plt','w');
fprintf(fid,'VARIABLES= "x"   ,   "y"   ,   "T" \n');
fprintf(fid,'ZONE I= %1.0f,  J= %g \n',X_number,Y_number);
for i=1: size(output,1)
    for j=1:size(output,2)
        fprintf(fid,'%f ',output(i,j));
    end 
        fprintf(fid,'\n');
end
fclose(fid);
% Export output to Exel
title4={'points','for','Gambit'};
title3={'temp','y at x=16'};
title2={'temp','x at y=32'};
title1={'x'   ,   'y'   ,   'T'};
title0={'ZONE','i=',X_number,'j=',Y_number};
recycle on
delete('temp_of_surface.xls');
xlswrite('temp_of_surface',[x_point(:),y_point(:),temp(:)],'Sheet1','A3');
xlswrite('temp_of_surface',title1,'Sheet1','A2');
xlswrite('temp_of_surface',title0,'Sheet1','A1');
recycle on
delete('temp_at_x16.xls');
xlswrite('temp_at_x16',[tempx16],'Sheet1','A3');
xlswrite('temp_at_x16',title3,'Sheet1','A2');
xlswrite('temp_at_x16',title0,'Sheet1','A1');
recycle on
delete('temp_at_y32.xls');
xlswrite('temp_at_y32',[tempy32],'Sheet1','A3');
xlswrite('temp_at_y32',title2,'Sheet1','A2');
xlswrite('temp_at_y32',title0,'Sheet1','A1');
% Export points to Exel for Gambit
recycle on
delete('points_for_gambit.xls');
xlswrite('points_for_gambit',title4,'Sheet1','A1');
xlswrite('points_for_gambit',[xp(:),yp(:),zp(:)],'Sheet1','A2');
format short