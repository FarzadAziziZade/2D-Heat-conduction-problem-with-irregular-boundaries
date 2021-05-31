clc;
clear all;
close all;
X_number=input('X number =');
b2=input('X Beta (1 < X Beta < inf ) =');
Y_number=input('Y number =');
b=input('Y Beta (1 < Y Beta < inf ) =');
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
teta(i)=pi+acot((-x_down(i)+x_up(i))/(((((16-64)/(-16*tan((pi/180)*50)-32+64*cot((pi/180)*70)))*(x_up(i)+16*tan((pi/180)*50))+16)-sqrt(x_down(i)))));
end
for i=1:X_number
    for j=1:Y_number
        x_point(i,j)=x_down(i)+y(i,j)*cos(teta(i));
        y_point(i,j)=sqrt(x_down(i))+y(i,j)*sin(teta(i));
    end
end
subplot(1,2,1)
line(x_point,y_point,'color','green')
hold on
for i=1:X_number
    line(x_point(i,:),y_point(i,:),'color','green')
    hold on
end
title(' Mesh ')
xlabel(' x ')
ylabel(' y ')
% boundary functions
subplot(1,2,2)
x_down2=linspace(0,32,100000);
y_down=sqrt(x_down2);
x_up2=linspace(-16*tan((pi/180)*50),32-64*cot((pi/180)*70),3);
y_up=((16-64)/(-16*tan((pi/180)*50)-32+64*cot((pi/180)*70)))*(x_up2+16*tan((pi/180)*50))+16;
x_left=linspace(-16*tan((pi/180)*50),0,3);
y_left=((0-16)/(0+16*tan((pi/180)*50)))*(x_left);
x_right=linspace(32-64*cot((pi/180)*70),32,3);
y_right=((sqrt(32)-64)/(64*cot((pi/180)*70)))*(x_right-32)+sqrt(32);
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