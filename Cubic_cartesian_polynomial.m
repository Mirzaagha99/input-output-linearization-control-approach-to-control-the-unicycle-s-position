s=0:0.001:1;
k=2;
theta_i=0;
theta_f=pi/2;
b=0.04;
x_i=0;
x_f=4;

y_i=0;
y_f=1;

alfa_x=k*cos(theta_f)-3*x_f;
alfa_y=k*sin(theta_f)-3*y_f;
beta_x=k*cos(theta_i)-3*x_i;
beta_y=k*sin(theta_i)-3*y_i;

tf=10;
time=linspace(0,tf,length(s));


a0=0;
a1=0;
a2=3/tf^2;
a3=-2/tf^3;
s_t=a3*time.^3+a2*time.^2+a1*time+a0;
s_t_dot=3*a3*time.^2+2*a2*time+a1;




for i = 1:length(s)
    si=s_t(i);
    
    x(i)=si^3*x_f-(si-1)^3*x_i+alfa_x*si^2*(si-1)+beta_x*si*(si-1)^2;
    y(i)=si^3*y_f-(si-1)^3*y_i+alfa_y*si^2*(si-1)+beta_y*si*(si-1)^2;

    xp(i)=alfa_x*si^2 + 3*si^2*x_f + beta_x*(si - 1)^2 + 2*alfa_x*si*(si - 1) + beta_x*si*(2*si - 2);
    yp(i)=alfa_y*si^2 + 3*si^2*y_f + beta_y*(si - 1)^2 + 2*alfa_y*si*(si - 1) + beta_y*si*(2*si - 2);
    
    xpp(i)=4*alfa_x*si + 2*beta_x*si + 6*si*x_f + 2*alfa_x*(si - 1) + 2*beta_x*(2*si - 2);
    ypp(i)=4*alfa_y*si + 2*beta_y*si + 6*si*y_f + 2*alfa_y*(si - 1) + 2*beta_y*(2*si - 2);
    
    theta(i)=atan2(yp(i),xp(i));
    
    v(i)=sqrt(xp(i)^2+yp(i)^2);
    v_t(i)=v(i)*s_t_dot(i);
   

    o(i)=(ypp(i)*xp(i)-xpp(i)*yp(i))/(xp(i)^2+yp(i)^2);
    o_t(i)=o(i)*s_t_dot(i); 
    

end
for i = 1:length(s)
    y_1_d(i,1) = [time(i)];
    y_1_d(i,2) = x(i)+b*cos(theta(i));
    y_2_d(i,1) = [time(i)];
    y_2_d(i,2) = y(i)+b*sin(theta(i));
    x_before(i,1)=[time(i)];
    x_before(i,2)=x(i);
    y_before(i,1)=[time(i)];
    y_before(i,2)=y(i);
    theta_before(i,1)=[time(i)];
    theta_before(i,2)=theta(i);
    v_before(i,1)=[time(i)];
    v_before(i,2)=v_t(i);
    o_before(i,1)=[time(i)];
    o_before(i,2)=o_t(i);
end
figure(1)
plot(x,y,'LineWidth',2) 
xlabel('x')
ylabel('y')

figure(2)
plot(time,theta,'LineWidth',2)
xlabel('t')
ylabel('theta')

figure(3)
plot(time,v_t,'LineWidth',2)
xlabel('t')
ylabel('v')

figure(4)
plot(time,o_t,"LineWidth",2)
xlabel('t')
ylabel('omega')

figure(5)
plot(time,s_t,'LineWidth',2)
xlabel('t')
ylabel('s')

figure(6)
plot(time,x,'LineWidth',2)
xlabel('t')
ylabel('x')

figure(7)
plot(time,y,'LineWidth',2)
xlabel('t')
ylabel('y')

