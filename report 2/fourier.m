clear all

close all

x = [-pi:pi/1000:pi]; %[a:b:c]で[a,a+b,a+2b,a+3b,...,c]と��?��?ベクトルを表現
size = 50;
error = zeros(size,1);

% n=1から50までの誤差を計算
for m = 1:size
  n = m + 100;
  b = zeros(n,1);
  te = 0;
  for k = 1:n
    b(k,1) = 1/pi * quad(@(x) x * sin(k*x),-pi,pi);
  end
  
  for i = 1:2001
    sn = 0;
    tmpx = x(1,i);
    for j = 1:n
      sn = sn + b(j,1) * sin(j*tmpx);
    end
    te = te + ((sn-x(1,i))*(sn-x(1,i)));
  end
  
  error(m,1) = te;
end
    
m = [101:1:100+size];
 
plot(m,error)