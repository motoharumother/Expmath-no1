#include <stdio.h>
#include <math.h>
//関数定義
double f(double x){
  double sigma = 13;
  double y= exp(-(x*x)/(2*sigma))/sqrt(2*M_PI*sigma);
  return y;
}

int main() {
  double dx = 1.0e-7;
  double _rlt = 0;
  double rlt_ = 0;
  double rlt;

  double _x;
  printf("Enter startpoint: ");
  scanf("%lf", &_x);

  double end;
  printf("Enter endpoint: ");
  scanf("%lf", &end);

  double x_ = _x;

  while (_x<=end) {
    _x += dx;
    _rlt += f(_x)*dx;
  }

  while (x_<=end) {
    rlt_ += f(x_)*dx;
    x_ += dx;
  }

  rlt = (_rlt+rlt_)/2;
  printf("%f\n", rlt);
}
