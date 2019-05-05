#include <stdio.h>
#include <math.h>
//関数定義
double f(double x){
  double y = exp(x)+(x*x) - 5;
  return y;
}

double df(double x){
  double h = 1.0e-5;
  double fprime =  (f(x+h)-f(x-h))/(2*h);
  return fprime;
}

int main() {
  double eps = 1.0e-8;

  double x_n, x_n_1;
  printf("startpoint:  ");
  scanf("%lf", &x_n);

  while (1) {
    x_n_1 = x_n - (f(x_n)/df(x_n));
    if (fabs(x_n-x_n_1)<eps) {
      break;
    }
    x_n = x_n_1;
  }
  printf("収束点は%lf\n", x_n_1);
  printf("f(%lf)=%lf\n",  x_n_1, f(x_n_1));
}
