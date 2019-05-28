#include <stdio.h>

//数列の計算関数を定義
double x(int n){
  //初項の準備
  double x_zero, x_one;
  double x_two;

  x_zero = 2;
  x_one = (double)2/3;

  for (int i=0; i<=n-2 ; i++) {
    x_two = ((double)7/3)*x_one - ((double)2/3)*x_zero;
    x_zero = x_one;
    x_one = x_two;
  }
  return x_two;
}

int main(){
  double ans;

  //数字を入手
  int n;
  printf("0を含めて自然数を入力してください：　");
  scanf("%d", &n);

  //n=0 or 1 の時は例外として処理
  if (n==0) {
    ans = (int)2;
    printf("a_%d = %lf\n", n, ans);
  }else if(n==1){
    ans = (double)2/3;
    printf("a_%d = %lf\n", n, ans);
  }else{
    ans = x(n);
    printf("X_%d = %lf\n", n, ans);
  }

  return 0;
}
