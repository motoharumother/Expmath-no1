#include <stdio.h>
#include <math.h>

int main(){
  int n = 1000;
  int a[n-2][2];
  //a には　（数字 (2~1000)、 素数 or 合成数 を表すラベル（1 or 0）） を表す。
  for(int i=0; i<n-1; i++){
    a[i][0] = i+2;
    //まずは全部素数表すラベルをつける
    a[i][1] = 1;
  }

  int l = 0;
  while(1) {
    // n の平方根を超えたら break
    if(a[l][1]>=sqrt(n)){
      break;
    }
    //素数でなかったら  ans[l][1] の値を０にする
    else if(a[l][1]==1){
      for(int m=l+1; m<n-1; m++){
        if(a[m][1]==1){
          //a[m][0] が a[l][0] で割り切れたら　a[m][1] のラベルは 0
          if(a[m][0]%a[l][0]==0){
            a[m][1] = 0;
          }
        }
      }
    }
    l += 1;
  }
  //print
  for(int k=0; k<n-1; k++){
    if(a[k][1]==1){
      printf("%d,", a[k][0]);
    }

  }
return 0;
}
