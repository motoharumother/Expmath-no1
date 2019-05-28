//グラムシュミットの正規直交化法

#include <stdio.h>
#include <math.h>

#define	MaxDim	4
// 簡単のため，ベクトルのサイズも数も4とするが，当然サイズと数は違ってよい
// その場合プログラムは微修正が必要


/*
ここに先週の VecIP から VecRNormalize までの全ての関数を書く
*/
double	VecIP(int dim, const double u[], const double v[])
{
	int i;
	double sum=0.0;
	// 0 に各成分の積を加えていく
	for (i = 0; i < dim; ++i) {
		sum = sum + u[i] * v[i];
	}
	return sum;
}

// ベクトルの表示
void	VecPrint(int dim, const double v[])
{
	int i;
	for (i = 0; i < dim; ++i) {
		printf("\t%f", v[i]); /* \t はタブ．*/
	}
	printf("\n");
}

// ベクトルを複写: src -> trg
void	VecCopy(int dim, const double src[], double trg[])
{
	int i;
	for (i = 0; i < dim; ++i){
		trg[i] = src[i];
	}
}

// ベクトルの和：u + v -> trg
void	VecPlus(int dim, const double u[], const double v[], double trg[])
{
	int i;
	for (i = 0; i < dim; ++i){
		trg[i] = u[i] + v[i];
	}
}

// ベクトルの和：trg + v -> trg （データ書き換え版）
void	VecRPlus(int dim, double trg[], const double v[])
{
	int i;
	for (i = 0; i < dim; ++i){
		trg[i] = trg[i] + v[i];
	}
}

// ベクトルの差: u - v -> trg
void	VecMinus(int dim, const double u[], const double v[], double trg[])
{
	//ここを埋める
	for(int i=0; i<dim; i++){
		trg[i] = u[i] - v[i];
	}
}

// ベクトルの差: trg - v -> trg （データ書き換え版）
void	VecRMinus(int dim, double trg[], const double v[])
{
	//ここを埋める
	for(int i=0; i<dim; i++){
		trg[i] -= v[i];
	}
}

// ベクトルのスカラー倍:  scale * src -> trg
void	VecScalarMultiply(int dim, double scale, const double src[], double trg[])
{
	//ここを埋める
	for(int i=0; i<dim; i++){
		trg[i] = scale * src[i];
	}
}

// ベクトルのスカラー倍:  scale * trg -> trg  （データ書き換え版）
void	VecRScalarMultiply(int dim, double scale, double trg[])
{
	//ここを埋める
	for(int i=0; i<dim; i++){
		trg[i] = scale * trg[i];
	}
}

// ベクトルのノルム
double VecNorm(int dim, const double u[])
{
	return sqrt(VecIP(dim, u, u));
}

/* ベクトル trg をスカラー倍して，長さ 1 のものに書き換える．
   関数値はベクトルの長さ（0ベクトルだと0.0，それ以外は1.0）
*/
double VecRNormalize(int dim, double trg[])
{
	//ここを埋める
	double norm = VecNorm(dim, trg);
	if(norm<=0){
		printf("this is 0 vector.\n");
		return 0.0;
	}else{
		double renm = 1/norm;
		VecRScalarMultiply(dim, renm, trg);
		return 1.0;
	}
}



/*
グラム-シュミットの直交化法は，ベクトルの列の長さに関して帰納的である．
そこで，d 項行ベクトル u[0], u[1], ..., u[k-1] が正規直交していて，
さらに，u[0], u[1], ..., u[k-1], u[k] が，線形独立であるとき，trg を書き換えて，
a) u[0], u[1], ..., u[k-1], trg は，正規直交
b) u[0], u[1], ..., u[k-1], trg で張られる部分空間は，
u[0], u[1], ..., u[k-1], u[k] で張られる部分空間と一致する
ようにする関数 VecOrthogolize(k, d, u[][],trg[]) を定義せよ．
ただし，元の u[0], u[1], ..., u[k] が線形独立ではないときは，関数の値は 0,
線形独立のときは，関数の値は，0 以外の double とする．
*/

double VecOrthogonalize(int k, const double u[][MaxDim], double trg[MaxDim])
{
	int i;
	double s;
	double temp[MaxDim];

	VecCopy(MaxDim, u[k], trg);
	if (VecNorm(MaxDim, trg) == 0.0) {
		return 0.0;
	}else{
		for (i = 0; i < k; ++i) {
			s = VecIP(MaxDim, u[i], trg);
			VecScalarMultiply(MaxDim, s, u[i], temp);
			VecRMinus(MaxDim, trg, temp);
			if (VecNorm(MaxDim, trg) == 0.0){
				return 0.0;
			}
		}
		return VecRNormalize(MaxDim, trg);
	}
}

//グラムシュミット
// 戻り値は独立なベクトルの本数

int VecGramSchmidt(int n, double v[][MaxDim], double trg[][MaxDim])
{
	int i;
	int k=0;
	double l;

	for (i = 0; i < n; ++i) {
		VecCopy(MaxDim, v[i], trg[i]);
		l=VecOrthogonalize(i, v, trg[i]);
		VecCopy(MaxDim, trg[i], v[i]);
		if (l == 0.0)
		k = k;
		else
		k=k+1;
	}
	return k;
}



int main(){

	double v[MaxDim][MaxDim] =  {	{9.4, 0.0, 0.0, 0.0},
	//{3.1, 1.0, 0.0, 0.0},//
	{4.7, 0.0, 0.0, 0.0},  
	{1.0, 6.0, 7.8, 0.0},
	{-4.4, 2.1, 0.7, - 5.7}};
	double trg[MaxDim][MaxDim];
	int dim = 4;
	int i, n;

	//先週の関数をコピーしたら次のコメントを外す
	
	for (i = 0; i < dim; ++i) {
	printf("v[%d] = ", i);
	VecPrint(MaxDim, v[i]);
	}


//グラムシュミットが完成したら次のコメントを外す

printf("\n");
n = VecGramSchmidt(dim, v, trg);
for (i = 0; i < dim; ++i) {
printf("trg[%d] = ",i);
VecPrint(MaxDim, trg[i]);
}
printf("\nthe number of independent vectors = %d\n", n);

return 0;
}
