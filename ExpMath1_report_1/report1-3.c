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
		double renm = 1.0/norm;
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


//Rを作る時のベクトルのノルムが欲しいので関数を定義する
double Get_Vec_Norm(int k, const double u[][MaxDim], double trg[MaxDim])
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
		double tp = VecNorm(MaxDim, trg);
		VecRNormalize(MaxDim, trg);

		return tp;
	}
}

int Get_Norm(int n, double v[][MaxDim],  double list[MaxDim])
{
	double trg[MaxDim][MaxDim];

	for (int i = 0; i < n; ++i) {
		double l=0;
		VecCopy(MaxDim, v[i], trg[i]);
		l=Get_Vec_Norm(i, v, trg[i]);
		VecCopy(MaxDim, trg[i], v[i]);
		list[i] = l;
	}
	return 0;
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
		if (l != 0.0){
			k += 1;
		}

	}
	return k;
}

//行列の転置 v -> T(v)
void Transpose(int n, double v[][MaxDim]){
	double list[n][n];
	for (int i=0; i<n; i++){
		VecCopy(MaxDim, v[i], list[i]);
	}
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			v[i][j] = list[j][i];
		}
	}
}

/*やること
１、Aを転置してA’　を作る
２、A’に直交化法を適用　→　Bを作る
３、Bを転置してＱを作る（こうすれば、Aの列ベクトルに対して直交化したことになる）
４、Rを証明どうりに作る
５、RxＱを返す
*/
double RtrnRQ(int n, double a[][MaxDim], double result[][MaxDim]){
	double a_prm[MaxDim][MaxDim];
	double b[MaxDim][MaxDim];
	double q[MaxDim][MaxDim];
	double r[MaxDim][MaxDim];
	double list[MaxDim];

	for(int i=0; i<n; i++){
		for(int j =0; j<n; j++){
			r[i][j] = 0.0;
		}
	}
	//step 1
	for(int i=0; i<n; i++){
		VecCopy(n, a[i], a_prm[i]);
	}

	Transpose(n, a_prm);
	//step 2
	VecGramSchmidt(n, a_prm, b);
	//step 3
	for(int i=0; i<n; i++){
		VecCopy(n, b[i], q[i]);
	}
	Transpose(n, q);
	//step 4

	//rの対角成分の取得
	for(int i=0; i<n; i++){
		VecCopy(n, a[i], a_prm[i]);
	}
	Transpose(n, a_prm);
	Get_Norm(n, a_prm, list);

	for(int i=0; i<n; i++){
		VecCopy(n, a[i], a_prm[i]);
	}
	Transpose(n, a_prm);

	for(int i=0; i<n; i++){
		for(int j=i; j<n; j++){
			if(i==j){
				r[i][j] = list[j];
			}else{
				r[i][j] = VecIP(n, a_prm[j], b[i]);
			}
		}
	}

	//step 5
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			result[i][j] = 0.0;
			for(int k=0; k<n; k++){
				result[i][j] += r[i][k]*q[k][j];
			}
		}
	}


/*テスト
	//test list q r の確認

	printf("list =");

		VecPrint(n, list);

		for (int i = 0; i < n; ++i) {
		printf("r[%d] = ",i);
		VecPrint(MaxDim, r[i]);
		}

		printf("\n");
		for (int i = 0; i < n; ++i) {
		printf("q[%d] = ",i);
		VecPrint(MaxDim, q[i]);
		}
		printf("\n");
*/
/*
	//test
	//しっかり　A=QR　か・・

	double test[MaxDim][MaxDim];
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			test[i][j] = 0.0;
			for(int k=0; k<n; k++){
				test[i][j] += q[i][k]*r[k][j];
			}
		}
	}
	for (int i = 0; i < n; ++i) {
	printf("Qr[%d] = ",i);
	VecPrint(MaxDim, test[i]);
	}
	*/
/*
	//test
	//しっかり　E=QxT(Q)　か・・
	//double test[MaxDim][MaxDim];
	double test_2[MaxDim][MaxDim];
	for(int i=0; i<n; i++){
		VecCopy(n, q[i], test_2[i]);
	}
	Transpose(n,test_2);
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			test[i][j] = 0.0;
			for(int k=0; k<n; k++){
				test[i][j] += q[i][k]*test_2[k][j];
			}
		}
	}
	for (int i = 0; i < n; ++i) {
	printf("QxT(Q)[%d] = ",i);
	VecPrint(MaxDim, test[i]);
	}
//ここまでテスト
*/
	return 0;

}



int main(){
	//double _v[MaxDim][MaxDim];
	double a[MaxDim][MaxDim]={{2,0,0,-1},{0,2,-1,0},{0,-1,2,0},{-1,0,0,2}};
	// ={{0,-1,0},{-1,0,0},{0,3,2}};
	//double b[MaxDim][MaxDim] = {{1,0,0},{0,2,0},{0,0,2}};
	//double p[MaxDim][MaxDim] = {{1,1,1},{1,2,2},{1,2,3}};
	//double in_p[MaxDim][MaxDim] = {{2,-1,0},{-1,2,-1},{0,-1,1}};
	int dim = 4;
/*
	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			_v[i][j] = 0.0;
			for(int k=0; k<dim; k++){
				_v[i][j] += in_p[i][k]*b[k][j];
			}
		}
	}
	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			a[i][j] = 0.0;
			for(int k=0; k<dim; k++){
						a[i][j] += _v[i][k]*p[k][j];
				}
			}
		}
*/

	printf("A0\n" );
	for (int i = 0; i < dim; ++i) {
	printf("result[%d] = ",i);
	VecPrint(MaxDim, a[i]);
	}
	printf("\n");

	for (int k=1; k<700; k++){
	double result[MaxDim][MaxDim];
	RtrnRQ(dim, a, result);
	printf("A%d\n", k);
	for (int i = 0; i < dim; ++i) {
	printf("result[%d] = ",i);
	VecPrint(MaxDim, result[i]);
	}
	//AにRxQを代入して次に進む
	for(int i=0; i<dim; i++){
	VecCopy(dim, result[i], a[i]);
	}
	printf("\n");
	}

return 0;
}
