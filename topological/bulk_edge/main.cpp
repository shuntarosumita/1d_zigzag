//-*- coding: utf-8 -*-
// Time-stamp: <2016-02-03 14:57:25 shunta>
// ----------------------------------------
// 1次元ジグザグチェーン構造の超伝導計算
// ----------------------------------------
#include "zigzag.h"
#include "topological.h"
typedef std::complex<double> Complex;
using namespace std;
using namespace Zigzag;

// メインプログラム
int main(int argc, char* argv[]) {
  // 物理量
  double haf;                            // 反強磁性的な分子場
  double mu;                             // 化学ポテンシャル
  Cvector order(2);                      // 秩序パラメータ
  Cmatrix H(DIM * N, DIM * N);           // ハミルトニアン
  Cmatrix U(DIM * N, DIM * N);           // H の固有ベクトル
  Dvector eigen(DIM * N);                // 固有値
  Zmatrix zmatrix;                       // 複素数行列
  Cmatrix H_cood(DIM, DIM);              // 座標にFourier変換されたハミルトニアン
  double k;                              // 波数
  int i, j, l;                           // ループ用変数
  size_t s1, s2;

  // 秩序パラメータの初期化
  Cvector order_init(2);
  order_init(0) = Complex(0.1, 0.0);
  order_init(1) = Complex(-0.1, 0.0);
  order = order_init;

  // ファイル入出力関連
  stringstream filename;                 // ファイル名
  ofstream output_file[1];               // 出力ファイル

  // 入力されたパラメータを設定
  if(argc > 3) {
    t[0] = atof(argv[1]);
    t[1] = 1.0;
    if(t[0] < 0.0 || t[0] > 1.0) {
      cerr << "t_1 range: [0.0, 1.0]" << endl;
      return 0;
    }

    mu = atof(argv[2]);
    if(mu <= -3.5 || mu > 2.5) {
      cerr << "mu range: (-3.5, 2.5]" << endl;
      return 0;
    }

    haf = atof(argv[3]);
    if(haf < 0.0 || haf > 4.0) {
      cerr << "h^AF range: [0.0, 4.0]" << endl;
      return 0;
    }
  }
  else {
    cerr << "Usage: " << argv[0] << " t_1 mu h^AF" << endl;
    return 0;
  }

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/t" << t[0] << "_mu" << mu << "_h" << haf << ".d";
  output_file[0].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 最初の行のコメントを挿入
  output_file[0] << "# i eigen[i]" << endl;

  // 座標のループ
  for(i = 0; i < N; i++) {
    for(j = 0; j < N; j++) {

      // 座標Hamiltonianの初期化
      for(s1 = 0; s1 < H_cood.size1(); s1++){
        for(s2 = 0; s2 < H_cood.size2(); s2++){
          H_cood(s1, s2) = 0.0;
        }
      }

      // PBCを取り除く
      // if(abs(i - j) < (double)N * 0.5) {

        // Fourier変換
        for(l = 0; l < N; l++) {
          k = 2.0 * N_inv * M_PI * (double)l - M_PI;

          // Hamiltonianを定義する
          zmatrix.init_zmatrix(order, haf, mu, 0.0, k);

          for(s1 = 0; s1 < H_cood.size1(); s1++){
            for(s2 = 0; s2 < H_cood.size2(); s2++){
              H_cood(s1, s2) += N_inv * exp(imag_unit * k * (Complex)(i - j)) * zmatrix.h(s1, s2);
            }
          }
        }

      // }

      for(s1 = 0; s1 < H_cood.size1(); s1++){
        for(s2 = 0; s2 < H_cood.size2(); s2++){
          H((size_t)(DIM * i) + s1, (size_t)(DIM * j) + s2) = H_cood(s1, s2);
        }
      }

    }
  }

  Topological::eigen(H, U, eigen);

  // 固有エネルギーを出力する
  for(i = 0; i < DIM * N; i++) {
    output_file[0] << fixed << setprecision(8) << i << "  " << eigen[i] << endl;
  }

  output_file[0].close();

  return 0;
}

//
// 行列 A の固有値と固有ベクトルを求める
//
void Topological::eigen(Cmatrix &A, Cmatrix &U, Dvector &eigen) {
  Dvector eigen_copy(A.size1());
  matrix<Complex, column_major> A_copy(A.size1(), A.size2());
  long info;

  for(size_t i = 0; i < A.size1(); i++) {
    for(size_t j = i; j < A.size2(); j++) {
      A_copy(i, j) = A(i, j);
    }
  }

  info = lapack::heev('V', 'U', A_copy, eigen_copy, lapack::optimal_workspace());
  BOOST_UBLAS_CHECK(info == 0, internal_logic());

  U = A_copy;
  eigen.swap(eigen_copy);
}
