//-*- coding: utf-8 -*-
// Time-stamp: <2016-02-03 14:26:12 shunta>
// ----------------------------------------
// 1次元ジグザグチェーン構造の超伝導計算
// ----------------------------------------
#include "zigzag.h"
typedef std::complex<double> Complex;
using namespace std;
using namespace Zigzag;

// メインプログラム
int main(int argc, char* argv[]) {
  // 物理量
  const double haf[4] = {0.00, 0.10, 0.20, 0.30}; // 反強磁性的な分子場
  double k;            // 波数
  static Zigzag::Zmatrix zmatrix;        // 複素数行列
  long info;                             // 対角化zheevの関数値

  // ファイル入出力関連
  stringstream filename;                 // ファイル名
  ofstream eigen_values[sizeof(haf) / sizeof(haf[0])]; // 出力ファイル

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  int i, band_num; // ループで必要な変数

  for(int haf_num = 0; haf_num < sizeof(haf) / sizeof(haf[0]); haf_num++) {

    // 出力ファイルを作成
    filename << fixed << setprecision(2) << "./data/h" << haf[haf_num] << ".d";
    eigen_values[haf_num].open(filename.str().c_str());
    filename.str("");                      // バッファのクリア
    filename.clear(stringstream::goodbit); // フラグのクリア

    // 最初の行のコメントを挿入
    eigen_values[haf_num] << "# wave_num  eigen_v[0]  eigen_v[1]" << endl;

    // -pi から pi まで波数を変えて各点で物理量を計算する(正常状態)
    for(i = 0; i <= N; i++) {
      k = 2.0 * N_inv * M_PI * (double)i - M_PI;

      // Hamiltonianを定義する
      zmatrix.init_zmatrix(0.0, haf[haf_num], 1.0, 0.0, k);

      // 対角化する
      info = zmatrix.zheev();

      // 温度と波数を出力する
      eigen_values[haf_num] << fixed << setprecision(8) << k;

      for(band_num = 0; band_num < DIM; band_num++) {
        // 固有エネルギーを出力する
        eigen_values[haf_num] << fixed << setprecision(8) << "  " << zmatrix.e[band_num];
      }
      eigen_values[haf_num] << endl;

    }

    eigen_values[haf_num].close();

  }

  return 0;
}
