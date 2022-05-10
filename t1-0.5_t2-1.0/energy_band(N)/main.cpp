//-*- coding: utf-8 -*-
// Time-stamp: <2016-01-06 13:41:56 shunta>
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
  const double haf[1] = {0.20}; // 反強磁性的な分子場
  double k, k_fermi[2][6][2];            // 波数
  double mu_fermi;                       // 化学ポテンシャル
  static Zigzag::Zmatrix zmatrix;        // 複素数行列
  long info;                             // 対角化zheevの関数値
  int max_enecount = 6000;
  double dos[max_enecount][2];           // 状態密度
  double epsilon = 0.01;                 // Lorentzianのパラメータ

  // ファイル入出力関連
  stringstream filename;                 // ファイル名
  ofstream eigen_values[sizeof(haf) / sizeof(haf[0])]; // 出力ファイル
  ofstream fermi[2];                                   // 出力ファイル
  bool in_fermi_surface[2][6];           // フェルミ面の内側かどうかを判断する変数
  ofstream dos_file[sizeof(haf) / sizeof(haf[0])];

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);
  mkdir("./data/eigen_values", 0755);
  mkdir("./data/fermi", 0755);
  mkdir("./data/dos", 0755);

  int i, band_num, mu_num, enecount; // ループで必要な変数

  // 出力ファイルを作成
  fermi[0].open("./data/fermi/lower_band.d");
  fermi[1].open("./data/fermi/upper_band.d");

  // 最初の行のコメントを挿入
  fermi[0] << "# haf  mu  k_fermi[0]  k_fermi[1]  k_mid" << endl;
  fermi[1] << "# haf  mu  k_fermi[0]  k_fermi[1]  k_mid" << endl;

  for(int haf_num = 0; haf_num < sizeof(haf) / sizeof(haf[0]); haf_num++) {
    // 変数の初期化
    for(mu_num = 0; mu_num < 6; mu_num++) {
      in_fermi_surface[0][mu_num] = false;
      in_fermi_surface[1][mu_num] = false;

      for(band_num = 0; band_num <= DIM - 1; band_num++) {
        k_fermi[band_num][mu_num][0] = 10.0;
        k_fermi[band_num][mu_num][1] = 10.0;
      }
    }

    // 出力ファイルを作成
    filename << fixed << setprecision(2) << "./data/eigen_values/h" << haf[haf_num] << ".d";
    eigen_values[haf_num].open(filename.str().c_str());
    filename.str("");                      // バッファのクリア
    filename.clear(stringstream::goodbit); // フラグのクリア

    filename << fixed << setprecision(2) << "./data/dos/h" << haf[haf_num] << ".d";
    dos_file[haf_num].open(filename.str().c_str());
    filename.str("");                      // バッファのクリア
    filename.clear(stringstream::goodbit); // フラグのクリア

    // 最初の行のコメントを挿入
    eigen_values[haf_num] << "# wave_num  eigen_v[0]  eigen_v[1]" << endl;
    dos_file[haf_num] << "# energy DOS" << endl;

    for(enecount = 0; enecount < max_enecount; enecount++) {
      dos[enecount][0] = 6.0 * (double)enecount / (double)max_enecount - 3.0;
    }

    // -pi から pi まで波数を変えて各点で物理量を計算する(正常状態)
    for(i = 0; i <= N; i++) {
      k = 2.0 * N_inv * M_PI * (double)i - M_PI;

      // Hamiltonianを定義する
      zmatrix.init_zmatrix(0.0, haf[haf_num], 1.0, 0.0, k);

      // 対角化する
      info = zmatrix.zheev();

      // 温度と波数を出力する
      eigen_values[haf_num] << fixed << setprecision(8) << k;

      for(band_num = 0; band_num <= DIM - 1; band_num++) {
        // 固有エネルギーを出力する
        eigen_values[haf_num] << fixed << setprecision(8) << "  " << zmatrix.e[band_num];

        // フェルミ波数を求める
        for(mu_num = 0; mu_num < 6; mu_num++) {
          mu_fermi = (double)mu_num - 3.0;
          if(zmatrix.e[band_num] - mu_fermi <= 0.0 && in_fermi_surface[band_num][mu_num] == false) {
            k_fermi[band_num][mu_num][0] = k;
            in_fermi_surface[band_num][mu_num] = true;
          }
          else if(zmatrix.e[band_num] - mu_fermi > 0.0 && in_fermi_surface[band_num][mu_num] == true) {
            k_fermi[band_num][mu_num][1] = k - 2.0 * N_inv * M_PI;
            in_fermi_surface[band_num][mu_num] = false;
          }
        }
      }
      eigen_values[haf_num] << endl;

      // 状態密度を計算する
      for(int j = 0; j <= DIM / 2 - 1; j++) {
        for(enecount = 0; enecount < max_enecount; enecount++) {
          dos[enecount][1] += N_inv / M_PI * epsilon / (pow(dos[enecount][0] - zmatrix.e[2 * j], 2.0) + pow(epsilon, 2.0));
        }
      }

    }

    // フェルミ波数を出力する
    for(band_num = 0; band_num <= DIM - 1; band_num++) {
      for(mu_num = 0; mu_num < 6; mu_num++) {
        mu_fermi = (double)mu_num - 3.0;
        if(k_fermi[band_num][mu_num][0] != 10.0 && k_fermi[band_num][mu_num][1] != 10.0) {
          fermi[band_num] << fixed << setprecision(8) << haf[haf_num] << "  " << mu_fermi << "  " << k_fermi[band_num][mu_num][0] << "  " << k_fermi[band_num][mu_num][1] << "  " << k_fermi[band_num][mu_num][0] + k_fermi[band_num][mu_num][1] << endl;
        }
      }
    }

    // 状態密度を出力する
    for(enecount = 0; enecount < max_enecount; enecount++) {
      dos_file[haf_num] << fixed << setprecision(8) << dos[enecount][0] << "  " << dos[enecount][1] << endl;
    }

    eigen_values[haf_num].close();
    dos_file[haf_num].close();

  }


  fermi[0].close();
  fermi[1].close();

  return 0;
}
