//-*- coding: utf-8; mode: c++ -*-
// Time-stamp: <2015-11-09 20:10:37 shunta>

//
// ジグザグチェーンの行列関数ライブラリ
//

#include "zigzag.h"

//
// フェルミ分布関数
//
inline double Zigzag::fermi(double e, double beta) {
  return 1.0 / (exp(beta * e) + 1.0);
}

//
// 秩序パラメータを計算する
//   order_o: 秩序パラメータ
//   haf: 反強磁性的な分子場
//   mu: 化学ポテンシャル
//   q: Cooper対の重心運動量
//   T: 温度
//   beta: 逆温度
//
vector<Complex> Zigzag::calculate_order_parameter(vector<Complex> order_o, double haf, double mu, double q, double T, double beta) {
#ifdef _OPENMP
  Zigzag::Zmatrix zmatrix;                   // 複素数行列
#else
  static Zigzag::Zmatrix zmatrix;            // 複素数行列
#endif
  vector<Complex> order_n(2, 0.0);           // 秩序パラメータ(固有ベクトルから計算する)
  int i, j;                                  // ループ用変数
  long info;                                 // 対角化zheevの関数値

  // 電子数のループ
  for(i = 0; i <= N; i++) {

    // Hamiltonianを定義する
    zmatrix.init_zmatrix(order_o, haf, mu, q, 2.0 * N_inv * M_PI * (double)i - M_PI);
    
    // 対角化する
    info = zmatrix.zheev();

    // 固有ベクトルをもとに秩序パラメータを計算する
    for(j = 0; j < DIM; j++) {
      order_n[0] += - 0.5 * N_inv * V[0] * (conj(zmatrix.h[j][7]) * zmatrix.h[j][0] - conj(zmatrix.h[j][4]) * zmatrix.h[j][3]) * fermi(zmatrix.e[j], beta);
      order_n[1] += - 0.5 * N_inv * V[1] * (conj(zmatrix.h[j][6]) * zmatrix.h[j][1] - conj(zmatrix.h[j][5]) * zmatrix.h[j][2]) * fermi(zmatrix.e[j], beta);
    }
  }

  // 計算した固有ベクトルを返す。
  return order_n;
}

//
// 自由エネルギーを計算する
//   order: 秩序パラメータ
//   haf: 反強磁性的な分子場
//   mu: 化学ポテンシャル
//   q: Cooper対の重心運動量
//   T: 温度
//   beta: 逆温度
//
double Zigzag::calculate_free_energy(vector<Complex> order, double haf, double mu, double q, double T, double beta) {
  double k;                                  // 波数
#ifdef _OPENMP
  Zigzag::Zmatrix zmatrix;                   // 複素数行列
#else
  static Zigzag::Zmatrix zmatrix;            // 複素数行列
#endif
  double omega;                              // 自由エネルギー
  int i, j;                                  // ループ用変数
  long info;                                 // 対角化zheevの関数値

  omega = (double)N / V[0] * pow(abs(order[0]), 2.0) + (double)N / V[1] * pow(abs(order[1]), 2.0);

  // 電子数のループ
  for(i = 0; i <= N; i++) {
    k = 2.0 * N_inv * M_PI * (double)i - M_PI;

    // Hamiltonianを定義する
    zmatrix.init_zmatrix(order, haf, mu, q, k);

    // 対角化する
    info = zmatrix.zheev();

    // 自由エネルギーを計算する
    for(j = 4; j <= DIM - 1; j++) {
      omega += - T * log(1.0 + exp(- beta * zmatrix.e[j])) - 0.5 * zmatrix.e[j];
    }
    omega += - 2.0 * zmatrix.next_nearest_neighbour_hopping(- k + q * 0.5);
  }

  // 計算した自由エネルギーの1サイトあたりの物理量を返す。
  return omega * N_inv;
}

//
// 物理量をファイルに出力する
//   order: 秩序パラメータ
//   haf: 反強磁性的な分子場
//   mu: 化学ポテンシャル
//   q: Cooper対の重心運動量
//   T: 温度
//   beta: 逆温度
//
void Zigzag::output(ofstream& eigen_values, ofstream& physical_quantities, vector<Complex> order, double haf, double mu, double q, double T, double beta) {
  double k, wave_num;                                 // 波数(wave_numは実際に出力する波数)
#ifdef _OPENMP
  Zigzag::Zmatrix zmatrix;                            // 複素数行列
#else
  static Zigzag::Zmatrix zmatrix;                     // 複素数行列
#endif
  double E;                                           // エネルギーの期待値
  double omega;                                       // 自由エネルギー
  double omega_n = 0.0;                               // 自由エネルギー(正常状態)
  double S = 0.0;                                     // エントロピー
  int i, j;                                           // ループ用変数
  long info;                                          // 対角化zheevの関数値

  E = (double)N / V[0] * pow(abs(order[0]), 2.0) + (double)N / V[1] * pow(abs(order[1]), 2.0);
  omega = (double)N / V[0] * pow(abs(order[0]), 2.0) + (double)N / V[1] * pow(abs(order[1]), 2.0);

  // -pi から pi まで波数を変えて各点で物理量を計算する(超伝導状態)
  for(i = 0; i <= N; i++) {
    k = 2.0 * N_inv * M_PI * (double)i - M_PI;

    // Hamiltonianを定義する
    zmatrix.init_zmatrix(order, haf, mu, q, k);

    // 対角化する
    info = zmatrix.zheev();

    // 波数を設定する(piを超えたときは折り返す)
    if(k + q * 0.5 <= M_PI) {
      wave_num = k + q * 0.5;
    }
    else {
      wave_num = k + q * 0.5 - 2.0 * M_PI;
    }

    // 温度と波数を出力する
    eigen_values << fixed << setprecision(8) << T << "  " << wave_num;

    // 固有エネルギーを出力する
    for(j = 0; j <= DIM - 1; j++) {
      eigen_values << fixed << setprecision(8) << "  " << zmatrix.e[j];
    }
    eigen_values << endl;

    // エネルギーの期待値と自由エネルギーを計算する
    for(j = 4; j <= DIM - 1; j++) {
      E += zmatrix.e[j] * (fermi(zmatrix.e[j], beta) - 0.5);
      omega += - T * log(1.0 + exp(- beta * zmatrix.e[j])) - 0.5 * zmatrix.e[j];
    }
    E += - 2.0 * zmatrix.next_nearest_neighbour_hopping(- k + q * 0.5);
    omega += - 2.0 * zmatrix.next_nearest_neighbour_hopping(- k + q * 0.5);

  }

  // -pi から pi まで波数を変えて各点で物理量を計算する(正常状態)
  for(i = 0; i <= N; i++) {
    k = 2.0 * N_inv * M_PI * (double)i - M_PI;

    // Hamiltonianを定義する(秩序パラメータを0とする)
    zmatrix.init_zmatrix(order_zero, haf, mu, q, k);

    // 対角化する
    info = zmatrix.zheev();

    // 自由エネルギーを計算する
    for(j = 4; j <= DIM - 1; j++) {
      omega_n += - T * log(1.0 + exp(- beta * zmatrix.e[j])) - 0.5 * zmatrix.e[j];
    }
    omega_n += - 2.0 * zmatrix.next_nearest_neighbour_hopping(- k + q * 0.5);
  }

  // 1サイトあたりの物理量にする
  E *= N_inv;
  omega *= N_inv;
  omega_n *= N_inv;

  // エントロピーを定義する
  S = beta * (E - omega);

  // 各種物理量を出力する
  physical_quantities << fixed << setprecision(8) << T << "  " << q << "  "<< real(order[0]) << "  " << imag(order[0]) << "  "
                      << real(order[1]) << "  " << imag(order[0]) << "  " << E << "  " << S << "  " << omega << "  " << omega_n << endl;

}
