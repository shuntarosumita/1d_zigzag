//-*- coding: utf-8; mode: c++ -*-
// Time-stamp: <2015-11-02 13:07:13 shunta>

//
// ジグザグチェーンの行列関数ライブラリ
//

#include "zigzag.h"

//
// dInをunitで丸める
//
double Zigzag::Round(double dIn, double unit) {
  double dOut;

  dOut = dIn / unit;
  dOut = (double)(int)(dOut + 0.5);

  return dOut * unit;
}

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
Complex Zigzag::calculate_order_parameter(Complex order_o, double haf, double mu, double q, double T, double beta) {
  static Zigzag::Zmatrix zmatrix;            // 複素数行列
  Complex order_n = 0.0;                     // 秩序パラメータ(固有ベクトルから計算する)
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
      order_n += - 0.5 * N_inv * V[0] * (conj(zmatrix.h[j][7]) * zmatrix.h[j][0] - conj(zmatrix.h[j][4]) * zmatrix.h[j][3]) * fermi(zmatrix.e[j], beta);
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
double Zigzag::calculate_free_energy(Complex order, double haf, double mu, double q, double T, double beta) {
  double k;                                  // 波数
  static Zigzag::Zmatrix zmatrix;            // 複素数行列
  double omega;                              // 自由エネルギー
  int i, j;                                  // ループ用変数
  long info;                                 // 対角化zheevの関数値

  omega = 2 * N / V[0] * pow(abs(order), 2.0);

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
void Zigzag::output(ofstream& eigen_values, ofstream& physical_quantities, Complex order, double haf, double mu, double q, double T, double beta) {
  double k, wave_num;                                 // 波数(wave_numは実際に出力する波数)
  static Zigzag::Zmatrix zmatrix;                     // 複素数行列
  double E = 2 * N / V[0] * pow(abs(order), 2.0);     // エネルギーの期待値
  double omega = 2 * N / V[0] * pow(abs(order), 2.0); // 自由エネルギー
  double omega_n = 0.0;                               // 自由エネルギー(正常状態)
  double S = 0.0;                                     // エントロピー
  int i, j;                                           // ループ用変数
  long info;                                          // 対角化zheevの関数値

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
    zmatrix.init_zmatrix(0.0, haf, mu, q, k);

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
  physical_quantities << fixed << setprecision(8) << T << "  " << q << "  "<< abs(order)
                      << "  " << E << "  " << S << "  " << omega << "  " << omega_n << endl;

}
