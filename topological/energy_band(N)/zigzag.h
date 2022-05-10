//-*- coding: utf-8; mode: c++ -*-
// Time-stamp: <2016-02-03 14:21:26 shunta>

//
// ジグザグチェーンの行列関数ライブラリヘッダファイル
//

// 二重インクルードを防止する。
#ifndef ZIGZAG_H
#define ZIGZAG_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex>
#include <vector>
#include <sys/stat.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef std::complex<double> Complex;
using namespace std;

namespace Zigzag {

  // lapack
  extern "C" {
    void zheev_(const char& jobz, const char& uplo, const long& n, Complex** a, const long& lda, double* w, Complex* work, const long& lwork, double* rwork, long& info, long jobzlen, long uplolen);
  };

  // 各種固定パラメータ
  const long N = 2000;                          // 電子数
  const double N_inv = 0.0005;                  // 電子数の逆数
  const long DIM = 4;                           // Hamiltonianの次元
  const double V[2] = { 1.5, 1.5 };             // 電子間相互作用の強さ(副格子a, b)
  const double t[2] = { 0.1, 1.0 };             // 飛び移り積分(hoppingのパラメータ)
  const double alpha = 0.4;                     // 反対称スピン軌道相互作用の結合定数
  const double delta = 1.0e-2;                  // 収束条件

  //
  // 複素数行列クラス
  //
  class Zmatrix {
  public:
    Complex h[DIM][DIM];                        // Hamiltonian
    double e[DIM];                              // 固有値を入れる配列

    //
    // コンストラクタ
    //
    Zmatrix() {

      // すべての要素を0にする
      for(int i = 0; i < DIM; i++) {
        for(int j = 0; j < DIM; j++) {
          h[i][j] = 0;
        }
      }
    }

    //
    // コンストラクタ
    //   order: 秩序パラメータ
    //   haf: 反強磁性的な分子場
    //   mu: 化学ポテンシャル
    //   q: Cooper対の重心運動量
    //   k: 波数
    //
    Zmatrix(Complex order, double haf, double mu, double q, double k) {
      init_zmatrix(order, haf, mu, q, k);
    }

    //
    // 初期化
    //
    void init_zmatrix(Complex order, double haf, double mu, double q, double k) {
      // すべての要素を0にする
      for(int i = 0; i < DIM; i++) {
        for(int j = 0; j < DIM; j++) {
          h[i][j] = 0.0;
        }
      }

      h[0][0] = next_nearest_neighbour_hopping(k) + alpha * sin(k);
      h[1][0] = nearest_neighbour_hopping(k);
      h[1][1] = next_nearest_neighbour_hopping(k) - alpha * sin(k);
      h[2][1] = - haf;
      h[2][2] = next_nearest_neighbour_hopping(k) + alpha * sin(k);
      h[3][0] = - haf;
      h[3][2] = nearest_neighbour_hopping(k);
      h[3][3] = next_nearest_neighbour_hopping(k) - alpha * sin(k);
    }

    //
    // Hamiltonianを対角化する
    //
    long zheev(void) {
      long info;
#ifdef _OPENMP
      Complex work[4 * DIM];
      double rwork[4 * DIM];
#else
      static Complex work[4 * DIM];
      static double rwork[4 * DIM];
#endif

      zheev_('V', 'U', DIM, (Complex**)h, DIM, e, work, 4 * DIM, rwork, info, 1, 1);

      return info;
    }

    //
    // 異なる副格子間の最近接電子のホッピング
    //   k: 波数
    //
    double nearest_neighbour_hopping(double k) {
      return - 2.0 * t[0] * cos(k * 0.5);
    }

    //
    // 同じ副格子内の最近接電子のホッピング
    //   k: 波数
    //
    double next_nearest_neighbour_hopping(double k) {
      return - 2.0 * t[1] * cos(k);
    }
  };

  // フェルミ分布関数
  inline double fermi(double e, double beta);

  // 秩序パラメータを計算する
  Complex calculate_order_parameter(Complex order_o, double haf, double mu, double q, double T, double beta);

  // 自由エネルギーを計算する
  double calculate_free_energy(Complex order, double haf, double mu, double q, double T, double beta);

  // 物理量をファイルに出力する
  void output(ofstream& eigen_values, ofstream& physical_quantities, Complex order, double haf, double mu, double q, double T, double beta);
}

#endif // ZIGZAG_H
