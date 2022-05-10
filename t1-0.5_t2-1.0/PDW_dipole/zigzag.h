//-*- coding: utf-8; mode: c++ -*-
// Time-stamp: <2015-12-14 14:31:12 shunta>

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
  const long DIM = 8;                           // Hamiltonianの次元
  const double V[2] = { 1.5, 1.5 };             // 電子間相互作用の強さ(副格子a, b)
  const vector<Complex> order_zero(2, Complex(0.0, 0.0)); // ゼロギャップ
  const double t[2] = { 0.5, 1.0 };             // 飛び移り積分(hoppingのパラメータ)
  const double alpha = 0.4;                     // 反対称スピン軌道相互作用の結合定数
  const double delta = 1.0e-2;                  // 収束条件
  const Complex imag_unit(0.0, 1.0);            // 虚数単位

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
    Zmatrix(vector<Complex> order, double haf, double mu, double q, double k) {
      init_zmatrix(order, haf, mu, q, k);
    }

    //
    // 初期化
    //
    void init_zmatrix(vector<Complex> order, double haf, double mu, double q, double k) {
      double nn_hopping[2] = { next_nearest_neighbour_hopping(k + q * 0.5), next_nearest_neighbour_hopping(- k + q * 0.5) };
      double as[2] = { alpha * sin(k + q * 0.5), alpha * sin(- k + q * 0.5) };
      double off_diag[2] = { nearest_neighbour_hopping(k + q * 0.5), nearest_neighbour_hopping(- k + q * 0.5) };

      // すべての要素を0にする
      for(int i = 0; i < DIM; i++) {
        for(int j = 0; j < DIM; j++) {
          h[i][j] = 0;
        }
      }

      // 1列目
      h[0][0] = nn_hopping[0] + as[0] - mu;

      // 2列目
      h[1][0] = off_diag[0];
      h[1][1] = nn_hopping[0] - as[0] - mu;

      // 3列目
      h[2][1] = - haf;
      // h[2][1] = imag_unit * haf;
      h[2][2] = nn_hopping[0] + as[0] - mu;

      // 4列目
      h[3][0] = - haf;
      // h[3][0] = imag_unit * haf;
      h[3][2] = off_diag[0];
      h[3][3] = nn_hopping[0] - as[0] - mu;

      // 5列目
      h[4][3] = - order[0];
      h[4][4] = - nn_hopping[1] - as[1] + mu;

      // 6列目
      h[5][2] = - order[1];
      h[5][4] = - off_diag[1];
      h[5][5] = - nn_hopping[1] + as[1] + mu;

      // 7列目
      h[6][1] = order[1];
      h[6][5] = haf;
      // h[6][5] = imag_unit * haf;
      h[6][6] = - nn_hopping[1] - as[1] + mu;

      // 8列目
      h[7][0] = order[0];
      h[7][4] = haf;
      // h[7][4] = imag_unit * haf;
      h[7][6] = - off_diag[1];
      h[7][7] = - nn_hopping[1] + as[1] + mu;
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
  vector<Complex> calculate_order_parameter(vector<Complex> order_o, double haf, double mu, double q, double T, double beta);

  // 自由エネルギーを計算する
  double calculate_free_energy(vector<Complex> order, double haf, double mu, double q, double T, double beta);

  // 物理量をファイルに出力する
  void output(ofstream& eigen_values, ofstream& physical_quantities, vector<Complex> order, double haf, double mu, double q, double T, double beta);
}

#endif // ZIGZAG_H
