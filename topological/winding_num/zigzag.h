//-*- coding: utf-8; mode: c++ -*-
// Time-stamp: <2016-02-18 17:35:46 shunta>

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
//#include <vector>
#include <sys/stat.h>

// Boost
#include <boost/numeric/ublas/vector.hpp>     // ベクトル用のヘッダ
#include <boost/numeric/ublas/matrix.hpp>     // 普通の行列用のヘッダ
#include <boost/numeric/ublas/triangular.hpp> // 三角行列用のヘッダ．前進消去，後退代入に必要
#include <boost/numeric/ublas/lu.hpp>         // LU分解，前進消去，後退代入用のヘッダ
#include <boost/numeric/bindings/lapack/heev.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/io.hpp>         // ストリーム入出力用のヘッダ

// OPENMP
#ifdef _OPENMP
#include <omp.h>
#endif

typedef std::complex<double> Complex;
typedef boost::numeric::ublas::matrix<Complex> Cmatrix;
typedef boost::numeric::ublas::vector<Complex> Cvector;
typedef std::vector<double> Dvector;
using namespace std;
using namespace boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

namespace Zigzag {

  // 各種固定パラメータ
  const long N = 200;                          // 電子数
  const double N_inv = 0.005;                  // 電子数の逆数
  const long DIM = 8;                           // Hamiltonianの次元
  const double V[2] = { 1.5, 1.5 };             // 電子間相互作用の強さ(副格子a, b)
  static Dvector t(2);                          // 飛び移り積分(hoppingのパラメータ)
  const double alpha = 0.4;                     // 反対称スピン軌道相互作用の結合定数
  const double delta = 1.0e-2;                  // 収束条件
  const Complex imag_unit(0.0, 1.0);            // 虚数単位

  //
  // 複素数行列クラス
  //
  class Zmatrix {
  public:
    Cmatrix h;                                  // Hamiltonian
    Dvector e;                                  // 固有値を入れる配列

    //
    // コンストラクタ
    //
    Zmatrix() {
      h = zero_matrix<Complex>(DIM, DIM);
      e = Dvector(DIM, 0.0);
    }

    //
    // コンストラクタ
    //   order: 秩序パラメータ
    //   haf: 反強磁性的な分子場
    //   mu: 化学ポテンシャル
    //   q: Cooper対の重心運動量
    //   k: 波数
    //
    Zmatrix(Cvector &order, double haf, double mu, double q, double k) {
      init_zmatrix(order, haf, mu, q, k);
    }

    //
    // 初期化
    //
    void init_zmatrix(Cvector &order, double haf, double mu, double q, double k) {
      double nn_hopping[2] = { next_nearest_neighbour_hopping(k + q * 0.5), next_nearest_neighbour_hopping(- k + q * 0.5) };
      double as[2] = { alpha * sin(k + q * 0.5), alpha * sin(- k + q * 0.5) };
      Complex off_diag[2] = { nearest_neighbour_hopping(k + q * 0.5), conj(nearest_neighbour_hopping(- k + q * 0.5)) };

      h = zero_matrix<Complex>(DIM, DIM);
      e = Dvector(DIM, 0.0);

      // 1行目
      h(0, 0) = nn_hopping[0]  + as[0] - mu;
      h(0, 1) = off_diag[0];
      h(0, 3) = - haf;
      //h(0, 3) = imag_unit * haf;
      h(0, 7) = order(0);

      // 2行目
      h(1, 1) = nn_hopping[0] - as[0] - mu;
      h(1, 2) = - haf;
      //h(1, 2) = imag_unit * haf;
      h(1, 6) = order(1);

      // 3行目
      h(2, 2) = nn_hopping[0] + as[0] - mu;
      h(2, 3) = off_diag[0];
      h(2, 5) = - order(1);

      // 4行目
      h(3, 3) = nn_hopping[0] - as[0] - mu;
      h(3, 4) = - order(0);

      // 5行目
      h(4, 4) = - nn_hopping[1] - as[1] + mu;
      h(4, 5) = - off_diag[1];
      h(4, 7) = haf;
      //h(4, 7) = imag_unit * haf;

      // 6行目
      h(5, 5) = - nn_hopping[1] + as[1] + mu;
      h(5, 6) = haf;
      //h(5, 6) = imag_unit * haf;

      // 7行目
      h(6, 6) = - nn_hopping[1] - as[1] + mu;
      h(6, 7) = - off_diag[1];

      // 8行目
      h(7, 7) = - nn_hopping[1] + as[1] + mu;

      for(int i = 0; i < h.size1(); i++) {
        for(int j = 0; j < i; j++) {
          h(i, j) = conj(h(j, i));
        }
      }
    }

    //
    // Hamiltonianを対角化する
    //
    long zheev(void) {
      Dvector e_copy(h.size1());
      matrix<Complex, column_major> h_copy(h.size1(), h.size2());
      long info;

      for(size_t i = 0; i < h.size1(); i++) {
        for(size_t j = i; j < h.size2(); j++) {
          h_copy(i, j) = h(i, j);
        }
      }

      info = lapack::heev('V', 'U', h_copy, e_copy, lapack::optimal_workspace());
      BOOST_UBLAS_CHECK(info == 0, internal_logic());

      h = h_copy;
      e.swap(e_copy);

      return info;
    }

    //
    // 異なる副格子間の最近接電子のホッピング
    //   k: 波数
    //
    Complex nearest_neighbour_hopping(double k) {
      // return - 2.0 * t[0] * cos(k * 0.5);
      return - 2.0 * t[0] * cos(k * 0.5) * exp(imag_unit * k * 0.5);
    }

    //
    // 同じ副格子内の最近接電子のホッピング
    //   k: 波数
    //
    double next_nearest_neighbour_hopping(double k) {
      return - 2.0 * t[1] * cos(k);
    }
  };

  // // フェルミ分布関数
  // inline double fermi(double e, double beta);

  // // 秩序パラメータを計算する
  // Cvector calculate_order_parameter(Cvector order_o, double haf, double mu, double q, double T, double beta);

  // // 自由エネルギーを計算する
  // double calculate_free_energy(Cvector order, double haf, double mu, double q, double T, double beta);

  // // 物理量をファイルに出力する
  // void output(ofstream& eigen_values, ofstream& physical_quantities, Cvector order, double haf, double mu, double q, double T, double beta);
}

#endif // ZIGZAG_H
