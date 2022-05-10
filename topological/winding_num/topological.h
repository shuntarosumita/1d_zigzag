//-*- coding: utf-8; mode: c++ -*-
// Time-stamp: <2016-02-25 19:34:39 shunta>

//
// トポロジカル計算の行列関数ライブラリヘッダファイル
//

// 二重インクルードを防止する。
#ifndef TOPOLOGICAL_H
#define TOPOLOGICAL_H

#include "zigzag.h"

typedef std::complex<double> Complex;
typedef boost::numeric::ublas::matrix<Complex> Cmatrix;
typedef boost::numeric::ublas::vector<Complex> Cvector;
using namespace std;
using namespace Zigzag;
using namespace boost::numeric::ublas;

namespace Topological {

  // 各種固定パラメータ
  const long k_num = 10000;                      // 波数点の数

  // M * M 行列の逆行列を計算する
  Cmatrix invert(Cmatrix &A);

  // 行列 A の固有値と固有ベクトルを求める
  void eigen(Cmatrix &A, Cmatrix &U, Dvector &eigen);

  // ハミルトニアンを変換して非対角成分を得る
  Cmatrix get_offdiagonal_Hamiltonian(Cmatrix &Hamiltonian, Cmatrix &U_Gamma);

  // winding numberを計算する
  double winding_num(Cvector &order, Dvector &t_input, double haf, double mu);

  // winding numberを計算する(その2)
  double winding_num2(Cvector &order, Dvector &t_input, double haf, double mu);

  // ハミルトニアンをバンド表示するユニタリ行列を定義
  Cmatrix get_U_band(double haf, double k);

  // ハミルトニアンをバンド表示して波動関数を射影
  void project_wavefunc(Cvector &order, Dvector &t_input, double haf, double mu);

  // ハミルトニアンをバンド表示して各バンドのwinding numberを計算する
  boost::numeric::ublas::vector<double> winding_num_band(Cvector &order, Dvector &t_input, double haf, double mu);

  // Berry 位相を計算する
  double Berry_phase(Cvector &order, Dvector &t_input, double haf, double mu);
}

#endif // TOPOLOGICAL_H
