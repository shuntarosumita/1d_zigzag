//-*- coding: utf-8; mode: c++ -*-
// Time-stamp: <2016-02-25 19:49:21 shunta>

//
// トポロジカル計算の行列関数ライブラリ
//

#include "zigzag.h"
#include "topological.h"
using namespace Zigzag;

//
// 行列 A の逆行列を計算する
//
Cmatrix Topological::invert(Cmatrix &A) {
  BOOST_UBLAS_CHECK(A.size1() == A.size2(), external_logic());

  Cmatrix A_copy(A);                       // A をコピーする
  Cmatrix inverse
    = identity_matrix<Complex>(A.size1()); // 出力する逆行列(単位行列で初期化)
  permutation_matrix<> pm(A.size1());      // LU分解の際の枢軸選択によって起こる行交換の情報が入る行列

  // LU分解を行う
  lu_factorize(A_copy, pm);

  // 行交換を考慮した前進消去と後退代入
  lu_substitute(A_copy, pm, inverse);

  BOOST_UBLAS_CHECK(inverse.size1() == A.size1(), internal_logic());
  BOOST_UBLAS_CHECK(inverse.size2() == A.size2(), internal_logic());

#if BOOST_UBLAS_TYPE_CHECK
  BOOST_UBLAS_CHECK(detail::expression_type_check(prod(A, inverse), identity_matrix<Cmatrix::value_type>(A.size1())), internal_logic());
#endif

  return inverse;

}

//
// 行列 A の固有値と固有ベクトルを求める
//
void Topological::eigen(Cmatrix &A, Cmatrix &U, Dvector &eigen) {
  Dvector eigen_copy(A.size1());  // 固有ベクトル
  matrix<Complex, column_major>
    A_copy(A.size1(), A.size2()); // A をコピーする
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

//
// ハミルトニアンを変換して非対角成分を得る
//
Cmatrix Topological::get_offdiagonal_Hamiltonian(Cmatrix &Hamiltonian, Cmatrix &U_Gamma) {
  size_t M = Hamiltonian.size1();
  Cmatrix H_offdiagonal(M, M);           // 非対角な形に変換されたハミルトニアン
  Cmatrix offdiagonal(M / 2, M / 2);     // 非対角化されたハミルトニアンの成分
  int i, j;                              // ループ用変数

  // カイラル演算子を対角化する基底を用いてハミルトニアンを非対角な形に変換
  // 行列の積 H_offdiagonal = U_Gamma * Hamiltonian
  H_offdiagonal = prod(U_Gamma, Hamiltonian);

  // 行列の積 H_offdiagonal = H_offdiagonal * U_Gamma^dagger
  H_offdiagonal = prod(H_offdiagonal, herm(U_Gamma));

  // 非対角成分だけを取り出す
  for(i = 0; i < M / 2; i++) {
    for(j = 0; j < M / 2; j++) {
      offdiagonal(i, j) = H_offdiagonal(i, j + M / 2);
    }
  }

  return offdiagonal;
}

//
// winding numberを計算する
//
double Topological::winding_num(Cvector &order, Dvector &t_input, double haf, double mu) {
  double k;                                   // 波数
  //// double dk = M_PI / (double)k_num;           // 波数の差
  double dk = 2.0 * M_PI / (double)k_num;     // 波数の差
  Zmatrix zmatrix;                            // 複素数行列
  Cmatrix U_Gamma
    = zero_matrix<Complex>(DIM, DIM);         // カイラル演算子を対角化するユニタリ行列
  Cmatrix prev_offdiagonal(DIM / 2, DIM / 2); // 非対角化されたハミルトニアンの成分(1つ前)
  Cmatrix offdiagonal(DIM / 2, DIM / 2);      // 非対角化されたハミルトニアンの成分
  Cmatrix calc(DIM / 2, DIM / 2);             // 計算過程で用いる行列
  int i, j;                                   // ループ用変数
  Complex wn = 0.0;                           // winding number

  // static変数に値を渡しておく
  t = t_input;

  // カイラル演算子を対角化するユニタリ行列 U_Gamma を定義
  U_Gamma(0, 0) = 1.0;
  U_Gamma(0, 5) = - imag_unit;
  U_Gamma(1, 1) = 1.0;
  U_Gamma(1, 4) = - imag_unit;
  U_Gamma(2, 2) = 1.0;
  U_Gamma(2, 7) = - imag_unit;
  U_Gamma(3, 3) = 1.0;
  U_Gamma(3, 6) = - imag_unit;
  U_Gamma(4, 1) = imag_unit;
  U_Gamma(4, 4) = - 1.0;
  U_Gamma(5, 0) = imag_unit;
  U_Gamma(5, 5) = - 1.0;
  U_Gamma(6, 3) = imag_unit;
  U_Gamma(6, 6) = - 1.0;
  U_Gamma(7, 2) = imag_unit;
  U_Gamma(7, 7) = - 1.0;

  U_Gamma /= sqrt(2.0);

  // ハミルトニアンを変換して非対角成分 q(-pi) を得る
  //// k = 0;
  k = - M_PI;
  zmatrix.init_zmatrix(order, haf, mu, 0.0, k);
  prev_offdiagonal = Topological::get_offdiagonal_Hamiltonian(zmatrix.h, U_Gamma);

  // 積分を計算するための波数のループ
  for(i = 1; i <= k_num; i++) {
    // 波数を定義
    //// k = dk * (double)i;
    k = dk * (double)i - M_PI;

    // ハミルトニアンを変換して非対角成分 q(k) を得る
    zmatrix.init_zmatrix(order, haf, mu, 0.0, k);
    offdiagonal = Topological::get_offdiagonal_Hamiltonian(zmatrix.h, U_Gamma);

    // 行列の積 q(k)^{-1} * ((q(k) - q(k - dk)) / dk) を計算
    //// calc = prod(Topological::invert(offdiagonal), (offdiagonal - prev_offdiagonal) / dk);
    calc = prod(Topological::invert(offdiagonal), (offdiagonal - prev_offdiagonal) / dk)
      - prod(herm(Topological::invert(offdiagonal)), herm(offdiagonal - prev_offdiagonal) / dk);

    // calc のトレースを計算
    for(j = 0; j < DIM / 2; j++) {
      wn += calc(j, j);
    }

    // 次のループのために1つずらした非対角成分を保存しておく
    prev_offdiagonal = offdiagonal;
  }

  //// wn /= imag_unit * (Complex)k_num;
  wn /= 2.0 * imag_unit * (Complex)k_num;

  return real(wn);

}

//
// winding numberを計算する(その2)
//
double Topological::winding_num2(Cvector &order, Dvector &t_input, double haf, double mu) {
  double k;                                   // 波数
  double dk = 2.0 * M_PI / (double)k_num;     // 波数の差
  Zmatrix zmatrix;                            // 複素数行列
  Cmatrix Gamma
    = zero_matrix<Complex>(DIM, DIM);         // カイラル演算子
  Cmatrix calc(DIM, DIM);                     // 計算過程で用いる行列
  Cmatrix Hamiltonian(DIM, DIM);              // 計算過程で用いる行列
  Cmatrix prev_Hamiltonian(DIM, DIM);         // 計算過程で用いる行列
  int i, j;                                   // ループ用変数
  Complex wn = 0.0;                           // winding number

  // static変数に値を渡しておく
  t = t_input;

  // カイラル演算子 Gamma を定義
  Gamma(0, 5) = - imag_unit;
  Gamma(1, 4) = - imag_unit;
  Gamma(2, 7) = - imag_unit;
  Gamma(3, 6) = - imag_unit;
  Gamma(4, 1) = imag_unit;
  Gamma(5, 0) = imag_unit;
  Gamma(6, 3) = imag_unit;
  Gamma(7, 2) = imag_unit;

  // ハミルトニアンを変換して非対角成分 q(-pi) を得る
  k = - M_PI;
  zmatrix.init_zmatrix(order, haf, mu, 0.0, k);
  prev_Hamiltonian = zmatrix.h;

  // 積分を計算するための波数のループ
  for(i = 1; i <= k_num; i++) {
    // 波数を定義
    k = dk * (double)i - M_PI;

    // ハミルトニアンを変換して非対角成分 q(k) を得る
    zmatrix.init_zmatrix(order, haf, mu, 0.0, k);
    Hamiltonian = zmatrix.h;

    // 行列の積 Gamma * H(k)^{-1} * (H(k) - H(k - dk)) を計算
    calc = prod(Gamma, invert(Hamiltonian));
    calc = prod(calc, Hamiltonian - prev_Hamiltonian);

    // calc のトレースを計算
    for(j = 0; j < DIM; j++) {
      wn += calc(j, j);
    }

    // 次のループのために1つずらした非対角成分を保存しておく
    prev_Hamiltonian = Hamiltonian;
  }

  wn /= - 4.0 * M_PI * imag_unit;

  return real(wn);

}

//
// ハミルトニアンをバンド表示するユニタリ行列を定義
//
Cmatrix Topological::get_U_band(double haf, double k) {
  Zmatrix zmatrix;
  Cmatrix U_band
    = zero_matrix<Complex>(DIM, DIM);
  double calc;
  double delta = 1.0e-7;
  int i, j;

  // 1列目
  U_band(0, 0) = sqrt(pow(zmatrix.nearest_neighbour_hopping(k) - haf, 2.0) + pow(alpha * sin(k), 2.0)) + alpha * sin(k);
  U_band(1, 0) = zmatrix.nearest_neighbour_hopping(k) - haf;
  U_band(2, 0) = U_band(0, 0);
  U_band(3, 0) = U_band(1, 0);
    
  // 2列目
  U_band(0, 1) = sqrt(pow(zmatrix.nearest_neighbour_hopping(k) + haf, 2.0) + pow(alpha * sin(k), 2.0)) + alpha * sin(k);
  U_band(1, 1) = zmatrix.nearest_neighbour_hopping(k) + haf;
  U_band(2, 1) = - U_band(0, 1);
  U_band(3, 1) = - U_band(1, 1);
    
  // 3列目
  U_band(0, 2) = sqrt(pow(zmatrix.nearest_neighbour_hopping(k) + haf, 2.0) + pow(alpha * sin(k), 2.0)) - alpha * sin(k);
  U_band(1, 2) = - (zmatrix.nearest_neighbour_hopping(k) + haf);
  U_band(2, 2) = - U_band(0, 2);
  U_band(3, 2) = - U_band(1, 2);
    
  // 4列目
  U_band(0, 3) = sqrt(pow(zmatrix.nearest_neighbour_hopping(k) - haf, 2.0) + pow(alpha * sin(k), 2.0)) - alpha * sin(k);
  U_band(1, 3) = - (zmatrix.nearest_neighbour_hopping(k) - haf);
  U_band(2, 3) = U_band(0, 3);
  U_band(3, 3) = U_band(1, 3);
    
  // 5列目
  U_band(4, 4) = - sqrt(pow(zmatrix.nearest_neighbour_hopping(- k) - haf, 2.0) + pow(alpha * sin(- k), 2.0)) - alpha * sin(- k);
  U_band(5, 4) = - zmatrix.nearest_neighbour_hopping(- k) + haf;
  U_band(6, 4) = U_band(4, 4);
  U_band(7, 4) = U_band(5, 4);

  // 6列目
  U_band(4, 5) = sqrt(pow(zmatrix.nearest_neighbour_hopping(- k) + haf, 2.0) + pow(alpha * sin(- k), 2.0)) + alpha * sin(- k);
  U_band(5, 5) = zmatrix.nearest_neighbour_hopping(- k) + haf;
  U_band(6, 5) = - U_band(4, 5);
  U_band(7, 5) = - U_band(5, 5);
    
  // 7列目
  U_band(4, 6) = sqrt(pow(zmatrix.nearest_neighbour_hopping(- k) + haf, 2.0) + pow(alpha * sin(- k), 2.0)) - alpha * sin(- k);
  U_band(5, 6) = - (zmatrix.nearest_neighbour_hopping(- k) + haf);
  U_band(6, 6) = - U_band(4, 6);
  U_band(7, 6) = - U_band(5, 6);
    
  // 8列目
  U_band(4, 7) = sqrt(pow(zmatrix.nearest_neighbour_hopping(- k) - haf, 2.0) + pow(alpha * sin(- k), 2.0)) - alpha * sin(- k);
  U_band(5, 7) = - (zmatrix.nearest_neighbour_hopping(- k) - haf);
  U_band(6, 7) = U_band(4, 7);
  U_band(7, 7) = U_band(5, 7);

  // 対角成分が正になるように規格化
  for(i = 0; i < DIM; i++) {
    calc = 0.0;
    for(j = 0; j < DIM; j++) {
      calc += pow(abs(U_band(j, i)), 2.0);
    }
    calc = sqrt(calc);

    if(calc >= delta) {
      for(j = 0; j < DIM; j++) {
        U_band(j, i) = U_band(j, i) / calc;
      }
    }

    if(real(U_band(i, i)) < 0.0) {
      for(j = 0; j < DIM; j++) {
        U_band(j, i) = - U_band(j, i);
      }
    }
  }

  return U_band;
}

//
// ハミルトニアンをバンド表示して波動関数を射影
//
void Topological::project_wavefunc(Cvector &order, Dvector &t_input, double haf, double mu) {
  double k;                                   // 波数
  double dk = 2.0 * M_PI / (double)k_num;     // 波数の差
  Zmatrix zmatrix;                            // 複素数行列
  Cmatrix zmatrix_half(DIM / 2, DIM / 2);     // ハミルトニアンの左上および右上
  Cmatrix eigen_vectors(DIM / 2, DIM / 2);    // 固有ベクトル
  Dvector eigen_values(DIM / 2);              // 固有値
  Cmatrix U_band(DIM, DIM);                   // バンド表示に変換するユニタリ行列
  int i, l;                                   // ループ用変数

  // ファイル入出力関連
  stringstream filename;                      // ファイル名
  ofstream output_file;                       // 出力ファイル

  // static変数に値を渡しておく
  t = t_input;

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/t" << t[0] << "_mu" << mu << "_h" << haf << ".d";
  output_file.open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 積分を計算するための波数のループ
  for(l = 0; l <= k_num; l++) {
    // 波数を定義
    k = dk * (double)l - M_PI;
    output_file << fixed << setprecision(8) << k;

    U_band = Topological::get_U_band(haf, k);

    // 波数 k のハミルトニアンを定義する
    zmatrix.init_zmatrix(order, haf, mu, 0.0, k);

    // バンド表示ハミルトニアンを用いて波動関数を射影
    zmatrix.h = prod(herm(U_band), zmatrix.h);
    zmatrix.h = prod(zmatrix.h, U_band);

    // テスト用出力
    // for(i = 0; i < DIM; i++) {
    //   for(j = 0; j < DIM; j++) {
    //     cout << fixed << setprecision(2) << zmatrix.h(i, j) << "  ";
    //   }
    //   cout << endl;
    // }
    // cout << endl;

    // ファイル出力
    for(i = 0; i < DIM / 2; i++) {
      output_file << fixed << setprecision(8) << "  " << real(zmatrix.h(i, i)) << "  " << real(zmatrix.h(i, i + DIM / 2)) 
                  << "  " << sqrt(pow(real(zmatrix.h(i, i)), 2.0) + pow(real(zmatrix.h(i, i + DIM / 2)), 2.0));
    }
    output_file << endl;
  }

  output_file.close();
}

//
// ハミルトニアンをバンド表示して各バンドのwinding numberを計算する
//
boost::numeric::ublas::vector<double> Topological::winding_num_band(Cvector &order, Dvector &t_input, double haf, double mu) {
  double k;                                   // 波数
  double dk = 2.0 * M_PI / (double)k_num;     // 波数の差
  Zmatrix zmatrix;                            // 複素数行列
  Cmatrix U_band(DIM, DIM);                   // バンド表示に変換するユニタリ行列
  Cvector prev_m1(DIM / 2);
  Cvector prev_m2(DIM / 2);
  Cvector m1(DIM / 2);
  Cvector m2(DIM / 2);
  Cvector wn_band
    = zero_vector<Complex>(DIM / 2);          // バンド表示したwinding number
  int i, l;                                   // ループ用変数

  // static変数に値を渡しておく
  t = t_input;

  // 波数 k = -pi のハミルトニアンを定義する
  k = - M_PI;
  zmatrix.init_zmatrix(order, haf, mu, 0.0, k);

  // ハミルトニアンをバンド表示するユニタリ行列を定義
  U_band = Topological::get_U_band(haf, k);

  // バンド表示ハミルトニアンを変換して非対角成分 q(-pi) を得る
  zmatrix.h = prod(herm(U_band), zmatrix.h);
  zmatrix.h = prod(zmatrix.h, U_band);
  for(i = 0; i < DIM / 2; i++) {
    prev_m1(i) = zmatrix.h(i, i) / sqrt(pow(abs(zmatrix.h(i, i)), 2.0) + pow(abs(zmatrix.h(i, i + DIM / 2)), 2.0));
    prev_m2(i) = zmatrix.h(i, i + DIM / 2) / sqrt(pow(abs(zmatrix.h(i, i)), 2.0) + pow(abs(zmatrix.h(i, i + DIM / 2)), 2.0));
  }

  // 積分を計算するための波数のループ
  for(l = 1; l <= k_num; l++) {
    // 波数を定義
    //// k = dk * (double)l;
    k = dk * (double)l - M_PI;

    // 波数 k のハミルトニアンを定義する
    zmatrix.init_zmatrix(order, haf, mu, 0.0, k);

    // ハミルトニアンをバンド表示するユニタリ行列を定義
    U_band = Topological::get_U_band(haf, k);

    // バンド表示ハミルトニアンを変換して非対角成分 q(k) を得る
    zmatrix.h = prod(herm(U_band), zmatrix.h);
    zmatrix.h = prod(zmatrix.h, U_band);
    for(i = 0; i < DIM / 2; i++) {
      m1(i) = zmatrix.h(i, i) / sqrt(pow(abs(zmatrix.h(i, i)), 2.0) + pow(abs(zmatrix.h(i, i + DIM / 2)), 2.0));
      m2(i) = zmatrix.h(i, i + DIM / 2) / sqrt(pow(abs(zmatrix.h(i, i)), 2.0) + pow(abs(zmatrix.h(i, i + DIM / 2)), 2.0));
    }

    // 積 q(k)^{-1} * ((q(k) - q(k - dk)) / dk) を計算
    for(i = 0; i < DIM / 2; i++) {
      wn_band(i) += m1(i) * (m2(i) - prev_m2(i)) / dk - m2(i) * (m1(i) - prev_m1(i)) / dk;
    }

    // 次のループのために1つずらした非対角成分を保存しておく
    prev_m1 = m1;
    prev_m2 = m2;
  }

  wn_band /= (Complex)k_num;

  boost::numeric::ublas::vector<double> wn_band_output(DIM / 2);
  for(i = 0; i < DIM / 2; i++) {
    wn_band_output(i) = real(wn_band(i));
  }
  
  return wn_band_output;

}

//
// Berry 位相を計算する
//
double Topological::Berry_phase(Cvector &order, Dvector &t_input, double haf, double mu) {
  double k;                                   // 波数
  double dk = 2.0 * M_PI / (double)k_num;     // 波数の差
  Zmatrix zmatrix;                            // 複素数行列
  Cmatrix prev_eigen_vectors(DIM, DIM);       // 固有ベクトル(1つ前)
  Cmatrix eigen_vectors(DIM, DIM);            // 固有ベクトル
  int i, j;                                   // ループ用変数
  Complex calc;                               // 計算用
  double bp = 0.0;                            // Berry phase

  // static変数に値を渡しておく
  t = t_input;

  // ハミルトニアンを対角化して固有ベクトル |u(-pi)> を得る
  zmatrix.init_zmatrix(order, haf, mu, 0.0, - M_PI);
  zmatrix.zheev();
  prev_eigen_vectors = zmatrix.h;

  // 積分を計算するための波数のループ
  for(i = 1; i <= k_num; i++) {
    // 波数を定義
    k = dk * (double)i - M_PI;

    // ハミルトニアンを対角化して固有ベクトル |u(k)> を得る
    zmatrix.init_zmatrix(order, haf, mu, 0.0, k);
    zmatrix.zheev();
    eigen_vectors = zmatrix.h;

    // 
    Cmatrix diag = identity_matrix<Complex>(DIM);
    for(j = DIM / 2; j < DIM; j++) {
      calc = inner_prod(column(eigen_vectors, j), column(prev_eigen_vectors, j));
      if (real(calc) < 0.0) {
        diag(j, j) = - 1.0;
      }
    }
    eigen_vectors = prod(eigen_vectors, diag);

    for(j = DIM / 2; j < DIM; j++) {
      bp += imag(inner_prod(column(eigen_vectors, j), column(prev_eigen_vectors, j)));
    }

    // 次のループのために1つずらした固有ベクトルを保存しておく
    prev_eigen_vectors = eigen_vectors;
  }

  bp /= M_PI;

  return bp;

}

// //
// // Berry 位相を計算する
// //
// double Topological::Berry_phase(Cvector &order, Dvector &t_input, double haf, double mu) {
//   double k;                                   // 波数
//   double dk = 2.0 * M_PI / (double)k_num;     // 波数の差
//   Zmatrix zmatrix;                            // 複素数行列
//   Cmatrix prev_eigen_vectors(DIM, DIM);       // 固有ベクトル(1つ前)
//   Cmatrix eigen_vectors(DIM, DIM);            // 固有ベクトル
//   int i, j;                                   // ループ用変数
//   Complex calc;                               // 計算用
//   Complex bp = 0.0;                           // Berry phase

//   // static変数に値を渡しておく
//   t = t_input;

//   // ハミルトニアンを対角化して固有ベクトル |u(-pi)> を得る
//   zmatrix.init_zmatrix(order, haf, mu, 0.0, - M_PI);
//   zmatrix.zheev();
//   prev_eigen_vectors = zmatrix.h;

//   // 積分を計算するための波数のループ
//   for(i = 1; i <= k_num; i++) {
//     // 波数を定義
//     k = dk * (double)i - M_PI;

//     // ハミルトニアンを対角化して固有ベクトル |u(k)> を得る
//     zmatrix.init_zmatrix(order, haf, mu, 0.0, k);
//     zmatrix.zheev();
//     eigen_vectors = zmatrix.h;

//     // 
//     Cmatrix diag = identity_matrix<Complex>(DIM);

//     for(j = DIM / 2; j < DIM; j++) {
//       calc = inner_prod(column(eigen_vectors, j), column(prev_eigen_vectors, j));
//       if (real(calc) < 0.0) {
//         diag(j, j) = - 1.0;
//       }
//     }
//     eigen_vectors = prod(eigen_vectors, diag);

//     // 次のループのために1つずらした固有ベクトルを保存しておく
//     prev_eigen_vectors = eigen_vectors;
//   }

//   // ハミルトニアンを対角化して固有ベクトル |u(-pi)> を得る
//   zmatrix.init_zmatrix(order, haf, mu, 0.0, - M_PI);
//   zmatrix.zheev();
//   prev_eigen_vectors = zmatrix.h;

//   for(j = DIM / 2; j < DIM; j++) {
//     cout << inner_prod(column(eigen_vectors, j), column(prev_eigen_vectors, j)) << endl;
//     bp += inner_prod(column(eigen_vectors, j), column(prev_eigen_vectors, j));
//   }

//   return real(bp);

// }
