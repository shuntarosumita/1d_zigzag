//-*- coding: utf-8 -*-
// Time-stamp: <2016-02-19 16:31:04 shunta>
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
  double haf;                            // 反強磁性的な分子場
  double mu;                             // 化学ポテンシャル
  vector< vector<Complex> > order(3, vector<Complex>(2)); // 秩序パラメータ(Steffensenの反復法を使うため3つ定義)
  double q, q_candidate, q_opti;         // Cooper対の重心運動量
  double T;                              // 温度
  double beta;                           // 逆温度
  double omega, omega_n, cond_opti;      // 自由エネルギー

  // 秩序パラメータの初期化
  vector<Complex> order_init;
  order_init.push_back(Complex(1.0, 0.0));
  order_init.push_back(Complex(1.0, 0.0));

  // 最小2乗法に必要なパラメータ
  double one;
  double x, x2, x3, x4;
  double y, xy, x2y;

  // ファイル入出力関連
  stringstream filename;                 // ファイル名
  ofstream output_file[3];               // 出力ファイル

  // 入力されたパラメータを設定
  if(argc > 2) {
    mu = atof(argv[1]);
    if(mu <= -3.5 || mu > 2.5) {
      cerr << "mu range: (-3.5, 2.5]" << endl;
      return 0;
    }

    haf = atof(argv[2]);
    if(haf < 0.0 || haf > 4.0) {
      cerr << "h^AF range: [0.0, 4.0]" << endl;
      return 0;
    }
  }
  else {
    cerr << "Usage: " << argv[0] << " mu h^AF" << endl;
    return 0;
  }

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);
  filename << fixed << setprecision(2) << "./data/mu" << mu;
  mkdir(filename.str().c_str(), 0755);
  mkdir((filename.str() + "/eigen_values").c_str(), 0755);
  mkdir((filename.str() + "/physical_quantities").c_str(), 0755);
  mkdir((filename.str() + "/phase_diagram").c_str(), 0755);
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/mu" << mu << "/eigen_values/FFLO_h" << haf << ".d";
  output_file[0].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  filename << fixed << setprecision(2) << "./data/mu" << mu << "/physical_quantities/FFLO_h" << haf << ".d";
  output_file[1].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  filename << fixed << setprecision(2) << "./data/mu" << mu << "/phase_diagram/FFLO_h" << haf << ".d";
  output_file[2].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 最初の行のコメントを挿入
  output_file[0] << "# T wave_num eigen_v[0] eigen_v[1] eigen_v[2] eigen_v[3] eigen_v[4] eigen_v[5] eigen_v[6] eigen_v[7]" << endl;
  output_file[1] << "# T q order_a_ABS order_a_ARG order_b_ABS order_b_ARG free_energy(S) free_energy(N)" << endl;
  output_file[2] << "# T haf" << endl;

  int loopcount, i, j;                   // ループで必要な変数

  // 温度を変化させながら計算・出力する
  cout << "#T haf q order_a_R order_a_I order_b_R order_b_I omega omega_n" << endl;
#pragma omp parallel for private(j, loopcount, q, q_candidate, q_opti, T, beta, omega, omega_n, cond_opti, one, x, x2, x3, x4, y, xy, x2y) firstprivate(order) schedule(guided, 2)
  for(i = 1; i <= 100; i++) {
    T = 0.001 * (double)i;
    beta = 1.0 / T;
    cond_opti = 0.0;
    q_candidate = 0.0;

    // 重心運動量の候補を大雑把に見つける
    for(j = -10; j <= 10; j++) {
      q = 0.01 * (double)j;

      // 秩序パラメータを計算
      // Steffensenの反復法を用いた
      loopcount = 0;
      order[0] = order_init;
      order[1] = calculate_order_parameter(order[0], haf, mu, q, T, beta);
      order[2] = calculate_order_parameter(order[1], haf, mu, q, T, beta);
      while(abs((order[0][0] - order[1][0]) / order[0][0]) > delta || abs((order[0][1] - order[1][1]) / order[0][1]) > delta) {
// #ifndef _OPENMP
//         cout << fixed << setprecision(8) << loopcount << "  " << order[0][0] << "  " << order[0][1] << endl;
// #endif
        loopcount += 1;
        order[0][0] = abs(order[0][0] - pow(order[1][0] - order[0][0], 2.0) / (order[0][0] - 2.0 * order[1][0] + order[2][0]));
        order[0][1] = abs(order[0][1] - pow(order[1][1] - order[0][1], 2.0) / (order[0][1] - 2.0 * order[1][1] + order[2][1]));

        if(abs(arg(order[0][0]) - arg(order[0][1])) > 0.05 * M_PI) {
          order[0][1] = polar(abs(order[0][1]), arg(order[0][0]));
        }

        order[1] = calculate_order_parameter(order[0], haf, mu, q, T, beta);
        order[2] = calculate_order_parameter(order[1], haf, mu, q, T, beta);

        // 200回試行しても収束しなければ秩序パラメータは0とみなしてループを抜ける
        if(loopcount >= 200) {
          order[0][0] = Complex(0.0, 0.0);
          order[0][1] = Complex(0.0, 0.0);
          break;
        }
      }

      // 求まった秩序パラメータを用いて自由エネルギーを計算
      omega = calculate_free_energy(order[0], haf, mu, q, T, beta);
      omega_n = calculate_free_energy(order_zero, haf, mu, q, T, beta);

// #ifdef _OPENMP
// #pragma omp critical
//       {
//         cout << omp_get_thread_num() << fixed << setprecision(8) << "  " << T << "  " << haf << "  " << q << "  " << real(order[0][0]) << "  " << imag(order[0][0]) << "  " << real(order[0][1]) << "  " << imag(order[0][1]) << "  " << omega << "  " << omega_n << endl;
//       }
// #else
//       cout << fixed << setprecision(8) << T << "  " << haf << "  " << q << "  " << real(order[0][0]) << "  " << imag(order[0][0]) << "  " << real(order[0][1]) << "  " << imag(order[0][1]) << "  " << omega << "  " << omega_n << endl;
// #endif

#pragma omp critical
      {
        // 自由エネルギー最小の点を最適解として保存
        if(omega - omega_n < cond_opti) {
          cond_opti = omega - omega_n;
          q_candidate = q;
        }
      }
    }

// #ifndef _OPENMP
//     cout << fixed << setprecision(8) << "q_candidate = " << q_candidate << endl;
// #endif

    // 重心運動量を細かく調べ，最小2乗法で2次関数をフィッティングして最適なものを見つける
    one = 0.0;
    x = 0.0, x2 = 0.0, x3 = 0.0, x4 = 0.0;
    y = 0.0, xy = 0.0, x2y = 0.0;
    q_opti = 0.0;

    for(j = -5; j <= 5; j++) {
      q = q_candidate + 0.001 * (double)j;

      // 秩序パラメータを計算
      // Steffensenの反復法を用いた
      loopcount = 0;
      order[0] = order_init;
      order[1] = calculate_order_parameter(order[0], haf, mu, q, T, beta);
      order[2] = calculate_order_parameter(order[1], haf, mu, q, T, beta);
      while(abs((order[0][0] - order[1][0]) / order[0][0]) > delta || abs((order[0][1] - order[1][1]) / order[0][1]) > delta) {
// #ifndef _OPENMP
//         cout << fixed << setprecision(8) << loopcount << "  " << order[0][0] << "  " << order[0][1] << endl;
// #endif
        loopcount += 1;
        order[0][0] = abs(order[0][0] - pow(order[1][0] - order[0][0], 2.0) / (order[0][0] - 2.0 * order[1][0] + order[2][0]));
        order[0][1] = abs(order[0][1] - pow(order[1][1] - order[0][1], 2.0) / (order[0][1] - 2.0 * order[1][1] + order[2][1]));

        if(abs(arg(order[0][0]) - arg(order[0][1])) > 0.05 * M_PI) {
          order[0][1] = polar(abs(order[0][1]), arg(order[0][0]));
        }

        order[1] = calculate_order_parameter(order[0], haf, mu, q, T, beta);
        order[2] = calculate_order_parameter(order[1], haf, mu, q, T, beta);

        // 200回試行しても収束しなければ秩序パラメータは0とみなしてループを抜ける
        if(loopcount >= 200) {
          order[0][0] = Complex(0.0, 0.0);
          order[0][1] = Complex(0.0, 0.0);
          break;
        }
      }

      // 求まった秩序パラメータを用いて自由エネルギーを計算
      omega = calculate_free_energy(order[0], haf, mu, q, T, beta);
      omega_n = calculate_free_energy(order_zero, haf, mu, q, T, beta);

// #ifdef _OPENMP
// #pragma omp critical
//       {
//         cout << omp_get_thread_num() << fixed << setprecision(8) << "  " << T << "  " << haf << "  " << q << "  " << real(order[0][0]) << "  " << imag(order[0][0]) << "  " << real(order[0][1]) << "  " << imag(order[0][1]) << "  " << omega << "  " << omega_n << endl;
//       }
// #else
//       cout << fixed << setprecision(8) << T << "  " << haf << "  " << q << "  " << real(order[0][0]) << "  " << imag(order[0][0]) << "  " << real(order[0][1]) << "  " << imag(order[0][1]) << "  " << omega << "  " << omega_n << endl;
// #endif

      // 最小2乗法で計算するための値
      one += 1.0;
      x += q;
      x2 += pow(q, 2.0);
      x3 += pow(q, 3.0);
      x4 += pow(q, 4.0);
      y += omega - omega_n;
      xy += q * (omega - omega_n);
      x2y += pow(q, 2.0) * (omega - omega_n);
    }

    // 最小2乗法を使って最適な重心運動量を求める
    q_opti = - 0.5 * (x * x2 * x2y - one * x3 * x2y + one * x4 * xy - pow(x2, 2.0) * xy + x2 * x3 * y - x * x4 * y) / (one * x2 * x2y - pow(x, 2.0) * x2y + x * x2 * xy - one * x3 * xy + x * x3 * y - pow(x2, 2.0) * y);

    // 最適な重心運動量での秩序パラメータを計算
    // Steffensenの反復法を用いた
    loopcount = 0;
    order[0] = order_init;
    order[1] = calculate_order_parameter(order[0], haf, mu, q_opti, T, beta);
    order[2] = calculate_order_parameter(order[1], haf, mu, q_opti, T, beta);
    while(abs((order[0][0] - order[1][0]) / order[0][0]) > delta || abs((order[0][1] - order[1][1]) / order[0][1]) > delta) {
// #ifndef _OPENMP
//       cout << fixed << setprecision(8) << loopcount << "  " << order[0][0] << "  " << order[0][1] << endl;
// #endif
      loopcount += 1;
      order[0][0] = abs(order[0][0] - pow(order[1][0] - order[0][0], 2.0) / (order[0][0] - 2.0 * order[1][0] + order[2][0]));
      order[0][1] = abs(order[0][1] - pow(order[1][1] - order[0][1], 2.0) / (order[0][1] - 2.0 * order[1][1] + order[2][1]));

      if(abs(arg(order[0][0]) - arg(order[0][1])) > 0.05 * M_PI) {
        order[0][1] = polar(abs(order[0][1]), arg(order[0][0]));
      }

      order[1] = calculate_order_parameter(order[0], haf, mu, q_opti, T, beta);
      order[2] = calculate_order_parameter(order[1], haf, mu, q_opti, T, beta);

      // 200回試行しても収束しなければ秩序パラメータは0とみなしてループを抜ける
      if(loopcount >= 200) {
        order[0][0] = Complex(0.0, 0.0);
        order[0][1] = Complex(0.0, 0.0);
        break;
      }
    }

    // 求まった秩序パラメータを用いて自由エネルギーを計算
    omega = calculate_free_energy(order[0], haf, mu, q_opti, T, beta);
    omega_n = calculate_free_energy(order_zero, haf, mu, q_opti, T, beta);

#ifndef _OPENMP
    cout << fixed << setprecision(8) << T << "  " << haf << "  " << q_opti << "  " << real(order[0][0]) << "  " << imag(order[0][0]) << "  " << real(order[0][1]) << "  " << imag(order[0][1]) << "  " << omega << "  " << omega_n << endl;
#endif

    // 最適な重心運動量で物理量を出力
#pragma omp critical
    {
      output(output_file[0], output_file[1], order[0], haf, mu, q_opti, T, beta);
      if( ( abs(order[0][0]) >= 1.0e-5 || abs(order[0][1]) >= 1.0e-5 ) && omega - omega_n < 0.0 ) {
        output_file[2] << fixed << setprecision(8) << T << "  " << haf << endl;
      }
    }

  }

  output_file[0].close();
  output_file[1].close();
  output_file[2].close();

  return 0;
}
