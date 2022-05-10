//-*- coding: utf-8 -*-
// Time-stamp: <2015-12-31 11:12:02 shunta>
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
  Complex order[3], order_opti;          // 秩序パラメータ(Steffensenの反復法を使うため3つ定義)
  double q, q_candidate, q_opti;         // Cooper対の重心運動量
  double T;                              // 温度
  double beta;                           // 逆温度
  double omega, omega_n, cond_opti;      // 自由エネルギー

  // 最小2乗法に必要なパラメータ
  double one;
  double x, x2, x3, x4;
  double y, xy, x2y;

  // ファイル入出力関連
  stringstream filename;                 // ファイル名
  ofstream output_file[1]; // 出力ファイル

  // 入力されたパラメータを設定
  if(argc > 2) {
    T = atof(argv[1]);
    if(T <= 0 || T > 0.5) {
      cerr << "Temperature range: (0, 0.5]" << endl;
      return 0;
    }
    beta = 1.0 / T;

    haf = atof(argv[2]);
    if(haf <= 0 || haf > 4.0) {
      cerr << "h^AF range: (0, 4.0]" << endl;
      return 0;
    }
  }
  else {
    cerr << "Usage: " << argv[0] << " Temperature h^AF" << endl;
    return 0;
  }

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  int loopcount, i, j;                   // ループで必要な変数

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/h" << haf << setprecision(3) << "_T" << T << ".d";
  output_file[0].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 最初の行のコメントを挿入
  output_file[0] << "# T mu q" << endl;

  // 化学ポテンシャルの値を変えながら，重心運動量としてよさそうなものを見つける
#pragma omp parallel for private(loopcount, j, mu, order, q, q_candidate, q_opti, omega, omega_n, cond_opti, one, x, x2, x3, x4, y, xy, x2y) schedule(guided, 2)
  for(i = -20; i <= 40; i++) {
    mu = 0.05 * (double)i;
    cond_opti = 0.0;
    q_candidate = 0.0;

    // 重心運動量の候補を大雑把に見つける
    for(j = -10; j <= 10; j++) {
      q = 0.01 * (double)j;

      // ギャップ方程式を満たす(収束条件に入っている)秩序パラメータを見つける
      // Steffensenの反復法を用いた
      loopcount = 0;
      order[0] = 1.0;
      order[1] = calculate_order_parameter(order[0], haf, mu, q, T, beta);
      order[2] = calculate_order_parameter(order[1], haf, mu, q, T, beta);
      while(abs((order[0] - order[1]) / order[0]) > delta) {
        loopcount += 1;
        order[0] = order[0] - pow(order[1] - order[0], 2.0) / (order[0] - 2.0 * order[1] + order[2]);
        order[1] = calculate_order_parameter(order[0], haf, mu, q, T, beta);
        order[2] = calculate_order_parameter(order[1], haf, mu, q, T, beta);

        // 200回試行しても収束しなければデータ点を捨ててループを抜ける
        if(loopcount >= 200) {
          break;
        }
      }

      if(loopcount < 200) {
        // 求まった秩序パラメータを用いて自由エネルギーを計算
        omega = calculate_free_energy(order[0], haf, mu, q, T, beta);
        omega_n = calculate_free_energy(0.0, haf, mu, q, T, beta);
#ifdef _OPENMP
#pragma omp critical
        {
          cout << omp_get_thread_num() << fixed << setprecision(8) << ", T = " << T << ", haf = " << haf << ", mu = " << mu << ", q = " << q << ", order = " << abs(order[0]) << ", omega = " << omega << ", cond_energy = " << omega - omega_n << endl;
        }
#else
        cout << fixed << setprecision(8) << "T = " << T << ", haf = " << haf << ", mu = " << mu << ", q = " << q << ", order = " << abs(order[0]) << ", omega = " << omega << ", cond_energy = " << omega - omega_n << endl;
#endif

#pragma omp critical
        {
          // 凝縮エネルギー最小の点を最適解として保存
          if(omega - omega_n < cond_opti) {
            cond_opti = omega - omega_n;
            q_candidate = q;
          }
        }
      } else {
#ifdef _OPENMP
#pragma omp critical
        {
          cout << omp_get_thread_num() << fixed << setprecision(8) << ", T = " << T << ", haf = " << haf << ", mu = " << mu << ", q = " << q << ", missed data" << endl;
        }
#else
        cout << fixed << setprecision(8) << "T = " << T << ", haf = " << haf << ", mu = " << mu << ", q = " << q << ", missed data" << endl;
#endif
      }
    }

#ifndef _OPENMP
    cout << fixed << setprecision(8) << "q_candidate = " << q_candidate << endl;
#endif

    // 重心運動量を細かく調べ，最小2乗法で2次関数をフィッティングして最適なものを見つける
    one = 0.0;
    x = 0.0, x2 = 0.0, x3 = 0.0, x4 = 0.0;
    y = 0.0, xy = 0.0, x2y = 0.0;
    q_opti = 0.0;

    for(j = -5; j <= 5; j++) {
      q = q_candidate + 0.001 * (double)j;

      // ギャップ方程式を満たす(収束条件に入っている)秩序パラメータを見つける
      // Steffensenの反復法を用いた
      loopcount = 0;
      order[0] = 1.0;
      order[1] = calculate_order_parameter(order[0], haf, mu, q, T, beta);
      order[2] = calculate_order_parameter(order[1], haf, mu, q, T, beta);
      while(abs((order[0] - order[1]) / order[0]) > delta) {
        loopcount += 1;
        order[0] = order[0] - pow(order[1] - order[0], 2.0) / (order[0] - 2.0 * order[1] + order[2]);
        order[1] = calculate_order_parameter(order[0], haf, mu, q, T, beta);
        order[2] = calculate_order_parameter(order[1], haf, mu, q, T, beta);

        // 200回試行しても収束しなければデータ点を捨ててループを抜ける
        if(loopcount >= 200) {
          break;
        }
      }

      if(loopcount < 200) {
        // 求まった秩序パラメータを用いて自由エネルギーを計算
        omega = calculate_free_energy(order[0], haf, mu, q, T, beta);
        omega_n = calculate_free_energy(0.0, haf, mu, q, T, beta);
#ifdef _OPENMP
#pragma omp critical
        {
          cout << omp_get_thread_num() << fixed << setprecision(8) << ", T = " << T << ", haf = " << haf << ", mu = " << mu << ", q = " << q << ", order = " << abs(order[0]) << ", omega = " << omega << ", cond_energy = " << omega - omega_n << endl;
        }
#else
        cout << fixed << setprecision(8) << "T = " << T << ", haf = " << haf << ", mu = " << mu << ", q = " << q << ", order = " << abs(order[0]) << ", omega = " << omega << ", cond_energy = " << omega - omega_n << endl;
#endif
        one += 1.0;
        x += q;
        x2 += pow(q, 2.0);
        x3 += pow(q, 3.0);
        x4 += pow(q, 4.0);
        y += omega - omega_n;
        xy += q * (omega - omega_n);
        x2y += pow(q, 2.0) * (omega - omega_n);
      } else {
#ifdef _OPENMP
#pragma omp critical
        {
          cout << omp_get_thread_num() << fixed << setprecision(8) << ", T = " << T << ", haf = " << haf << ", mu = " << mu << ", q = " << q << ", missed data" << endl;
        }
#else
        cout << fixed << setprecision(8) << "T = " << T << ", haf = " << haf << ", mu = " << mu << ", q = " << q << ", missed data" << endl;
#endif
      }
    }

    // 最小2乗法を使って最適な重心運動量を求める
    q_opti = - 0.5 * (x * x2 * x2y - one * x3 * x2y + one * x4 * xy - pow(x2, 2.0) * xy + x2 * x3 * y - x * x4 * y) / (one * x2 * x2y - pow(x, 2.0) * x2y + x * x2 * xy - one * x3 * xy + x * x3 * y - pow(x2, 2.0) * y);

#ifndef _OPENMP
    cout << fixed << setprecision(8) << "q_opti = " << q_opti << endl << endl;
#endif

    // ファイル出力
#pragma omp critical
    {
      output_file[0] << fixed << setprecision(8) << T << "  " << mu << "  " << q_opti << endl;
    }

  }

  cout << endl;

  output_file[0].close();

  return 0;
}
