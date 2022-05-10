//-*- coding: utf-8 -*-
// Time-stamp: <2016-01-03 14:26:34 shunta>
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
  Complex order[3];                      // 秩序パラメータ(Steffensenの反復法を使うため3つ定義)
  double q;                              // Cooper対の重心運動量
  double T;                              // 温度
  double beta;                           // 逆温度
  double omega, omega_n;                 // 自由エネルギー

  // ファイル入出力関連
  stringstream filename;                 // ファイル名
  ofstream output_file[1]; // 出力ファイル

  // 入力されたパラメータを設定
  if(argc > 3) {
    T = atof(argv[1]);
    if(T <= 0 || T > 0.5) {
      cerr << "Temperature range: (0, 0.5]" << endl;
      return 0;
    }
    beta = 1.0 / T;

    mu = atof(argv[2]);
    if(mu <= -3.5 || mu > 2.5) {
      cerr << "mu range: (-3.5, 2.5]" << endl;
      return 0;
    }

    haf = atof(argv[3]);
    if(haf <= 0 || haf > 4.0) {
      cerr << "h^AF range: (0, 4.0]" << endl;
      return 0;
    }
  }
  else {
    cerr << "Usage: " << argv[0] << " Temperature mu h^AF" << endl;
    return 0;
  }

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);
  filename << fixed << setprecision(2) << "./data/mu" << mu;
  mkdir(filename.str().c_str(), 0755);
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/mu" << mu << "/h" << haf << setprecision(3) << "_T" << T << ".d";
  output_file[0].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 最初の行のコメントを挿入
  output_file[0] << "# T q order cond_energy" << endl;

  int loopcount, i;                      // ループで必要な変数

  // 温度を変化させながら計算・出力する
  cout << "#T haf q order omega omega_n" << endl;

#pragma omp parallel for private(loopcount, order, q, omega, omega_n) schedule(guided, 2)
  // 重心運動量を変えながら凝縮エネルギーを計算
  for(i = -200; i <= 200; i++) {
    q = 0.001 * (double)i;

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

      // 100000回試行しても収束しなければデータ点を捨ててループを抜ける
      if(loopcount >= 100000) {
        break;
      }
    }

    if(loopcount < 100000) {
      // 求まった秩序パラメータを用いて自由エネルギーを計算
      omega = calculate_free_energy(order[0], haf, mu, q, T, beta);
      omega_n = calculate_free_energy(0.0, haf, mu, q, T, beta);
#ifdef _OPENMP
#pragma omp critical
      {
        cout << omp_get_thread_num() << fixed << setprecision(8) << "  " << T << "  " << haf << "  " << q << "  " << abs(order[0]) << "  " << omega << "  " << omega_n << endl;
      }
#else
        cout << fixed << setprecision(8) <<  T << "  " << haf << "  " << q << "  " << abs(order[0]) << "  " << omega << "  " << omega_n << endl;
#endif
#pragma omp critical
      {
        output_file[0] << fixed << setprecision(8) << T << "  " << q << "  " << abs(order[0]) << "  " << omega - omega_n << endl;
      }
    }
  }

  output_file[0].close();

  return 0;
}
