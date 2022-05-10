//-*- coding: utf-8 -*-
// Time-stamp: <2016-03-11 16:45:14 shunta>
// ----------------------------------------
// 1次元ジグザグチェーン構造の超伝導計算
// ----------------------------------------
#include "zigzag.h"
typedef std::complex<double> Complex;
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
  int loopcount, enecount;               // ループ用変数
  int max_enecount = 6000;
  double dosS[max_enecount][2];          // 超伝導状態密度
  double dosN[max_enecount][2];          // 正常状態密度
  double epsilon = 0.003;                 // Lorentzianのパラメータ

  // ファイル入出力関連
  stringstream filename;                 // ファイル名
  ofstream output_file[2];               // 出力ファイル

  // 入力されたパラメータを設定
  if(argc > 4) {
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

    mu = atof(argv[3]);
    if(mu <= -3.5 || mu > 2.5) {
      cerr << "mu range: (-3.5, 2.5]" << endl;
      return 0;
    }

    q = atof(argv[4]);
    if(q < - M_PI || q > M_PI) {
      cerr << "q range: [-pi, pi]" << endl;
      return 0;
    }
  }
  else {
    cerr << "Usage: " << argv[0] << " Temperature h^AF mu q" << endl;
    return 0;
  }

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);
  mkdir("./data/energy_band", 0755);
  mkdir("./data/dos", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(3) << "./data/energy_band/T" << T << "_haf" << haf << "_mu" << mu << "_q" << q << ".d";
  output_file[0].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  filename << fixed << setprecision(3) << "./data/dos/T" << T << "_haf" << haf << "_mu" << mu << "_q" << q << ".d";
  output_file[1].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 最初の行のコメントを挿入
  output_file[0] << "# wave_num eigen_v[0] eigen_v[2] eigen_v[4] eigen_v[6]" << endl;
  output_file[1] << "# energy DOS(S) DOS(N)" << endl;

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
    for(enecount = 0; enecount < max_enecount; enecount++) {
      dosS[enecount][0] = 6.0 * (double)enecount / (double)max_enecount - 3.0;
      dosN[enecount][0] = 6.0 * (double)enecount / (double)max_enecount - 3.0;
    }

    // 求まった秩序パラメータを用いて自由エネルギーを計算
    omega = calculate_free_energy(order[0], haf, mu, q, T, beta);
    omega_n = calculate_free_energy(0.0, haf, mu, q, T, beta);

    cout << "# order cond_energy" << endl;
    cout << fixed << setprecision(8) << abs(order[0]) << "  " << omega - omega_n << endl;

    double k, wave_num;                                 // 波数(wave_numは実際に出力する波数)
    static Zigzag::Zmatrix zmatrix;                     // 複素数行列
    int i, j;                                           // ループ用変数
    long info;                                          // 対角化zheevの関数値

    // -pi から pi まで波数を変えて各点で物理量を計算する(超伝導状態)
    for(i = 0; i <= N; i++) {
      k = 2.0 * N_inv * M_PI * (double)i - M_PI;

      // Hamiltonianを定義する
      zmatrix.init_zmatrix(order[0], haf, mu, q, k);

      // 対角化する
      info = zmatrix.zheev();

      // 波数を設定する(piを超えたときは折り返す)
      if(k + q * 0.5 <= M_PI) {
        wave_num = k + q * 0.5;
      }
      else {
        wave_num = k + q * 0.5 - 2.0 * M_PI;
      }

      // 波数を出力する
      output_file[0] << fixed << setprecision(8) << wave_num;

      // 固有エネルギーを出力する
      for(j = 0; j <= DIM / 2 - 1; j++) {
        output_file[0] << fixed << setprecision(8) << "  " << zmatrix.e[2 * j];
      }
      output_file[0] << endl;

      // 状態密度を計算する
      for(j = 0; j <= DIM / 2 - 1; j++) {
        for(enecount = 0; enecount < max_enecount; enecount++) {
          dosS[enecount][1] += N_inv / M_PI * epsilon / (pow(dosS[enecount][0] - zmatrix.e[2 * j], 2.0) + pow(epsilon, 2.0));
        }
      }

    }

    // -pi から pi まで波数を変えて各点で物理量を計算する(超伝導状態)
    for(i = 0; i <= N; i++) {
      k = 2.0 * N_inv * M_PI * (double)i - M_PI;

      // Hamiltonianを定義する
      zmatrix.init_zmatrix(Complex(0.0, 0.0), haf, mu, q, k);

      // 対角化する
      info = zmatrix.zheev();

      // 状態密度を計算する
      for(j = 0; j <= DIM / 2 - 1; j++) {
        for(enecount = 0; enecount < max_enecount; enecount++) {
          dosN[enecount][1] += N_inv / M_PI * epsilon / (pow(dosS[enecount][0] - zmatrix.e[2 * j], 2.0) + pow(epsilon, 2.0));
        }
      }

    }


  }
  else {
    cout << "missed." << endl;
  }

  // 状態密度を出力する
  for(enecount = 0; enecount < max_enecount; enecount++) {
    output_file[1] << fixed << setprecision(8) << dosS[enecount][0] << "  " << dosS[enecount][1] << "  " << dosN[enecount][1] << endl;
  }

  output_file[0].close();
  output_file[1].close();

  return 0;
}
