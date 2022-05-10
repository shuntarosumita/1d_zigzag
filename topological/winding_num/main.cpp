//-*- coding: utf-8 -*-
// Time-stamp: <2016-02-25 19:53:41 shunta>
// ----------------------------------------
// 1次元ジグザグチェーン構造の超伝導計算
// ----------------------------------------
#include "zigzag.h"
#include "topological.h"
typedef std::complex<double> Complex;
using namespace std;
using namespace Topological;

// メインプログラム
int main(int argc, char* argv[]) {
  // 物理量
  double haf;                            // 反強磁性的な分子場
  double mu;                             // 化学ポテンシャル
  Cvector order(2);                      // 秩序パラメータ
  double wn;                             // winding number
  boost::numeric::ublas::vector<double>
    wn_band(DIM / 2);                    // バンド表示したwinding number
  double bp;                             // Berry phase

  // 秩序パラメータの初期化
  Cvector order_init(2);
  double Delta = 0.01;
  order_init(0) = Complex(Delta, 0.0);
  order_init(1) = Complex(-Delta, 0.0);
  order = order_init;

  // 入力されたパラメータを設定
  if(argc > 3) {
    t[0] = atof(argv[1]);
    t[1] = 1.0;
    if(t[0] < 0.0 || t[0] > 1.0) {
      cerr << "t_1 range: [0.0, 1.0]" << endl;
      return 0;
    }

    mu = atof(argv[2]);
    if(mu <= -3.5 || mu > 2.5) {
      cerr << "mu range: (-3.5, 2.5]" << endl;
      return 0;
    }

    haf = atof(argv[3]);
    if(haf < 0.0 || haf > 4.0) {
      cerr << "h^AF range: [0.0, 4.0]" << endl;
      return 0;
    }
  }
  else {
    cerr << "Usage: " << argv[0] << " t_1 mu h^AF" << endl;
    return 0;
  }

  wn = winding_num(order, t, haf, mu);
  cout << fixed << setprecision(8) << "Total winding number = " << wn << endl;

  // wn = winding_num2(order, t, haf, mu);
  // cout << fixed << setprecision(8) << "Total winding number2 = " << wn << endl;

  // bp = Berry_phase(order, t, haf, mu);
  // cout << fixed << setprecision(8) << "Berry phase = " << bp << endl;

  // project_wavefunc(order, t, haf, mu);

  // wn_band = winding_num_band(order, t, haf, mu);
  // cout << fixed << setprecision(8) << "Band winding number = " << wn_band << endl;

  return 0;
}
