// 三角不等式を満たす距離空間を作るプログラム
// グリッド上に適当に点をばらまいてそれのマンハッタン距離を求める
// O(N^2 log N)
#include<iostream>
#include<random>
#include<vector>
#include<algorithm>
#include<cmath>

using namespace std;

random_device rnd;
mt19937 mt(rnd());

int main() {
    // グラフの頂点数
    int N;
    cin >> N;
    // 返す距離空間は行列で作る
    vector<vector<int> > mat(N, vector<int>(N));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int tmp = abs(i-j);
            mat[i][j] = min(tmp, N-tmp);
        }
    }
    // 出力
    cout << N << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
    return 0;
}
