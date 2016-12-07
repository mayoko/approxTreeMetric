#include<iostream>
#include<random>
#include<vector>
#include<algorithm>
#include<unordered_map>
#include<cassert>
#include<cstring>
#include<queue>

using namespace std;

random_device rnd;
mt19937 mt(rnd());
// 距離の最大値
const int MAX = 114514;
// 頂点数の最大値
const int MAXV = 3333;
// 頂点数(入力)
int N;

struct Edge {
    int to;
    int cost;
    int cap;
    Edge() {}
    Edge(int t, int c) : to(t), cost(c) {}
};

// 距離空間
int dist[MAXV][MAXV];
// requirement
int R[MAXV][MAXV];

// 木(出力すべきもの)
vector<Edge> T[2*MAXV];
// 木の頂点に対応する頂点群がシングルトンになっているとき, 対応する頂点をメモしておく(これも出力すべきもの)
int corr[2*MAXV];
// 今のところの木の頂点数
int treeNum;

// 距離空間が与えられるので, 近似的に木距離にする
// O(N^2)
namespace approxTreeMetric {
    // 木が何段組みになるかってヤツ
    int delta;
    // 例の乱数
    double beta;
    // 例のランダム順列
    int perm[MAXV];
    // 距離空間で最大の距離と最小の距離
    int mini = MAX, maxi = 0;
    // ある頂点が, 今順列の何番目まで進んでいるか
    int cur[MAXV];

    // 深さ depth のところで頂点群 vs が残っているので適当に分割して木を作る
    // 集合 vs に対応する木の頂点は node とラベリングされてる
    void dfs(const vector<int>& vs, int depth, int node) {
        // 新しい頂点についたので
        ++treeNum;
        if (vs.size() == 1) {
            corr[node] = vs[0];
            return;
        }
        // vs を最も近い頂点でグループ分け
        unordered_map<int, vector<int> > mp;
        double radius = beta * pow(2, delta-depth-2)*mini;
        int cost = pow(2, delta-depth)*mini;
        for (int el : vs) {
            int& now = cur[el];
            while (now < N) {
                if (dist[el][perm[now]] <= radius) break;
                ++now;
            }
            assert(now < N);
            mp[perm[now]].push_back(el);
        }
        for (auto p : mp) {
            T[node].emplace_back(treeNum, cost);
            dfs(p.second, depth+1, treeNum);
        }
    }

    // 与えられた入力を用いて解く
    void solve() {
        // 最小値, 最大値を持ってくる
        mini = MAX, maxi = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i != j) mini = min(mini, dist[i][j]);
                maxi = max(maxi, dist[i][j]);
            }
        }
        // 木が何段必要なのかを求める
        {
            int now = mini;
            while (1) {
                if (now >= maxi) break;
                now *= 2;
                delta++;
            }
        }
        cerr << delta << endl;
        // ランダム順列と乱数を用意しておく
        for (int i = 0; i < N; i++)
            perm[i] = i;
        shuffle(perm, perm+N, mt19937());
        // とりあえず一様乱数で良いや
        // TODO: 1/(x ln 2) の形の乱数を用意する
        uniform_real_distribution<double> uni(1.0, 2.0);
        beta = uni(mt);
        // 木の作成
        vector<int> vs(N);
        for (int i = 0; i < N; i++)
            vs[i] = i;
        memset(corr, -1, sizeof(corr));
        dfs(vs, 0, 0);
        // 木の出力
        cerr << treeNum << endl;
        for (int i = 0; i < treeNum; i++) {
            cerr << i << " " << corr[i] << endl;
            for (Edge e : T[i]) {
                cerr << e.to << " " << e.cost << endl;
            }
            cerr << endl;
        }
    }
}

// 木のどこまでが実際の頂点と対応してるか的な
int maxID;
// 木グラフから距離に変換する
namespace treeToMetric {
    // 出力すべき距離
    int ans[MAXV][MAXV];
    // 対応する頂点
    int rcorr[2*MAXV];
    // 頂点の親
    int par[2*MAXV];
    int lca[2*MAXV][2*MAXV];
    // 深さ(辺のコストは考慮してない)
    int depth[2*MAXV];
    // 辺のコストを考慮した深さ
    int length[2*MAXV];

    void dfs(int v, int d, int l) {
        depth[v] = d;
        length[v] = l;
        for (Edge e : T[v]) {
            dfs(e.to, d+1, l+e.cost);
        }
    }

    void solve() {
        int N = treeNum;
        for (int i = 0; i < N; i++) {
            maxID = max(maxID, corr[i]);
            if (corr[i] != -1) rcorr[corr[i]] = i;
            int edgeNum = T[i].size();
            for (int j = 0; j < edgeNum; j++) {
                par[T[i][j].to] = i;
            }
        }
        dfs(0, 0, 0);
        // lca を求める
        {
            vector<int> order(N);
            for (int i = 0; i < N; i++)
                order[i] = i;
            sort(order.begin(), order.end(), [&](int i, int j) {return depth[i] < depth[j];});
            for (int i = 0; i < N; i++) {
                int v = order[i];
                for (int j = 0; j < N; j++) {
                    int u = order[j];
                    if (v == u) lca[v][u] = v;
                    else if (depth[v] > depth[u]) lca[v][u] = lca[par[v]][u];
                    else if (depth[v] < depth[u]) lca[v][u] = lca[v][par[u]];
                    else if (depth[v] > 0) lca[v][u] = lca[par[v]][u];
                    else lca[v][u] = 0;
                }
            }
        }
        for (int i = 0; i <= maxID; i++) {
            for (int j = 0; j <= maxID; j++) {
                int v = rcorr[i], u = rcorr[j];
                int l = lca[v][u];
                ans[i][j] = length[v]-length[l] + length[u]-length[l];
                //cerr << ans[i][j] << " " ;
            }
            //cerr << endl;
        }
        vector<double> muls;
        double avg = 0, mini = MAX, maxi = 0, sum = 0, total = 0;
        for (int i = 0; i <= maxID; i++) {
            for (int j = i+1; j <= maxID; j++) {
                total += ans[i][j];
                double tmp = 1. * ans[i][j] / dist[i][j];
                assert(tmp >= 0.9999);
                mini = min(mini, tmp);
                maxi = max(maxi, tmp);
                muls.push_back(tmp);
                avg += tmp;
                sum += dist[i][j];
            }
        }
        sort(muls.begin(), muls.end());
        cerr << "total ratio: " << total/sum << endl;
        int num = (maxID+1)*maxID/2;
        avg /= num;
        cerr << "avg ratio: " << avg << endl;
        cerr << "min ratio: " << mini << endl;
        cerr << "med ratio: " << muls[num/2] << endl;
        cerr << "max ratio: " << maxi << endl;
    }

    bool exist[2*MAXV];
    void getVertex(int v, vector<int>& vs) {
        for (const Edge& e : T[v]) {
            getVertex(e.to, vs);
        }
        if (T[v].size() == 0) vs.push_back(v);
    }
    void efs(int v, vector<int>& vs) {
        // 一番上の頂点から一番辺のコストを削れる場所を探す
        int best = -1, bestValue = -1;
        for (int i = 0; i < T[v].size(); ++i) {
            const Edge& e = T[v][i];
            vector<int> ws;
            getVertex(e.to, ws);
            memset(exist, false, sizeof(exist));
            for (int el : ws)
                exist[el] = true;
            int mini = 1e9;
            for (int el : ws) {
                for (int i = 0; i < treeNum; ++i) if (!exist[i] && corr[i] != -1) {
                    int v = corr[el], u = corr[i];
                    assert(ans[v][u] - dist[v][u] >= 0);
                    mini = min(mini, ans[v][u] - dist[v][u]);
                }
            }
            mini = min(mini, e.cost);
            if (mini > bestValue) {
                bestValue = mini;
                best = i;
            }
        }
        if (T[v].size() == 0) vs.push_back(v);
        else {
            T[v][best].cost -= bestValue;
            memset(exist, false, sizeof(exist));
            vector<int> ws;
            getVertex(T[v][best].to, ws);
            for (int el : ws)
                exist[el] = true;
            for (int el : ws) {
                for (int i = 0; i < treeNum; ++i) if (!exist[i] && corr[i] != -1) {
                    int v = corr[el], u = corr[i];
                    ans[v][u] -= bestValue;
                    ans[u][v] -= bestValue;
                }
            }
        }
        for (int i = 0; i < T[v].size(); ++i) {
            const Edge& e = T[v][i];
            vector<int> ws;
            efs(e.to, ws);
            vs.insert(vs.end(), ws.begin(), ws.end());
        }
    }
    //void efs(int v, vector<int>& vs) {
    //    bool leaf = true;
    //    // 一番上の頂点から一番辺のコストを削れる場所を探す
    //    int best = -1, bestValue = -1;
    //    vector<vector<int> > wss;
    //    for (int i = 0; i < T[v].size(); ++i) {
    //        const Edge& e = T[v][i];
    //        leaf = false;
    //        vector<int> ws;
    //        efs(e.to, ws);
    //        vs.insert(vs.end(), ws.begin(), ws.end());
    //        wss.push_back(ws);
    //    }
    //    for (int i = 0; i < T[v].size(); ++i) {
    //        const Edge& e = T[v][i];
    //        vector<int> ws = wss[i];
    //        memset(exist, false, sizeof(exist));
    //        for (int el : ws)
    //            exist[el] = true;
    //        int mini = 1e9;
    //        for (int el : ws) {
    //            for (int i = 0; i < treeNum; ++i) if (!exist[i] && corr[i] != -1) {
    //                int v = corr[el], u = corr[i];
    //                assert(ans[v][u] - dist[v][u] >= 0);
    //                mini = min(mini, ans[v][u] - dist[v][u]);
    //            }
    //        }
    //        mini = min(mini, e.cost);
    //        if (mini > bestValue) {
    //            bestValue = mini;
    //            best = i;
    //        }
    //    }
    //    if (leaf) vs.push_back(v);
    //    else {
    //        T[v][best].cost -= bestValue;
    //        memset(exist, false, sizeof(exist));
    //        vector<int> ws;
    //        efs(T[v][best].to, ws);
    //        for (int el : ws)
    //            exist[el] = true;
    //        for (int el : ws) {
    //            for (int i = 0; i < treeNum; ++i) if (!exist[i] && corr[i] != -1) {
    //                int v = corr[el], u = corr[i];
    //                ans[v][u] -= bestValue;
    //                ans[u][v] -= bestValue;
    //            }
    //        }
    //    }
    //}


    void update() {
        // まず一番上の頂点を見つける
        int root = 0;
        while (T[root].size() == 1) {
            root = T[root][0].to;
        }
        if (T[root].size() == 0) {
            cerr << "only one vertex" << endl;
            return;
        }
        vector<int> vs;
        efs(root, vs);
    }
}

int main() {
    cin >> N;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cin >> dist[i][j];
        }
    }
    // 近似する
    cerr << "############################" << endl;
    cerr << "approx Tree Metric" << endl;
    approxTreeMetric::solve();
    // 木距離を求める
    cerr << "############################" << endl;
    cerr << "tree to metric" << endl;
    treeToMetric::solve();
    cerr << "############################" << endl;
    cerr << "tree metric update" << endl;
    treeToMetric::update();
    cerr << "############################" << endl;
    cerr << "new metric" << endl;
    treeToMetric::solve();
}
