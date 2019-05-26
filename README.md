# scgra
This reposity provides MATLAB implementation of SCGRA (Sequential Conjugate Gradient-Restoration Algorithm), which solves optimal control problem.
I'm afraid that documents are written in Japanese.

最適制御の数値解法の一つである、SCGRA (Sequential Conjugate Gradient-Restoration Algorithm)のMATLAB実装を公開します。
アルゴリズムの詳細は[こちら](./doc/doc.pdf)で解説してあります。

# Example
[RunCar](./src/RunCar.m)がExampleです。
最適制御の教科書の例題に載っているような、車の急加速・急停車問題を解くものです。
入力の上限値に拘束があるため、Bang-Bang解が得られるはずです。


# Ref
1. A.K.Wu and A.Miele, "Sequential conjugate gradient-restoration algorithm for optimal control problems with non-differential constraints and general boundary conditions, part I", Optimal Control Applications and Methods, Volume 1, Issue 1, January/March 1980, pp.69--88
2. 麥谷高志、"状態量不等式拘束を伴う最適問題の数値解法に関する研究"、東京大学博士論文、1996			
3. 原田正範、"最適制御問題による地面付近における滑空機の距離最大飛行の解析に関する研究"、東海大学博士論文、1996
