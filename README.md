# README

## 概要
-----

1次元、2次元、3次元の自由電子ガスに対してLindhard関数を計算するプログラム。

## 実行方法
-----
コマンド
```
$ g++ jellium_Lindhard.cpp
```
でコンパイルする。
実行は、実行するフォルダ内にresultという名前のフォルダを作った上で
```
$ ./a.out
```
で実行する。計算が終わると、

* result/q_chi_w=0.txt : 周波数omega=0のLindhard関数のq依存性
* result/q_w_chi.txt : Lindhard関数のq, w依存性

に計算結果が書き込まれる。コマンド
```
$ gnuplot "q_chi_w=0.plt"
$ gnuplot "q_w_Im_chi.plt"
```
を実行すると、resultフォルダ内に計算結果のプロット画像が生成される。

