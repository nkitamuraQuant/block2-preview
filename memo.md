# Translational (K point) / Lz symmetry memo

対象ファイル:

- `README.md`
- `docs/source/index.rst`
- `docs/source/user/advanced.rst`
- `docs/source/user/keywords.rst`

## 結論

`Translational (K point) / Lz symmetry` は、block2 の DMRG/MPS/MPO 計算に追加できる状態対称性として扱われている。
ドキュメント上で明示されている主な利用対象は、次の2系統。

1. 並進対称性を持つモデルハミルトニアン
   - 例: 1D Hubbard model in momentum space
   - `block2main` 入力では `model hubbard_kspace ...` と `k_symmetry` を組み合わせる。

2. Lz 対称性を使える分子系
   - 例: diatomic molecule
   - FCIDUMP を K/LZ 対称性に合うよう特殊に生成する必要がある。
   - `advanced.rst` には `C2 x Lz` と `C2 x Z4` の入力例がある。

## 使える手法・モード

### DMRG ground-state 計算

- デフォルトの `block2main` 計算タイプは DMRG。
- `k_symmetry` を入力に追加すると、通常の DMRG 計算を K/LZ 対称性付きで実行できる。
- `advanced.rst` の LZ Symmetry 節では、`C2 x Lz` および `C2 x Z4` 対称性を使った DMRG 入力例が示されている。

最小イメージ:

```text
sym d2h
orbitals FCIDUMP
k_symmetry
k_irrep 0

nelec ...
spin ...
irrep ...
```

### Spin-adapted / non-spin-adapted DMRG

- `keywords.rst` の `k_symmetry` 説明では、K symmetry は spin-adapted mode と non-spin-adapted mode の両方で使えるとされている。
- つまり、通常の空間軌道ベースの `SU2` mode と `SZ` mode が対象。
- 一方で `use_general_spin` は現状 `k_symmetry` と併用不可。

整理:

| モード | K/LZ symmetry |
| --- | --- |
| spin-adapted `SU2` | 利用可 |
| non-spin-adapted `SZ` | 利用可 |
| general spin orbital (`use_general_spin`) | 併用不可 |

### Periodic Hubbard model in momentum space

- `keywords.rst` の `model` では、`model hubbard_kspace 16 1 2` が periodic Hubbard model の momentum-space 計算として説明されている。
- この場合、`k_symmetry` を付けると translational symmetry を利用できる。
- `k_symmetry` を付けなければ、同じ momentum-space model でも translational symmetry は使わない。

例:

```text
model hubbard_kspace 16 1 2
k_symmetry
```

補足:

- `model hubbard_periodic ...` は periodic Hubbard model の実空間側の指定として説明されている。
- k 点対称性を利用する対象として明示されているのは `hubbard_kspace`。

### Diatomic molecule / Lz symmetry

- `advanced.rst` の LZ Symmetry 節では、diatomic molecule で LZ symmetry を利用できると説明されている。
- 分子系では、FCIDUMP をそのまま使うのではなく、K/LZ symmetry を使える形に生成する必要がある。
- 例として `C2 x Lz` symmetry 用の FCIDUMP 生成スクリプトと、対応する DMRG 入力が示されている。

使う主なキーワード:

- `k_symmetry`
- `k_irrep`
- `k_mod`
- `sym`
- `irrep`

`k_irrep` は target state の LZ/K 既約表現番号を指定する。

### 有限群に制限した Lz: `Z4` / `Z2`

- LZ は有限群ではないため、軌道数が多い場合に初期 MPS の bond dimension が大きくなり、最初の sweep が重くなることがある。
- `k_mod` を使うと LZ を modular arithmetic で有限群に制限できる。
- `advanced.rst` では `k_mod 4` による `C2 x Z4` の例と、`k_mod 2` の結果が示されている。

整理:

| 指定 | 意味 |
| --- | --- |
| `k_mod 0` | 元の infinite LZ group |
| `k_mod 4` | `Z4` に制限 |
| `k_mod 2` | `Z2` に制限 |

### Occupation number 初期値を使う K symmetry DMRG

- 高い対称性では Davidson が局所解に詰まりやすい、と `advanced.rst` に注意がある。
- また、LZ symmetry では初期 MPS が大きくなりやすい。
- 対策として、`warmup occ` と `occ` による occupation number 初期値が紹介されている。
- occupation number は同じ integral を使った K symmetry あり/なし DMRG、または CCSD/MP2 などから得ることができる。
- `cbias` を併用して occupation を少し均す例も示されている。

例:

```text
k_symmetry
k_irrep 0
warmup occ
occ ...
cbias 0.2
```

## 併用・制約

### point group symmetry との併用

- `advanced.rst` では、point group symmetry は k symmetry と併用できると明記されている。
- 例では `sym d2h` と `k_symmetry` を同時に使っている。

### SG / general spin との関係

- `keywords.rst` では `use_general_spin` が `k_symmetry` と現在併用不可とされている。
- したがって、`-DUSE_SG` が必要な general spin orbital 計算は `k_symmetry` の対象外。

### ビルド要件

- K symmetry を使うには `-DUSE_KSYMM` が必要。
- `advanced.rst` では `-DUSE_KSYMM=ON` が default と説明されている。

### MPS/MPO の互換性

- `k_symmetry` ありで生成した MPS/MPO と、`k_symmetry` なしで生成した MPS/MPO は相互に reload できない。
- restart 系の計算では、作成時と同じ `k_symmetry` の有無を入力側で揃える必要がある。

### integral の対称性チェック

- `symmetrize_ints` は point group または K symmetry を破る integral element の許容値を設定する。
- 対称性を破る integral element は計算では捨てられる。
- このキーワードは、計算を続行するかエラーにするかの閾値制御。

## ドキュメントから読める範囲での「手法」一覧

| 手法・対象 | K/LZ symmetry の扱い | 根拠 |
| --- | --- | --- |
| 標準 DMRG | `k_symmetry` で追加利用可 | `advanced.rst` LZ Symmetry 入力例、`keywords.rst` `k_symmetry` |
| spin-adapted DMRG (`SU2`) | 利用可 | `keywords.rst` `k_symmetry` |
| non-spin-adapted DMRG (`SZ`) | 利用可 | `keywords.rst` `k_symmetry` |
| general spin orbital 計算 | 併用不可 | `keywords.rst` `use_general_spin` |
| periodic Hubbard model in momentum space | `model hubbard_kspace ...` + `k_symmetry` で利用 | `keywords.rst` `model` |
| translational symmetry を持つ model Hamiltonian | 利用可と説明あり | `advanced.rst` LZ Symmetry |
| diatomic molecule の LZ symmetry | 特殊生成 FCIDUMP + `k_symmetry` で利用 | `advanced.rst` LZ Symmetry |
| `C2 x Lz` | 入力例あり | `advanced.rst` LZ Symmetry |
| `C2 x Z4` / `Z2` | `k_mod` で有限群化 | `advanced.rst` LZ Symmetry |
| point group + K/LZ | 併用可 | `advanced.rst` LZ Symmetry |
| occupation number 初期値付き DMRG | K symmetry 計算の収束補助として利用例あり | `advanced.rst` Initial Guess with Occupation Numbers |

## 実務メモ

- k 点対称性を使うなら、まずビルドが `-DUSE_KSYMM=ON` か確認する。
- model 系なら `model hubbard_kspace ...` に `k_symmetry` を足すのが明示ルート。
- 分子系なら FCIDUMP の作り方が重要。K/LZ symmetry 用に `k_sym`, `k_mod`, `orb_sym` を正しく入れる。
- target の K/LZ セクターは `k_irrep` で指定する。
- LZ で初期 MPS が重い場合は `k_mod 4` や `k_mod 2`、または `warmup occ` を検討する。
- 高対称性では Davidson が局所解に詰まりやすいので、noise を大きめにし、Davidson threshold を小さめにした custom schedule が推奨されている。
