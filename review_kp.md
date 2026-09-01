# Review: block2 K point / k-mesh support

## Summary

`run.py` のように PySCF PBC の `kpts` から block2 の `SU2K` DMRG 入力を作る経路では、現状の `block2` の K symmetry は多次元 `kmesh` をそのまま表現できない。

理由は、block2 側の K symmetry が

- 各軌道の `KSYM`: 1 個の整数ラベル
- 全体の `KMOD`: 1 個の整数 modulus
- target の `KISYM` / `k_irrep`: 1 個の整数セクター

として実装されており、数学的には単一の巡回群 `Z_KMOD`、または `KMOD=0` の無限整数群的な Lz ラベルを扱う設計になっているためである。

したがって、`kmesh = [4, 1, 1]` のような 1D mesh は自然に `Z4` として扱えるが、`kmesh = [2, 2, 2]` のような mesh は本来 `Z2 x Z2 x Z2` であり、単純に `Z8` として扱うことは一般に正しくない。

## run.py conversion path

`run.py` では、PySCF PBC の KRHF 計算結果を block2 の K symmetry 付き FCIDUMP に変換している。

重要箇所:

- `run.py:8`: `_block2_k_labels_1d(cell, kpts, kmesh)` を定義。
- `run.py:10`: `kmesh > 1` の軸を調べる。
- `run.py:11-12`: 非自明な軸が 1 個でなければ例外を投げる。
- `run.py:15-17`: 非自明軸の scaled k-point を使って、スカラーの `labels` を作る。
- `run.py:19-24`: PySCF の `kconserv` とスカラー `labels` の保存則が一致するか確認する。
- `run.py:58`: FCIDUMP 生成時に `_block2_k_labels_1d` を必ず使う。
- `run.py:111`: 各軌道へ `fd.k_sym` としてスカラー k label を割り当てる。
- `run.py:112`: `fd.k_mod = nk` として単一 modulus を設定する。
- `run.py:113`: `fd.k_isym = int(target_k) % nk` として target も単一整数にする。
- `run.py:209`: 実例では `kmesh = [4, 1, 1]`。

この関数名と実装から、少なくとも `run.py` の converter は意図的に 1D kmesh 専用である。

## FCIDUMP.K format observed

現在の `FCIDUMP.K` の先頭は次の形だった。

```text
&FCI NORB=   4,NELEC=   8,MS2=   0,
  ORBSYM=1,1,1,1,
  KSYM=0,1,2,3,
  KMOD=4,
  ISYM=   1,
  IGENERAL=1,
  ITGENERAL=1,
&END
```

ここで、`KSYM` は軌道ごとの配列だが、`KMOD` は 1 個の整数である。

重要箇所:

- `FCIDUMP.K:3`: `KSYM=0,1,2,3,`
- `FCIDUMP.K:4`: `KMOD=4,`

現状仕様で `KMOD=2,2,2,` のような複数値を持つ設計にはなっていない。

## Core FCIDUMP implementation

`src/core/integral.hpp` の `FCIDUMP` 実装でも、`kmod` は単一整数として扱われている。

重要箇所:

- `src/core/integral.hpp:710-711`: `KMOD=` を書き出すときに `k_mod()` を 1 個だけ出力する。
- `src/core/integral.hpp:1064`: `symmetrize(const vector<int> &ksym, int kmod)` の signature が `int kmod`。
- `src/core/integral.hpp:1075-1076`: 1 電子積分の k selection rule が `(ksym[i] + kmod - ksym[j]) % kmod`。
- `src/core/integral.hpp:1083-1087`: 2 電子積分の k selection rule が `(ksym[i] - ksym[j] + ksym[k] - ksym[l]) % kmod`。
- `src/core/integral.hpp:1167`: `set_k_mod(int kmod)`。
- `src/core/integral.hpp:1169-1170`: `int k_mod() const` で、`params["kmod"]` を `Parsing::to_int` して返す。
- `src/core/integral.hpp:1173-1177`: `k_isym` も単一整数として保存・取得する。

このため、FCIDUMP レベルで多成分 k label や多成分 modulus を表す余地は現在ない。

## Symmetry quantum number implementation

K symmetry は `SU2K` / `SZK` などの symmetry quantum number に、point group と組み合わされた形で埋め込まれる。

重要箇所:

- `src/core/hamiltonian.hpp:125-131`: `Hamiltonian::combine_orb_sym` が `orb_sym`, `k_sym`, `k_mod` を受け取り、各軌道の symmetry を `S::pg_combine(orb_sym[i], k_sym[i], k_mod)` に変換する。
- `src/core/symmetry.hpp:763`: `SZKLong::pg_k_mod()` は 1 個の k modulus を返す。
- `src/core/symmetry.hpp:798-806`: `SZKLong::operator+` は 1 個の k 成分について mod 加算する。
- `src/core/symmetry.hpp:830`: `SZKLong::pg_combine(int pg, int k = 0, int kmod = 0)` は `pg`, `k`, `kmod` の 3 スカラーを pack する。
- `src/core/symmetry.hpp:1347`: `SU2KLong::pg_k_mod()` も 1 個の k modulus を返す。
- `src/core/symmetry.hpp:1385-1393`: `SU2KLong::operator+` は 1 個の k 成分について mod 加算する。
- `src/core/symmetry.hpp:1444`: `SU2KLong::pg_combine(int pg, int k = 0, int kmod = 0)` も `pg`, `k`, `kmod` の 3 スカラーを pack する。

ここでも、`kx, ky, kz` を別々の成分として持つ構造ではない。

## block2main input path

`block2main` の入力でも、K symmetry は単一の `k_irrep` / `k_mod` として扱われる。

重要箇所:

- `pyblock2/driver/parser.py:42`: known keys に `k_symmetry`, `k_irrep`, `k_mod` がある。
- `pyblock2/driver/block2main:832`: `k_irrep` を整数列として読むが、各値は単一の k sector。
- `pyblock2/driver/block2main:1020-1027`: `k_mod` が指定された場合、`fcidump.k_mod` を単一整数にし、`fcidump.k_sym` と `fcidump.k_isym` をその modulus で丸める。
- `pyblock2/driver/block2main:1028-1030`: target は `SX.pg_combine(..., fcidump.k_isym, fcidump.k_mod)` で作る。
- `pyblock2/driver/block2main:1036-1040`: 複数 target sector を列挙する場合も、各 `iiksym` は単一スカラー。
- `pyblock2/driver/block2main:1061-1068`: `fcidump.symmetrize(k_sym, k_mod)` を呼び、`HamiltonianQC.combine_orb_sym(orb_sym, k_sym, k_mod)` で orbitals の対称性へ統合する。

`k_irrep` に複数値を書けるとしても、それは複数の候補 target sector を列挙するためであり、1 個の target が多成分 k vector を持つという意味ではない。

## Hubbard k-space model

`model hubbard_kspace 16 1 2` は、周期 Hubbard model の momentum-space 計算として説明されているが、実装上は 1D Hubbard chain の k-space 表現である。

重要箇所:

- `docs/source/user/keywords.rst:37`: `model hubbard 16 1 2` は 1-dimensional non-periodic Hubbard model と説明される。
- `docs/source/user/keywords.rst:38`: `model hubbard_periodic 16 1 2` は periodic Hubbard model。
- `docs/source/user/keywords.rst:39`: `model hubbard_kspace 16 1 2` は periodic Hubbard model in momentum space。
- `pyblock2/driver/block2main:849-862`: `hubbard_kspace` は `HubbardKSpaceFCIDUMP(n_sites, const_t, const_u)` を作る。
- `pyblock2/driver/block2main:856`: ログ表示は `1D %s model : L = ...`。
- `src/core/hubbard.hpp:84`: `HubbardKSpaceFCIDUMP(uint16_t n_sites, double t, double u)`。
- `src/core/hubbard.hpp:96-106`: `ksym` を `0,1,...,n_sites-1` の 1D ラベルとして設定する。
- `src/core/hubbard.hpp:111`: `kmod = n_sites`。
- `src/core/hubbard.hpp:115-119`: one-body dispersion は `-2 t cos(2 pi i / n_sites + pi)`。これは 1D tight-binding chain の分散。
- `src/core/hubbard.hpp:121-130`: two-body selection rule は `(i - j + k - l) mod n = 0`。

したがって、ここでいう `hubbard_kspace` は 1D periodic Hubbard model の momentum-space representation と読むべきである。

## Documentation context

`docs/source/user/advanced.rst` の LZ Symmetry 節も、K/LZ symmetry の主な対象を 1D 的な追加量子数として説明している。

重要箇所:

- `docs/source/user/advanced.rst:561`: `LZ Symmetry` 節。
- `docs/source/user/advanced.rst:564`: diatomic molecules または translational symmetry を持つ model Hamiltonian、例として `1D Hubbard model in momentum space`。
- `docs/source/user/advanced.rst:566`: K space symmetry には `-DUSE_KSYMM=ON` が必要。
- `docs/source/user/advanced.rst:568-572`: `k_symmetry` keyword による追加対称性。
- `docs/source/user/advanced.rst:576`: 分子では K/LZ symmetry 用に特殊な FCIDUMP 生成が必要。
- `docs/source/user/advanced.rst:718`: `k_irrep` は target state の LZ eigenvalue を指定。
- `docs/source/user/advanced.rst:736`: `k_mod` は modulus。`k_mod = 0` は infinite LZ group。
- `docs/source/user/advanced.rst:742-744`: `k_mod 4` による `C2 x Z4` 例。
- `docs/source/user/advanced.rst:766`: `k_mod 2` 例。

このドキュメントも、複数方向の k vector ではなく、単一の Lz/K ラベルを想定している。

## Mathematical issue: `2 x 2 x 2` kmesh

`kmesh = [2, 2, 2]` の場合、本来の有限並進群は

```text
Z2 x Z2 x Z2
```

である。

一方、現在の block2 K symmetry に無理に流し込むと、取り得る表現は単一 modulus の

```text
Z8
```

になる。

この 2 つは同じ 8 個の元を持つが、群として明確に異なる。

違い:

- `Z8` は巡回群で、位数 8 の生成元を持つ。
- `Z2 x Z2 x Z2` は全ての非単位元の位数が 2。
- `Z8` の位数 2 の元は 1 個だけ。
- `Z2 x Z2 x Z2` の非単位元 7 個はすべて位数 2。

したがって同型ではない。

保存則も異なる。

`Z8` では、整数ラベル `a` に対して

```text
(a_i - a_j + a_k - a_l) mod 8 = 0
```

を判定する。

`Z2 x Z2 x Z2` では、各成分ごとに

```text
(kx_i - kx_j + kx_k - kx_l) mod 2 = 0
(ky_i - ky_j + ky_k - ky_l) mod 2 = 0
(kz_i - kz_j + kz_k - kz_l) mod 2 = 0
```

を判定する必要がある。

`Z2` では引き算と足し算は同じなので、bit 表現を使うなら XOR 的な成分別演算になる。これは一般に integer `mod 8` 加算とは一致しない。

## When can multi-dimensional mesh be encoded into one cyclic group?

有限 Abel 群

```text
Z_nx x Z_ny x Z_nz
```

が単一の巡回群 `Z_N`, `N = nx * ny * nz`, と同型になるのは、基本的には `nx`, `ny`, `nz` が互いに素な場合である。

例えば:

- `Z2 x Z3` は `Z6` と同型。
- `Z2 x Z2` は `Z4` と同型ではない。
- `Z2 x Z2 x Z2` は `Z8` と同型ではない。
- `Z4 x Z4` は `Z16` と同型ではない。

そのため、一般の 2D/3D kmesh を単一 `KMOD = prod(kmesh)` に潰す方法は、一般解ではない。

特に電子構造の k-mesh では `2x2x2`, `4x4x1`, `4x4x4` のように互いに素でない mesh が普通に出るため、単一 `Z_N` への置き換えは危険である。

## Practical conclusion

現状仕様で安全に使えるのは次のケースである。

1. `kmesh = [N, 1, 1]`, `[1, N, 1]`, `[1, 1, N]` のような 1D mesh。
2. Lz のような単一整数量子数。
3. 特殊に `Z_nx x Z_ny x ...` が `Z_N` と同型になる場合。ただし実装側の k label mapping と PySCF の `kconserv` がその群演算と一致する必要がある。

安全ではないケース:

1. `kmesh = [2, 2, 2]` を `KMOD=8` として扱う。
2. `kmesh = [4, 4, 1]` を `KMOD=16` として扱う。
3. `kmesh = [N_x, N_y, N_z]` を単純な flatten index で 1D 化し、`mod prod(kmesh)` で保存則を判定する。

## What would be needed for true multi-dimensional kmesh support

本当に 2D/3D kmesh を block2 の K symmetry として扱うには、少なくとも次の設計変更が必要になる。

1. FCIDUMP metadata
   - `KSYM` を 1 軌道 1 整数ではなく、1 軌道あたり多成分 k vector として表す。
   - `KMOD` も 1 整数ではなく、例えば `KMOD=2,2,2` のような vector にする。

2. Integral symmetrization
   - `src/core/integral.hpp` の `symmetrize(ksym, kmod)` を多成分対応にする。
   - selection rule を成分ごとの modulo 判定にする。

3. Symmetry quantum numbers
   - `SU2KLong` / `SZKLong` の内部 pack 形式を拡張する。
   - `pg_combine(pg, k, kmod)` を多成分 k に対応させる。
   - `operator+`, `operator-`, `pg_mul`, `pg_inv`, `pg_equal`, `hash`, `to_str` なども多成分 k の群演算に合わせる。

4. Python driver / input parser
   - `k_mod` / `k_irrep` の vector 指定を許す。
   - target sector も多成分 k vector として扱う。

5. PySCF converter
   - `run.py` の `_block2_k_labels_1d` を一般化する。
   - `cell.get_scaled_kpts(kpts)` から `(kx, ky, kz)` label を作り、PySCF `kconserv` と成分別保存則を照合する。

ただし、これは単なる converter 修正では足りない。block2 本体の symmetry 型と FCIDUMP 仕様を触る比較的大きな変更になる。

## References

### Source files

- `run.py:8-24`
- `run.py:58`
- `run.py:111-113`
- `run.py:209-210`
- `FCIDUMP.K:1-8`
- `src/core/integral.hpp:710-711`
- `src/core/integral.hpp:1064-1110`
- `src/core/integral.hpp:1167-1177`
- `src/core/hamiltonian.hpp:125-131`
- `src/core/symmetry.hpp:763`
- `src/core/symmetry.hpp:798-806`
- `src/core/symmetry.hpp:830`
- `src/core/symmetry.hpp:1347`
- `src/core/symmetry.hpp:1385-1393`
- `src/core/symmetry.hpp:1444`
- `pyblock2/driver/parser.py:42`
- `pyblock2/driver/block2main:832`
- `pyblock2/driver/block2main:849-862`
- `pyblock2/driver/block2main:1020-1040`
- `pyblock2/driver/block2main:1061-1068`
- `src/core/hubbard.hpp:84-130`

### Documentation files

- `docs/source/user/advanced.rst:561-576`
- `docs/source/user/advanced.rst:718`
- `docs/source/user/advanced.rst:736-744`
- `docs/source/user/advanced.rst:766`
- `docs/source/user/keywords.rst:35-40`
