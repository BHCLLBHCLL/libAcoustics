# libAcoustics Farassat 1A 分析与修改记录

> 分支：`v2412`（`origin/v2412`）
> 基准：最初 `v2412` 为 `04f8a7a`；本文写作时的分支头为 `72eea13`，本次在其上增加 patch 表面法向修复。
> 源码：`/home/ubuntu/programs/libAcoustics`
> 日期：2026-11-16（2026-11-17 补充实测对比、逐面 OAP/SPL/相干性分析与信息融合）

本文把四类内容分开：

1. **性能分析**：最初 `v2412` 在 Windows 8 rank 算例上的平方分配、并行和内存问题；
2. **性能与正确性修改**：当前分支已经合入的导数分配、采样、Farassat 观测网格、面速度、多面求和等；
3. **patch 法向修复**：`laptop_wall` 的 7 个 patch 是闭合壁面，原质心启发式只翻转了 14.5% 的面，导致壁面压力载荷的符号和相消关系错误。本文给出源码路径、几何验证、修复方法和与 Cradle scFLOW 的差异；
4. **2026-11-17 实测复核与信息融合**：新库 `_m1` 与旧库 `_all` 的逐面 OAP、时域漂移、SPL、1/3 倍频程、A 计权、20–50 kHz 全频带、相关系数和能量预算对比。结论：法向修复解决了直流/低频（wall 20–200 Hz −20 dB），但 5–20 kHz 的壁面宽带几乎不变（+0.05 dB），总级仍比叶轮单独高约 19 dB。第 9 节把这些实测与 §5 的几何/源码分析逐项对应。

> 注意：性能节的绝对时间来自最初的 Windows 8 rank 运行；法向节的 OAP 数据来自 Linux 48 rank 的 `openfoam_fwh_m4_all`（旧库）和 `openfoam_fwh_m4_all_m1`（新库）运行；第 9 节另用 `openfoam_fwh_m4`（叶轮单独）作基线。几套数据的网格、进程数和时间步不同，**不要直接混用绝对时间**。第 9 节的比较统一用 `0.03–0.0529 s` 公共窗口；旧数据已停止在 `0.05291 s`，新数据仍在运行。

---

## 目录

- [0. 结论摘要](#0-结论摘要)
- [1. 算例、版本与记号](#1-算例版本与记号)
- [2. 性能分析（最初 v2412）](#2-性能分析最初-v2412)
- [3. 性能修改（已合入）](#3-性能修改已合入)
- [4. 正确性修改（已合入）](#4-正确性修改已合入)
- [5. patch 表面法向定向（本次修复，重点）](#5-patch-表面法向定向本次修复重点)
- [6. 修改清单](#6-修改清单)
- [7. 验证与复现](#7-验证与复现)
- [8. 阶段性结论](#8-阶段性结论)
- [9. 实测复核与信息融合（2026-11-17）](#9-实测复核与信息融合2026-11-17)

---

## 0. 结论摘要

### 0.1 性能

最初 `v2412` 每步约 **98 分钟**，其中流动约 2–3 分钟，其余在 Farassat 1A 探针。主因是 `AcousticAnalogy::tdPrimitive::dot` 在每个面上调用 `List::resize`，而 `List` 没有容量预留，长度每变 1 就重新分配并搬移整张表。一个面上的代价是面数的平方；rank 3 的 `laptop_wall` 面数最多，临界路径约 95 分钟。按当时步长跑到 `endTime = 0.1` 约需 **36 年**。

修改后：导数历史一次分配；每个源面只采样一次；Farassat 改为每麦克风一条观测时间网格；`gatherAndMerge` 移出积分路径；`p'` 写各面之和。声学开销降到流动以下，内存从按面历史（rank 3 约 280 GB/对象）降到每麦克风一条信号（数 MB 量级）。

### 0.2 正确性

已合入的修正包括：

- `Cf0` 每个采样步都更新，面速度是**一个 `deltaT` 的位移**，不再是 `t = 0` 以来的总位移；
- 多个积分面时，`acousticData` 中的 `p'` 是**各面之和**，不再只保留最后一个面；
- `sampledPatch` 表面法向直接用 `-Sf`，不再使用质心启发式（本次新增）。

### 0.3 本次法向修复

`laptop_wall` 由 7 个 wall patch 组成，合计是**闭合面**（`Σ Sf = 0`）。正确的 FW-H 法向是 patch 法向的反方向 `-Sf`（patch 的 `Sf` 从流体指向固体，FW-H 的 `n` 从固体指向流体）。

原代码用

```cpp
Cs = gSum(Cf)/N;
if (((Cf - Cs) & n) > 0)  n 不翻;
else                      n 翻;
```

对 `laptop_wall`：**85.5% 的面没有被翻转**，仍指向固体；只有 14.5% 被翻。结果

```text
Σ n_heur dS = 0.0224 m²      （应为 0）
Σ (-Sf) dS  = 0
```

其中 `air_to_Cover` 单独贡献了 `ΣSf_y = +0.01759 m²`，而启发式几乎原样保留了这个符号（正确应为 `-0.01759`）。cover 占壁面面积的 **73%**，是误差的主要来源。

运行数据（`openfoam_fwh_m4_all`，旧库，0.014–0.0164 s 窗口）：

| 量 | wall | impeller1 | impeller2 | 叶轮之和 |
|---|---:|---:|---:|---:|
| OAP RMS / Pa | **0.197** | 0.022 | 0.019 | 0.059 |
| 5–20 kHz RMS / Pa | **0.193** | 0.010 | 0.011 | 0.014 |
| 均值 / Pa | -0.051 | — | — | +0.044 |

总 `p'` RMS 从叶轮单独的 **0.0595 Pa** 增到 **0.194 Pa**（约 **+10.3 dB**）；5–20 kHz 约 **+23 dB**（相对叶轮之和；相对单个叶轮约 +25 dB）。wall 与叶轮的相关系数只有 **-0.007**，说明它不是叶轮力的物理反作用，而是壁面压力被错误法向积分出来的独立信号。

修复：`fwhFormulation::updateNormals()` 对 `isA<sampledPatch>` 的表面直接使用 `-Sf/|Sf|`；其他表面（STL/sampledSurface）保留质心启发式。叶轮是封闭凸体，原启发式本来就翻对了 99.8%，修复后仍用 `-Sf`，结果不变。

> **2026-11-17 实测修正（要点）**：上表的 0.197/0.193 Pa 和“独立信号”描述来自旧库的 **0.014–0.0164 s 早期窗口**，不能代表稳态。公共窗口 `0.03–0.0529 s` 下，旧 wall OAP 均值 `−0.1915 Pa`、AC RMS `0.1868 Pa`；新库 `−0.0765 / 0.1835 Pa`。法向修复使 wall 20–200 Hz 从 `0.0133` 降到 `0.0013 Pa`（−20.2 dB），但 5–20 kHz 从 `0.1440` 到 `0.1449 Pa`（+0.05 dB）几乎不变；总 OASPL `77.46→76.88 dB`，5–20 kHz `77.40→76.87 dB`，仍比叶轮单独高约 19 dB。旧/新总声压相关系数 `−0.58`，说明壁面项主要是“反号”而不是“消失”。完整数据、能量预算和机理见第 9 节。结论：法向错误只解释低频/直流；5–20 kHz 的 wall 项是壁面压力载荷本身。

---

## 1. 算例、版本与记号

### 1.1 版本

| 项目 | 值 |
|---|---|
| 分支 | `v2412` |
| 最初基准 | `04f8a7a` |
| 当前分支头（修改前） | `72eea13` |
| 求解器 | OpenFOAM v2412 `pimpleFoam -parallel` |
| 公式 | `Farassat1AFormulation` |

### 1.2 算例

**性能算例（Windows，8 rank）**

- 流动：`openfoam`（无声学）
- 流动加声学：`openfoam_fwh`
- 8 进程，`maxCo 0.6`，稳定后 `deltaT ≈ 5.13×10⁻⁷ s`
- 声学对象：`FwhImpellers`（`fixedResponseDelay false`）、`FwhImpellersFixedDelay`（`fixedResponseDelay true`）

**法向算例（Linux，48 rank）**

- 路径：`/home/ubuntu/shared/openfoam_fwh_m4_all`
- 网格：约 6.11 M cells、20.13 M faces
- 声源面：`laptop_wall` + 两个叶轮，共 804,532 面
- 声学对象：`FwhImpellers`、`FwhImpellersFixedDelay`、`FwhImpellersFixedDelay2`
- 观测点：`mic`，位置 `(0.0062 0.1 0.2)`
- 数据快照：旧库停止在 `0.05291 s`；新库 `_m1` 仍在运行，2026-11-17 23:52 已到 `0.06247 s`。第 9 节用两者公共的 `0.03–0.0529 s` 稳态窗口；§5.1 的 `0.014–0.0164 s` 只是早期快照

### 1.3 源面规模

| 面 | 面数 | 面积 / m² |
|---|---:|---:|
| `laptop_wall`（7 个 patch） | 576,666 | 0.19942 |
| `PartSurface_impeller1` | 113,933 | 0.003397 |
| `PartSurface_impeller2` | 113,933 | 0.003397 |
| 合计 | 804,532 | 0.20622 |

`laptop_wall` 的面数占 72%，面积是两个叶轮合计的约 29 倍。AMI 与 `open` 不在声源面上。

---

## 2. 性能分析（最初 v2412）

### 2.1 计时

- 第一步流动结束时 `ExecutionTime = 200.7 s`，此时还没有声学输出。
- 其后每步 `ExecutionTime` 间隔为 5792–5956 s，13 步平均 **5877 s**。
- 日志末尾停在第二个对象的 `Starting acoustics probe`，墙钟约 76598 s，约 13 步、21 小时。
- `endTime = 0.1`，步数约 `0.1 / 5.13e-7 ≈ 1.95×10⁵`；按 5877 s/步，总时间约 `1.15×10⁹ s`，约 36 年。

v2412 在每个时间步开头调用声学：`Time::run()` → `functionObjectList::execute()` 对每个 function object 先 `execute()` 再 `write()`。`FfowcsWilliamsHawkings::execute()` 为空，采样在 `AcousticAnalogy::write()` → `correct()`。因此日志顺序是：上一步流动的 `ExecutionTime`，然后 `writing force and moment files` 与 `Starting acoustics probe`，然后才是下一步的 Courant 数和 `Time =`。一步的 5877 s 等于本步声学探针加上本步流动；流动约 150–200 s，探针约 5700 s。

两个函数对象各跑一遍。每遍 18 行 `OAP`，即 6 个麦克风 × 3 个源面。日志里 `p' = 0` 是因为物理时间仍在约 `10⁻⁵ s`，而面到麦克风的传播时间 `R/c0` 为数毫秒；`valueAt` 在到达时间之前返回 0，但这个零结果仍然走完全部分配。

### 2.2 核心瓶颈：每个面重新分配整张导数表

时间导数在 `AcousticAnalogy::tdPrimitive::dot`（`Sources/lib/analogies/AcousticAnalogy/AcousticAnalogy.H`）。`old_`、`oldOld_` 是 `List<T>`，没有容量。`controlDict` 中 `adjustTimeStep yes`，`Farassat1AFormulation::initialize` 对每个麦克风、每个面调用 `setTimeStepVariable(true)`。于是每个面都走进：

```cpp
if ((oldOld_.size() <= idx) || (varTimeStep_))
{
    dvdt = (new_[idx] - old_[idx]) / deltat_;
    oldOld_.resize(idx + 1);
    ...
}
```

面循环下标为 `0, 1, 2, …, N-1`。某步第一次调用时表长为 `N`，`resize(1)` 把它缩成 1，之后每个面再长 1。OpenFOAM v2412 的 `List::doResize`（`src/OpenFOAM/containers/Lists/List/List.C`）在长度变化时始终：

```cpp
T* old = this->v_;
this->size_ = len;
this->v_ = new T[len];
std::move(old, (old + overlap), this->v_);
delete[] old;
```

`new T[len]` 会默认构造全部元素。一个面上的分配次数是 `1+2+…+N ≈ N²/2`。

`calculateAcousticPressure` 对每个面调用三次 `dot`：`Un_`（标量）、`L_`（矢量）、`M_`（矢量）。这三份量不依赖麦克风，存储却是 `List<List<tdPrimitive>>`，下标为 `[iObs][iSurf]`，6 个麦克风各算各存一遍。两个函数对象再乘 2。

第一步走的是 `old_.resize` 的同样路径（`old_` 从空表逐面长到 `N`）。之后每步因为 `varTimeStep_` 为真，走 `oldOld_.resize`。所以从第 2 步起步长就平坦，与历史长度无关；这和日志里 5792–5956 s 的稳定间隔一致。

数量级（rank 3 上 `laptop_wall` 的 189609 个面，一个麦克风，一个 24 字节矢量）：

```text
189609 × 189610 / 2 × 24 ≈ 400 GB 的分配、搬移和释放
```

两个矢量场约 800 GB，标量 `Un_` 再约 140 GB。6 个麦克风 × 两个函数对象之后，临界路径上是数 TB 的 `malloc` / 搬移 / `free`。Farassat 1A 每个面只有几十次浮点运算，被分配淹没。有效内存流量若按约 2 GB/s（大量小块分配，到不了带宽峰值），这一项就是大约一小时量级，与实测约 95 分钟的探针时间相符。

### 2.3 并行与负载

声压积分使用各 rank 的局部面，再对两个标量求和：

```cpp
reduce(intDotQdS_.value(iObs), sumOp<scalar>());
reduce(intFdS_.value(iObs), sumOp<scalar>());
```

调用在 `FfowcsWilliamsHawkings::correct` 里，顺序是每个麦克风、每个源面各一次。面与面之间全体同步。归约本身只有一个 `scalar`，时间可忽略；等待的是最慢的 rank。

scotch 平衡的是单元。源面分布（各 `processor*/constant/polyMesh/boundary` 的 `nFaces`，8 rank）：

| rank | laptop_wall | impeller1 | impeller2 |
|---:|---:|---:|---:|
| 0 | 81305 | 0 | 0 |
| 1 | 158701 | 0 | 0 |
| 2 | 95367 | 0 | 0 |
| 3 | 189609 | 0 | 0 |
| 4 | 32630 | 0 | 46711 |
| 5 | 0 | 0 | 67222 |
| 6 | 19054 | 49434 | 0 |
| 7 | 0 | 64499 | 0 |

`laptop_wall` 的平方工作量正比于 `Σ Nᵢ² ≈ 7.83×10¹⁰`。rank 3 一项 `189609² ≈ 3.60×10¹⁰`，占 46%。单核做完全部 576666 个面是 `3.33×10¹¹`，8 核在这张面上的加速比约 **4.3**，效率约一半。rank 5 和 rank 7 在这张面上没有面，整个壳体阶段都在等。

叶轮只落在两个 rank 上（impeller1 在 6 和 7，impeller2 在 4 和 5），而且排在 `laptop_wall` 之前。三个面之间有 `reduce`，叶轮时间和壳体时间相加，不能重叠。临界路径上，叶轮的 `N²` 约为壳体的四分之一，壳体约占八成。

没有按麦克风或按面的线程并行；麦克风循环在 rank 内部串行。

`nonUniformSurfaceMotion true` 使每步 `expire()` 再 `update()`，并对三个面调用 `PatchTools::gatherAndMerge`。合并结果只写入 `mergeList_`；声压使用的是局部 `surf.Sf()`、`surf.Cf()`，这次聚集不进入积分，但每步所有 rank 都参加通信。`update()` 里还有每个麦克风一次 `reduce(rMax)` / `reduce(rMin)`，以及每个面一次 `gSum(Cf)`，相对平方分配是次要开销。

`correct()` 把 `surfaceVelocity`、`surfaceDensity`、`surfacePressure` 放在麦克风循环内部。采样结果与麦克风无关，每个面仍采样 6 次，并且每次 `interpolation<Type>::New("cell", ...)`。`interpolate false`，采样只是从边界取值。这一项是数十秒量级，不是 5700 秒的来源。

### 2.4 平方项去掉之后的内存墙

每个麦克风、每个面、厚度项 `qds_` 和载荷项 `fds_` 各存一对 `DynamicList<scalar>`（时间与数值），类型为 `List<List<List<Pair<DynamicList<scalar>>>>>`，下标 `[iObs][iSurf][iFace]`。

`fixedResponseDelay false` 时，过期时间为 `ct - tauMax`，`tauMax = R_max / c0`。麦克风在壳体包络外约 0.5 m，最远距离大约 0.5–1.5 m，`c0 = 343`，窗口约 3–4 ms。除以 `5.13×10⁻⁷` s，大约 6000–8000 个样本。

`fixedResponseDelay true` 时，过期时间为 `ct - (tauMax - tauMin + responseDelay)`，`responseDelay = 0.01`。窗口至少再增加 `0.01 / 5.13e-7 ≈ 2×10⁴` 个样本。

`cleanFreq 100` 每 100 步调用 `clearExpiredData` → `findExpiredIndex` → `getNewPointData`。它只删窗口之外的样本，并把剩余历史整表拷贝出来；窗口长度不变。

rank 3、189609 个面、6 个麦克风、4 条双精度序列：

```text
189609 × 6 × 4 × 8 × K ≈ 3.64×10⁷ × K 字节
```

`K ≈ 7800` 时一个函数对象约 280 GB。两个对象再翻倍，而且这只是一个 rank。窗口填满之前进程会先耗尽内存。`valueAt` 对每个面从头二分，历史变长后是 `N log K`，相对平方分配仍然很小。

在信号到达之前（约 `R_min/c0 ≈ 1.5` ms，约 3000 步），`valueAt` 一直返回 0。按当时速度，这 3000 步的零结果就要大约半年。

### 2.5 调用链

```text
Time::run()
  functionObjectList::execute()
    AcousticAnalogy::write()          // 每步都进
      forces::write()                 // "writing force and moment files"
      correct()                       // "Starting acoustics probe"
        expire(); update()            // 动网格：重算 Cf，gatherAndMerge
        vS = (Cf - Cf0) / dt
        fwhFormulation::update()      // 重算 robs、法向、tobs
        for iObs
          for iSurf
            sample U, rho, p          // 与 iObs 无关，仍在循环内
            Farassat1AFormulation::observerAcousticPressure
              calculateAcousticPressure   // dot() → List::resize，O(N²)
              valueAt(qds), valueAt(fds)  // 二分历史
              reduce 两个 scalar
        clearExpiredData()            // 每 cleanFreq 步
```

相关文件：

| 文件 | 作用 |
|---|---|
| `Sources/lib/analogies/AcousticAnalogy/AcousticAnalogy.H` | `tdPrimitive::dot`，平方分配 |
| `Sources/lib/analogies/AcousticAnalogy/AcousticAnalogy.C` | `write()`，每步进入 `correct()` |
| `Sources/lib/analogies/FfowcsWilliamsHawkings/FfowcsWilliamsHawkings.C` | 采样循环、`gatherAndMerge` |
| `Sources/lib/analogies/FfowcsWilliamsHawkings/Farassat1AFormulation.C` | 1A 积分与 `reduce` |
| `Sources/lib/analogies/FfowcsWilliamsHawkings/fwhFormulation.C` | `valueAt`、历史裁剪、距离与法向 |
| `src/OpenFOAM/containers/Lists/List/List.C` | `doResize` 每次整表分配 |

---

## 3. 性能修改（已合入）

### 3.1 导数历史一次分配

`adjustTimeStep yes` 让每个面都执行 `oldOld_.resize(idx+1)`。`List::doResize` 没有容量，长度从 `N` 缩到 1 再逐面长回去，一个面的代价是 `N²/2` 次分配和搬移。rank 3 的 `laptop_wall` 有 189609 个面，再乘 6 个麦克风和两个函数对象，就是现在这约 95 分钟。

`resize()` 里一次配齐三张表，`dot()` 里只改元素；变步长仍用一阶差商，固定步长仍用原来的三点公式，数值与慢路径一致。

```cpp
void resize(const label n)
{
    new_.resize(n);
    old_.resize(n);
    oldOld_.resize(n);
}

T dot(const scalar newt, const label idx)
{
    if (newt > prevt_)
    {
        deltat_ = newt - prevt_;
        prevt_ = newt;
    }

    T dvdt;
    if (!haveOld_)
    {
        dvdt = new_[idx] / deltat_;
    }
    else if (varTimeStep_ || !haveOldOld_)
    {
        dvdt = (new_[idx] - old_[idx]) / deltat_;
    }
    else
    {
        dvdt = (3.0*new_[idx] - 4.0*old_[idx] + oldOld_[idx])
              / (2.0*deltat_);
    }

    oldOld_[idx] = old_[idx];
    old_[idx] = new_[idx];

    if (idx == new_.size() - 1)
    {
        if (haveOld_) haveOldOld_ = true;
        haveOld_ = true;
    }
    return dvdt;
}
```

`haveOld_`、`haveOldOld_` 初始为 false；面数变化时在 `resize()` 里清掉这两个标志。改动只在 `AcousticAnalogy.H` 的 `tdPrimitive`，重编 `libAcoustics.dll`/`libAcoustics.so` 后两个函数对象一起生效。

这一步之后，导数循环从数 TB 的分配变成约 `10⁸` 次浮点运算，rank 3 上不到 1 秒；一步时间回到流动的约 150–200 秒，再加上仍在的线性声学开销。

### 3.2 采样、归约与 `gatherAndMerge`

平方项消失后，这些是声学里还能看见的部分：

- `L_`、`M_`、`Un_` 与麦克风无关，却按 6 个麦克风各存各算。现在每个源面只留一份 `tdPrimitive`。
- `surfaceVelocity`、`surfaceDensity`、`surfacePressure` 原来放在麦克风循环里面，每个面采样 6 次。现在提到 `iObs` 循环外，每个面采样一次。
- `nonUniformSurfaceMotion` 原来每步对三个面做 `PatchTools::gatherAndMerge`，结果只进 `mergeList_`，积分用的是本地 `Sf()`、`Cf()`。现在积分路径上不再做这次聚集。
- Farassat 的 `p'` 写各面之和；日志里每个面的 `OAP` 仍单独打印。

### 3.3 Farassat 改为每麦克风一条观测时间网格

每个麦克风、每个面、厚度项和载荷项原来各有一条 `DynamicList`。保留窗口约 `R/c0`，以 `deltaT ≈ 5.13×10⁻⁷` 计是 6000–8000 个样本；`fixedResponseDelay` 还要加上 `responseDelay 0.01`。rank 3 上一个函数对象窗口填满时约 280 GB，两个对象再翻倍；只修 `resize` 的话，跑几千步后会先耗尽内存。

延迟时刻的面积分对时间插值是线性运算，可以把每个面的贡献直接散到该麦克风的到达时刻上，不必保存每个面的历史。每个麦克风一条公共时间轴，长度约 `(endTime + R_max/c0) / dt_ac`。`endTime = 0.1`、声学步长 `10⁻⁵ s` 时，6 个麦克风一共几 MB。各 rank 先在本地累加，再对这条短数组做一次求和。`valueAt` 的逐面二分和 `getNewPointData` 的整表拷贝随之消失。

到达时刻 `t + R/c0` 落在两个相邻声学时刻之间时，按线性权重分到这两个点上；这与原来「先对每个面插值再求和」是同一种线性插值，结果不会逐位相同。

### 3.4 预期

| 措施 | 一步时间 | 能否跑到 0.1 s |
|---|---|---|
| 现状（最初 v2412） | 约 5900 s | 不能，约 36 年，且内存会先满 |
| 关闭 `adjustTimeStep` | 前两步仍约 5900 s，之后降到数分钟 | 不能，历史内存仍在 |
| 修好 `resize`，并去掉重复采样和 `gatherAndMerge` | 约 3–4 分钟，由流动决定 | 不能，数千步后内存满 |
| 再改成麦克风信号累加 | 仍约 3–4 分钟 | 能，声学内存为数 MB |
| 再把 `probeFrequency` 提到 20、只留一个函数对象 | 流动时间不变，声学按倍数下降 | 能 |

算例侧的 `probeFrequency`、函数对象个数和 `laptop_wall` 仍由算例字典决定，库没有改这些默认值。可选关键字 `acousticDeltaT` 为 0 时，Farassat 输出网格采用第一次声学采样的 `deltaT`。

---

## 4. 正确性修改（已合入）

### 4.1 面速度：一个 `deltaT` 的位移

`nonUniformSurfaceMotion` 使用

```text
vS = (当前面心 - Cf0) / deltaT
```

最初 `Cf0` 只在跳过声学采样的那些步里更新（`probeI % probeFrequency != 0`）。`probeFrequency` 为 1 时没有跳过的步，`Cf0` 停在初始化时的位置，也就是 `t = 0`。第一步的位移还是一个 `deltaT`；从第二步起，分子变成从初始位置到现在的总位移，除数仍是一个 `deltaT`。这个速度不是表面在该步的运动速度。

叶轮转速约 `628 rad/s`。总位移随时间增大后，虚假马赫数会穿过 1；Farassat 的 `1/(1-Mr)` 随之发散。`openfoam_fwh_m3`（`adjustTimeStep no`，`deltaT = 1e-5`）里可以看到：到 `t = 0.005` 最大 Courant 数仍约 5，静止的 `laptop_wall` 声压停留在数帕到约 20 Pa，两个叶轮从 `t ≈ 0.0019` 的 1 Pa 量级在 `t = 0.00219` 升到 `1e7`–`1e10`，随后到 `1e17`。

修复：每次采样之后把 `Cf0` 写成当前面心；采样窗口之外网格仍在运动时也更新 `Cf0`。这样面速度就是这一个时间步的位移，不再是 `t = 0` 以来的总位移。

### 4.2 多个积分面：`p'` 写各面之和

`SoundObserver::apressure(scalar)` 会覆盖该麦克风当前的 `p'`。最初的 `correct()` 在每个源面里都调用它，因此 `acousticData` 里的 `p'` 只剩循环中的最后一个面。日志里的逐面 `OAP` 仍是各面自己的值；FW-H 要求的是全部积分面上的和。

修复：Farassat 1A 在 `correct()` 结束时把各面的 `OAP` 求和后写入 `apressure`；GT 仍按原来的逐面 `valueAt` 求和。日志里每个面的 `OAP` 仍单独打印。

### 4.3 Windows/GCC 7 兼容

Windows 上的 OpenFOAM 自带 GCC 7 的 `libstdc++`，没有 GCC 13 MinGW 从 `<iostream>` 引用的 `std::ios_base_library_init`。`mingwStdCompat.C` 在本库内提供该符号，只为让 DLL 能被加载，不参与声压计算。

### 4.4 与最初算法的差异

- 厚度项和载荷项仍是原来的 Farassat 1A，查询时刻仍是 `ct + deltaT×1e-6`。
- 导数仍是第一步 `新值/deltaT`，变步长一阶差分，定步长从第三步起 BDF2；每个面只推进一次。
- 最初每个面保存整段 `(时间, q)` 和 `(时间, f)`，读声压时在查询时刻对每个面线性插值再求和。现在把 `q+f` 涂到间距为 `acousticDeltaT` 的均匀时间网格上再读取。网格步长等于流动步长、查询落在节点上时，节点值等于各面线性插值之和；查询比节点晚 `deltaT×1e-6`，带进相邻节点的比例是 `1e-6`。网格比流动步长粗时，输出信号的时间分辨率低于最初的逐面历史。
- 结果文件中的 `p'` 现在是各面之和；日志里每个面的 `OAP` 仍单独打印。
- 每次采样之后把 `Cf0` 写成当前面心；采样窗口之外也更新，避免窗口开始后的第一步把多步位移除以一个 `deltaT`。
- `PatchTools::gatherAndMerge` 已从积分路径去掉；积分本来就用本进程上的面，再对声压标量做 `reduce`。
- 体积四极子项最初的库就没有，这次也没有加。

---

## 5. patch 表面法向定向（本次修复，重点）

### 5.1 现象：加入 `laptop_wall` 后总声压增大很多

`openfoam_fwh_m4_all` 在 `system/fwhControl` 中把 7 个静态壁面合成 `laptop_wall`，与两个叶轮一起作为 Farassat 1A 的源面：

```foam
laptop_wall
{
    type            patch;
    patches
    (
        air_to_case2
        air_to_case1
        air_to_Cover
        air_to_fin_1
        air_to_fin_2
        air_to_CPU
        air_to_Cu_block
    );
    interpolate     false;
}
```

从数据上看，wall 的贡献远大于理论预期。下表是 **0.014–0.0164 s 的早期快照**；wall 的均值随时间漂移，稳态窗口（0.03–0.0529 s）的数值不同（旧库 wall RMS `0.2675 Pa`、均值 `−0.1915 Pa`）。新旧对比必须以同一窗口为准，见 §9.1。

| 量 | wall | impeller1 | impeller2 | 叶轮之和 |
|---|---:|---:|---:|---:|
| 面积 / m² | 0.19942 | 0.003397 | 0.003397 | 0.006794 |
| OAP RMS（0.014–0.0164 s）/ Pa | **0.197** | 0.022 | 0.019 | 0.059 |
| 5–20 kHz RMS / Pa | **0.193** | 0.010 | 0.011 | 0.014 |
| OAP 均值 / Pa | -0.051 | — | — | +0.044 |

- 总 `p'` RMS：加 wall **0.194 Pa**，叶轮单独 **0.0595 Pa**，约 **+10.3 dB**。
- 5–20 kHz：约 **+23 dB**（相对叶轮之和；相对单个叶轮约 +25 dB）。
- wall 与叶轮 OAP 的相关系数 **−0.007**（早期窗口），说明 wall 不是叶轮力的反作用，而是一个独立的宽频信号；稳态窗口仍是 −0.04 量级（§9.4.2）。
- 在相同时间窗内，两个算例的叶轮 OAP 完全一致（0.0590 vs 0.0595，均值都是 0.0442），所以增量全部来自 wall。

此外，wall 的均值随 `responseDelay` 变化很大：

| 对象 | wall 均值 / Pa |
|---|---:|
| `FwhImpellers`（无延迟） | -0.051 |
| `FwhImpellersFixedDelay`（0.01 s） | **-0.576** |
| `FwhImpellersFixedDelay2`（0.002 s） | -0.066 |

同一个壁面、同一个流场，只改延迟参数，早期直流量就差 10 倍；这说明绝对压力/稳态载荷被积分了，不是物理声压。到 0.03 s 后，三者的 wall 均值收敛到 −0.14…−0.19 Pa（旧）/ −0.075 Pa（新），早期 10 倍差异主要来自启动段和 `tauMin` 时基（§9.4）。

### 5.2 源码路径：质心启发式

法向在 `fwhFormulation::initialize()` 和 `fwhFormulation::update()` 中计算（修改前）：

```cpp
//calculate normals
List<vector> Cs(fwh_.controlSurfaces_.size());
ni_.resize(fwh_.controlSurfaces_.size());
nl_.resize(fwh_.controlSurfaces_.size());
forAll(ni_, iSurf)
{
    const vectorField& Sf = fwh_.controlSurfaces_[iSurf].Sf();
    const vectorField& Cf = fwh_.controlSurfaces_[iSurf].Cf();
    ni_[iSurf].resize(Cf.size());
    nl_[iSurf].resize(Cf.size());
    Cs[iSurf] = gSum(Cf);
    scalar surfSize = scalar(Cf.size());
    reduce(surfSize, sumOp<scalar>());
    Cs[iSurf] /= surfSize;

    scalar magSf = 0.0;
    forAll(ni_[iSurf], iFace)
    {
        magSf = mag(Sf[iFace]);
        ni_[iSurf].value(iFace) = Sf[iFace]/magSf;

        if ( ((Cf[iFace] - Cs[iSurf]) & ni_[iSurf].value(iFace)) > 0 )
        {
            nl_[iSurf][iFace] = 1.0;
        }
        else
        {
            nl_[iSurf][iFace] = -1.0;
        }
        ni_[iSurf].value(iFace) *= nl_[iSurf][iFace];
    }
}
```

这个启发式假设：**表面围住固体，且 `Cs` 在固体内部**。对 `type patch` 的壁面，`Sf` 已经从流体指向固体，正确的 FW-H 法向就是：

```text
n = -Sf/|Sf|
```

根本不需要质心。

`Farassat1AFormulation::calculateAcousticPressure()` 直接使用这个 `n`：

```cpp
const vector n = ni_[iSurf].value(iFace);
...
L = Pf & n;
Un = U & n;
```

所以 85.5% 的面 `n` 反了，`L = p n` 和 `Un = U n` 也就反了；闭合壁面本该相消的压力积分被破坏。

### 5.3 几何验证

用 `constant/polyMesh` 的 `points`、`faces`、`boundary` 直接计算每个 patch 的 `Cf`、`Sf`、面积和启发式翻转比例（`Sf` 为从流体指向固体的 patch 法向）：

| patch | 面积 / m² | `ΣSf_y` | `Σn_heur_y` | 正确 `-ΣSf_y` | 未翻转（`flip+`） |
|---|---:|---:|---:|---:|---:|
| `air_to_Cover` | **0.14545** | +0.01759 | **+0.01759** | **-0.01759** | 44.0% |
| `air_to_case1` | 0.01651 | -0.00640 | +0.00187 | +0.00640 | 95.8% |
| `air_to_case2` | 0.01651 | -0.00640 | +0.00204 | +0.00640 | 95.8% |
| `air_to_fin_1` | 0.00930 | -0.00135 | +0.00135 | +0.00135 | 52.9% |
| `air_to_fin_2` | 0.00930 | -0.00135 | +0.00135 | +0.00135 | 52.8% |
| `air_to_CPU` | 0.00105 | 小 | 小 | 小 | 89.5% |
| `air_to_Cu_block` | 0.00131 | 小 | 小 | 小 | 81.4% |
| **wall 合计** | **0.19942** | **0**（闭合） | **+0.0224** | **0** | **85.5%** |

关键结论：

1. **7 个 patch 合计是闭合面**：`Σ Sf = 0`（数值上 ~1e-14）。
2. **正确法向 `-Sf` 的和严格为 0**，闭合面的一致性成立。
3. **启发式的 `Σ n dS = 0.0224 m²`**，是 wall 面积的 11%；闭合面的正确法向必须满足 `Σ n dS = 0`，所以这个法向场**数学上不闭合**。
4. `air_to_Cover` 的 `ΣSf_y = +0.01759 m²`，启发式几乎原样保留（`+0.01759`），正确应为 `-0.01759`；cover 占 wall 面积的 **73%**，是误差的主要来源。
5. 叶轮：`flip+` 只有 0.17%，`Σ n_heur dS ≈ 1.5e-6 m²`，是自洽的；所以叶轮结果看起来正常。

为什么质心会落在流体侧：`Cs` 用的是**未加权面心平均**。wall 的 576,666 个面中，两个 case patch 占 448,344 个（78%），但面积只有 0.033 m²；cover 只有 54,391 个面，却占 0.145 m²（73%）。面数加权使 `Cs` 落在两扇之间、靠近 cover 的流体侧（约 `(0.0062, -0.0039, 0.0731)`）。对大多数壁面，`(Cf - Cs)` 与 `Sf` 同向（都指向固体），启发式就不翻，法向指向固体内部；对凹面/背风面，符号又反过来，只有 14.5% 被翻。于是同一张闭合面上的法向一半对、一半错，`Σ n dS` 不为零。

### 5.4 为什么扇叶没有这个问题

两个叶轮各自是**封闭凸体**：`Sf` 从流体指向叶轮固体，`Cs` 在叶轮内部；`(Cf - Cs)` 从固体指向流体，与 `Sf` 反向，启发式翻成 `-Sf`，即从固体指向流体，正好是 FW-H 需要的外法向。几何检查也显示叶轮 `flip+` 只有 0.17%，`Σ n dS` 约 1e-6 m²。

所以问题不是“Farassat 公式错了”，而是**同一个启发式对“围住固体的凸面”有效，对“围住流体/非凸的壁面组”失效**。`laptop_wall` 是后者。

### 5.5 修复：patch 表面直接用 `-Sf`

修改 `Sources/lib/analogies/FfowcsWilliamsHawkings/fwhFormulation.C` 和 `.H`：

- 把法向计算抽成 `updateNormals()`；
- `initialize()` 和 `update()` 都调用它；
- 对 `isA<sampledPatch>(surf)` 的表面：

```cpp
const scalar magSf = mag(Sf[iFace]);
ni_[iSurf].value(iFace) = -Sf[iFace]/magSf;
nl_[iSurf][iFace] = -1.0;
```

- 对非 patch 表面（STL/sampledSurface）保留原来的质心启发式：

```cpp
vector Cs = gSum(Cf);
scalar surfSize = scalar(Cf.size());
reduce(surfSize, sumOp<scalar>());
Cs /= surfSize;

forAll(Cf, iFace)
{
    const scalar magSf = mag(Sf[iFace]);
    ni_[iSurf].value(iFace) = Sf[iFace]/magSf;
    nl_[iSurf][iFace] =
        (((Cf[iFace] - Cs) & ni_[iSurf].value(iFace)) > 0) ? 1.0 : -1.0;
    ni_[iSurf].value(iFace) *= nl_[iSurf][iFace];
}
```

这样：

- `laptop_wall` 的 7 个 patch 全部用 `-Sf`，`Σ n dS = 0`，闭合壁面的稳态/均匀压力分量相消；
- 两个叶轮也用 `-Sf`，和原启发式的结果一致（原来就翻对了 99.8%）；
- 只有非 patch 的 STL 表面继续用质心启发式。

**注意**：这里假设 `type patch` 选的是**固体壁面**。如果以后把 patch 当作可渗透面（例如风扇入口平面），`-Sf` 的方向会相反，需要另外加 `normalOrientation`/`permeable` 选项区分。当前库的 FW-H 用法（圆柱、叶轮、laptop_wall）都是固体壁面，因此先按固体壁面处理。

### 5.6 与 Cradle scFLOW 的差异

scFLOW 的壁面增量小，可能来自以下一项或多项，不能直接等同于“OpenFOAM 物理算错了”：

1. **法向处理**：scFLOW 用物面外法向（从固体指向流体），壁面载荷与叶片载荷反相抵消；修复前的库有 85.5% 壁面反了，变成相加（本次已修，见 §5.5；实测只改变了低频/直流，见 §9）。
2. **表面集**：scFLOW 的 FW-H 面是否包含 cover、内部 CPU/Cu block、fins？`laptop_wall` 的 cover 占 73% 面积；如果 scFLOW 只取风扇/壳体局部，增量自然小。
3. **封闭/可渗透面**：scFLOW 可能用包围风扇的封闭可渗透面，而不是所有固体壁面；自由场 1/r 积分对内部壁面也没有遮挡/反射，会高估。
4. **网格/湍流**：cover 的平均面尺度约 1.6 mm，case 约 0.27 mm，叶轮约 0.17 mm；cover 的网格最粗，5–20 kHz 的壁面压力脉动可能包含未解析的湍流/数值分量。当前算例的 max Courant 实测约 4.5–11（`adjustTimeStep no`，`deltaT = 1e-5`），也会放大高频噪声。

因此，先修法向，再和 scFLOW 对齐表面集、网格和 CFL，才有可比性。

### 5.7 仍建议的后续

法向修复是必要条件，但不一定能把 wall 降到“理论上的小量”，还需要：

> **验证后的状态（2026-11-17）**：下面第 1–4 条仍未做；第 5 条已完成——`_all_m1` 已用新库重跑。实测表明第 3、4 条是当前高频问题的首要候选，见 §9.7–9.8。

- **壁面压力去均值/用脉动量**：`surfacePressure` 只乘 `rho`、减 `pInf = 0`（`FfowcsWilliamsHawkings.C`），没有减时均。对闭合面且法向正确时，均匀分量会相消；非均匀稳态载荷仍可能产生直流。建议按面减时均，或直接采 `p'`。
- **`tauMin` 全局耦合**：`fwhFormulation::update()` 的 `rMin_/tauMin_` 对所有 `controlSurfaces_` 取最小值，`tobs = ct + R/c0 - (tauMin - responseDelay)`。加入 wall 后，wall 上离麦克风更近的面把 `tauMin` 拉小，所有叶轮面的 `tobs` 都变。这解释了 `FwhImpellersFixedDelay` 的叶轮 OAP 在早期从 1.06 变成 0.40——是时基变了，不是物理。建议按面/按 patch 取 `tauMin`，或先用 `FwhImpellers`（`fixedResponseDelay false`）做对比。
- **cover 网格/CFL**：法向修完后 5–20 kHz 壁面仍大（本次已确认，几乎没有下降），把 `deltaT` 降到约 `2e-6`（或开 `adjustTimeStep`、`maxCo 0.6`）复算一小段，确认高频分量是不是数值。
- **表面集对齐**：确认 scFLOW 是否包含 cover/内部壁面、是否用封闭可渗透面；否则“增量大小”不是同一口径。
- **旧 `_all` 运行**：它加载的是旧库，wall 项是错的；已用新库从 0 重跑为 `_all_m1`（2026-11-17，48 rank），第 9 节即其公共窗口结果。旧 `_all` 只保留作对照，不再作为最终结论。

---

## 6. 修改清单

### 6.1 按主题

**性能**

| 修改 | 位置 |
|---|---|
| `tdPrimitive` 一次分配导数历史，差分公式不变 | `AcousticAnalogy.H` |
| 每个源面只采样一次，`L`/`M`/`Un` 每步只推进一次 | `FfowcsWilliamsHawkings.C`、`Farassat1AFormulation.C` |
| 去掉积分路径上的 `gatherAndMerge` | `FfowcsWilliamsHawkings.C` |
| Farassat 改为观测时间网格，不再保存逐面 `qds_`/`fds_` | `Farassat1AFormulation.C`、`fwhFormulation.C` |
| 可选 `acousticDeltaT`，0 表示采用第一次采样的 `deltaT` | `FfowcsWilliamsHawkings.C` |
| Farassat 的 `p'` 写各面之和 | `FfowcsWilliamsHawkings.C` |
| 采样后更新 `Cf0`；采样窗口外也更新 | `FfowcsWilliamsHawkings.C`、`AcousticAnalogy.C` |
| 提供 `ios_base_library_init`，供 GCC 7 运行库加载 | `mingwStdCompat.C` |

**正确性（本次新增）**

| 修改 | 位置 |
|---|---|
| patch 表面法向直接用 `-Sf`，不再用质心启发式 | `fwhFormulation.C`、`fwhFormulation.H` |

### 6.2 按文件

| 文件 | 本次改动 |
|---|---|
| `Sources/lib/analogies/FfowcsWilliamsHawkings/fwhFormulation.H` | 新增 `updateNormals()` 声明 |
| `Sources/lib/analogies/FfowcsWilliamsHawkings/fwhFormulation.C` | 新增 `updateNormals()`；`initialize()`/`update()` 改为调用它；`sampledPatch` 用 `-Sf`；非 patch 保留质心 |
| `libAcoustics_analysis.md` | 本文：补充法向分析、几何验证、逐面 OAP、SPL/A 计权/相干性实测与信息融合，并按章节重排 |

### 6.3 未改动的部分

- Farassat 1A 的厚度项、载荷项、`1/(1-Mr)`、到达时间 `ct + R/c0`、除以 `4π` 不变。
- `fwhFormulation` 的 `tobs`/`tauMin` 逻辑不变（后续可单独改）。
- 算例字典不变；`laptop_wall` 仍是 7 个 patch 的 `type patch`。
- 体积四极子项仍未加。

---

## 7. 验证与复现

### 7.1 编译

```bash
source /home/ubuntu/programs/OpenFOAM-v2412/etc/bashrc
cd /home/ubuntu/programs/libAcoustics/Sources/lib
wmake -j 8
```

本次为了不打断正在运行的 48 rank 算例（它仍映射着旧的 `libAcoustics.so`），把 `Make/files` 临时改成

```text
LIB = $(FOAM_USER_LIBBIN)/libAcoustics_fix
```

编译通过，产物为 `libAcoustics_fix.so`；提交前已把 `Make/files` 改回 `libAcoustics`。2026-11-16 23:54 已用 `makeLibrary.sh` 正式生成并安装 `libAcoustics.so`，第 9 节的 `_m1` 运行加载的就是它（`/proc/<pid>/maps` 确认）。

另用 `Tutorials/vortexShed` 的 `cylinder` patch 做了 2 步单进程冒烟测试（`libs ("libAcoustics_fix.so")`）：FWH 对象正常加载、采样并输出 `OAP`，无崩溃。该测试只验证 `sampledPatch` 分支的运行路径；圆柱是凸体，新旧法向本来就一致，不能用来验证 `laptop_wall` 的符号。`laptop_wall` 的符号由 7.2 的几何检查和 9.2 的运行数据说明。

### 7.2 几何检查

检查对象：`constant/polyMesh` 的 `boundary`、`points`、`faces`。对每个 patch：

- `Cf`：面顶点平均；
- `Sf`：按面顶点顺序的 Newell 面积矢量；
- `ΣSf`、`Σn_heur dS`、`flip+` 比例。

关键结果：

```text
Σ Sf (wall)           = 0                 // 闭合
Σ (-Sf) dS (wall)     = 0                 // 正确法向
Σ n_heur dS (wall)    = 0.0224 m²         // 当前启发式
flip+ (wall)          = 85.5%             // 未翻转，仍指向固体
air_to_Cover ΣSf_y    = +0.01759 m²
air_to_Cover n_heur_y = +0.01759 m²       // 应翻成 -0.01759
impeller flip+        = 0.17%
impeller Σ n_heur dS  ≈ 1.5e-6 m²
```

### 7.3 运行对比

从 `log.pimpleFoam` 解析每个 `Starting acoustics probe` 后的 `OAP` 行。三个对象 × 三个源面，每个对象每个时间步 3 行 `OAP`，顺序是 `impeller1`、`impeller2`、`laptop_wall`。

`openfoam_fwh_m4_all`（旧库）在 0.014–0.0164 s：

| 量 | wall | impeller1 | impeller2 |
|---|---:|---:|---:|
| OAP RMS / Pa | 0.197 | 0.022 | 0.019 |
| 5–20 kHz RMS / Pa | 0.193 | 0.010 | 0.011 |

`openfoam_fwh_m4`（叶轮单独）同窗口叶轮之和 0.0595 Pa；all-wall 叶轮之和 0.0590 Pa，均值都是 0.0442。总 `p'` 0.194 vs 0.0595，约 +10.3 dB。

修复后的实测见第 9 节。先给结论：公共窗口 0.03–0.0529 s 下，wall 20–200 Hz 从 0.0133 降到 0.0013 Pa（−20.2 dB），但 5–20 kHz 从 0.1440 到 0.1449 Pa（+0.05 dB）；总 OASPL 77.46→76.88 dB，L_A 73.62→74.00 dBA。wall 仍约占总波动能量 97.5%，总级仍比叶轮单独高约 19 dB。第 9 节给出逐面 OAP、时域漂移、相关系数和能量预算。

### 7.4 复现步骤

1. 停掉加载旧库的运行，`wmake` 生成新的 `libAcoustics.so`（已完成）。
2. 从 0 或保存的时间重跑 `pimpleFoam -parallel`，保持 `system/fwhControl` 不变。本次已从 0 重跑为 `openfoam_fwh_m4_all_m1`（48 rank）。
3. 用 `parse OAP` 脚本比较 wall/impeller 的 OAP；先看 `FwhImpellers`（无延迟），再看 `FixedDelay`。本次结果见 §9。
4. 如果 wall 仍大（本次已确认，5–20 kHz 几乎不变）：把 `laptop_wall` 拆成 `cover`、`case`、`fins` 等独立 surface，定位是哪一块；同时检查 cover 的面尺度和 max Courant，并按 §9.8 做远场 observer 与 dt/CFL 复算。

---

## 8. 阶段性结论

1. **性能问题**的根因是 `tdPrimitive::dot` 的平方分配；当前分支已用一次分配、单次采样、观测时间网格解决，声学开销降到流动以下，内存从按面历史降到每麦克风一条信号。
2. **`laptop_wall` 的质心启发式确实错了**：7 个 wall patch 合计是闭合面，但启发式只翻了 14.5% 的面，`Σ n dS = 0.0224 m²`，cover 的净面积矢量没有被翻。修复为 `sampledPatch` 直接用 `-Sf/|Sf|` 后，几何上 `Σ n dS = 0`。
3. **2026-11-17 实测（公共窗口 0.03–0.0529 s）：法向修复只解决了一部分问题**：
   - wall 均值 `−0.1915 → −0.0765 Pa`（−60%）；20–200 Hz `0.0133 → 0.0013 Pa`（−20.2 dB），200 Hz–1 kHz −12.8 dB；
   - 5–20 kHz `0.1440 → 0.1449 Pa`（**+0.05 dB**）；20–50 kHz `0.1130 → 0.1095 Pa`（−0.27 dB）；
   - 总 OASPL `77.46 → 76.88 dB`（−0.58），`L_A` `73.62 → 74.00 dBA`（+0.38）；
   - wall AC RMS `0.1868 → 0.1835 Pa`，仍占总方差 **97.5%**，比叶轮两面之和高 **16 dB**；
   - 旧/新总声压相关系数 **−0.58**（wall −0.65），说明壁面项主要是反号，不是消失。
4. **机理**：低频空间相干，`Σ n dS` 是否为零决定相消/相加，所以修复有效；高频空间不相干，翻转法向只改符号图案、不改方差，所以 SPL 不变。5–20 kHz 的 wall 项是壁面压力载荷本身（cover 网格/CFL、近场、表面集、未解析高频）。
5. **结论**：法向 bug 已修复（几何闭合性 + 低频/直流 + 叶轮不变性）；但“wall 高频贡献过大、总声压高约 19 dB”的问题未解决。验收标准见 §9.8。
6. **仍待做**：拆分 wall、远场 observer、dt/CFL、按面去均值、`tauMin` 按面、与 scFLOW 对齐；旧 `_all` 的 wall 绝对值不能作为最终结论，叶轮单独稳态 AC 约 `0.04 Pa` 仍是可信基线（早期 OAP `0.0595 Pa` 含均值，口径不同）。

---

## 9. 实测复核与信息融合（2026-11-17）

> 本节把 §5 的几何/源码分析、§7.3 的早期 OAP 快照，与 2026-11-17 新库 `_m1` 的逐面 OAP、SPL、1/3 倍频程、A 计权、20–50 kHz、时域漂移和相干性结果合并。除特别注明外，所有新旧对比都用 **0.03–0.0529 s** 的公共稳态窗口，避免启动瞬态、不同结束时间和不同 `responseDelay` 时基造成的口径差。

### 9.1 对照实验设计与有效性

| 项目 | 旧库 `_all` | 新库 `_all_m1` | 叶轮单独基线 |
|---|---|---|---|
| 算例 | `openfoam_fwh_m4_all` | `openfoam_fwh_m4_all_m1` | `openfoam_fwh_m4` |
| 法向 | 质心启发式（`Σ n dS = 0.0224 m²`） | `sampledPatch` 用 `-Sf`（`Σ n dS = 0`） | 同旧库，但无 `laptop_wall` 源面 |
| 库 | 旧 `libAcoustics.so` | 2026-11-16 23:54 构建的 `libAcoustics.so` | 旧库 |
| 数据终点 | 0.05291 s（已停止） | ≥0.06247 s（仍在运行） | 0.10000 s |
| 声源面 | 2 叶轮 + `laptop_wall`（7 patch） | 同左 | 2 叶轮 |
| 对象 | `FwhImpellers`、`FixedDelay`、`FixedDelay2` | 同左 | 同左 |

**对照有效性检查：**

1. **新库确实生效**：`_m1` 进程的 `/proc/<pid>/maps` 指向 2026-11-16 23:54 的 `libAcoustics.so`；该库由 §5.5 的 `updateNormals()` 代码编译而来。
2. **流动与时间基一致**：同一窗口内两个叶轮面的 OAP 新旧逐位相同（`FwhImpellers`：impeller1 RMS 0.0331 Pa、均值 0.0249 Pa；impeller2 RMS 0.0320 Pa、均值 0.0256 Pa；`FixedDelay`/`FixedDelay2` 同样一致到 1e-4 Pa）。因此所有差异都可归因于 `laptop_wall`。
3. **窗口公共且稳态**：0–0.01 s 有启动瞬态（峰值可达 57.9 Pa），0.01–0.03 s 仍在建立；0.03 s 后 wall 均值和叶轮 OAP 都进入稳定段。旧数据只到 0.0529 s，所以公共窗口取 0.03–0.0529 s（2291 点，`dt = 1e-5 s`，`df ≈ 43.65 Hz`）。
4. **基线只用于量级**：`openfoam_fwh_m4` 是叶轮单独算例，网格/流场与 `_all` 不完全相同，不能做逐点相减；但它给出了“没有 wall 时总级应在什么量级”的独立参照。

**本文的量的定义**（避免第 9.2 节出现“RMS”口径混淆）：

- **时域 AC RMS**：`sqrt(mean((p−mean(p))²))`，含 0–50 kHz 全部频率；
- **PSD 带内 RMS**：`sqrt(∫_band PSD df)`，单边 Hann 窗；
- **窄带 SPL**：单根谱线幅值 `A/√2` 对 `p_ref`；
- **1/3 倍频程 / OASPL / L_A**：由 PSD 带内积分，A 计权按 IEC 61672-1、1 kHz 归一化；
- **OAP**：`log.pimpleFoam` 中每个源面单独输出的观测点声压贡献，不是三个面的总 `p'`；
- 参考声压 `p_ref = 20 µPa`。

### 9.2 预期签名与实测结果总表

把 §5 的几何结论、§5.7 的预期和本次实测放在同一张表里，可以直接看出“哪一部分被修好了、哪一部分没有”：

| 检查项 | §5 几何/源码预期 | 旧库实测 | 新库实测 | 判定 |
|---|---|---:|---:|---|
| wall 法向闭合性 `Σ n dS` | 0.0224 → 0 m² | 0.0224 m² | **0** | ✅ 已修 |
| wall OAP 均值 | 稳态/均匀载荷应相消 | −0.1915 Pa | **−0.0765 Pa**（−60%） | ⚠️ 部分 |
| wall 20–200 Hz RMS | 低频相干分量应大幅下降 | 0.0133 Pa | **0.0013 Pa（−20.2 dB）** | ✅ 已修 |
| wall 200 Hz–1 kHz RMS | 同上 | 0.0061 Pa | **0.0014 Pa（−12.8 dB）** | ✅ 已修 |
| wall 1–5 kHz RMS | 应小幅下降 | 0.0062 Pa | 0.0050 Pa（−1.9 dB） | ◐ 有限 |
| wall 5–20 kHz RMS | 若为法向问题则应下降 | 0.1440 Pa | **0.1449 Pa（+0.05 dB）** | ❌ 未变 |
| wall 20–50 kHz RMS | 同上 | 0.1130 Pa | 0.1095 Pa（−0.27 dB） | ❌ 未变 |
| 叶轮 OAP | 修复不应改变 | 0.0331 / 0.0320 Pa | **0.0331 / 0.0320 Pa** | ✅ 对照通过 |
| 总 `p'` 均值漂移 | 旧漂移应消失 | −0.10 → −0.18 Pa | **稳定在 −0.026 Pa** | ✅ 已修 |
| 总时域 AC RMS | 预期接近叶轮单独 ~0.04 Pa | 0.1916 Pa | **0.1863 Pa（−0.24 dB）** | ❌ 未达 |
| 总 OASPL 20–20 kHz | 预期接近 57.7 dB | 77.46 dB | **76.88 dB（−0.58 dB）** | ❌ 未达 |
| 总 5–20 kHz | 预期显著下降 | 77.40 dB | **76.87 dB（−0.53 dB）** | ❌ 未变 |
| 总 L_A 20–20 kHz | 预期下降 | 73.62 dBA | **74.00 dBA（+0.38 dB）** | ❌ 略升 |

这张表是全文的核心：**法向修复在低频/直流上完全符合几何预期，在高频上几乎无效**。下面给出每一行的原始数据和解释。

### 9.3 逐面 OAP：wall 的均值、波动与频带分解

#### 9.3.1 `FwhImpellers` 的 wall 频带分解

对 wall OAP 时程做与总声压相同的 Hann/PSD 处理（0.03–0.0529 s）：

| 量 | 旧库 `_all` | 新库 `_m1` | 变化 |
|---|---:|---:|---:|
| 均值 / Pa | −0.1915 | **−0.0765** | 直流降 60% |
| 时域 AC RMS / Pa | 0.1868 | **0.1835** | **−0.15 dB** |
| 0–20 Hz RMS / Pa | 0.0033 | 0.0000 | 漂移消失 |
| 20–200 Hz RMS / Pa | 0.0133 | **0.0013** | **−20.2 dB** |
| 200 Hz–1 kHz RMS / Pa | 0.0061 | **0.0014** | **−12.8 dB** |
| 1–5 kHz RMS / Pa | 0.0062 | 0.0050 | −1.9 dB |
| 5–20 kHz RMS / Pa | 0.1440 | **0.1449** | **+0.05 dB** |
| 20–50 kHz RMS / Pa | 0.1130 | 0.1095 | −0.27 dB |
| 20–20 kHz RMS / Pa | 0.1449 | 0.1450 | 0.00 dB |
| 20–50 kHz 全带 RMS / Pa | 0.1837 | 0.1818 | −0.09 dB |

> 注意：wall 的能量并不只在 5–20 kHz。20–50 kHz 的 RMS 也有 0.11 Pa 量级；也就是说，`dt = 1e-5 s`（`fs = 100 kHz`）下 20 kHz 以上的“不可听”高频占了 wall 波动的很大一部分。第 9.5 节的总级会看到这一点。

#### 9.3.2 三个对象的 wall 均值/波动

| 对象 | wall 均值 / Pa（旧→新） | wall RMS / Pa（旧→新） | wall **AC** RMS / Pa（旧→新） | AC 变化 |
|---|---:|---:|---:|---:|
| `FwhImpellers` | −0.1915 → −0.0765 | 0.2675 → 0.1989 | **0.1868 → 0.1835** | −0.15 dB |
| `FwhImpellersFixedDelay` | −0.1411 → −0.0747 | 0.2360 → 0.1920 | **0.1892 → 0.1768** | −0.59 dB |
| `FwhImpellersFixedDelay2` | −0.1872 → −0.0761 | 0.2606 → 0.1898 | **0.1813 → 0.1739** | −0.36 dB |

三个对象的规律一致：**均值大幅下降，AC 波动几乎不变**。RMS 的下降几乎全部来自均值（RMS² = mean² + AC²），不是来自脉动幅值。

#### 9.3.3 叶轮面的不变性

`FwhImpellers` 的逐面 OAP（Pa）：

| 面 | 旧：均值 / RMS | 新：均值 / RMS | 5–20 kHz RMS（旧=新） | 20–50 kHz RMS（旧=新） |
|---|---:|---:|---:|---:|
| impeller1 | 0.0249 / 0.0331 | 0.0249 / 0.0331 | 0.0075 | 0.0216 |
| impeller2 | 0.0256 / 0.0320 | 0.0256 / 0.0320 | 0.0081 | 0.0189 |
| 两面之和（方均根） | — | — | 0.0110 | 0.0287 |

`FixedDelay` 与 `FixedDelay2` 的叶轮面同样逐位一致。这是本次对比最重要的控制变量：**法向修复只动了 wall，没有动叶轮，也没有改流动**。

#### 9.3.4 能量预算

新库 `FwhImpellers`、0.03–0.0529 s：

```text
wall AC RMS                 = 0.1835 Pa
impeller1 AC RMS           = 0.0218 Pa
impeller2 AC RMS           = 0.0193 Pa
impeller 两面方均根        = 0.0291 Pa
总 AC RMS（实测）          = 0.1863 Pa
wall 占总方差比例          = 0.1835² / (0.1835² + 0.0291²) ≈ 97.5%
wall / impeller 幅值比     = 20 lg(0.1835 / 0.0291) ≈ +16.0 dB
wall 5–20k / impeller 5–20k = 20 lg(0.1449 / 0.0110) ≈ +22.4 dB
```

也就是说，**总声压的 97.5% 方差来自 wall**；把 wall 完全消掉，总 AC RMS 应落到 ~0.03 Pa（约 64 dB）量级，而不是现在的 0.186 Pa。这是判断“问题是否解决”的最直接能量指标。

### 9.4 时域漂移、相关性与“反号而非消失”

#### 9.4.1 总 `p'` 的分段均值

每 5 ms 一段的均值（Pa），`FwhImpellers`：

| 区间 / s | 旧库 `_all` | 新库 `_m1` | 叶轮单独 |
|---|---:|---:|---:|
| 0.030–0.035 | −0.097 | −0.024 | +0.051 |
| 0.035–0.040 | −0.132 | −0.027 | +0.052 |
| 0.040–0.045 | −0.146 | −0.026 | +0.049 |
| 0.045–0.050 | −0.176 | −0.027 | +0.049 |
| 0.050–0.0529 | −0.164 | −0.027 | +0.051 |

旧库的均值在 23 ms 内漂了约 −0.08 Pa，这是一个低频斜坡；它的能量落在 FFT 的第一根谱线（43.65 Hz）附近，表现为 40–130 Hz 的“虚假低频”。新库的均值稳定在 −0.026 Pa（残差是壁面非均匀稳态载荷，见 9.7 的“去均值”建议）。

#### 9.4.2 相关系数

| 比较 | `FwhImpellers` | `FixedDelay` | `FixedDelay2` |
|---|---:|---:|---:|
| 总 `p'`：旧 vs 新 | **−0.584** | −0.589 | −0.568 |
| wall OAP：旧 vs 新 | **−0.649** | −0.615 | −0.602 |
| wall vs 叶轮（新库） | −0.041 | −0.079 | −0.077 |
| wall vs 总 `p'`（新库） | 0.976 | 0.987 | 0.986 |

三个结论：

1. **旧/新负相关**：说明法向修复把 wall 信号中相当大的一部分**反了号**，而不是把它消掉；`−0.65` 而不是 `−1`，是因为旧法向只错了 85.5% 的面（cover 等），仍有 14.5% 原本就是对的。
2. **wall 与叶轮近似不相关**：wall 不是叶轮力的物理反作用，而是一个独立的壁面压力载荷信号。早期窗口的 −0.007 在稳态窗口仍然是 −0.04 量级。
3. **wall 与总 `p'` 相关 0.98**：总声压基本就是 wall 信号，再次说明“总级没降”是因为 wall 的波动没降。

#### 9.4.3 方差守恒的直观解释

如果 wall 的脉动在空间上完全不相干，那么把 85.5% 的面反号只会改变各面贡献的**符号图案**，不会改变总方差：

```text
E[ |Σ_i p_i n_i dS_i|² ] ≈ Σ_i E[p_i²] dS_i²     （交叉项因不相干而平均掉）
```

方差只依赖 `p_i²` 和面积，不依赖 `n_i` 的正负。所以高频 SPL 对法向符号不敏感——这正是 9.3.1 中 5–20 kHz 只变 +0.05 dB 的原因。相反，低频压力在壁面上空间相干，交叉项不可忽略，`Σ n dS` 是否为零决定了它们相消还是相加，所以法向修复对低频有效。**同一套法向修复在低频和高频表现出完全不同的效果，不是修复不稳定，而是两个频段的物理相干性不同。**

### 9.5 总级、A 计权与全频带

#### 9.5.1 三个对象的总级与四段频带

`FwhImpellers`：

| 指标 | 旧库 | 新库 | 叶轮单独 | Δ(新−旧) |
|---|---:|---:|---:|---:|
| OASPL 20 Hz–20 kHz / dB | 77.46 | 76.88 | 57.73 | **−0.58** |
| L_A 20 Hz–20 kHz / dBA | 73.62 | 74.00 | 56.83 | **+0.38** |
| 时域 AC RMS / Pa | 0.1916 | 0.1863 | 0.0401 | −0.24 dB |
| 20–100 Hz / dB | 56.71 | 28.82 | 36.18 | **−27.89** |
| 100 Hz–1 kHz / dB | 52.45 | 34.68 | 36.85 | **−17.78** |
| 1–5 kHz / dB | 51.87 | 47.93 | 44.31 | −3.94 |
| 5–20 kHz / dB | 77.40 | 76.87 | 57.46 | **−0.53** |

`FwhImpellersFixedDelay`：

| 指标 | 旧库 | 新库 | 叶轮单独 | Δ |
|---|---:|---:|---:|---:|
| OASPL / dB | 77.92 | 77.22 | 57.98 | −0.70 |
| L_A / dBA | 74.08 | 74.30 | 56.93 | +0.22 |
| 时域 AC RMS / Pa | 0.1921 | 0.1768 | 0.0415 | −0.72 dB |
| 20–100 Hz / dB | 58.80 | 33.43 | 41.69 | −25.37 |
| 100 Hz–1 kHz / dB | 52.70 | 34.97 | 38.34 | −17.72 |
| 1–5 kHz / dB | 50.80 | 47.27 | 41.47 | −3.53 |
| 5–20 kHz / dB | 77.85 | 77.22 | 57.73 | −0.63 |

`FwhImpellersFixedDelay2`：

| 指标 | 旧库 | 新库 | 叶轮单独 | Δ |
|---|---:|---:|---:|---:|
| OASPL / dB | 77.55 | 76.92 | 57.71 | −0.63 |
| L_A / dBA | 73.68 | 74.01 | 56.78 | +0.33 |
| 时域 AC RMS / Pa | 0.1850 | 0.1741 | 0.0418 | −0.53 dB |
| 20–100 Hz / dB | 57.00 | 27.81 | 37.21 | −29.18 |
| 100 Hz–1 kHz / dB | 51.87 | 34.49 | 36.92 | −17.38 |
| 1–5 kHz / dB | 51.44 | 47.66 | 43.99 | −3.78 |
| 5–20 kHz / dB | 77.49 | 76.92 | 57.44 | −0.58 |

三个对象完全同构：**低频 −18…−29 dB，高频 −0.5…−0.6 dB，总级 −0.6…−0.7 dB**。`FixedDelay` 的绝对值和 `FwhImpellers` 略有差别，是因为 `tauMin`/延迟改变了观测时间基（§5.7），但“低频降、高频不降”的结论不受影响。

#### 9.5.2 到 Nyquist 的全频带（`FwhImpellers`）

| 频带 | 旧 OASPL | 新 OASPL | 叶轮单独 | 旧 L_A | 新 L_A | 叶轮单独 L_A |
|---|---:|---:|---:|---:|---:|---:|
| 20–20 kHz | 77.46 | 76.88 | 57.73 | 73.62 | **74.00** | 56.83 |
| 20–50 kHz | 79.48 | 79.30 | 65.93 | 73.83 | **74.15** | 57.19 |
| 仅 20–50 kHz | 75.19 | 75.60 | 65.22 | 60.66 | 59.42 | 46.17 |

- 未计权 OASPL 算到 50 kHz 会再高约 2 dB；20–50 kHz 单独也有 75 dB 量级，叶轮单独的 20–50 kHz 也有 65.2 dB。这部分很可能是**未解析的时间/网格高频**，不能当作可听声。
- A 计权在 20 kHz 以上衰减很快，所以 `L_A` 对 20–50 kHz 不敏感（20–20k 与 20–50k 只差约 0.2 dB）；A 计权对比的结论是稳的。
- 新库 `L_A` 反而略高（+0.38 dBA），原因见 9.5.4：高频能量在 5–10 kHz 增加、在 16–20 kHz 减少。

#### 9.5.3 1/3 倍频程关键点（未计权，`FwhImpellers`）

| fc / Hz | 旧库 / dB | 新库 / dB | 叶轮单独 / dB |
|---:|---:|---:|---:|
| 40 | 56.64 | **26.85** | 35.29 |
| 80 | 38.69 | **24.43** | 28.84 |
| 125 | 46.98 | **22.54** | 27.08 |
| 200 | 32.93 | **14.18** | 25.47 |
| 400 | 43.82 | **14.92** | 28.17 |
| 800 | 44.55 | **30.12** | 31.68 |
| 1000 | 45.49 | **28.78** | 34.37 |
| 2000 | 46.73 | 43.31 | 43.19 |
| 5000 | 46.61 | **49.94** | 30.63 |
| 8000 | 71.15 | **72.82** | 57.01 |
| 12500 | 71.32 | **69.45** | 41.98 |
| 20000 | 70.08 | **68.24** | 50.32 |

- 40–1000 Hz：新库比旧库低 14–32 dB，部分带甚至低于叶轮单独基线（说明 wall 与叶轮在低频发生相消，或基线算例的低频本身不同；见 9.1 的基线口径说明）。
- 2 kHz 以上：新旧几乎重合；5–10 kHz 新库反而高 1–7 dB，12.5–20 kHz 低 1–2 dB。这是 `L_A` 略升的直接原因。

#### 9.5.4 A 计权 1/3 倍频程（dBA，`FwhImpellers`）

| fc / Hz | 旧库 | 新库 | 叶轮单独 |
|---:|---:|---:|---:|
| 200 | 22.1 | **3.3** | 14.6 |
| 500 | 41.4 | **25.7** | 21.7 |
| 1000 | 45.5 | **28.8** | 34.4 |
| 2000 | 47.9 | 44.5 | 44.4 |
| 4000 | 42.6 | 40.1 | 22.2 |
| 5000 | 47.2 | **50.5** | 31.2 |
| 8000 | 70.0 | **71.7** | 55.9 |
| 10000 | 56.9 | **59.5** | 37.6 |
| 16000 | 67.1 | **65.2** | 36.0 |
| 20000 | 60.7 | **58.9** | 41.0 |

`L_A` 是“能量加权 + A 计权”的结果：5–10 kHz 的上升（+1.7…+3.3 dBA）超过了 16–20 kHz 的下降（−1.8…−1.9 dBA），所以总 `L_A` 从 73.62 升到 74.00 dBA。这不是“高频变差”，而是**频谱重心向 5–10 kHz 移动**；未计权总能量仍然略降（−0.58 dB）。

#### 9.5.5 窄带第一根谱线与低频

旧库窄带谱在 43.65 Hz（`df` 的第一根非零谱线）为 **58.4 dB**，新库 **28.6 dB**，叶轮单独 37.1 dB；130.9 Hz 为 50.3 / 14.4 / 26.9 dB。43.65 Hz 正是 22.9 ms 窗口的频率分辨率，旧库的慢漂移在这里泄漏成“低频音”；新库均值稳定后，这根谱线随之消失。这个细节与 9.4.1 的分段均值表互相印证。

### 9.6 为什么低频有效、高频无效：机理汇总

1. **低频：空间相干 → `Σ n dS` 决定相消。** 闭合壁面上的均匀/低频压力近似同相，积分 `∫ p n dS ≈ p Σ n dS`。旧法向 `Σ n dS = 0.0224 m²`，把本应相消的稳态载荷变成一个净源；新法向 `Σ n dS = 0`，均匀分量严格相消，慢漂移和 20–200 Hz 虚假声随之消失（−20.2 dB）。
2. **高频：空间不相干 → 方差与符号无关。** 5–20 kHz 的壁面压力相关长度远小于 patch 尺度，各面贡献的交叉项平均为零，`E|Σ p_i n_i dS_i|² ≈ Σ E[p_i²] dS_i²`。把 85.5% 的面反号只改变符号图案，不改变方差，所以 SPL 只变 +0.05 dB。
3. **20–50 kHz 是另一块大头。** wall 在 20–50 kHz 的 RMS 有 0.11 Pa（约 75 dB），与 5–20 kHz 的 0.145 Pa 同量级；这部分很可能是 `dt = 1e-5 s`、max Courant ~4.5–11 和 cover 粗网格下的未解析/数值高频。法向修复对它的预期本来就是“无效”。
4. **近场项可能放大 wall。** 麦克风在 `(0.0062, 0.1, 0.2)`，位于空气网格包围盒内、离壁面较近；Farassat 1A 含 `1/r²` 近场项，wall 的高频不一定全是传播声。需要用远场 observer 验证 `1/r` 衰减。
5. **`tauMin` 只影响延迟对象。** `FwhImpellers`（`fixedResponseDelay false`）没有延迟，仍然出现同样的高频；所以高频残留不是 `tauMin` 全局耦合造成的。`tauMin` 只解释 `FixedDelay` 与 `FwhImpellers` 的绝对差异。

> 一句话：**法向修复解决的是“闭合面稳态载荷不闭合”这个数学错误；5–20 kHz 的 wall 项是壁面压力载荷本身的物理/数值内容，不在这个修复的作用域内。**

### 9.7 与 scFLOW 差异的再评估（融合 §5.6）

按“现有证据强度”排序，而不是按猜测排序：

| 假设 | 现有证据 | 判定 | 下一步 |
|---|---|---|---|
| **表面集不同**：scFLOW 可能不含 cover/内部壁面，或用封闭可渗透面 | cover 占 wall 面积 73%，新库 5–20 kHz 几乎不变；wall 能量集中在 cover 量级 | **最强** | 拆分 wall，逐 patch 输出 OAP；对齐 scFLOW 表面集 |
| **cover 网格/CFL**：cover 平均面尺度 ~1.6 mm，case ~0.27 mm，叶轮 ~0.17 mm；max Courant ~4.5–11 | 20–50 kHz 有 0.11 Pa；高频对法向不敏感，符合未解析高频特征 | **强** | `deltaT ≈ 2e-6`（maxCo≈1）或加密 cover，短段复算 |
| **近场 `1/r²`**：mic 距壁面近，wall 高频可能不是传播声 | wall 5–20k 0.145 Pa 远大于叶轮 0.011 Pa；无法从单点判断 | **中** | 加 0.5/1/2 m 远场 observer，检查是否按 `1/r` 衰减 |
| **固体壁面 vs 可渗透面**：自由场 1/r 积分不含遮挡/反射，内部壁面会被高估 | `laptop_wall` 含 CPU/Cu block 等内部面；scFLOW 口径未知 | **中** | 与 scFLOW 对齐表面集后再比绝对增量 |
| **法向符号**：旧法向 `Σ n dS ≠ 0` | 低频 −20 dB、均值 −60%，几何闭合性已验证；高频不变 | **低频已证实；高频已排除** | 无需再改法向 |
| **`tauMin` 全局耦合** | `FwhImpellers` 无延迟仍有高频；仅 `FixedDelay` 绝对值不同 | **高频已排除；时基问题仍在** | 按面/按 patch 取 `tauMin` |

### 9.8 验收标准与下一步（按优先级）

| 优先级 | 动作 | 预期/判据 | 验收标准 |
|---|---|---|---|
| P0 | 把 `laptop_wall` 拆成 `Cover`、`case1/2`、`fins`、`CPU/Cu_block` 独立 surface | 找出 5–20 kHz 的主要贡献块 | 单块贡献 >80% 即可锁定 |
| P0 | 加 0.5 / 1 / 2 m 远场 observer | 传播声应按 `1/r` 衰减；近场衰减更快 | 若 1 m 处 wall 高频降 >6 dB，说明当前以近场为主 |
| P1 | `deltaT ≈ 2e-6`（或 `adjustTimeStep` + `maxCo 0.6`）短段复算 | 若 5–20 kHz 降 >3 dB，说明是数值高频 | 与基准段比较同带宽 RMS |
| P1 | 壁面压力按面减时均 / 直接采 `p'` | 消除非均匀稳态载荷的直流 | wall 均值 <0.01 Pa；20–200 Hz <0.001 Pa |
| P2 | `tauMin` 按面/按 patch | `FixedDelay` 的叶轮 OAP 应与 `FwhImpellers` 一致 | 叶轮 OAP 差 <5% |
| P2 | 与 scFLOW 对齐表面集、observer 距离、网格 | 口径一致后再比绝对增量 | 同一表面集、同一 observer、同量级网格 |

**“问题解决”的量化标准（建议）：**

- wall 5–20 kHz RMS ≤ 0.02 Pa（当前 0.145 Pa，需降约 17 dB）；
- wall AC RMS ≤ 0.05 Pa（当前 0.184 Pa）；
- 总 OASPL 20–20 kHz 落在叶轮单独基线 +3 dB 以内，即约 **58–61 dB**（当前 76.9 dB）；
- `L_A` ≤ 60 dBA（当前 74.0 dBA）；
- wall 均值 ≤ 0.01 Pa（当前 −0.0765 Pa）。

在这些标准达到之前，只能说“法向 bug 已修复”，不能说“wall 声压异常已解决”。

### 9.9 图、数据与复现

**本次新增（2026-11-17，三对象、时域+1/3 倍频程+窄带）：**

- [诊断图：wall OAP 时程 + 频带柱状图](/home/ubuntu/shared/spl_compare/m1_vs_all_diagnosis.png)
- [FwhImpellers 三面板图](/home/ubuntu/shared/spl_compare/m1_vs_all_FwhImpellers.png)
- [FwhImpellersFixedDelay 三面板图](/home/ubuntu/shared/spl_compare/m1_vs_all_FwhImpellersFixedDelay.png)
- [FwhImpellersFixedDelay2 三面板图](/home/ubuntu/shared/spl_compare/m1_vs_all_FwhImpellersFixedDelay2.png)
- [完整报告](/home/ubuntu/shared/spl_compare/m1_vs_all_report.md)
- 数值：`m1_vs_all_<对象>_1_3oct.csv`、`m1_vs_all_<对象>_narrowband.csv`
- 脚本：[m1_vs_all_analysis.py](/home/ubuntu/shared/spl_compare/m1_vs_all_analysis.py)、[m1_vs_all_diagnosis.py](/home/ubuntu/shared/spl_compare/m1_vs_all_diagnosis.py)

**上一轮补充（A 计权、20–50 kHz）：**

- [A 计权对比图（FwhImpellers）](/home/ubuntu/shared/spl_compare/spl_compare_FwhImpellers_A.png)
- [A 计权对比图（FwhImpellersFixedDelay）](/home/ubuntu/shared/spl_compare/spl_compare_FwhImpellersFixedDelay_A.png)
- [A 计权与全频带汇总](/home/ubuntu/shared/spl_compare/spl_compare_summary.md)
- [未计权对比图（FwhImpellers）](/home/ubuntu/shared/spl_compare/spl_compare_FwhImpellers.png)
- [上一轮绘图脚本](/home/ubuntu/shared/spl_compare/make_spl_compare.py)

**复现命令：**

```bash
python3 /home/ubuntu/shared/spl_compare/m1_vs_all_analysis.py
python3 /home/ubuntu/shared/spl_compare/m1_vs_all_diagnosis.py
```

脚本只依赖 Python 标准库 + PIL，不需要 numpy；使用精确 DFT（Bluestein），因此 `df` 与记录长度一致（43.65 Hz），不是补零后的近似值。

> 新库 `_m1` 仍在运行。等它跑到更晚的时间后，重新运行上述脚本即可刷新到更长的窗口；结论所依赖的 0.03–0.0529 s 公共窗口不受影响。若要把窗口延长到 0.06 s，需要同时确认旧库数据已无法延长，否则只能用新库自己的时间演化做趋势检查，不能做新旧同窗对比。

### 9.10 数据量与频率分辨率可用性

这一节回答“现有数据够不够做 20 Hz–20 kHz / 5–20 kHz 的 SPL”，只讨论记录长度和频率分辨率，不改变 §9.5–9.7 的物理/数值结论。

#### 9.10.1 当前记录

| 项目 | 值 |
|---|---|
| 新库 `_m1` 全记录 | `t = 1e-5 … ≥0.06247 s`，约 6250 点 |
| 新旧公共窗口 | `0.03–0.0529 s`，2291 点 |
| 采样 | `fs = 100 kHz`，Nyquist **50 kHz** |
| 频率分辨率 | 全记录 `df ≈ 16 Hz`；公共窗口 `df ≈ 43.65 Hz` |
| 20 kHz 采样点数 | 5 点/周期；SPL 够用，波形峰值不一定够 |
| 20 Hz 周期 | 0.05 s；全记录约 **1.25 个周期**，公共窗口约 **0.46 个周期** |

#### 9.10.2 分目标判定

| 分析目标 | 当前数据 | 判定 |
|---|---|---|
| **20–20 kHz 总 OASPL / L_A** | 5–20 kHz 主导，公共窗口 bin 数足够 | **够用**（低频带缺失，但对总级影响小） |
| **5–20 kHz 1/3 倍频程** | 每带 26–106 个 bin（公共）/72–289 个（全记录） | **够用** |
| **5–20 kHz 窄带/音调** | `df ≈ 43.65 Hz`（公共）/16 Hz（全记录） | **仅能分开 >44 Hz / >16 Hz 的峰**；高分辨率音调不够 |
| **≥200–315 Hz 1/3 倍频程** | 公共窗 1–2 bin，全记录 3–5 bin | **边际/勉强够**；公共窗误差约 1–2 dB |
| **20–100 Hz 1/3 倍频程** | 公共窗 <1 bin，全记录 0.3–1.4 bin | **不够** |
| **20 Hz 窄带 SPL** | 全记录 1.25 个周期，`df ≈ 16 Hz` | **不够**；一个 bin 不能代表 20 Hz |
| **20 kHz 上限** | `fs = 100 kHz`，Nyquist 50 kHz | **采样够**；物理精度仍受 CFL/网格限制 |

#### 9.10.3 5–20 kHz 的 bin 数与统计不确定度

1/3 倍频程带宽为 `fc × 0.231`；对宽带稳态信号，带级标准不确定度粗略为 `1σ ≈ 4.34/√(B·T)` dB：

| fc | 带宽 | 公共窗 bin 数 | 全记录 bin 数 | 公共窗 1σ | 全记录 1σ |
|---:|---:|---:|---:|---:|---:|
| 5 kHz | 1.16 kHz | ~26 | ~72 | ~0.84 dB | ~0.51 dB |
| 8 kHz | 1.85 kHz | ~42 | ~115 | ~0.67 dB | ~0.40 dB |
| 10 kHz | 2.31 kHz | ~53 | ~144 | ~0.60 dB | ~0.36 dB |
| 16 kHz | 3.70 kHz | ~85 | ~230 | ~0.47 dB | ~0.29 dB |
| 20 kHz | 4.62 kHz | ~106 | ~289 | ~0.42 dB | ~0.26 dB |

- 5–20 kHz 的 1/3 倍频程在当前公共窗口已有 **26–106 个 bin**，全记录 **72–289 个 bin**；数据量不是 5–20 kHz 的限制。
- 实际观测到的旧/新 wall 5–20 kHz 差只有 **+0.05 dB**，`L_A` 差 +0.38 dB，都在上表 1σ 量级内；能稳妥下的结论是“**没有明显改善**”，不能把 +0.38 dB 当作确定的高频恶化。
- 旧库数据只到 0.0529 s，所以新旧对比只能用 0.0229 s 公共窗；新库自己的 5–20 kHz 谱可以用 0.0625 s 全记录，结论不会变。

#### 9.10.4 20–100 Hz 为什么不够

1/3 倍频程 20 Hz 的带宽约 **4.6 Hz**（17.8–22.4 Hz）。要每带有 3–5 个独立 bin：

```text
df ≤ 1–1.5 Hz  →  T ≥ 0.7–1 s
```

同时 20 Hz 周期 0.05 s，至少需要约 10 个周期（0.5 s）才比较稳。当前：

- 公共窗 `df ≈ 43.65 Hz`：20 Hz 带只有 **0.1 个 bin**，25/31.5/50/63 Hz 也基本为空；
- 全记录 `df ≈ 16 Hz`：20 Hz 带 **0.29 个 bin**，25/31.5/50/63 Hz 为 0.36/0.45/0.72/0.91 个 bin；
- 20 Hz 全记录只有 **1.25 个周期**，公共窗 **0.46 个周期**。

所以图里 20–63 Hz 曲线断开是记录长度限制，不是零声压；这些带不能作为定量结论。

- 100 Hz 带：带宽约 23 Hz，`T ≥ 0.2 s` 可得到约 4–5 个 bin；
- 200 Hz 带：`T ≥ 0.1 s` 左右；
- 20 Hz 带：`T ≥ 0.5 s`（最好 1 s）。

#### 9.10.5 建议

1. **5–20 kHz 的 1/3 倍频程/宽带 SPL**：当前数据够用；为了降低方差，新库单独分析可用 0.0625 s 全记录，新旧对比仍用 0.03–0.0529 s 公共窗。
2. **20–20 kHz 的 OASPL/L_A**：当前数据够用；低频带缺失对总级影响小，但应在图上标注 20–100 Hz 不可用。
3. **要完整 20 Hz–20 kHz 的 1/3 倍频程曲线**：20–100 Hz 需要把记录延长到 **0.5–1 s**；按当前约 13.7 s/步，从 0.0625 s 到 0.5 s 约 **7 天**，到 1 s 约 **15 天**。
4. **不要用零填充、Welch 平均或降采样“补”20 Hz**：它们不能增加物理记录长度；20 Hz 需要更长的物理时间。
5. **5–20 kHz 的窄带峰**：如果只关心宽带 SPL，当前 `df` 够；如果要分辨窄带音调，至少用 0.0625 s 全记录（`df ≈ 16 Hz`），并确认峰间距 >16 Hz。

> 总结：**5–20 kHz 的数据量够；20–20 kHz 总级够；20–100 Hz 的 1/3 倍频程/20 Hz 窄带不够。** 当前 5–20 kHz 的结论不是被数据量限制，而是壁面高频本身（cover 网格/CFL、近场、表面集）需要按 §9.8 继续排查。
