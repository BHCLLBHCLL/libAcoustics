# libAcoustics Farassat 1A 分析与修改记录

> 分支：`v2412`（`origin/v2412`）
> 基准：最初 `v2412` 为 `04f8a7a`；本文写作时的分支头为 `72eea13`，本次在其上增加 patch 表面法向修复。
> 源码：`/home/ubuntu/programs/libAcoustics`
> 日期：2026-11-16

本文把三类内容分开：

1. **性能分析**：最初 `v2412` 在 Windows 8 rank 算例上的平方分配、并行和内存问题；
2. **性能与正确性修改**：当前分支已经合入的导数分配、采样、Farassat 观测网格、面速度、多面求和等；
3. **本次新增的 patch 法向修复**：`laptop_wall` 的 7 个 patch 是闭合壁面，原质心启发式只翻转了 14.5% 的面，导致壁面压力载荷的符号和相消关系错误。本文给出源码路径、几何验证、运行数据、修复方法和与 Cradle scFLOW 的差异。

> 注意：性能节的绝对时间来自最初的 Windows 8 rank 运行；法向节的 OAP 数据来自 Linux 48 rank 的 `openfoam_fwh_m4_all` 运行。两套数据的网格、进程数和时间步不同，**不要直接混用绝对时间**。

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
- [8. 结论与后续](#8-结论与后续)

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
- 数据快照：旧库运行到 `t ≈ 0.0164 s`；本文的法向表取自 0.014–0.0164 s 的晚期窗口

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

从数据上看，wall 的贡献远大于理论预期：

| 量 | wall | impeller1 | impeller2 | 叶轮之和 |
|---|---:|---:|---:|---:|
| 面积 / m² | 0.19942 | 0.003397 | 0.003397 | 0.006794 |
| OAP RMS（0.014–0.0164 s）/ Pa | **0.197** | 0.022 | 0.019 | 0.059 |
| 5–20 kHz RMS / Pa | **0.193** | 0.010 | 0.011 | 0.014 |
| OAP 均值 / Pa | -0.051 | — | — | +0.044 |

- 总 `p'` RMS：加 wall **0.194 Pa**，叶轮单独 **0.0595 Pa**，约 **+10.3 dB**。
- 5–20 kHz：约 **+23 dB**（相对叶轮之和；相对单个叶轮约 +25 dB）。
- wall 与叶轮 OAP 的相关系数 **-0.007**，说明 wall 不是叶轮力的反作用，而是一个独立的宽频信号。
- 在相同时间窗内，两个算例的叶轮 OAP 完全一致（0.0590 vs 0.0595，均值都是 0.0442），所以增量全部来自 wall。

此外，wall 的均值随 `responseDelay` 变化很大：

| 对象 | wall 均值 / Pa |
|---|---:|
| `FwhImpellers`（无延迟） | -0.051 |
| `FwhImpellersFixedDelay`（0.01 s） | **-0.576** |
| `FwhImpellersFixedDelay2`（0.002 s） | -0.066 |

同一个壁面、同一个流场，只改延迟参数，直流量就差 10 倍；这说明绝对压力/稳态载荷被积分了，不是物理声压。

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

1. **法向处理**：scFLOW 用物面外法向（从固体指向流体），壁面载荷与叶片载荷反相抵消；当前库的 85.5% 壁面反了，变成相加。
2. **表面集**：scFLOW 的 FW-H 面是否包含 cover、内部 CPU/Cu block、fins？`laptop_wall` 的 cover 占 73% 面积；如果 scFLOW 只取风扇/壳体局部，增量自然小。
3. **封闭/可渗透面**：scFLOW 可能用包围风扇的封闭可渗透面，而不是所有固体壁面；自由场 1/r 积分对内部壁面也没有遮挡/反射，会高估。
4. **网格/湍流**：cover 的平均面尺度约 1.6 mm，case 约 0.27 mm，叶轮约 0.17 mm；cover 的网格最粗，5–20 kHz 的壁面压力脉动可能包含未解析的湍流/数值分量。当前算例的 max Courant 约 4.5–10（`adjustTimeStep no`，`deltaT = 1e-5`），也会放大高频噪声。

因此，先修法向，再和 scFLOW 对齐表面集、网格和 CFL，才有可比性。

### 5.7 仍建议的后续

法向修复是必要条件，但不一定能把 wall 降到“理论上的小量”，还需要：

- **壁面压力去均值/用脉动量**：`surfacePressure` 只乘 `rho`、减 `pInf = 0`（`FfowcsWilliamsHawkings.C`），没有减时均。对闭合面且法向正确时，均匀分量会相消；非均匀稳态载荷仍可能产生直流。建议按面减时均，或直接采 `p'`。
- **`tauMin` 全局耦合**：`fwhFormulation::update()` 的 `rMin_/tauMin_` 对所有 `controlSurfaces_` 取最小值，`tobs = ct + R/c0 - (tauMin - responseDelay)`。加入 wall 后，wall 上离麦克风更近的面把 `tauMin` 拉小，所有叶轮面的 `tobs` 都变。这解释了 `FwhImpellersFixedDelay` 的叶轮 OAP 在早期从 1.06 变成 0.40——是时基变了，不是物理。建议按面/按 patch 取 `tauMin`，或先用 `FwhImpellers`（`fixedResponseDelay false`）做对比。
- **cover 网格/CFL**：法向修完后如果 5–20 kHz 壁面仍大，把 `deltaT` 降到约 `2e-6`（或开 `adjustTimeStep`、`maxCo 0.6`）复算一小段，确认高频分量是不是数值。
- **表面集对齐**：确认 scFLOW 是否包含 cover/内部壁面、是否用封闭可渗透面；否则“增量大小”不是同一口径。
- **当前 `_all` 运行**：它加载的是旧库，wall 项目前是错的；修复后需要停掉、重编 `libAcoustics.so`、从 0（或保存的时间）重跑。

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
| `libAcoustics_analysis.md` | 本文：补充法向分析、数据、几何验证和修改清单，并按章节重排 |

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

编译通过，产物为 `libAcoustics_fix.so`；提交前已把 `Make/files` 改回 `libAcoustics`。正式使用前需要在算例停止后重新 `wmake`，生成/安装 `libAcoustics.so`，或把算例的 `libs` 指向 `libAcoustics_fix.so` 做短段验证。

另用 `Tutorials/vortexShed` 的 `cylinder` patch 做了 2 步单进程冒烟测试（`libs ("libAcoustics_fix.so")`）：FWH 对象正常加载、采样并输出 `OAP`，无崩溃。该测试只验证 `sampledPatch` 分支的运行路径；圆柱是凸体，新旧法向本来就一致，不能用来验证 `laptop_wall` 的符号。`laptop_wall` 的符号由 7.2 的几何检查和 7.3 的运行数据说明。

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

修复后的预期：

- wall 的 `Σ n dS` 从 0.0224 m² 变为 0；
- wall OAP 明显下降，总 `p'` 接近叶轮单独的 0.06 Pa；
- 5–20 kHz 壁面分量应显著下降；如果仍大，再查 cover 网格/CFL。

### 7.4 复现步骤

1. 停掉当前 `_all` 运行（它加载旧库）。
2. `wmake` 生成新的 `libAcoustics.so`。
3. 从 0 或保存的时间重跑 `pimpleFoam -parallel`，保持 `system/fwhControl` 不变。
4. 用 `parse OAP` 脚本比较 wall/impeller 的 OAP；先看 `FwhImpellers`（无延迟），再看 `FixedDelay`。
5. 如果 wall 仍大：把 `laptop_wall` 拆成 `cover`、`case`、`fins` 等独立 surface，定位是哪一块；同时检查 cover 的面尺度和 max Courant。

---

## 8. 结论与后续

1. **性能问题**的根因是 `tdPrimitive::dot` 的平方分配；当前分支已用一次分配、单次采样、观测时间网格解决，声学开销降到流动以下，内存从按面历史降到每麦克风一条信号。
2. **`laptop_wall` 声压异常增大**的根因是 patch 法向的质心启发式：7 个 wall patch 合计是闭合面，但启发式只翻了 14.5% 的面，`Σ n dS = 0.0224 m²`，cover 的净面积矢量没有被翻。壁面压力载荷因此不闭合，wall OAP 0.197 Pa（叶轮 0.059 Pa），5–20 kHz 0.193 Pa（叶轮 0.014 Pa）。
3. **修复**：`sampledPatch` 表面直接用 `-Sf/|Sf|`；叶轮和 wall 都回到正确外法向；非 patch 表面保留质心启发式。编译已通过，几何检查表明修复后 `Σ n dS = 0`。
4. **仍待做**：壁面压力去均值、`tauMin` 按面/按 patch、cover 网格/CFL 检查、与 scFLOW 表面集对齐；当前 `_all` 运行需要停掉并用新库重跑。
5. **不要用旧 `_all` 的 wall 数据下结论**：它是在错误法向下积分的，wall 的绝对值不可信；叶轮单独的 0.059 Pa 仍是可信的基线。
