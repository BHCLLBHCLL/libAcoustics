# libAcoustics Farassat 1A 并行性与核心效率分析

对照算例：

- 流动：`D:\training\cradle\laptop\laptop_thermal_steady_scaled_v3_transient\openfoam`
- 流动加声学：`D:\training\cradle\laptop\laptop_thermal_steady_scaled_v3_transient\openfoam_fwh`

库源码（与正在运行的 `libAcoustics.dll` 同一套）：

- `D:\training\openfoam\libAcoustics-2412`
- WSL：`/home/sdcll221/work/tries/of_mingw/libAcoustics-2412`

求解器为 OpenFOAM v2412 `pimpleFoam -parallel`，8 进程。声学对象为 `FwhImpellers`（`fixedResponseDelay false`）和 `FwhImpellersFixedDelay`（`fixedResponseDelay true`），公式均为 `Farassat1AFormulation`。

## 结论

每步大约 98 分钟。流动求解约占 2–3 分钟，其余时间在 libAcoustics 的 Farassat 1A 探针里。

热循环是 `AcousticAnalogy::tdPrimitive::dot`。`adjustTimeStep yes` 使每个面都调用 `List::resize`。`List` 没有容量预留，长度每变 1 就重新分配并搬移整张表。一个面上的代价是面数的平方。8 个进程都参加求和归约，但平方工作量集中在壳体面最多的 rank 3 上，其余进程堵在 `reduce`。

按当前步长跑到 `endTime = 0.1` 大约要三十多年。去掉平方分配之后，按麦克风、按面保存的延迟时间历史会成为下一堵内存墙。对应的修改顺序、预期时间和算例侧取舍见文末「效率提升方案」。

## 计时

两边网格相同，都是 8 进程、`maxCo 0.6`，稳定后 `deltaT ≈ 5.13×10⁻⁷` s。

`openfoam`（无声学）流动一步大约 110–155 秒。

`openfoam_fwh` 的 `log.pimpleFoam`：

- 第一步流动结束时 `ExecutionTime = 200.7 s`，此时还没有声学输出。
- 其后每步 `ExecutionTime` 间隔为 5792–5956 秒，13 步平均 **5877 秒**。
- 日志末尾停在第二个对象的 `Starting acoustics probe`，墙钟约 76598 秒，约 13 步、21 小时。

`endTime = 0.1`，步数约 `0.1 / 5.13e-7 ≈ 1.95×10⁵`。按 5877 秒/步，总时间约 `1.15×10⁹` 秒，约 36 年。

v2412 在每个时间步开头调用声学。`Time::run()` → `functionObjectList::execute()` 对每个 function object 先 `execute()` 再 `write()`。`FfowcsWilliamsHawkings::execute()` 为空，采样在 `AcousticAnalogy::write()` → `correct()`。因此日志顺序是：上一步流动的 `ExecutionTime`，然后 `writing force and moment files` 与 `Starting acoustics probe`，然后才是下一步的 Courant 数和 `Time =`。一步的 5877 秒等于本步声学探针加上本步流动。流动与无声学算例同量级（约 150–200 秒），探针约 5700 秒。

两个函数对象各跑一遍。每遍 18 行 `OAP`，即 6 个麦克风 × 3 个源面。日志里 `p' = 0` 是因为物理时间仍在约 `10⁻⁵` s，而面到麦克风的传播时间 `R/c0` 为数毫秒。`valueAt` 在到达时间之前返回 0。这个零结果仍然走完全部分配。

## 源面规模

全局面数（`constant/polyMesh/boundary`）：

| 面 | 面数 |
|---|---:|
| `laptop_wall`（`air_to_case1/2`、`air_to_Cover`、`air_to_fin_1/2`、`air_to_CPU`、`air_to_Cu_block`） | 576666 |
| `PartSurface_impeller1` | 113933 |
| `PartSurface_impeller2` | 113933 |
| 合计 | 804532 |

`laptop_wall` 占源面的 72%。两个叶轮各是一张 sampledSurface，壳体上的 7 个 wall 合成一张面。AMI 与开口不在源面上。

## 核心瓶颈：每个面重新分配整张导数表

时间导数在 `AcousticAnalogy::tdPrimitive::dot`（`Sources/lib/analogies/AcousticAnalogy/AcousticAnalogy.H`）。`old_`、`oldOld_` 是 `List<T>`，没有容量。`controlDict` 中 `adjustTimeStep yes`，`Farassat1AFormulation::initialize` 对每个麦克风、每个面调用 `setTimeStepVariable(true)`。于是每个面都走进：

```cpp
if ((oldOld_.size() <= idx) || (varTimeStep_))
{
    dvdt = (new_[idx] - old_[idx]) / deltat_;
    oldOld_.resize(idx + 1);
    ...
}
```

面循环下标为 `0, 1, 2, …, N-1`。某步第一次调用时表长为 `N`，`resize(1)` 把它缩成 1，之后每个面再长 1。

OpenFOAM v2412 `List::doResize`（`src/OpenFOAM/containers/Lists/List/List.C`）在长度变化时始终：

```cpp
T* old = this->v_;
this->size_ = len;
this->v_ = new T[len];
std::move(old, (old + overlap), this->v_);
delete[] old;
```

`new T[len]` 会默认构造全部元素。一个面上的分配次数是 `1+2+…+N ≈ N²/2`。

`calculateAcousticPressure` 对每个面调用三次 `dot`：`Un_`（标量）、`L_`（矢量）、`M_`（矢量）。这三份量不依赖麦克风，存储却是 `List<List<tdPrimitive>>`，下标为 `[iObs][iSurf]`，6 个麦克风各算各存一遍。两个函数对象再乘 2。

第一步走的是 `old_.resize` 的同样路径（`old_` 从空表逐面长到 `N`）。之后每步因为 `varTimeStep_` 为真，走 `oldOld_.resize`。所以从第 2 步起步长就平坦，与历史长度无关。这和日志里 5792–5956 秒的稳定间隔一致。

数量级（rank 3 上 `laptop_wall` 的 189609 个面，一个麦克风，一个 24 字节矢量）：

`189609 × 189610 / 2 × 24 ≈ 400 GB` 的分配、搬移和释放。

两个矢量场约 800 GB，标量 `Un_` 再约 140 GB。6 个麦克风 × 两个函数对象之后，临界路径上是数 TB 的 `malloc` / 搬移 / `free`。Farassat 1A 每个面只有几十次浮点运算，被分配淹没。有效内存流量若按约 2 GB/s（大量小块分配，到不了带宽峰值），这一项就是大约一小时量级，与实测约 95 分钟的探针时间相符。

## 并行

声压积分使用各 rank 的局部面，再对两个标量求和：

```cpp
reduce(intDotQdS_.value(iObs), sumOp<scalar>());
reduce(intFdS_.value(iObs), sumOp<scalar>());
```

调用在 `FfowcsWilliamsHawkings::correct` 里，顺序是每个麦克风、每个源面各一次。面与面之间全体同步。归约本身只有一个 `scalar`，时间可忽略。等待的是最慢的 rank。

scotch 平衡的是单元。源面分布（各 `processor*/constant/polyMesh/boundary` 的 `nFaces`）：

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

没有按麦克风或按面的线程并行。麦克风循环在 rank 内部串行。

`nonUniformSurfaceMotion true` 使每步 `expire()` 再 `update()`，并对三个面调用 `PatchTools::gatherAndMerge`。合并结果只写入 `mergeList_`。声压使用的是局部 `surf.Sf()`、`surf.Cf()`，这次聚集不进入积分，但每步所有 rank 都参加通信。`update()` 里还有每个麦克风一次 `reduce(rMax)` / `reduce(rMin)`，以及每个面一次 `gSum(Cf)`，相对平方分配是次要开销。

`correct()` 把 `surfaceVelocity`、`surfaceDensity`、`surfacePressure` 放在麦克风循环内部。采样结果与麦克风无关，每个面仍采样 6 次，并且每次 `interpolation<Type>::New("cell", ...)`。`interpolate false`，采样只是从边界取值。这一项是数十秒量级，不是 5700 秒的来源。

## 平方项去掉之后的内存墙

每个麦克风、每个面、厚度项 `qds_` 和载荷项 `fds_` 各存一对 `DynamicList<scalar>`（时间与数值），类型为 `List<List<List<Pair<DynamicList<scalar>>>>>`，下标 `[iObs][iSurf][iFace]`。

`fixedResponseDelay false` 时，过期时间为 `ct - tauMax`，`tauMax = R_max / c0`。麦克风在壳体包络外约 0.5 m，最远距离大约 0.5–1.5 m，`c0 = 343`，窗口约 3–4 ms。除以 `5.13×10⁻⁷` s，大约 6000–8000 个样本。

`fixedResponseDelay true` 时，过期时间为 `ct - (tauMax - tauMin + responseDelay)`，`responseDelay = 0.01`。窗口至少再增加 `0.01 / 5.13e-7 ≈ 2×10⁴` 个样本。

`cleanFreq 100` 每 100 步调用 `clearExpiredData` → `findExpiredIndex` → `getNewPointData`。它只删窗口之外的样本，并把剩余历史整表拷贝出来。窗口长度不变。

rank 3、189609 个面、6 个麦克风、4 条双精度序列：

`189609 × 6 × 4 × 8 × K ≈ 3.64×10⁷ × K` 字节。

`K ≈ 7800` 时一个函数对象约 280 GB。两个对象再翻倍，而且这只是一个 rank。窗口填满之前进程会先耗尽内存。`valueAt` 对每个面从头二分，历史变长后是 `N log K`，相对现在的平方分配仍然很小。

在信号到达之前（约 `R_min/c0 ≈ 1.5` ms，约 3000 步），`valueAt` 一直返回 0。按当前速度，这 3000 步的零结果就要大约半年。

## 调用链

```
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

## 效率提升方案

每步大约 98 分钟里，约 95 分钟是 `tdPrimitive::dot` 对每个面做一次 `List::resize`。先去掉这层平方分配，声学就会降到流动求解以下。要跑满 `endTime = 0.1`，还得把按面保存的延迟历史改成每个麦克风一条信号。

### 1. 去掉平方分配

`adjustTimeStep yes` 让每个面都执行 `oldOld_.resize(idx+1)`。`List::doResize` 没有容量，长度从 `N` 缩到 1 再逐面长回去，一个面的代价是 `N²/2` 次分配和搬移。rank 3 的 `laptop_wall` 有 189609 个面，再乘 6 个麦克风和两个函数对象，就是现在这约 95 分钟。

`resize()` 里一次配齐三张表，`dot()` 里只改元素。变步长仍用一阶差商，固定步长仍用原来的三点公式，数值与现在这条慢路径一致。

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

`haveOld_`、`haveOldOld_` 初始为 false。面数变化时在 `resize()` 里清掉这两个标志。改动只在 `AcousticAnalogy.H` 的 `tdPrimitive`，重编 `libAcoustics.dll` 后两个函数对象一起生效。

这一步之后，导数循环从数 TB 的分配变成约 `10⁸` 次浮点运算，rank 3 上不到 1 秒。一步时间回到流动的约 150–200 秒，再加上仍在的线性声学开销。

重编之前可以把 `adjustTimeStep` 关掉，`deltaT` 取当前库朗数允许的值（约 `5×10⁻⁷`）。固定步长从第三步起走原来的三点公式，不再每步 `resize`。前两步仍是平方复杂度。导数从现在的一阶变为二阶。历史内存不变。

### 2. 去掉重复工作和无用通信

平方项消失后，这些会成为声学里还能看见的部分。

- `L_`、`M_`、`Un_` 与麦克风无关，却按 6 个麦克风各存各算。每个源面只留一份 `tdPrimitive`。
- `surfaceVelocity`、`surfaceDensity`、`surfacePressure` 现在放在麦克风循环里面，每个面采样 6 次。提到 `iObs` 循环外，每个面采样一次。
- `nonUniformSurfaceMotion` 每步对三个面做 `PatchTools::gatherAndMerge`，结果只进 `mergeList_`，积分用的是本地 `Sf()`、`Cf()`。积分路径上不再做这次聚集。

### 3. 历史改成每个麦克风一条信号

每个麦克风、每个面、厚度项和载荷项各有一条 `DynamicList`。保留窗口约 `R/c0`，以 `deltaT ≈ 5.13×10⁻⁷` 计是 6000–8000 个样本；`fixedResponseDelay` 还要加上 `responseDelay 0.01`。rank 3 上一个函数对象窗口填满时约 280 GB，两个对象再翻倍。`cleanFreq 100` 只删窗口外的样本，窗口长度不变。只修 `resize` 的话，跑几千步后会先耗尽内存。

延迟时刻的面积分对时间插值是线性运算，可以把每个面的贡献直接散到该麦克风的到达时刻上，不必保存每个面的历史。每个麦克风一条公共时间轴，长度约 `(endTime + R_max/c0) / dt_ac`。`endTime = 0.1`、声学步长 `10⁻⁵` s 时，6 个麦克风一共几 MB。各 rank 先在本地累加，再对这条短数组做一次求和。`valueAt` 的逐面二分和 `getNewPointData` 的整表拷贝随之消失。

到达时刻 `t + R/c0` 落在两个相邻声学时刻之间时，按线性权重分到这两个点上。这与现在「先对每个面插值再求和」是同一种线性插值，结果不会逐位相同。

### 4. 算例侧

库改完之后，流动仍是每步约 150–200 秒，因为 `maxCo 0.6` 把 `deltaT` 压在 `5×10⁻⁷`。声学采样不必跟着这个步长。叶轮转速约 100 Hz，`probeFrequency 20` 对应大约 100 kHz，仍高于叶片通过频率。采样步数和历史长度都按这个倍数下降。`maxDeltaT 2×10⁻⁴` 若真的被用上，这个间隔会穿过叶片通过频带，那时再把 `probeFrequency` 收回 1。

两个函数对象各做一遍完整积分。比较延迟公式时保留两个；只需要一条声压曲线时留一个，时间和内存减半。

`laptop_wall` 有 576666 个面，占源面 72%，也是 rank 3 的负载来源。只要风扇表面噪声时可以先不算这张面。它是固体壁面载荷，不是叶轮表面的替代。

### 预期

| 措施 | 一步时间 | 能否跑到 0.1 s |
|---|---|---|
| 现状 | 约 5900 s | 不能，约 36 年，且内存会先满 |
| 关闭 `adjustTimeStep` | 前两步仍约 5900 s，之后降到数分钟 | 不能，历史内存仍在 |
| 修好 `resize`，并去掉重复采样和 `gatherAndMerge` | 约 3–4 分钟，由流动决定 | 不能，数千步后内存满 |
| 再改成麦克风信号累加 | 仍约 3–4 分钟 | 能，声学内存为数 MB |
| 再把 `probeFrequency` 提到 20、只留一个函数对象 | 流动时间不变，声学按倍数下降 | 能 |
