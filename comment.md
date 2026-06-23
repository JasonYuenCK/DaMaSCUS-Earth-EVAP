# DaMaSCUS-SUN-EVAP 代码修改报告

修改范围：

* 捕获模式及“固定捕获数”逻辑保持不变。
* **彻底删除全部蒸发时间分峰分析代码。**
* **彻底删除全部性能诊断统计代码，不新增性能诊断开关。**
* 保留快照调度、进度显示和单轨迹 wall-time 安全截断所必需的计时。
* 保留蒸发物理诊断功能，但与 `evaporation_times.txt` 分离。

---

## 1. 删除蒸发时间分峰分析

### 1.1 `include/Parameter_Scan.hpp`

从 `Configuration` 删除：

```cpp
bool evaporation_mode_bincount_enabled;
std::vector<double> evaporation_mode_boundaries_log10_s;
std::vector<std::string> evaporation_mode_labels;
bool evaporation_mode_include_truncated;
```

这些字段目前位于配置类中。

---

### 1.2 `src/Parameter_Scan.cpp`

从 `Configuration::Import_Parameter_Scan_Parameter()` 中完整删除以下配置解析：

```ini
evaporation_mode_bincount_enabled
evaporation_mode_boundaries_log10_s
evaporation_mode_labels
evaporation_mode_include_truncated
```

包括：

* 默认边界 `{4.5, 11.1}`；
* 默认标签 `P1_fast`、`P2_theory`、`P3_tail`；
* 边界排序检查；
* 标签数量检查；
* `include_truncated` 兼容逻辑。

当前完整配置解析块位于 `Import_Parameter_Scan_Parameter()` 中。

同时从 `Configuration::Print_Summary()` 删除这些参数的终端打印。

---

### 1.3 `src/main.cpp`

删除：

```cpp
data_set.Configure_Evaporation_Mode_Bincount(
    cfg.evaporation_mode_bincount_enabled,
    cfg.evaporation_mode_boundaries_log10_s,
    cfg.evaporation_mode_labels,
    cfg.evaporation_mode_include_truncated
);
```

当前调用位于数据生成前。

---

### 1.4 `include/Data_Generation.hpp`

完整删除：

```cpp
struct EvaporationModeBincount
{
    unsigned long int count;
    unsigned long int truncated_count;
    std::array<double, NUM_BINS> dt_hist;
    std::array<double, NUM_BINS> v2dt_hist;
    std::array<double, NUM_BINS> dt_sq_hist;
    std::array<double, NUM_BINS> v2dt_sq_hist;
};
```

同时从 `Simulation_Data` 删除：

```cpp
bool evaporation_mode_bincount_enabled;
std::vector<double> evaporation_mode_boundaries_log10_s;
std::vector<std::string> evaporation_mode_labels;
bool evaporation_mode_include_truncated;
std::vector<EvaporationModeBincount> evaporation_mode_bincounts;
```

删除公开方法：

```cpp
void Configure_Evaporation_Mode_Bincount(
    bool enabled,
    const std::vector<double>& boundaries_log10_s,
    const std::vector<std::string>& labels,
    bool include_truncated
);
```

当前结构、成员和接口集中在该头文件。

---

### 1.5 `src/Data_Generation.cpp`

删除辅助函数：

```cpp
int Classify_Evaporation_Mode(...);
void Accumulate_Mode_Bincount(...);
```

当前函数负责按人为边界分类和累计分峰 bincount。

删除成员配置函数：

```cpp
Simulation_Data::Configure_Evaporation_Mode_Bincount(...)
```

删除轨迹完成后的分峰分类：

```cpp
if(evaporation_mode_bincount_enabled
   && rec.survival_valid
   && rec.event_observed
   && Has_Positive_Evaporation_Time(rec.lifetime_unbinding))
{
    int mode = Classify_Evaporation_Mode(...);
    Accumulate_Mode_Bincount(...);
}
```

删除 `Perform_MPI_Reductions()` 中：

```cpp
for(auto& mode_bincount : evaporation_mode_bincounts)
{
    MPI_Allreduce(...);
}
```

删除全部输出生成：

```text
evaporation_mode_summary.txt
evaporation_mode_bincount.txt
```

以及与以下内容有关的所有代码：

* mode label；
* mode boundary；
* mode count；
* mode error；
* mode radial bincount。

当前输出代码位于 `Write_Output_Files()`。

完成后，源代码中不应再出现：

```text
EvaporationMode
evaporation_mode
P1_fast
P2_theory
P3_tail
mode_bincount
```

---

## 2. 删除全部性能诊断统计

不增加：

```cpp
performance_diagnostics_enabled
```

而是直接删除整套统计和输出。

---

### 2.1 `include/Data_Generation.hpp`

从 `Simulation_Data` 删除：

```cpp
double total_wall_time_captured;
double total_wall_time_not_captured;

static constexpr int WALL_TIME_BINS;
static constexpr double WALL_TIME_LOG_MIN;
static constexpr double WALL_TIME_LOG_MAX;

std::array<unsigned long int, WALL_TIME_BINS>
    wall_time_hist_captured;
std::array<unsigned long int, WALL_TIME_BINS>
    wall_time_hist_not_captured;

unsigned long int wall_time_overflow_captured;
unsigned long int wall_time_overflow_not_captured;

unsigned long int total_rk45_steps_captured;
unsigned long int total_rk45_steps_not_captured;

static constexpr int STEP_COUNT_BINS;
static constexpr double STEP_COUNT_LOG_MIN;
static constexpr double STEP_COUNT_LOG_MAX;

std::array<unsigned long int, STEP_COUNT_BINS>
    step_count_hist_captured;
std::array<unsigned long int, STEP_COUNT_BINS>
    step_count_hist_not_captured;

unsigned long int step_count_overflow_captured;
unsigned long int step_count_overflow_not_captured;
```

当前整套成员位于头文件的计算时间统计区域。

如果 `termination_reason_counts` 不再用于任何物理输出，也一并删除：

```cpp
std::array<unsigned long int,
           TRAJECTORY_TERMINATION_REASON_COUNT>
    termination_reason_counts;
```

保留已有的汇总计数即可：

```cpp
number_of_complete_evaporation_particles
number_of_censored_captured_particles
number_of_invalid_survival_captured_particles
```

---

### 2.2 `src/Data_Generation.cpp` 顶部辅助代码

删除：

```cpp
constexpr int SNAPSHOT_WALL_TIME_BINS = 120;
constexpr int SNAPSHOT_STEP_COUNT_BINS = 140;
```

删除仅为性能直方图服务的：

```cpp
template<typename T, size_t N>
bool Find_Nonzero_Bin_Range(...);
```

如果删除 `termination_reason_counts`，同时删除：

```cpp
const char* TerminationReason_Name(...);
int TerminationReason_Index(...);
unsigned long int Count_Incomplete_Physical_Terminations(...);
```

但保留：

```cpp
TrajectoryTerminationReason
TrajectoryTerminationInvalidatesSurvival(...)
```

因为它们仍用于轨迹物理状态判断。

---

### 2.3 简化 `Rank_Snapshot_State`

删除：

```cpp
double current_trajectory_wall_sec;
double current_trajectory_physical_sec;

double total_wall_time_captured;
double total_wall_time_not_captured;

uint64_t total_rk45_steps_captured;
uint64_t total_rk45_steps_not_captured;

std::array<uint64_t, SNAPSHOT_WALL_TIME_BINS>
    wall_time_hist_captured;
std::array<uint64_t, SNAPSHOT_WALL_TIME_BINS>
    wall_time_hist_not_captured;

uint64_t wall_time_overflow_captured;
uint64_t wall_time_overflow_not_captured;

std::array<uint64_t, SNAPSHOT_STEP_COUNT_BINS>
    step_count_hist_captured;
std::array<uint64_t, SNAPSHOT_STEP_COUNT_BINS>
    step_count_hist_not_captured;

uint64_t step_count_overflow_captured;
uint64_t step_count_overflow_not_captured;
```

当前这些字段直接进入每个 rank 的 checkpoint。

保留：

```cpp
double rank_elapsed_wall_sec;
```

因为它仍用于：

* 判断最终 rank state 是否属于某个快照时间；
* 快照调度；
* 运行完成状态。

---

### 2.4 简化 `Snapshot_Report_State`

删除所有：

```cpp
total_wall_time_*
total_rk45_steps_*
wall_time_hist_*
step_count_hist_*
*_overflow_*
```

保留：

```cpp
total_trajectories
captured_particles
snapshot_bincount_captured_samples
snapshot_bincount_not_captured_samples

captured_dt_hist
captured_v2dt_hist
captured_dt_sq_hist
captured_v2dt_sq_hist

not_captured_dt_hist
not_captured_v2dt_hist
not_captured_dt_sq_hist
not_captured_v2dt_sq_hist
```

---

### 2.5 修改 checkpoint 二进制读写

从以下函数删除性能字段的序列化和反序列化：

```cpp
Write_Rank_Snapshot_State(...)
Read_Rank_Snapshot_State(...)
```

删除当前写入或读取的：

```cpp
state.total_wall_time_captured
state.total_wall_time_not_captured
state.total_rk45_steps_captured
state.total_rk45_steps_not_captured
state.wall_time_hist_*
state.wall_time_overflow_*
state.step_count_hist_*
state.step_count_overflow_*
```

当前读取逻辑包含所有这些字段。

修改 checkpoint 格式版本号，避免新代码误读旧二进制文件：

```cpp
constexpr uint32_t SNAPSHOT_FORMAT_VERSION = 2;
```

旧版本 checkpoint 检测到后直接忽略，不尝试兼容解析。

---

### 2.6 简化快照累计

从：

```cpp
Accumulate_Snapshot_Report_State(...)
```

删除：

```cpp
report.total_wall_time_captured += ...
report.total_wall_time_not_captured += ...
report.total_rk45_steps_captured += ...
report.total_rk45_steps_not_captured += ...

wall_time_hist 循环
step_count_hist 循环
overflow 累计
```

当前快照累计函数同时处理物理 bincount 和性能统计。

---

### 2.7 简化快照文本输出

从 `Write_Snapshot_Report_File()` 删除：

```text
total_wall_time_captured_sec
total_wall_time_not_captured_sec
mean_wall_time_captured_sec
mean_wall_time_not_captured_sec
total_rk45_steps_captured
total_rk45_steps_not_captured
mean_rk45_steps_captured
mean_rk45_steps_not_captured
ranks_ready
per-rank status
current trajectory wall time
current trajectory physical time
```

当前这些字段写入快照头部。

快照文件只保留：

```text
snapshot index
snapshot target wall time
total trajectories
captured particles
captured bincount sample count
not-captured bincount sample count
bincount histogram
```

删除：

```cpp
Snapshot_Load_Diagnostics
Write_Snapshot_Diagnostics(...)
Write_Snapshot_Text_Status(...)
Format_Rank_Status(...)
```

当快照尚未收齐所有 rank 时，`Try_Write_Merged_Snapshot()` 直接返回 `false`，不再生成等待状态或 rank 性能状态文件。

---

### 2.8 删除主循环中的每轨迹计时

删除：

```cpp
auto traj_t0 = std::chrono::high_resolution_clock::now();

Trajectory_Result trajectory =
    simulator.Simulate(IC, DM, mpi_rank);

auto traj_t1 = std::chrono::high_resolution_clock::now();

double traj_wall_sec =
    std::chrono::duration_cast<...>(traj_t1 - traj_t0).count();
```

删除：

```cpp
wall_time_bin(...)
step_count_bin(...)
```

删除 captured 和 not-captured 分支中的：

```cpp
total_wall_time_* += traj_wall_sec;
total_rk45_steps_* += trajectory.total_rk45_steps;
wall_time_hist_*[...];
step_count_hist_*[...];
*_overflow_*++;
```

这些计算当前对每条轨迹无条件执行。

保留：

```cpp
const double trajectory_completion_wall_time_sec =
    elapsed_since_start();
```

因为快照蒸发日志仍需判断事件属于哪个快照时间段。

---

### 2.9 删除 MPI 性能 reduction

从 `Perform_MPI_Reductions()` 删除：

```cpp
MPI_Allreduce(... total_wall_time_captured ...);
MPI_Allreduce(... total_wall_time_not_captured ...);

MPI_Allreduce(... wall_time_hist_captured ...);
MPI_Allreduce(... wall_time_hist_not_captured ...);

MPI_Allreduce(... total_rk45_steps_captured ...);
MPI_Allreduce(... total_rk45_steps_not_captured ...);

MPI_Allreduce(... step_count_hist_captured ...);
MPI_Allreduce(... step_count_hist_not_captured ...);

MPI_Allreduce(... overflow ...);
```

当前对应代码集中在函数尾部。

保留：

```cpp
MPI_Allreduce(MPI_IN_PLACE,
              &computing_time,
              1,
              MPI_DOUBLE,
              MPI_MAX,
              MPI_COMM_WORLD);
```

它仅用于总运行时间显示，不属于逐轨迹性能诊断。

---

### 2.10 删除性能输出文件

完整删除生成：

```text
computation_time_summary.txt
```

包括其中：

* wall-time 总量与均值；
* RK45 步数总量与均值；
* wall-time histogram；
* RK45-step histogram；
* overflow；
* termination reason 表。

当前文件生成代码位于 `Write_Output_Files()` 尾部。

---

## 3. 删除 `Trajectory_Result` 中纯诊断 RK45 计数

### 3.1 `include/Simulation_Trajectory.hpp`

从 `Trajectory_Result` 删除：

```cpp
unsigned long int total_rk45_steps;
```

修改构造函数：

```cpp
Trajectory_Result(
    const Event& event_ini,
    const Event& event_final,
    unsigned long int nScat,
    TrajectoryBincount bc
);
```

当前字段和构造函数参数位于结果结构中。

从 `Trajectory_Simulator` 删除：

```cpp
unsigned long int total_rk45_steps_current_traj;
```

当前成员只用于生成性能统计。

---

### 3.2 `src/Simulation_Trajectory.cpp`

修改构造函数：

```cpp
Trajectory_Result::Trajectory_Result(
    const Event& event_ini,
    const Event& event_final,
    unsigned long int nScat,
    TrajectoryBincount bc
)
: initial_event(event_ini),
  final_event(event_final),
  number_of_scatterings(nScat),
  bincount(std::move(bc))
{}
```

删除所有：

```cpp
total_rk45_steps_current_traj = 0;
total_rk45_steps_current_traj += time_steps;
```

修改返回：

```cpp
return Trajectory_Result(
    initial_condition,
    current_event,
    number_of_scatterings,
    current_bincount
);
```

当前返回仍携带累计 RK45 步数。

保留 `Propagate_Freely()` 中的局部：

```cpp
unsigned long int time_steps;
```

它仍用于：

```cpp
while(time_steps < maximum_time_steps)
```

以及 `MaxFreeSteps` 安全截断，不能删除。

---

## 4. 统一蒸发输出

### 4.1 `evaporation_times.txt`

无论是否开启蒸发物理诊断，始终使用固定简单格式：

```text
# rank trajectory_id lifetime_unbinding_sec
0  123  1.2345678900e+08
```

只写：

```cpp
rec.survival_valid
&& rec.event_observed
&& std::isfinite(rec.lifetime_unbinding)
&& rec.lifetime_unbinding >= 0.0
```

删除当前根据诊断开关改变 `evaporation_times.txt` 列结构的行为。

---

### 4.2 修复快照日志硬编码

当前快照合并调用：

```cpp
Write_Evaporation_Time_Log_Block(
    ...,
    snapshot_records,
    true
);
```

这里的 `true` 强制启用诊断格式。

应改为简单蒸发事件追加接口，例如：

```cpp
Append_Completed_Evaporation_Events(
    evaporation_log_path,
    snapshot_records,
    log_state
);
```

该函数不接收 diagnostics 参数。

---

### 4.3 保留独立的蒸发物理诊断文件

当：

```ini
evaporation_diagnostics_enabled = true;
```

额外写：

```text
evaporation_diagnostics.txt
```

将当前：

```text
evaporation_summary.txt
```

重命名为该文件。

诊断文件可继续包含：

* 捕获时间；
* 最终解缚散射时间；
* 边界逃逸时间；
* termination reason；
* 能量漂移；
* 散射次数；
* censoring 和 invalid 标记。

它不包含：

* 单轨迹 wall time；
* RK45 步数；
* wall-time histogram；
* 性能统计。

当前诊断文件生成位置位于 `Write_Output_Files()`。

---

### 4.4 最终普通模式输出

修改后普通运行只产生：

```text
bincount.txt
evaporation_times.txt
```

开启蒸发物理诊断后增加：

```text
evaporation_diagnostics.txt
```

不再产生：

```text
evaporation_mode_summary.txt
evaporation_mode_bincount.txt
computation_time_summary.txt
```

---

## 5. 修改蒸发日志追加机制

### 5.1 删除每次追加前的全文件扫描

当前日志追加流程会重复读取已有日志并重建状态。应改为在 rank 0 内存中维护：

```cpp
struct EvaporationLogState
{
    int last_snapshot_index = 0;
    bool final_block_written = false;
    uint64_t records_written = 0;
    std::unordered_set<uint64_t> written_record_keys;
};
```

记录 key 可由：

```cpp
uint64_t key =
    (static_cast<uint64_t>(rank) << 48)
    ^ trajectory_id;
```

构造。

---

### 5.2 修改接口

将：

```cpp
Write_Evaporation_Time_Log_Block(...)
```

拆分为：

```cpp
bool Initialize_Evaporation_Log(
    const std::string& path,
    double mass_gev,
    double sigma_cm2,
    EvaporationLogState& state
);

bool Append_Completed_Evaporation_Events(
    const std::string& path,
    const std::vector<EvaporationRecord>& records,
    EvaporationLogState& state
);
```

运行期间只处理新记录。

仅在断点恢复时调用一次：

```cpp
Recover_Evaporation_Log_State(...)
```

删除正常追加路径中的：

```cpp
Read_Evaporation_Log_State(...)
Evaporation_Log_Has_Snapshot_Block(...)
```

---

## 6. 修改快照 checkpoint 存储

将当前：

```text
rank_<rank>_snapshot_<index>.bin
```

改为每个 rank 一个文件：

```text
rank_<rank>.bin
```

文件内使用固定大小 slot：

```cpp
offset =
    sizeof(FileHeader)
    + snapshot_index * sizeof(RankSnapshotRecord);
```

使用：

```cpp
pwrite(...)
pread(...)
```

按 snapshot index 定位。

每个 slot 包含：

```cpp
uint32_t format_version;
uint32_t snapshot_index;
uint64_t run_id;
uint8_t ready;
Rank_Snapshot_State state;
uint64_t checksum;
```

写入顺序：

1. `ready = 0`；
2. 写主体和 checksum；
3. 最后写 `ready = 1`。

这样不再为每个 rank、每个快照创建和删除独立文件。

---

## 7. 修复自由传播中的散射光深积分

当前代码先完成 RK45 步进，再在步末计算散射率，并将 (0.1/\Gamma) 仅应用于下一步。

应修改 `Trajectory_Simulator::Propagate_Freely()`。

### 7.1 步进前计算散射率

```cpp
double rate_before = 0.0;

if(r_before < rSun)
{
    rate_before =
        solar_model.Total_DM_Scattering_Rate(
            DM,
            r_before,
            v_before
        );

    if(std::isfinite(rate_before) && rate_before > 0.0)
    {
        constexpr double MAX_OPTICAL_DEPTH_STEP = 0.05;

        particle_propagator.time_step =
            std::min(
                particle_propagator.time_step,
                MAX_OPTICAL_DEPTH_STEP / rate_before
            );
    }
}
```

---

### 7.2 使用梯形光深增量

步进后：

```cpp
double rate_after = 0.0;

if(r_after < rSun)
{
    rate_after =
        solar_model.Total_DM_Scattering_Rate(
            DM,
            r_after,
            v_after
        );
}

double delta_tau =
    0.5 * (rate_before + rate_after) * actual_dt;
```

替换当前：

```cpp
minus_log_xi -= actual_dt * total_rate;
```

---

### 7.3 定位步内散射时刻

当：

```cpp
delta_tau >= minus_log_xi
```

计算：

```cpp
double fraction =
    minus_log_xi / delta_tau;

fraction = std::clamp(fraction, 0.0, 1.0);
```

然后将散射事件定位在该 RK 步内部，而不是直接使用步末状态。

至少需要对以下量插值：

```cpp
time
radius
radial velocity
phi
```

更稳妥的实现是：

* 保存步前 propagator；
* 二分缩小 RK 步；
* 直到累计光深误差小于指定容差；
* 使用该位置作为 `Scatter()` 输入。

---

## 8. 修复热靶速度方向奇点

当前代码计算：

```cpp
double aux =
    sqrt(1.0 - pow(unit_vector_DM[2], 2.0));
```

并使用：

```cpp
sin_theta / aux
```

当 DM 速度接近 (z) 轴时会出现除零。

替换为稳定正交基：

```cpp
libphysica::Vector u = vel_DM.Normalized();

libphysica::Vector reference =
    std::fabs(u[2]) < 0.9
    ? libphysica::Vector({0.0, 0.0, 1.0})
    : libphysica::Vector({1.0, 0.0, 0.0});

libphysica::Vector e1 =
    reference.Cross(u).Normalized();

libphysica::Vector e2 =
    u.Cross(e1).Normalized();

libphysica::Vector unit_vector_T =
      cos_theta * u
    + sin_theta * cos_phi * e1
    + sin_theta * sin_phi * e2;
```

并在 `vel_DM.Norm()` 非正或非有限时直接报数值错误。

---

## 9. 增加散射目标采样保护

当前 `Sample_Target()` 直接执行：

```cpp
double sum = rate_electron / total_rate;
```

没有检查 `total_rate`。

增加：

```cpp
if(!std::isfinite(total_rate) || total_rate <= 0.0)
{
    throw std::runtime_error(
        "Sample_Target(): total scattering rate "
        "is non-positive or non-finite."
    );
}
```

同时逐项检查：

```cpp
rate_nuclei_cache[i] >= 0
rate_electron >= 0
std::isfinite(rate_nuclei_cache[i])
std::isfinite(rate_electron)
```

出现非法散射率时，将轨迹标记为：

```cpp
TrajectoryTerminationReason::NumericalFailure
```

而不是继续采样或直接 `std::exit()` 终止整个 MPI 作业。

---

## 10. 修改测试代码

删除所有涉及以下内容的测试：

```text
evaporation mode boundaries
evaporation mode labels
mode classification
mode bincount
mode output files
wall-time histogram
RK45-step histogram
computation_time_summary.txt
```

新增以下测试。

### 输出 contract

默认模式必须存在：

```text
bincount.txt
evaporation_times.txt
```

默认模式不得存在：

```text
evaporation_mode_summary.txt
evaporation_mode_bincount.txt
computation_time_summary.txt
evaporation_diagnostics.txt
```

开启蒸发物理诊断后只额外存在：

```text
evaporation_diagnostics.txt
```

### 热靶速度

测试：

```cpp
vel_DM = {0.0, 0.0, 1.0};
vel_DM = {0.0, 0.0, -1.0};
vel_DM = {1e-15, 0.0, 1.0};
```

要求采样速度所有分量均有限。

### 散射率

测试：

```text
total_rate = 0
total_rate = NaN
one isotope rate < 0
```

要求返回数值失败，不出现除零。

### 光深收敛

比较：

```text
MAX_OPTICAL_DEPTH_STEP = 0.10
MAX_OPTICAL_DEPTH_STEP = 0.05
MAX_OPTICAL_DEPTH_STEP = 0.02
```

检查：

```text
capture result
first scattering radius
total scattering count
evaporation lifetime
radial bincount
```

### 捕获模式回归

保持当前：

```text
sample_size = target captured count
```

并验证上述删除和传播修改不会改变捕获模式的停止条件。
