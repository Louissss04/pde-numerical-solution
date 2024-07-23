# pde-numerical-solution

## 介绍
本仓库用于存储<偏微分方程数值解>暑期学校算法代码。

## 文件夹结构

- `FDMscheme/`：包含实现不同有限差分方法的函数。
    - `FDM3points.m`：三点方案
    - `LeapFrog.m`：Leap-Frog 方案
    - `Boxscheme.m`：Box 方案
    - `BeamWarming.m`：Beam-Warming 方案
    - `Implict_central.m`：隐式-中心差分方案
    - `NonlinearHyperbolicSolver.m`：非线性双曲方程求解方案

- `Example/`：演示上述方法的脚本。
    - `FDMAdvection3points1D.m`：调用 `FDM3points` 函数的示例
    - `FDMAdvectionLeapFrog1D.m`：调用 `LeapFrog` 函数的示例
    - `FDMAdvectionBox1D.m`：调用 `Boxscheme` 函数的示例
    - `FDMAdvectionBeamWarming1D.m`：调用 `BeamWarming` 函数的示例
    - `FDMAdvectionImplict1D.m`：调用 `Implict_central` 函数的示例
    - `FDMNonLinHyper1D.m`：调用 `NonlinearHyperbolicSolver` 函数的示例

## 文件说明

### FDMscheme/FDM3points.m
该函数使用不同的三点有限差分方法求解偏微分方程 `u_t + a u_x = 0`。

**语法：**
```matlab
u = FDM3points(a, delta_t, delta_x, x_start, x_end, t_start, t_end, scheme, phi, g)
```

**输入参数：**
- `a`：对流速度
- `delta_t`：时间步长
- `delta_x`：空间步长
- `x_start`：空间起点
- `x_end`：空间终点
- `t_start`：时间起点
- `t_end`：时间终点
- `scheme`：使用的方案（'LF', 'LW', 'CD', 'BW', 'FW'）
- `phi`：边界条件函数句柄，`u(t, 0)`
- `g`：初始条件函数句柄，`u(0, x)`

**输出参数：**
- `u`：解矩阵，其中每一行表示一个时间步长的状态向量。

### FDMscheme/LeapFrog.m
该函数使用 Leap-Frog 方法求解偏微分方程 `u_t + a u_x = 0`。

**语法：**
```matlab
u = LeapFrog(a, delta_t, delta_x, x_start, x_end, t_start, t_end, phi, g)
```

**输入参数：**
- `a`：对流速度
- `delta_t`：时间步长
- `delta_x`：空间步长
- `x_start`：空间起点
- `x_end`：空间终点
- `t_start`：时间起点
- `t_end`：时间终点
- `phi`：边界条件函数句柄，`u(t, 0)`
- `g`：初始条件函数句柄，`u(0, x)`

**输出参数：**
- `u`：解矩阵，其中每一行表示一个时间步长的状态向量。

### FDMscheme/Boxscheme.m
该函数使用 Box 方案求解偏微分方程 `u_t + a u_x = 0`。

**语法：**
```matlab
u = Boxscheme(a, delta_t, delta_x, x_start, x_end, t_start, t_end, phi, g)
```

**输入参数：**
- `a`：对流速度
- `delta_t`：时间步长
- `delta_x`：空间步长
- `x_start`：空间起点
- `x_end`：空间终点
- `t_start`：时间起点
- `t_end`：时间终点
- `phi`：边界条件函数句柄，`u(t, 0)`
- `g`：初始条件函数句柄，`u(0, x)`

**输出参数：**
- `u`：解矩阵，其中每一行表示一个时间步长的状态向量。

### FDMscheme/BeamWarming.m
该函数使用 Beam-Warming 方法求解偏微分方程 `u_t + a u_x = 0`。

**语法：**
```matlab
u = BeamWarming(a, delta_t, delta_x, x_start, x_end, t_start, t_end, phi, g)
```

**输入参数：**
- `a`：对流速度
- `delta_t`：时间步长
- `delta_x`：空间步长
- `x_start`：空间起点
- `x_end`：空间终点
- `t_start`：时间起点
- `t_end`：时间终点
- `phi`：边界条件函数句柄，`u(t, 0)`
- `g`：初始条件函数句柄，`u(0, x)`

**输出参数：**
- `u`：解矩阵，其中每一行表示一个时间步长的状态向量。

### FDMscheme/Implict_central.m
该函数使用隐式中心差分方案求解偏微分方程 `u_t + a u_x = 0`。

**语法：**
```matlab
u = Implict_central(a, delta_t, delta_x, x_start, x_end, t_start, t_end, phi, g)
```

**输入参数：**
- `a`：对流速度
- `delta_t`：时间步长
- `delta_x`：空间步长
- `x_start`：空间起点
- `x_end`：空间终点
- `t_start`：时间起点
- `t_end`：时间终点
- `phi`：边界条件函数句柄，`u(t, 0)`
- `g`：初始条件函数句柄，`u(0, x)`

**输出参数：**
- `u`：解矩阵，其中每一行表示一个时间步长的状态向量。

### FDMscheme/NonlinearHyperbolicSolver.m
该函数使用不同的数值格式求解非线性双曲方程 `u_t + f(u)_x = 0`。

**语法：**
```matlab
u = NonlinearHyperbolicSolver(f, delta_t, delta_x, x_start, x_end, t_start, t_end, scheme, phi, g)
```

**输入参数：**
- `f`：非线性通量函数
- `delta_t`：时间步长
- `delta_x`：空间步长
- `x_start`：空间起点
- `x_end`：空间终点
- `t_start`：时间起点
- `t_end`：时间终点
- `scheme`：使用的数值格式（'BW', 'CD', 'LF', 'LW', 'LWorigin', 'CIR', 'MCIR'）
- `phi`：边界条件函数句柄，`u(t, 0)`
- `g`：初始条件函数句柄，`u(0, x)`

**输出参数：**
- `u`：解矩阵，其中每一行表示一个时间步长的状态向量。

