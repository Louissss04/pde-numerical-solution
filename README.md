# pde-numerical-solution

# 介绍
本仓库用于存储<偏微分方程数值解>暑期学校算法代码

## 文件夹结构
- `FDMscheme/`：包含实现不同有限差分方法的函数。
    - `FDM3points.m`：三点方案
    - `LeapFrog.m`：Leap-Frog 方案
    - `Boxscheme.m`：Box 方案
    - `BeamWarming.m`：Beam-Warming 方案
    - `Implict_central.m`：隐式-中心差分方案
- `Example/`：演示上述方法的脚本。
    - `FDMAdvection3points1D.m`：调用 `FDM3points` 函数的示例
    - `Example_LeapFrog.m`：调用 `LeapFrog` 函数的示例
    - `Example_Boxscheme.m`：调用 `Boxscheme` 函数的示例
    - `Example_BeamWarming.m`：调用 `BeamWarming` 函数的示例
    - `Example_Implict_central.m`：调用 `Implict_central` 函数的示例

## 文件说明
### FDMscheme/FDM3points.m
该函数使用不同的三点有限差分方法求解偏微分方程 `u_t + a u_x = 0`。

**语法：**
```matlab
u = FDM3points(a, delta_t, delta_x, x_start, x_end, t_start, t_end, scheme, phi, g)
```



