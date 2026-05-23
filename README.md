# FDGB-Sim
Formate-Driven Glycine Biosynthesis Simulator


# 项目简介

本项目为2026年iGEM竞赛开发的一碳代谢与氮代谢集成动态仿真系统，模拟甲酸（C1源）+ 氨（N源）→ 甘氨酸/丝氨酸的细胞工厂代谢通路，支持参数调优、通量分析与敏感性评估，助力合成生物学路径设计与实验验证。  

# 功能特性

• 动态ODE仿真：基于刚性问题求解器（LSODA）模拟8种代谢物的时序变化；  

• 多酶动力学模型：覆盖5种核心酶（FtfL/Fch/MtdA/GCS/SHMT）的米氏-门滕动力学（含可逆反应与产物抑制）；  

• 敏感性分析：量化±10%参数扰动对甘氨酸产量的影响，输出排序表格与条形图；  

• 可视化界面：Streamlit交互式面板，实时展示代谢流、浓度变化与关键指标；  

• 参数灵活可调：支持Vmax、初始浓度、模拟时长的自定义配置。  

# 代谢模型（核心反应网络）

| 酶代号 | 酶全称 | 反应方向 | 关键特征 |
|:---|:---|:---|:---|
| `FtfL` | 甲酸-THF连接酶 | 不可逆：<br>Formate + THF + ATP → 10-Formyl-THF | 三元底物饱和动力学 |
| `Fch` | 10-甲酰-THF环水解酶 | 不可逆：<br>10-Formyl-THF → 5,10-Methenyl-THF | 体正向催化 |
| `MtdA` | 5,10-亚甲基-THF脱氢酶 | 不可逆：<br>5,10-Methylene-THF + NADPH → 5,10-Methenyl-THF + NADP⁺ | 体正向催化 |
| `GCS` | 甘氨酸合成酶系统 | 正向：5,10-Methenyl-THF + NH₃ + HCO₃⁻ + NADH → Glycine + THF<br><br>逆向：Glycine + THF → 5,10-Methenyl-THF + NH₃ + HCO₃⁻ + NADH | 甘氨酸产物抑制 |
| `SHMT` | 丝氨酸羟甲基转移酶 | 可逆：<br>Glycine + 5,10-Methylene-THF ⇌ Serine + THF | THF依赖型双向反应 |
  

# 酶动力学方程（米氏-门滕模型）

所有方程含底物耗尽截断：任意底物浓度<1e-8 mM时，速率强制归零（避免数值发散）。  


## 1. FtfL（三元不可逆反应）

![FtfL 反应](./eq1.svg)

## 2. Fch（不可逆）

![Fch 反应](./eq2.svg)

## 3. MtdA（不可逆）

![Mtda 反应](./eq3.svg)

## 4. GCS（主正向+弱逆向+产物抑制）

![GCS 反应](./eq4.svg)

## 5. SHMT（THF依赖型可逆反应）

![SHMT 反应](./eq5.svg)

## 其中：

<img src="https://latex.codecogs.com/svg.image?v_{\text{fwd}}=V_{\text{max,SHMT}}\cdot&space;\frac{[Gly][mTHF]}{K_{Gly}K_{mTHF}}" alt="Equation 1" />

<img src="https://latex.codecogs.com/svg.image?v_{\text{rev}}=\frac{V_{\text{max,SHMT}}}{K_{\text{eq}}}\cdot&space;\frac{[Ser][THF]}{K_{Ser}K_{THF}}K_{\text{eq}}=1.2&space;&space;&space;&space;&space;" alt="Equation 1" />


# ODE系统与数值求解

状态变量（浓度，单位：mM）

y=[F, THF, F10, MH4F, mTHF, Gly, Ser, NH_3]  

微分方程组（时间导数）

| 变量 | 含义 | 导数方程 |
|:---:|:---:|:---:|
| `[F]` | 甲酸 | ![](https://latex.codecogs.com/svg.image?\frac{d[F]}{dt}=-v_1) |
| `[THF]` | 四氢叶酸 | ![](https://latex.codecogs.com/svg.image?\frac{d[THF]}{dt}=-v_1+v_4+v_5) |
| `[F10]` | 10-甲酰-THF | ![](https://latex.codecogs.com/svg.image?\frac{d[F10]}{dt}=v_1-v_2) |
| `[MH4F]` | 5,10-亚甲基-THF | ![](https://latex.codecogs.com/svg.image?\frac{d[MH4F]}{dt}=v_2-v_3) |
| `[mTHF]` | 5,10-次甲基-THF | ![](https://latex.codecogs.com/svg.image?\frac{d[mTHF]}{dt}=v_3-v_4-v_5) |
| `[Gly]` | 甘氨酸 | ![](https://latex.codecogs.com/svg.image?\frac{d[Gly]}{dt}=v_4-v_5) |
| `[Ser]` | 丝氨酸 | ![](https://latex.codecogs.com/svg.image?\frac{d[Ser]}{dt}=v_5) |
| `[NH_3]` | 氨 | ![](https://latex.codecogs.com/svg.image?\frac{d[NH_3]}{dt}=-v_4) |


# 求解器配置

• 算法：LSODA（自适应切换刚性/非刚性模式）；  

• 时间步长：最大步长5.0秒，初始步长0.001秒；  

• 输出节点：3000个均匀时间点（覆盖全程动态）。  

# 软件操作与界面功能

## 4.1 侧边栏参数设置

| 模块 | 可调参数 | 默认值 | 物理意义 |
|:---|:---|:---|:---|
| Vmax | FtfL / Fch / MtdA / GCS / SHMT | `DEFAULT_VMAX` (见代码) | 最大催化速率 (mM/min) |
| 初始浓度 | Formate / THF / Gly / Ser / NH₃ / HCO₃⁻ / ATP / NADPH / NADH | `DEFAULT_INIT` (见代码) | 初始底物池 (mM) |
| 模拟时长 | Simulation Time | 120 分钟 | 总反应时间 |

## 4.2 主界面布局
| 区域 | 功能 | 说明 |
|:---|:---|:---|
| **左侧主面板** | 🚀 Run Simulation 按钮<br>📈 模拟结果图表（3子图）<br>📊 关键指标卡片 | 图表保留，不受敏感性分析干扰 |
| **右侧侧边栏** | 📊 Sensitivity Analysis 按钮<br>📝 敏感性分析结果（表格+条形图） | 需先运行模拟，自动读取左侧结果 |
  

## 4.3 敏感性分析方法

1. 扰动范围：各参数±10%（受生理上下限约束）；  
2. 灵敏度指数：<img src="https://latex.codecogs.com/svg.image?SI=\frac{\Delta&space;Gly&space;/&space;Gly_{\text{base}}}{\Delta&space;Param&space;/&space;Param_{\text{base}}}" alt="Equation 1" />
4. 输出内容：排序表格、正负效应条形图、实验优先级建议。  

# 本地部署与运行指南

## 5.1 环境准备

• 安装 Python 3.8及以上版本（推荐通过https://www.python.org/或https://www.anaconda.com/安装）。  

## 5.2 依赖安装

代码所需全部依赖如下，通过pip一键安装：  
pip install streamlit numpy scipy matplotlib pandas
  
注：若安装缓慢，可添加国内镜像源（如-i https://pypi.tuna.tsinghua.edu.cn/simple）。  

## 5.3 运行应用

1. 将代码文件（如9.py）保存到本地文件夹（例：D:\glycine_model）；  
2. 打开终端，进入代码所在目录：  
cd D:\glycine_model
  
3. 启动Streamlit应用：  
streamlit run 9.py
  
4. 浏览器会自动弹出访问地址（默认：http://localhost:8501），即可开始使用。  

## 5.4 常见问题排查

| 问题现象 | 解决方法 |
|:---|:---|
| 端口被占用<br>（提示 Address already in use） | 换端口运行：<br>`streamlit run 9.py --server.port 8502`<br>（数字可自定义） |
| 依赖安装失败<br>（如 scipy 编译错误） | 1. 升级 pip：`pip install --upgrade pip`<br>2. Windows 推荐使用 Conda：`conda install scipy` |
| 运行时提示 ModuleNotFoundError | 确认虚拟环境已激活，且依赖已正确安装：<br>重新执行 `pip install -r requirements.txt` |
  

# 许可证

本项目采用 MIT License，详情见LICENSE文件。  
版权声明：Copyright (c) 2026 [Air2Protein]。  

# 引用方式

若为学术研究使用，请按以下格式引用：  
@misc{FAFS_Sim_2026,
  title={FAFS-Sim: A Dynamic Kinetic Model for Formate-to-Amino Acid Conversion in Cell Factories},
  author={[Air2Protein]},
  year={2026},
  publisher={GitHub},
  url={https://github.com/[JoyceWang321]/FAFS-Sim}
}
  

# 致谢

感谢iGEM社区的开源精神支持，感谢指导老师与团队成员的共同努力！  

---  
祝2026 iGEM竞赛取得佳绩！ 🧪🏆
