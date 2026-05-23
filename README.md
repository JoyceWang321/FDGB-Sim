# FDGB-Sim
Formate-Driven Glycine Biosynthesis Simulator
🧬 FAFS-Sim: Formate-to-Amino Acid Factory Simulation

!https://img.shields.io/badge/iGEM-2024-yellowgreen.svg !https://img.shields.io/badge/License-MIT-blue.svg  

📌 项目简介

本项目为iGEM竞赛开发的一碳代谢与氮代谢集成动态仿真系统，模拟甲酸（C1源）+ 氨（N源）→ 甘氨酸/丝氨酸的细胞工厂代谢通路，支持参数调优、通量分析与敏感性评估，助力合成生物学路径设计与实验验证。  

✨ 功能特性

• 动态ODE仿真：基于刚性问题求解器（LSODA）模拟8种代谢物的时序变化；  

• 多酶动力学模型：覆盖5种核心酶（FtfL/Fch/MtdA/GCS/SHMT）的米氏-门滕动力学（含可逆反应与产物抑制）；  

• 敏感性分析：量化±10%参数扰动对甘氨酸产量的影响，输出排序表格与条形图；  

• 可视化界面：Streamlit交互式面板，实时展示代谢流、浓度变化与关键指标；  

• 参数灵活可调：支持Vmax、初始浓度、模拟时长的自定义配置。  

🔬 代谢模型（核心反应网络）

反应编号 酶（代码命名） 正向反应（主要方向） 逆向反应（次要方向） 关联物质变化（ODE导数项）
1 FtfL（甲酸-THF连接酶） 甲酸(F) + THF + ATP → 10-甲酰-THF(F10) 无（仅正向驱动） \frac{dF}{dt}=-v_1；\frac{dTHF}{dt}=-v_1；\frac{dF10}{dt}=v_1
2 Fch（10-甲酰-THF环水解酶） 10-甲酰-THF(F10) → 5,10-亚甲基-THF(MH4F) 5,10-亚甲基-THF → 10-甲酰-THF \frac{dF10}{dt}=-v_2；\frac{dMH4F}{dt}=v_2
3 MtdA（5,10-亚甲基-THF脱氢酶） 5,10-亚甲基-THF(MH4F) + NADPH → 5,10-次甲基-THF(mTHF) + NADP 5,10-次甲基-THF + NADP → 5,10-亚甲基-THF + NADPH \frac{dMH4F}{dt}=-v_3；\frac{dmTHF}{dt}=v_3
4 GCS（甘氨酸合成酶/裂解系统） 5,10-次甲基-THF(mTHF) + NH₃ + HCO₃⁻ + NADH → 甘氨酸(Gly) + THF 无（含甘氨酸产物抑制） \frac{dmTHF}{dt}=-v_4；\frac{dGly}{dt}=v_4；\frac{dNH_3}{dt}=-v_4；\frac{dTHF}{dt}=v_4
5 SHMT（丝氨酸羟甲基转移酶） 甘氨酸(Gly) + mTHF → 丝氨酸(Ser) + THF 丝氨酸 + THF → 甘氨酸 + mTHF \frac{dGly}{dt}=-v_5；\frac{dmTHF}{dt}=-v_5；\frac{dSer}{dt}=v_5；\frac{dTHF}{dt}=v_5
  

⚙️ 速率方程（酶动力学模型）

所有方程基于米氏-门滕动力学，含底物耗尽截断（浓度<1e-8 mM时速率归零）：  

1. FtfL（三元底物不可逆反应）

v_1=\frac{V_{\max,\text{FtfL}} \cdot [F] \cdot [THF] \cdot [ATP]}{K_{F}K_{THF}K_{ATP} + K_{THF}K_{ATP}[F] + K_{F}K_{ATP}[THF] + K_{F}K_{THF}[ATP] + [F][THF][ATP]}  

2. Fch（可逆双底物反应）

v_2=\frac{\frac{V_{\max,\text{Fch}}}{K_{F10}}[F10] - \frac{V_{\max,\text{Fch}}/0.54}{K_{MH4F}}[MH4F]}{1 + \frac{[F10]}{K_{F10}} + \frac{[MH4F]}{K_{MH4F}}}  
注：逆向反应隐含平衡常数0.54。  

3. MtdA（可逆四元反应）

v_3=\frac{\frac{V_{\max,\text{MtdA}}}{K_{MH4F}K_{NADPH}}[MH4F][NADPH] - \frac{V_{\max,\text{MtdA}}/4.0}{K_{mTHF}K_{NADP}}[mTHF][NADP]}{1 + \frac{[MH4F]}{K_{MH4F}} + \frac{[NADPH]}{K_{NADPH}} + \frac{[mTHF]}{K_{mTHF}} + \frac{[NADP]}{K_{NADP}} + \frac{[MH4F][NADPH]}{K_{MH4F}K_{NADPH}}}  

4. GCS（不可逆+甘氨酸产物抑制）

v_4=V_{\max,\text{GCS}} \cdot \frac{[mTHF]}{K_{mTHF}+[mTHF]} \cdot \frac{[NH_3]}{K_{NH_3}+[NH_3]} \cdot \frac{[HCO_3^-]}{K_{HCO_3}+[HCO_3^-]} \cdot \frac{1}{1 + \frac{[Gly]}{K_{\text{inh}}}}  
注：甘氨酸抑制常数K_{\text{inh}}=0.1 mM。  

5. SHMT（可逆双底物反应）

v_5=\frac{v_{\text{fwd}} - v_{\text{rev}}}{1 + \frac{[Gly]}{K_{Gly}} + \frac{[mTHF]}{K_{mTHF}} + \frac{[Ser]}{K_{Ser}} + \frac{[THF]}{K_{THF}} + \frac{[Gly][mTHF]}{K_{Gly}K_{mTHF}} + \frac{[Ser][THF]}{K_{Ser}K_{THF}}}  
其中：  
v_{\text{fwd}}=V_{\max,\text{SHMT}} \cdot \frac{[Gly][mTHF]}{K_{Gly}K_{mTHF}}（正向速率）；  
v_{\text{rev}}=\frac{V_{\max,\text{SHMT}}}{K_{\text{eq}}} \cdot \frac{[Ser][THF]}{K_{Ser}K_{THF}}（逆向速率，K_{\text{eq}}=1.2）。  

📊 ODE系统与数值求解

状态变量（浓度，单位：mM）

y=[F, THF, F10, MH4F, mTHF, Gly, Ser, NH_3]  

微分方程组（时间导数）

\begin{cases}  
\frac{dF}{dt}=-v_1 \\  
\frac{dTHF}{dt}=-v_1 + v_4 + v_5 \\  
\frac{dF10}{dt}=v_1 - v_2 \\  
\frac{dMH4F}{dt}=v_2 - v_3 \\  
\frac{dmTHF}{dt}=v_3 - v_4 - v_5 \\  
\frac{dGly}{dt}=v_4 - v_5 \\  
\frac{dSer}{dt}=v_5 \\  
\frac{dNH_3}{dt}=-v_4  
\end{cases}  

求解器配置

• 算法：LSODA（自适应处理刚性问题）；  

• 精度：相对误差rtol=1e-8、绝对误差atol=1e-10；  

• 步长：最大步长30秒、初始步长0.001秒；  

• 采样点：3000个均匀时间点（t_{\text{end}}=模拟时长\times60秒）。  

🛠️ 安装与运行

环境要求

• Python 3.8+；  

• 依赖库：streamlit, numpy, scipy, matplotlib, pandas。  

快速启动

# 克隆仓库
git clone https://github.com/[JoyceWang321]/FAFS-Sim.git
cd FAFS-Sim

# 安装依赖
pip install -r requirements.txt

# 运行仿真系统
streamlit run app.py
  

🎮 使用方法

1. 侧边栏参数配置

参数类型 可调范围（单位） 说明
Vmax FtfL(0.1-10)、Fch(0.1-10)、MtdA(0.1-20)、GCS(0.001-1)、SHMT(0.1-50)（mM/min） 酶最大催化速率
初始浓度 甲酸(0-200)、THF(0-5)、甘氨酸(0-10)、丝氨酸(0-100)、氨(0-200)、HCO₃⁻(0-100)、ATP(0-5)、NADPH(0-5)、NADH(0-5)（mM） 代谢物初始池
模拟时长 1-180分钟 总反应时间
  

2. 主界面操作

• 左侧面板：点击「🚀 Run Simulation」启动仿真，实时显示4张动态图表（氨基酸合成、叶酸循环、底物耗尽、GCS活性）；  

• 右侧面板：点击「📊 Sensitivity Analysis」执行±10%参数扰动分析，输出灵敏度排序与实验优先级建议。  

📜 核心假设与限制条件

核心假设

1. 辅因子（ATP/NADPH/NADH/HCO₃⁻）浓度恒定（无动态变化）；  
2. 米氏动力学适用，仅GCS含甘氨酸竞争性抑制；  
3. 底物浓度<1e-8 mM时速率归零（避免数值发散）；  
4. 无额外旁路反应，代谢流仅沿定义通路转化。  

限制条件

• 参数输入需符合上述范围，否则自动截断；  

• 求解器最大步长≤30秒，确保刚性阶段稳定性；  

• 所有浓度非负（强制y=\max(y,0)）。  

📄 许可证

本项目采用 MIT License，详情见LICENSE文件。  
版权声明：Copyright (c) 2026 [Air2Protein]。  

📚 引用方式

若为学术研究使用，请按以下格式引用：  
@misc{FAFS_Sim_2026,
  title={FAFS-Sim: A Dynamic Kinetic Model for Formate-to-Amino Acid Conversion in Cell Factories},
  author={[Air2Protein]},
  year={2026},
  publisher={GitHub},
  url={https://github.com/[JoyceWang321]/FAFS-Sim}
}
  

🙏 致谢

感谢iGEM社区的开源精神支持，感谢指导老师与团队成员的共同努力！  

---  
祝iGEM竞赛取得佳绩！ 🧪✨