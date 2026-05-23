# FDGB-Sim
Formate-Driven Glycine Biosynthesis Simulator
甲酸制氨基酸细胞工厂模拟仿真系统说明
一、模型架构与代谢逻辑 
1.1 核心反应网络（新版实现） 
酶代号	酶全称	反应方向	关键特征
FtfL
甲酸-THF连接酶	不可逆：Formate + THF + ATP → 10-Formyl-THF	三元底物饱和动力学
Fch
10-甲酰-THF环水解酶	不可逆：10-Formyl-THF → 5,10-Methenyl-THF	仅正向催化
MtdA
5,10-亚甲基-THF脱氢酶	不可逆：5,10-Methylene-THF + NADPH → 5,10-Methenyl-THF + NADP⁺	仅正向催化
GCS
甘氨酸合成酶系统	正向：5,10-Methenyl-THF + NH₃ + HCO₃⁻ + NADH → Glycine + THF
逆向：Glycine + THF → 5,10-Methenyl-THF + NH₃ + HCO₃⁻ + NADH	甘氨酸产物抑制
SHMT
丝氨酸羟甲基转移酶	可逆：Glycine + 5,10-Methylene-THF ⇌ Serine + THF	THF依赖型双向反应
1.2 物质守恒与边界条件 
1.	THF总库守恒：游离THF + 结合态叶酸（F10/MH4F/mTHF）总量动态平衡
2.	底物耗尽截断：任意底物浓度＜1e-8 mM时，对应反应速率强制归零（避免数值发散）
3.	初始条件约束：mTHF初始值 = max(0.5, Vmax_GCS × 20)（确保启动期叶酸供应）

二、酶动力学方程 
2.1 FtfL（三元不可逆反应） 
 
2.2 Fch（不可逆） 
 

2.3 MtdA（不可逆） 
 
2.4 GCS（主正向+弱逆向+产物抑制） 
 

2.5 SHMT（THF依赖型可逆反应） 
 
其中：
 

三、ODE系统与数值求解 
3.1 微分方程定义 
变量	含义	导数方程
 
甲酸	 

 
四氢叶酸	 

 
10-甲酰-THF	 

 
5,10-亚甲基-THF	 

 
5,10-次甲基-THF	 

 
甘氨酸	 

 
丝氨酸	 

 
氨	 

3.2 求解器配置 
1.	算法：LSODA（自动切换刚性/非刚性模式）
2.	时间步长：最大步长5.0秒，初始步长0.001秒
3.	误差控制：相对误差1e-6，绝对误差1e-8
4.	输出节点：3000个均匀时间点（覆盖全程动态）

四、软件操作与界面功能 
4.1 侧边栏参数设置 
模块	可调参数	默认值	物理意义
Vmax
FtfL/Fch/MtdA/GCS/SHMT	见代码DEFAULT_VMAX	最大催化速率（mM/min）
初始浓度
Formate/THF/Gly/Ser/NH₃/HCO₃⁻/ATP/NADPH/NADH	见代码DEFAULT_INIT	初始底物池（mM）
模拟时长
Simulation Time	120分钟	总反应时间
4.2 主界面布局 
区域	功能	说明
左侧主面板
🚀 Run Simulation按钮
📈 模拟结果图表（3子图）
📊 关键指标卡片	-
右侧侧边栏
📊 Sensitivity Analysis按钮
📝 敏感性分析结果（表格+条形图）	需先运行模拟，自动读取左侧结果
4.3 敏感性分析方法 
1.	扰动范围：各参数±10%（受生理上下限约束）
2.	灵敏度指数： 
3.	输出内容：排序表格、正负效应条形图、实验优先级建议

五、本地部署与运行指南 
5.1 环境准备 
1.	安装Python：需Python 3.8及以上版本（推荐通过Python官网或Anaconda安装）。
5.2 依赖安装 
代码所需全部依赖如下，通过pip一键安装：
pip install streamlit numpy scipy matplotlib pandas
注：若安装缓慢，可添加国内镜像源（如-i https://pypi.tuna.tsinghua.edu.cn/simple）。
5.3 运行应用 
1.	将代码文件（如9.py）保存到本地文件夹（例：D:\glycine_model）。
2.	打开终端，进入代码所在目录：
cd D:\glycine_model
3.	启动Streamlit应用：
streamlit run 9.py
4.	浏览器会自动弹出访问地址（默认：http://localhost:8501），即可开始使用。
5.4 常见问题排查 
问题现象	解决方法
端口被占用（提示Address already in use）	换端口运行：streamlit run 9.py --server.port 8502（数字可自定义）
依赖安装失败（如scipy编译错误）	升级pip后重试：pip install --upgrade pip，或安装预编译包（Windows推荐用conda install scipy）
运行时提示ModuleNotFoundError	确认依赖已正确安装（可重新执行pip install命令）

