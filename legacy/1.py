import streamlit as st
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# ==========================================
# 1. 侧边栏：参数输入区 (The Input Panel)
# ==========================================
st.set_page_config(layout="wide", page_title="甲酸→甘氨酸 代谢模拟器")

st.sidebar.header("🔬 模拟参数调节")
st.sidebar.markdown("拖动滑块修改参数，右侧会自动更新结果。")

# --- 酶动力学参数 ---
st.sidebar.subheader("1️⃣ 酶动力学参数 (Vmax)")
# 使用 slider 代替硬编码的数字
vmax_ftfl = st.sidebar.slider("FtfL (EC 6.3.4.3)", 0.1, 10.0, 1.6, 0.1)
vmax_fch = st.sidebar.slider("Fch (EC 3.5.4.9)", 0.1, 10.0, 1.9, 0.1)
vmax_mtda = st.sidebar.slider("MtdA (EC 1.5.1.5)", 0.1, 20.0, 6.8, 0.1)
vmax_shmt = st.sidebar.slider("SHMT (EC 2.1.2.1)", 0.1, 50.0, 20.0, 0.1)

# --- 初始浓度 ---
st.sidebar.subheader("2️⃣ 初始浓度 (mM)")
init_formate = st.sidebar.number_input("初始甲酸浓度", 0.0, 100.0, 20.0)
init_thf = st.sidebar.number_input("初始 THF 浓度", 0.0, 1.0, 0.05)

# --- 模拟时间 ---
st.sidebar.subheader("3️⃣ 模拟设置")
sim_time = st.sidebar.slider("模拟时长 (分钟)", 10, 600, 120)

# ==========================================
# 2. 参数库 (现在从侧边栏读取)
# ==========================================
params_FtfL = {'Vmax': vmax_ftfl, 'Km_formate': 22.0, 'Km_THF': 0.8, 'Km_ATP': 0.021}
params_Fch = {'Vmax_f': vmax_fch, 'Keq': 0.54, 'Km_F10': 0.08, 'Km_MH4F': 0.08}
params_MtdA = {'Vmax_f': vmax_mtda, 'Keq': 4.0, 'Km_MH4F': 0.03, 'Km_NADPH': 0.01, 'Km_mTHF': 0.03, 'Km_NADP': 0.01}
params_SHMT = {'Vmax': vmax_shmt, 'Km_mTHF': 0.16}

# ==========================================
# 3. 速率方程 (保持不变)
# ==========================================
def rate_FtfL(F, THF, ATP):
    p = params_FtfL
    num = p['Vmax'] * F * THF * ATP
    den = (p['Km_formate']*p['Km_THF']*p['Km_ATP'] + p['Km_THF']*p['Km_ATP']*F + p['Km_formate']*p['Km_ATP']*THF + p['Km_formate']*p['Km_THF']*ATP + F*THF*ATP)
    return num / den if den != 0 else 0

def rate_Fch(F10, MH4F):
    p = params_Fch
    Vmax_r = p['Vmax_f'] / p['Keq']
    num = p['Vmax_f'] * F10 / p['Km_F10'] - Vmax_r * MH4F / p['Km_MH4F']
    den = 1 + F10/p['Km_F10'] + MH4F/p['Km_MH4F']
    return num / den if den != 0 else 0

def rate_MtdA(MH4F, NADPH, mTHF, NADP):
    p = params_MtdA
    Vmax_r = p['Vmax_f'] / p['Keq']
    num = (p['Vmax_f'] * MH4F * NADPH / (p['Km_MH4F'] * p['Km_NADPH']) - Vmax_r * mTHF * NADP / (p['Km_mTHF'] * p['Km_NADP']))
    den = (1 + MH4F/p['Km_MH4F'] + NADPH/p['Km_NADPH'] + mTHF/p['Km_mTHF'] + NADP/p['Km_NADP'] + MH4F * NADPH / (p['Km_MH4F'] * p['Km_NADPH']))
    return num / den if den != 0 else 0

def rate_SHMT(mTHF):
    p = params_SHMT
    return p['Vmax'] * mTHF / (p['Km_mTHF'] + mTHF)

# ==========================================
# 4. 主ODE系统 (保持不变)
# ==========================================
def system(y, t):
    F, THF, F10, MH4F, mTHF, Gly, ATP, ADP, Pi, NADPH, NADP, NADH, NAD = y
    v1 = rate_FtfL(F, THF, ATP)
    v2 = rate_Fch(F10, MH4F)
    v3 = rate_MtdA(MH4F, NADPH, mTHF, NADP)
    v4 = rate_SHMT(mTHF)
    
    dF_dt = -v1
    dTHF_dt = -v1 + v4
    dF10_dt = v1 - v2
    dMH4F_dt = v2 - v3
    dmTHF_dt = v3 - v4
    dGly_dt = v4
    dATP_dt = -v1
    dADP_dt = v1
    dPi_dt = v1
    dNADPH_dt = -v3
    dNADP_dt = v3
    dNADH_dt = -v4
    dNAD_dt = v4
    return [dF_dt, dTHF_dt, dF10_dt, dMH4F_dt, dmTHF_dt, dGly_dt,
            dATP_dt, dADP_dt, dPi_dt, dNADPH_dt, dNADP_dt, dNADH_dt, dNAD_dt]

# ==========================================
# 5. 运行模拟
# ==========================================
# 只有当用户点击按钮时才运行（防止卡顿）
if st.sidebar.button("🚀 开始模拟", type="primary"):
    
    # 初始条件
    y0 = [init_formate, init_thf, 0, 0, 0, 0, 2.0, 0.1, 1.0, 0.5, 0.05, 0.02, 0.2]
    t = np.linspace(0, sim_time, sim_time * 5)
    
    with st.spinner('模拟计算中...'):
        solution = odeint(system, y0, t)
    
    st.success('模拟完成！')

    # ==========================================
    # 6. 结果展示 (前端页面)
    # ==========================================
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("📈 底物与产物")
        fig1, ax1 = plt.subplots()
        ax1.plot(t, solution[:, 0], 'g-', label='Formate')
        ax1.plot(t, solution[:, 1], 'b-', label='THF')
        ax1.plot(t, solution[:, 5], 'r-', label='Glycine')
        ax1.set_xlabel('Time (min)')
        ax1.set_ylabel('Concentration (mM)')
        ax1.legend()
        ax1.grid(True)
        st.pyplot(fig1)

    with col2:
        st.subheader("🔄 叶酸循环中间体")
        fig2, ax2 = plt.subplots()
        ax2.plot(t, solution[:, 2], 'c-', label='10-F-THF')
        ax2.plot(t, solution[:, 3], 'm-', label='Methenyl-THF')
        ax2.plot(t, solution[:, 4], 'y-', label='Methylene-THF')
        ax2.set_xlabel('Time (min)')
        ax2.legend()
        ax2.grid(True)
        st.pyplot(fig2)

    # 显示关键数据
    st.subheader("📊 最终结果")
    final_gly = solution[-1, 5]
    final_formate = solution[-1, 0]
    
    metric_col1, metric_col2, metric_col3 = st.columns(3)
    metric_col1.metric("最终甘氨酸浓度", f"{final_gly:.4f} mM")
    metric_col2.metric("剩余甲酸浓度", f"{final_formate:.4f} mM")
    metric_col3.metric("转化率", f"{(final_gly/init_formate)*100:.2f} %")

else:
    st.info("请在左侧调整参数，然后点击 '开始模拟'。")