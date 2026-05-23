import streamlit as st
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import time

st.set_page_config(layout="wide", page_title="Glycine Metabolism Model")
st.sidebar.header("🔬 Metabolic Parameters")

# ============================================================================
# 1. DEFAULT VALUES
# ============================================================================
DEFAULT_VMAX = {
    "FtfL": 1.6,
    "Fch": 1.9,
    "MtdA": 6.8,
    "GCS": 0.2,
    "SHMT": 15.0
}

DEFAULT_KM = {
    "FtfL_F": 22.0, "FtfL_THF": 0.8, "FtfL_ATP": 0.021,
    "Fch_F10": 0.08, "Fch_MH4F": 0.08,
    "MtdA_MH4F": 0.03, "MtdA_NADPH": 0.01, "MtdA_mTHF": 0.03, "MtdA_NADP": 0.01,
    "GCS_mTHF": 0.0677, "GCS_NH3": 65.4, "GCS_HCO3": 3.4, "GCS_NADH": 0.02,
    "GCS_Gly_inh": 0.1, "GCS_Gly_rev": 5.0, "GCS_THF_rev": 0.2,
    "SHMT_Gly": 5.0, "SHMT_mTHF": 0.05, "SHMT_Ser": 3.0, "SHMT_THF": 0.2
}

DEFAULT_INIT = {
    "Formate": 50.0, "THF": 1.0, "Glycine": 0.0, "Serine": 5.0,
    "Ammonia": 100.0, "HCO3": 50.0, "ATP": 5.0, "NADPH": 1.0, "NADH": 0.5
}

# ============================================================================
# 2. RESET BUTTON
# ============================================================================
if st.sidebar.button("🔄 Reset to Defaults"):
    st.session_state.vmax_ftfl = DEFAULT_VMAX["FtfL"]
    st.session_state.vmax_fch = DEFAULT_VMAX["Fch"]
    st.session_state.vmax_mtda = DEFAULT_VMAX["MtdA"]
    st.session_state.vmax_gcs = DEFAULT_VMAX["GCS"]
    st.session_state.vmax_shmt = DEFAULT_VMAX["SHMT"]

# ============================================================================
# 3. VMAX INPUTS
# ============================================================================
st.sidebar.subheader("Vmax (mM/min)")
vmax_ftfl = st.sidebar.number_input("FtfL", 0.01, 10.0, DEFAULT_VMAX["FtfL"], 0.01, key="vmax_ftfl")
vmax_fch = st.sidebar.number_input("Fch", 0.01, 10.0, DEFAULT_VMAX["Fch"], 0.01, key="vmax_fch")
vmax_mtda = st.sidebar.number_input("MtdA", 0.01, 20.0, DEFAULT_VMAX["MtdA"], 0.01, key="vmax_mtda")
vmax_gcs = st.sidebar.number_input("GCS", 0.01, 10.0, DEFAULT_VMAX["GCS"], 0.01, key="vmax_gcs")
vmax_shmt = st.sidebar.number_input("SHMT", 0.01, 50.0, DEFAULT_VMAX["SHMT"], 0.01, key="vmax_shmt")

# ============================================================================
# 4. INITIAL CONCENTRATIONS
# ============================================================================
st.sidebar.subheader("Initial Concentrations (mM)")
init_formate = st.sidebar.number_input("Formate", 0.0, 200.0, DEFAULT_INIT["Formate"])
init_thf = st.sidebar.number_input("THF", 0.0, 5.0, DEFAULT_INIT["THF"])
init_gly = st.sidebar.number_input("Glycine", 0.0, 10.0, DEFAULT_INIT["Glycine"])
init_ser = st.sidebar.number_input("Serine", 0.0, 100.0, DEFAULT_INIT["Serine"])
init_nh3 = st.sidebar.number_input("Ammonia", 0.0, 200.0, DEFAULT_INIT["Ammonia"])
init_hco3 = st.sidebar.number_input("HCO3-", 0.0, 100.0, DEFAULT_INIT["HCO3"])
init_atp = st.sidebar.number_input("ATP", 0.0, 5.0, DEFAULT_INIT["ATP"])
init_nadph = st.sidebar.number_input("NADPH", 0.0, 5.0, DEFAULT_INIT["NADPH"])
init_nadh = st.sidebar.number_input("NADH", 0.0, 5.0, DEFAULT_INIT["NADH"])

# ============================================================================
# 5. FIXED KINETIC CONSTANTS
# ============================================================================
Km_FtfL = {"F": DEFAULT_KM["FtfL_F"], "THF": DEFAULT_KM["FtfL_THF"], "ATP": DEFAULT_KM["FtfL_ATP"]}
Km_Fch = {"F10": DEFAULT_KM["Fch_F10"], "MH4F": DEFAULT_KM["Fch_MH4F"]}
Km_MtdA = {"MH4F": DEFAULT_KM["MtdA_MH4F"], "NADPH": DEFAULT_KM["MtdA_NADPH"], 
           "mTHF": DEFAULT_KM["MtdA_mTHF"], "NADP": DEFAULT_KM["MtdA_NADP"]}
Km_GCS = {"mTHF": DEFAULT_KM["GCS_mTHF"], "NH3": DEFAULT_KM["GCS_NH3"], 
          "HCO3": DEFAULT_KM["GCS_HCO3"], "NADH": DEFAULT_KM["GCS_NADH"],
          "GCS_Gly_inh": DEFAULT_KM["GCS_Gly_inh"],
          "GCS_Gly_rev": DEFAULT_KM["GCS_Gly_rev"],
          "GCS_THF_rev": DEFAULT_KM["GCS_THF_rev"]}
Km_SHMT = {"Gly": DEFAULT_KM["SHMT_Gly"], "mTHF": DEFAULT_KM["SHMT_mTHF"], 
           "Ser": DEFAULT_KM["SHMT_Ser"], "THF": DEFAULT_KM["SHMT_THF"]}
Keq_shmt = 1.2

# ============================================================================
# 6. RATE EQUATIONS
# ============================================================================
def safe_divide(num, den):
    return num / den if den > 1e-12 else 0

def rate_ftfl(F, THF, ATP):
    if F < 1e-8 or THF < 1e-8 or ATP < 1e-8:
        return 0
    p = Km_FtfL
    num = vmax_ftfl * F * THF * ATP
    den = (p["F"]*p["THF"]*p["ATP"] + p["THF"]*p["ATP"]*F +
           p["F"]*p["ATP"]*THF + p["F"]*p["THF"]*ATP + F*THF*ATP)
    return safe_divide(num, den)

def rate_fch(F10, MH4F):
    if F10 < 1e-8:
        return 0
    p = Km_Fch
    num = vmax_fch * F10 / p["F10"]
    den = 1 + F10/p["F10"] + MH4F/p["MH4F"]
    return safe_divide(num, den)

def rate_mtda(MH4F, NADPH, mTHF, NADP):
    if MH4F < 1e-8:
        return 0
    p = Km_MtdA
    num = vmax_mtda * MH4F * NADPH / (p["MH4F"]*p["NADPH"])
    den = (1 + MH4F/p["MH4F"] + NADPH/p["NADPH"] +
           mTHF/p["mTHF"] + NADP/p["NADP"] +
           MH4F*NADPH/(p["MH4F"]*p["NADPH"]))
    return safe_divide(num, den)

def rate_gcs(mTHF, NH3, NADH, HCO3, Gly, THF):
    if any(c < 1e-8 for c in [mTHF, NH3, HCO3, NADH, Gly, THF]):
        return 0
    
    p = Km_GCS
    
    # 正向（合成）
    v_fwd = vmax_gcs * (
        mTHF/(p["mTHF"]+mTHF) *
        NH3/(p["NH3"]+NH3) *
        HCO3/(p["HCO3"]+HCO3) *
        NADH/(p["NADH"]+NADH)
    ) * (1/(1+Gly/p["GCS_Gly_inh"]))
    
    # 逆向（裂解）- 进一步降低逆向强度
    v_rev = vmax_gcs * 0.02 * (  # 从 0.05 降到 0.02
        Gly/(p["GCS_Gly_rev"]+Gly) *
        THF/(p["GCS_THF_rev"]+THF)
    )
    
    return safe_divide(v_fwd - v_rev, 1)

def rate_shmt(Gly, mTHF, Ser, THF):
    if Gly < 1e-8 and Ser < 1e-8:
        return 0
    
    p = Km_SHMT
    
    if THF < 0.01:
        v_fwd = 0
        v_rev = (vmax_shmt / Keq_shmt) * (Ser * THF) / (p["Ser"] * p["THF"])
    else:
        v_fwd = vmax_shmt * (Gly * mTHF) / (p["Gly"] * p["mTHF"])
        v_rev = (vmax_shmt / Keq_shmt) * (Ser * THF) / (p["Ser"] * p["THF"])
    
    den = (1 + Gly/p["Gly"] + mTHF/p["mTHF"] + 
           Ser/p["Ser"] + THF/p["THF"] +
           (Gly*mTHF)/(p["Gly"]*p["mTHF"]) +
           (Ser*THF)/(p["Ser"]*p["THF"]))
    
    return safe_divide(v_fwd - v_rev, den)

# ============================================================================
# 7. ODE SYSTEM
# ============================================================================
def system(t, y):
    y = np.maximum(y, 0)
    F, THF, F10, MH4F, mTHF, Gly, Ser, NH3 = y
    
    ATP = init_atp
    NADPH = init_nadph
    NADH = init_nadh
    NADP = 0.05
    HCO3 = init_hco3
    
    # THF 保护：过低时暂停 FtfL
    if THF < 0.005:
        v1 = 0
    else:
        v1 = rate_ftfl(F, THF, ATP)
    
    v2 = rate_fch(F10, MH4F)
    v3 = rate_mtda(MH4F, NADPH, mTHF, NADP)
    v4 = rate_gcs(mTHF, NH3, NADH, HCO3, Gly, THF)
    v5 = rate_shmt(Gly, mTHF, Ser, THF)
    
    dF = -v1
    dTHF = -v1 + v4 + v5
    dF10 = v1 - v2
    dMH4F = v2 - v3
    dmTHF = v3 - v4 - v5
    dGly = v4 - v5
    dSer = v5
    dNH3 = -v4
    
    return [dF, dTHF, dF10, dMH4F, dmTHF, dGly, dSer, dNH3]

# ============================================================================
# 8. SIMULATION SETTINGS
# ============================================================================
st.sidebar.subheader("Simulation Settings")
sim_minutes = st.sidebar.slider("Simulation Time (minutes)", 1, 180, 120)

if st.sidebar.button("🚀 Run Simulation", type="primary"):
    start_time = time.time()
    
    # 根据 GCS Vmax 动态调整 mTHF 初始浓度
    mthf_init = max(0.5, vmax_gcs * 20)  # 关键：GCS 越高，mTHF 初始越多
    y0 = [init_formate, init_thf, 1e-6, 1e-6, mthf_init, init_gly, init_ser, init_nh3]
    
    t_end = sim_minutes * 60
    t_eval = np.linspace(0, t_end, 3000)
    
    try:
        with st.spinner(f'Simulating {sim_minutes} minutes (GCS Vmax={vmax_gcs})...'):
            sol = solve_ivp(
                system, 
                [0, t_end], 
                y0, 
                method="LSODA",
                t_eval=t_eval,
                rtol=1e-6,           # 放宽容差
                atol=1e-8,
                max_step=5.0,        # 限制最大步长
                first_step=0.001,
                max_num_steps=100000 # 限制最大步数，防止无限循环
            )
        
        elapsed = time.time() - start_time
        st.write(f"⏱️ 模拟耗时: {elapsed:.2f} 秒")
        
        if sol.success:
            st.success(f'✅ Simulation completed! Steps: {sol.nfev}')
            
            t = sol.t / 60
            F, THF, F10, MH4F, mTHF, Gly, Ser, NH3 = sol.y
            
            # 检查是否有问题
            if np.min(THF) < 0.001:
                st.warning("⚠️ THF 接近零！建议降低 FtfL 活性或提高 mTHF 初始值。")
            if np.any(np.isnan(sol.y)):
                st.error("❌ 检测到 NaN 值！数值不稳定。")
            
            # 绘图
            fig, axes = plt.subplots(2, 2, figsize=(14, 10))
            
            axes[0,0].plot(t, F, 'g-', lw=2, label='Formate')
            axes[0,0].plot(t, Gly, 'r-', lw=3, label='Glycine')
            axes[0,0].plot(t, Ser, 'b-', lw=2, label='Serine')
            axes[0,0].set_xlabel('Time (min)')
            axes[0,0].set_ylabel('Concentration (mM)')
            axes[0,0].set_title(f'Glycine Metabolism (GCS Vmax={vmax_gcs})')
            axes[0,0].legend()
            axes[0,0].grid(True, alpha=0.3)
            
            axes[0,1].plot(t, THF, 'c-', lw=2, label='THF')
            axes[0,1].plot(t, mTHF, 'k-', lw=2, label='mTHF')
            axes[0,1].plot(t, F10, 'm-', lw=2, label='10-F-THF')
            axes[0,1].plot(t, MH4F, 'y-', lw=2, label='5,10-M-THF')
            axes[0,1].set_xlabel('Time (min)')
            axes[0,1].set_ylabel('Concentration (mM)')
            axes[0,1].set_title('Folate Cycle')
            axes[0,1].legend()
            axes[0,1].grid(True, alpha=0.3)
            
            axes[1,0].plot(t, Gly, 'r-', lw=3, label='Glycine')
            axes[1,0].plot(t, Ser, 'b-', lw=2, label='Serine')
            axes[1,0].set_xlabel('Time (min)')
            axes[1,0].set_ylabel('Concentration (mM)')
            axes[1,0].set_title('Glycine & Serine')
            axes[1,0].legend()
            axes[1,0].grid(True, alpha=0.3)
            
            axes[1,1].plot(t, F, 'g-', lw=2, label='Formate')
            axes[1,1].plot(t, NH3, 'c-', lw=2, label='Ammonia')
            axes[1,1].set_xlabel('Time (min)')
            axes[1,1].set_ylabel('Concentration (mM)')
            axes[1,1].set_title('Substrate Consumption')
            axes[1,1].legend()
            axes[1,1].grid(True, alpha=0.3)
            
            plt.tight_layout()
            st.pyplot(fig)
            
            # 结果显示
            st.subheader("📊 Results")
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Final Glycine", f"{Gly[-1]:.4f} mM")
            col2.metric("Final Serine", f"{Ser[-1]:.4f} mM")
            col3.metric("Min THF", f"{np.min(THF):.6f} mM")
            col4.metric("Steps", f"{sol.nfev}")
                
        else:
            st.error(f"❌ Simulation failed: {sol.message}")
            st.info("建议：降低 GCS Vmax 或增加 mTHF 初始浓度")
            
    except Exception as e:
        st.error(f"❌ Error: {str(e)}")
        st.info("建议：降低 GCS Vmax 到 0.05 以下")

else:
    st.info("Adjust parameters and click 'Run Simulation'.")