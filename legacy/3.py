import streamlit as st
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ============================================================================
# 全局设置：解决中文方框问题 & 使用英文标签
# ============================================================================
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.unicode_minus'] = False

st.set_page_config(layout="wide", page_title="E. coli Glycine Synthesis Simulator")
st.title("E. coli Glycine Synthesis Simulator (Based on Crowther 2008)")

# ============================================================================
# 侧边栏参数 (严格基于文献，不做随意修改)
# ============================================================================
with st.sidebar:
    st.header("Enzyme Kinetics (Vmax, mM/s)")
    
    vmax_ftfl = st.slider("FtfL (EC 6.3.4.3)", 0.1, 5.0, 1.6, 0.1)
    vmax_fch = st.slider("Fch (EC 3.5.4.9)", 0.1, 5.0, 1.9, 0.1)
    vmax_mtda = st.slider("MtdA (EC 1.5.1.5)", 0.1, 10.0, 6.8, 0.1)
    vmax_gcs = st.slider("GCS Complex", 0.01, 5.0, 0.5, 0.01)  # 关键修正：从20降到0.5
    vmax_dhfr = st.slider("DHFR (THF Regeneration)", 0.1, 10.0, 3.0, 0.1)
    vmax_regen = st.slider("NADH Regeneration Rate", 0.01, 5.0, 0.5, 0.01)
    
    st.markdown("---")
    st.header("Initial Concentrations (mM)")
    
    # 底物
    init_formate = st.number_input("Formate", 0.0, 200.0, 20.0)
    init_thf = st.number_input("THF", 0.0, 2.0, 0.2)
    init_dhf = st.number_input("DHF", 0.0, 2.0, 0.05)
    init_nh3 = st.number_input("Ammonia", 0.0, 200.0, 100.0)
    init_co2 = st.number_input("CO2", 0.0, 50.0, 5.0)
    
    # 辅因子 (Crowther 2008 Table 1)
    init_atp = st.number_input("ATP", 0.0, 5.0, 2.88)
    init_adp = st.number_input("ADP", 0.0, 5.0, 0.60)
    init_pi = st.number_input("Phosphate", 0.0, 5.0, 3.0)
    init_nadph = st.number_input("NADPH", 0.0, 5.0, 0.42)
    init_nadp = st.number_input("NADP+", 0.0, 5.0, 0.62)
    init_nadh = st.number_input("NADH", 0.0, 5.0, 0.53)
    init_nad = st.number_input("NAD+", 0.0, 5.0, 1.69)
    
    st.markdown("---")
    sim_hours = st.slider("Simulation Time (hours)", 0.1, 10.0, 2.0)

# ============================================================================
# 固定动力学常数 (Km values from Crowther 2008)
# ============================================================================
Km_FtfL = {"F": 22.0, "THF": 0.8, "ATP": 0.021}
Km_Fch = {"F10": 0.08, "MH4F": 0.08}
Km_MtdA = {"mTHF": 0.03, "NADP": 0.01, "MH4F": 0.03, "NADPH": 0.01}
Km_GCS = {"mTHF": 0.0677, "NH3": 65.4, "CO2": 3.4, "NADH": 0.02}
Km_DHFR = {"DHF": 0.01, "NADPH": 0.01}

# ============================================================================
# 速率方程 (包含再生机制)
# ============================================================================
def rate_ftfl(F, THF, ATP):
    p = Km_FtfL
    num = vmax_ftfl * F * THF * ATP
    den = (p["F"]*p["THF"]*p["ATP"] + p["THF"]*p["ATP"]*F + 
           p["F"]*p["ATP"]*THF + p["F"]*p["THF"]*ATP + F*THF*ATP)
    return num / den if den > 0 else 0

def rate_fch(F10, MH4F):
    p = Km_Fch
    v_rev = vmax_fch / 0.54
    num = (vmax_fch * F10 / p["F10"]) - (v_rev * MH4F / p["MH4F"])
    den = 1 + F10/p["F10"] + MH4F/p["MH4F"]
    return num / den if den > 0 else 0

def rate_mtda(mTHF, NADP, MH4F, NADPH):
    """Corrected: mTHF + NADP+ -> MH4F + NADPH, Keq=4"""
    p = Km_MtdA
    v_rev = vmax_mtda / 4.0
    num = (vmax_mtda * mTHF * NADP / (p["mTHF"]*p["NADP"])) - \
          (v_rev * MH4F * NADPH / (p["MH4F"]*p["NADPH"]))
    den = (1 + mTHF/p["mTHF"] + NADP/p["NADP"] + MH4F/p["MH4F"] + 
           NADPH/p["NADPH"] + mTHF*NADP/(p["mTHF"]*p["NADP"]))
    return num / den if den > 0 else 0

def rate_gcs(mTHF, NH3, CO2, NADH):
    p = Km_GCS
    sat_mTHF = mTHF / (p["mTHF"] + mTHF)
    sat_NH3 = NH3 / (p["NH3"] + NH3)
    sat_CO2 = CO2 / (p["CO2"] + CO2)
    sat_NADH = NADH / (p["NADH"] + NADH)
    return vmax_gcs * sat_mTHF * sat_NH3 * sat_CO2 * sat_NADH

def rate_dhfr(DHF, NADPH):
    """DHFR: DHF + NADPH -> THF + NADP+"""
    p = Km_DHFR
    v = vmax_dhfr * (DHF / (p["DHF"] + DHF)) * (NADPH / (p["NADPH"] + NADPH))
    return v

def rate_nadh_regen(NAD):
    """NADH regeneration via respiratory chain"""
    return vmax_regen * NAD

# ============================================================================
# ODE 系统 (16个状态变量)
# ============================================================================
def system(t, y):
    F, THF, DHF, F10, MH4F, mTHF, Gly, ATP, ADP, Pi, NADPH, NADP, NADH, NAD, NH3, CO2 = y
    
    # 计算反应速率
    v1 = rate_ftfl(F, THF, ATP)
    v2 = rate_fch(F10, MH4F)
    v3 = rate_mtda(mTHF, NADP, MH4F, NADPH)
    v4 = rate_gcs(mTHF, NH3, CO2, NADH)
    v5 = rate_dhfr(DHF, NADPH)
    v6 = rate_nadh_regen(NAD)
    
    # 微分方程
    dF = -v1
    dTHF = -v1 + v4 + v5
    dDHF = -v5
    dF10 = v1 - v2
    dMH4F = v2 - v3
    dmTHF = v3 - v4
    dGly = v4
    dATP = -v1
    dADP = v1
    dPi = v1
    dNADPH = -v3 - v5
    dNADP = v3 + v5
    dNADH = -v4 + v6
    dNAD = v4 - v6
    dNH3 = -v4
    dCO2 = -v4
    
    return [dF, dTHF, dDHF, dF10, dMH4F, dmTHF, dGly, dATP, dADP, dPi,
            dNADPH, dNADP, dNADH, dNAD, dNH3, dCO2]

# ============================================================================
# 运行模拟
# ============================================================================
if st.button("🚀 Run Precise Simulation", type="primary"):
    # 初始条件
    y0 = [init_formate, init_thf, init_dhf, 0, 0, 0, 0, init_atp, init_adp, init_pi,
          init_nadph, init_nadp, init_nadh, init_nad, init_nh3, init_co2]
    
    t_end = int(sim_hours * 3600)
    t_eval = np.linspace(0, t_end, 1000)  # 减少点数以提高速度
    
    with st.spinner('Solving ODE system...'):
        sol = solve_ivp(system, [0, t_end], y0, method='LSODA',
                       t_eval=t_eval, rtol=1e-6, atol=1e-9, max_step=10)
    
    st.success(f'Simulation completed! Time: {sol.t[-1]/3600:.2f} hours')
    t_hours = sol.t / 3600
    
    # ========================================================================
    # 绘图 (全英文标签，避免方框问题)
    # ========================================================================
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))
    
    # Plot 1: Main Reaction
    ax1 = axes[0, 0]
    ax1.plot(t_hours, sol.y[0], 'g-', lw=2, label='Formate')
    ax1.plot(t_hours, sol.y[1], 'b--', lw=2, label='THF')
    ax1.plot(t_hours, sol.y[14], 'c-', lw=2, label='Ammonia')
    ax1.plot(t_hours, sol.y[6], 'r-', lw=3, label='Glycine')
    ax1.set_xlabel('Time (hours)')
    ax1.set_ylabel('Concentration (mM)')
    ax1.set_title('Main Reaction: Formate + NH3 → Glycine')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Folate Cycle
    ax2 = axes[0, 1]
    ax2.plot(t_hours, sol.y[1], 'c-', lw=2, label='THF')
    ax2.plot(t_hours, sol.y[2], 'm-', lw=2, label='DHF')
    ax2.plot(t_hours, sol.y[3], 'm--', lw=2, label='10-formyl-THF')
    ax2.plot(t_hours, sol.y[4], 'y-', lw=2, label='Methenyl-THF')
    ax2.plot(t_hours, sol.y[5], 'k-', lw=2, label='Methylene-THF')
    ax2.set_xlabel('Time (hours)')
    ax2.set_ylabel('Concentration (mM)')
    ax2.set_title('Folate Cycle (with THF Regeneration)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Redox Cofactors
    ax3 = axes[1, 0]
    ax3.plot(t_hours, sol.y[10], 'g-', lw=2, label='NADPH')
    ax3.plot(t_hours, sol.y[11], 'b-', lw=2, label='NADP+')
    ax3.plot(t_hours, sol.y[12], 'r-', lw=2, label='NADH')
    ax3.plot(t_hours, sol.y[13], 'c-', lw=2, label='NAD+')
    ax3.set_xlabel('Time (hours)')
    ax3.set_ylabel('Concentration (mM)')
    ax3.set_title('Redox Cofactors (with NADH Regeneration)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Energy Status
    ax4 = axes[1, 1]
    ax4.plot(t_hours, sol.y[7], 'g-', lw=2, label='ATP')
    ax4.plot(t_hours, sol.y[8], 'b-', lw=2, label='ADP')
    ax4.plot(t_hours, sol.y[9], 'r-', lw=2, label='Pi')
    ax4.set_xlabel('Time (hours)')
    ax4.set_ylabel('Concentration (mM)')
    ax4.set_title('Energy Status')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: THF Health Ratio
    ax5 = axes[2, 0]
    thf_ratio = sol.y[1] / (sol.y[2] + 1e-6)  # Avoid division by zero
    ax5.plot(t_hours, thf_ratio, 'purple', lw=2, label='THF/DHF Ratio')
    ax5.axhline(y=4, color='gray', linestyle='--', label='Healthy Threshold')
    ax5.set_xlabel('Time (hours)')
    ax5.set_ylabel('Ratio')
    ax5.set_title('THF Regeneration Health')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: NADH Health Ratio
    ax6 = axes[2, 1]
    nadh_ratio = sol.y[12] / (sol.y[13] + 1e-6)
    ax6.plot(t_hours, nadh_ratio, 'orange', lw=2, label='NADH/NAD+ Ratio')
    ax6.axhline(y=0.3, color='gray', linestyle='--', label='Physiological Range')
    ax6.set_xlabel('Time (hours)')
    ax6.set_ylabel('Ratio')
    ax6.set_title('NADH Regeneration Health')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    st.pyplot(fig)
    
    # ========================================================================
    # 关键指标
    # ========================================================================
    st.subheader("📊 Simulation Results")
    col1, col2, col3, col4 = st.columns(4)
    
    final_gly = sol.y[6, -1]
    formate_used = sol.y[0, 0] - sol.y[0, -1]
    nh3_used = sol.y[14, 0] - sol.y[14, -1]
    
    col1.metric("Final Glycine", f"{final_gly:.4f} mM")
    col2.metric("Formate Consumed", f"{formate_used:.4f} mM")
    col3.metric("Ammonia Consumed", f"{nh3_used:.4f} mM")
    col4.metric("Carbon Yield", f"{(final_gly/formate_used)*100:.1f} %" if formate_used > 0 else "0 %")
    
    # ========================================================================
    # 健康检查
    # ========================================================================
    st.subheader("🔍 System Health Check")
    
    final_thf_ratio = sol.y[1, -1] / (sol.y[2, -1] + 1e-6)
    final_nadh_ratio = sol.y[12, -1] / (sol.y[13, -1] + 1e-6)
    
    if final_thf_ratio > 4:
        st.success(f"✅ THF regeneration is healthy (THF/DHF = {final_thf_ratio:.1f})")
    else:
        st.warning(f"⚠️ THF regeneration insufficient (THF/DHF = {final_thf_ratio:.1f})")
        
    if 0.2 < final_nadh_ratio < 0.5:
        st.success(f"✅ NADH regeneration is healthy (NADH/NAD+ = {final_nadh_ratio:.2f})")
    else:
        st.warning(f"⚠️ NADH regeneration abnormal (NADH/NAD+ = {final_nadh_ratio:.2f})")

else:
    st.info("Adjust parameters in the sidebar and click 'Run Precise Simulation'.")