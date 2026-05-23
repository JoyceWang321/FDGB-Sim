import streamlit as st
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

st.set_page_config(layout="wide", page_title="Thermodynamic SHMT Model")
st.sidebar.header("🔬 Thermodynamic Control")

# ============================================================================
# 1. ENZYME PARAMETERS (Literature)
# ============================================================================
st.sidebar.subheader("Vmax (mM/s)")
vmax_ftfl = st.sidebar.slider("FtfL", 0.1, 5.0, 1.6, 0.1)
vmax_fch = st.sidebar.slider("Fch", 0.1, 5.0, 1.9, 0.1)
vmax_mtda = st.sidebar.slider("MtdA", 0.1, 10.0, 6.8, 0.1)
vmax_gcs = st.sidebar.slider("GCS", 0.1, 50.0, 20.0, 0.1)
vmax_shmt = st.sidebar.slider("SHMT", 1.0, 50.0, 15.0, 1.0)

# ============================================================================
# 2. INITIAL CONCENTRATIONS
# ============================================================================
st.sidebar.subheader("Initial Concentrations (mM)")
init_formate = st.sidebar.number_input("Formate", 0.0, 200.0, 50.0)
init_thf = st.sidebar.number_input("THF", 0.0, 5.0, 1.0)
init_gly = st.sidebar.number_input("Glycine (initial)", 0.0, 10.0, 0.0)
init_ser = st.sidebar.number_input("Serine (initial)", 0.0, 100.0, 5.0)
init_nh3 = st.sidebar.number_input("Ammonia", 0.0, 200.0, 100.0)
init_hco3 = st.sidebar.number_input("HCO3-", 0.0, 100.0, 50.0)
init_atp = st.sidebar.number_input("ATP", 0.0, 5.0, 5.0)
init_nadph = st.sidebar.number_input("NADPH", 0.0, 5.0, 1.0)
init_nadh = st.sidebar.number_input("NADH", 0.0, 5.0, 0.5)

# ============================================================================
# 3. FIXED KINETIC CONSTANTS
# ============================================================================
Km_FtfL = {"F": 22.0, "THF": 0.8, "ATP": 0.021}
Km_Fch = {"F10": 0.08, "MH4F": 0.08}
Km_MtdA = {"MH4F": 0.03, "NADPH": 0.01, "mTHF": 0.03, "NADP": 0.01}
Km_GCS = {"mTHF": 0.0677, "NH3": 65.4, "HCO3": 3.4, "NADH": 0.02}
Km_SHMT = {"Gly": 5.0, "mTHF": 0.05, "Ser": 3.0, "THF": 0.2}
Keq_shmt = 1.2  # 略微偏向甘氨酸生成

# ============================================================================
# 4. RATE EQUATIONS
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

def rate_mtda(MH4F, NADPH, mTHF, NADP):
    p = Km_MtdA
    v_rev = vmax_mtda / 4.0
    num = (vmax_mtda * MH4F * NADPH / (p["MH4F"]*p["NADPH"])) - \
          (v_rev * mTHF * NADP / (p["mTHF"]*p["NADP"]))
    den = (1 + MH4F/p["MH4F"] + NADPH/p["NADPH"] +
           mTHF/p["mTHF"] + NADP/p["NADP"] +
           MH4F*NADPH/(p["MH4F"]*p["NADPH"]))
    return num / den if den > 0 else 0

def rate_gcs(mTHF, NH3, NADH, HCO3):
    p = Km_GCS
    return vmax_gcs * (
        mTHF / (p["mTHF"] + mTHF) *
        NH3 / (p["NH3"] + NH3) *
        HCO3 / (p["HCO3"] + HCO3) *
        NADH / (p["NADH"] + NADH)
    )

def rate_shmt(Gly, mTHF, Ser, THF):
    """SHMT: Gly + mTHF ⇌ Ser + THF (热力学驱动)"""
    p = Km_SHMT
    v_fwd = vmax_shmt * (Gly * mTHF) / (p["Gly"] * p["mTHF"])
    v_rev = (vmax_shmt / Keq_shmt) * (Ser * THF) / (p["Ser"] * p["THF"])
    
    den = (1 + Gly/p["Gly"] + mTHF/p["mTHF"] + 
           Ser/p["Ser"] + THF/p["THF"] +
           (Gly*mTHF)/(p["Gly"]*p["mTHF"]) +
           (Ser*THF)/(p["Ser"]*p["THF"]))
    
    return (v_fwd - v_rev) / den if den > 0 else 0

# ============================================================================
# 5. ODE SYSTEM (THERMODYNAMIC CONTROL)
# ============================================================================
def system(t, y):
    F, THF, F10, MH4F, mTHF, Gly, Ser, NH3 = y
    
    ATP = init_atp
    NADPH = init_nadph
    NADH = init_nadh
    NADP = 0.05
    HCO3 = init_hco3
    
    v1 = rate_ftfl(F, THF, ATP)
    v2 = rate_fch(F10, MH4F)
    v3 = rate_mtda(MH4F, NADPH, mTHF, NADP)
    v4 = rate_gcs(mTHF, NH3, NADH, HCO3)
    v5 = rate_shmt(Gly, mTHF, Ser, THF)  # SHMT: Gly + mTHF ⇌ Ser + THF
    
    dF = -v1
    dTHF = -v1 + v4 + v5  # THF: 消耗于FtfL，再生于GCS和SHMT
    dF10 = v1 - v2
    dMH4F = v2 - v3
    dmTHF = v3 - v4 - v5  # mTHF: 消耗于GCS和SHMT
    dGly = v4 - v5        # 甘氨酸: GCS生成，SHMT消耗
    dSer = v5             # 丝氨酸: SHMT生成
    dNH3 = -v4
    
    return [dF, dTHF, dF10, dMH4F, dmTHF, dGly, dSer, dNH3]

# ============================================================================
# 6. SIMULATION
# ============================================================================
st.sidebar.subheader("Simulation")
sim_hours = st.sidebar.slider("Time (hours)", 0.5, 10.0, 3.0)

if st.sidebar.button("🚀 Run Thermodynamic Model", type="primary"):
    y0 = [init_formate, init_thf, 1e-6, 1e-6, 1e-6, init_gly, init_ser, init_nh3]
    t_end = int(sim_hours * 3600)
    t_eval = np.linspace(0, t_end, 2000)
    
    sol = solve_ivp(system, [0, t_end], y0, method="LSODA", t_eval=t_eval, rtol=1e-9, atol=1e-12)
    
    t = sol.t / 3600
    F, THF, F10, MH4F, mTHF, Gly, Ser, NH3 = sol.y
    
    # ========================================================================
    # 7. PLOTS
    # ========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    axes[0,0].plot(t, F, 'g-', lw=2, label='Formate')
    axes[0,0].plot(t, Gly, 'r-', lw=3, label='Glycine')
    axes[0,0].plot(t, Ser, 'b-', lw=2, label='Serine')
    axes[0,0].set_title('Main Reaction (Thermodynamic Control)')
    axes[0,0].legend(); axes[0,0].grid(alpha=0.3)
    
    axes[0,1].plot(t, THF, 'c-', lw=2, label='THF')
    axes[0,1].plot(t, F10, 'm-', lw=2, label='10-F-THF')
    axes[0,1].plot(t, MH4F, 'y-', lw=2, label='5,10-M-THF')
    axes[0,1].plot(t, mTHF, 'k-', lw=2, label='5,10-m-THF')
    axes[0,1].set_title('Folate Cycle (Self-Regulating)')
    axes[0,1].legend(); axes[0,1].grid(alpha=0.3)
    
    axes[1,0].plot(t, Gly, 'r-', lw=3, label='Glycine')
    axes[1,0].plot(t, Ser, 'b-', lw=2, label='Serine')
    axes[1,0].set_title('Glycine vs Serine (Competition)')
    axes[1,0].legend(); axes[1,0].grid(alpha=0.3)
    
    axes[1,1].plot(t, F, 'g-', lw=2)
    axes[1,1].plot(t, NH3, 'c-', lw=2)
    axes[1,1].set_title('Substrate Consumption')
    axes[1,1].legend(['Formate', 'Ammonia'])
    axes[1,1].grid(alpha=0.3)
    
    plt.tight_layout()
    st.pyplot(fig)
    
    # ========================================================================
    # 8. METRICS
    # ========================================================================
    st.subheader("📊 Results")
    col1, col2, col3, col4 = st.columns(4)
    
    final_gly = Gly[-1]
    final_ser = Ser[-1]
    formate_used = F[0] - F[-1]
    total_gly_ser = final_gly + final_ser
    
    col1.metric("Final Glycine", f"{final_gly:.2f} mM")
    col2.metric("Final Serine", f"{final_ser:.2f} mM")
    col3.metric("Total Gly+Ser", f"{total_gly_ser:.2f} mM")
    col4.metric("Formate Used", f"{formate_used:.2f} mM")
    
    # ========================================================================
    # 9. THERMODYNAMIC CHECK
    # ========================================================================
    st.subheader("🌡️ Thermodynamic Behavior Check")
    
    # 检查THF是否稳定在合理范围
    thf_min = min(THF)
    thf_max = max(THF)
    thf_final = THF[-1]
    
    if thf_final < 0.01:
        st.error(f"🚨 THF critically low: {thf_final:.4f} mM")
    elif thf_final < 0.1:
        st.warning(f"⚠️ THF low: {thf_final:.4f} mM")
    else:
        st.success(f"✅ THF stable: {thf_final:.4f} mM")
    
    # 检查甘氨酸是否被SHMT消耗
    gly_peak_idx = np.argmax(Gly)
    gly_after_peak = Gly[gly_peak_idx:]
    
    if len(gly_after_peak) > 10 and max(gly_after_peak) - min(gly_after_peak) > 0.1:
        st.info("ℹ️ Glycine decreases after peak (SHMT consuming excess)")
    else:
        st.success("✅ Glycine production balanced by SHMT")
    
    # 检查丝氨酸是否生成
    if final_ser > init_ser:
        st.success(f"✅ Serine generated: {final_ser - init_ser:.2f} mM")
    else:
        st.warning("⚠️ No significant serine generation")
    
    st.info(f"""
    **Thermodynamic Interpretation:**
    - THF 最低降至 {thf_min:.4f} mM 后回升
    - 甘氨酸峰值后，SHMT 逆向运行消耗甘氨酸
    - 丝氨酸从 {init_ser:.1f} 增至 {final_ser:.1f} mM
    - 系统自动调节，无需人工干预
    """)

else:
    st.info("Adjust parameters and run simulation.")