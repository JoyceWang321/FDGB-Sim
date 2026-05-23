import streamlit as st
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

st.set_page_config(layout="wide", page_title="Stable Glycine Model")
st.sidebar.header("🔬 Stable Parameters")

# ============================================================================
# 1. DEFAULT VALUES
# ============================================================================
DEFAULT_VMAX = {
    "FtfL": 1.6,
    "Fch": 1.9,
    "MtdA": 6.8,
    "GCS": 0.10,
    "SHMT": 15.0
}

DEFAULT_KM = {
    "FtfL_F": 22.0, "FtfL_THF": 0.8, "FtfL_ATP": 0.021,
    "Fch_F10": 0.08, "Fch_MH4F": 0.08,
    "MtdA_MH4F": 0.03, "MtdA_NADPH": 0.01, "MtdA_mTHF": 0.03, "MtdA_NADP": 0.01,
    "GCS_mTHF": 0.0677, "GCS_NH3": 65.4, "GCS_HCO3": 3.4, "GCS_NADH": 0.02, "GCS_Gly_inh": 0.1,
    "SHMT_Gly": 5.0, "SHMT_mTHF": 0.05, "SHMT_Ser": 3.0, "SHMT_THF": 0.2
}

DEFAULT_INIT = {
    "Formate": 50.0,
    "THF": 1.0,
    "Glycine": 0.0,
    "Serine": 5.0,
    "Ammonia": 100.0,
    "HCO3": 50.0,
    "ATP": 5.0,
    "NADPH": 1.0,
    "NADH": 0.5
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
vmax_ftfl = st.sidebar.number_input("FtfL", 0.1, 10.0, DEFAULT_VMAX["FtfL"], 0.1, key="vmax_ftfl")
vmax_fch = st.sidebar.number_input("Fch", 0.1, 10.0, DEFAULT_VMAX["Fch"], 0.1, key="vmax_fch")
vmax_mtda = st.sidebar.number_input("MtdA", 0.1, 20.0, DEFAULT_VMAX["MtdA"], 0.1, key="vmax_mtda")
vmax_gcs = st.sidebar.number_input("GCS", 0.001, 1.0, DEFAULT_VMAX["GCS"], 0.001, key="vmax_gcs")
vmax_shmt = st.sidebar.number_input("SHMT", 0.1, 50.0, DEFAULT_VMAX["SHMT"], 0.1, key="vmax_shmt")

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
          "HCO3": DEFAULT_KM["GCS_HCO3"], "NADH": DEFAULT_KM["GCS_NADH"]}
Km_SHMT = {"Gly": DEFAULT_KM["SHMT_Gly"], "mTHF": DEFAULT_KM["SHMT_mTHF"], 
           "Ser": DEFAULT_KM["SHMT_Ser"], "THF": DEFAULT_KM["SHMT_THF"]}
Keq_shmt = 1.2
Ki_gly_gcs = DEFAULT_KM["GCS_Gly_inh"]

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
    if F10 < 1e-8 and MH4F < 1e-8:
        return 0
    p = Km_Fch
    v_rev = vmax_fch / 0.54
    num = (vmax_fch * F10 / p["F10"]) - (v_rev * MH4F / p["MH4F"])
    den = 1 + F10/p["F10"] + MH4F/p["MH4F"]
    return safe_divide(num, den)

def rate_mtda(MH4F, NADPH, mTHF, NADP):
    if MH4F < 1e-8 and mTHF < 1e-8:
        return 0
    p = Km_MtdA
    v_rev = vmax_mtda / 4.0
    num = (vmax_mtda * MH4F * NADPH / (p["MH4F"]*p["NADPH"])) - \
          (v_rev * mTHF * NADP / (p["mTHF"]*p["NADP"]))
    den = (1 + MH4F/p["MH4F"] + NADPH/p["NADPH"] +
           mTHF/p["mTHF"] + NADP/p["NADP"] +
           MH4F*NADPH/(p["MH4F"]*p["NADPH"]))
    return safe_divide(num, den)

def rate_gcs(mTHF, NH3, NADH, HCO3, Gly):
    if mTHF < 1e-8 or NH3 < 1e-8 or HCO3 < 1e-8 or NADH < 1e-8:
        return 0
    p = Km_GCS
    inhibition = 1 / (1 + Gly/Ki_gly_gcs)
    base_rate = vmax_gcs * (
        mTHF / (p["mTHF"] + mTHF) *
        NH3 / (p["NH3"] + NH3) *
        HCO3 / (p["HCO3"] + HCO3) *
        NADH / (p["NADH"] + NADH)
    )
    return base_rate * inhibition

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
# 7. ODE SYSTEM (WITH NON-NEGATIVITY ENFORCEMENT)
# ============================================================================
def system(t, y):
    # 强制所有浓度非负
    y = np.maximum(y, 0)
    
    F, THF, F10, MH4F, mTHF, Gly, Ser, NH3 = y
    
    ATP = init_atp
    NADPH = init_nadph
    NADH = init_nadh
    NADP = 0.05
    HCO3 = init_hco3
    
    v1 = rate_ftfl(F, THF, ATP)
    v2 = rate_fch(F10, MH4F)
    v3 = rate_mtda(MH4F, NADPH, mTHF, NADP)
    v4 = rate_gcs(mTHF, NH3, NADH, HCO3, Gly)
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
sim_minutes = st.sidebar.slider("Simulation Time (minutes)", 1, 180, 120)  # 扩展到180分钟

if st.sidebar.button("🚀 Run Stable Simulation", type="primary"):
    y0 = [init_formate, init_thf, 1e-6, 1e-6, 1e-6, init_gly, init_ser, init_nh3]
    t_end = sim_minutes * 60
    t_eval = np.linspace(0, t_end, 3000)  # 增加采样点
    
    try:
        with st.spinner(f'Simulating {sim_minutes} minutes...'):
            sol = solve_ivp(
                system, 
                [0, t_end], 
                y0, 
                method="LSODA",
                t_eval=t_eval,
                rtol=1e-8, 
                atol=1e-10,
                max_step=30.0,  # 最大步长30秒
                first_step=0.001
            )
        
        if sol.success:
            st.success(f'✅ Simulation completed successfully! Time: {sol.t[-1]/60:.1f} minutes')
            
            t = sol.t / 60
            F, THF, F10, MH4F, mTHF, Gly, Ser, NH3 = sol.y
            
            # 检查是否有负值
            min_vals = [np.min(F), np.min(THF), np.min(Gly), np.min(Ser)]
            if any(val < -1e-6 for val in min_vals):
                st.warning("⚠️ Negative concentrations detected (numerical artifact)")
            
            # ========================================================================
            # 9. PLOTS
            # ========================================================================
            fig, axes = plt.subplots(2, 2, figsize=(14, 10))
            
            axes[0,0].plot(t, F, 'g-', lw=2, label='Formate')
            axes[0,0].plot(t, Gly, 'r-', lw=3, label='Glycine')
            axes[0,0].plot(t, Ser, 'b-', lw=2, label='Serine')
            axes[0,0].set_xlabel('Time (min)')
            axes[0,0].set_ylabel('Concentration (mM)')
            axes[0,0].set_title(f'Stable Simulation ({sim_minutes} min)')
            axes[0,0].legend()
            axes[0,0].grid(True, alpha=0.3)
            
            axes[0,1].plot(t, THF, 'c-', lw=2, label='THF')
            axes[0,1].plot(t, mTHF, 'k-', lw=2, label='mTHF')
            axes[0,1].plot(t, F10, 'm-', lw=2, label='10-F-THF')
            axes[0,1].plot(t, MH4F, 'y-', lw=2, label='5,10-M-THF')
            axes[0,1].set_xlabel('Time (min)')
            axes[0,1].set_ylabel('Concentration (mM)')
            axes[0,1].set_title('Folate Cycle (Stable)')
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
            
            # ========================================================================
            # 10. METRICS
            # ========================================================================
            st.subheader("📊 Stable Results")
            col1, col2, col3, col4 = st.columns(4)
            
            final_gly = Gly[-1]
            final_ser = Ser[-1]
            formate_used = F[0] - F[-1]
            
            col1.metric("Final Glycine", f"{final_gly:.4f} mM")
            col2.metric("Final Serine", f"{final_ser:.4f} mM")
            col3.metric("Formate Used", f"{formate_used:.2f} mM")
            col4.metric("GCS Vmax", f"{vmax_gcs:.3f} mM/min")
            
            # 稳定性检查
            st.subheader("🔍 Stability Check")
            if sol.nfev > 10000:
                st.info(f"Integration required {sol.nfev} function evaluations")
            if len(sol.t) < 100:
                st.warning("Few time steps - possible stiffness issues")
            else:
                st.success("Integration stable")
                
        else:
            st.error(f"❌ Simulation failed: {sol.message}")
            
    except Exception as e:
        st.error(f"❌ Error during simulation: {str(e)}")
        st.info("Try reducing simulation time or adjusting parameters.")

else:
    st.info("Adjust parameters and click 'Run Stable Simulation'.")