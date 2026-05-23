import streamlit as st
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

st.set_page_config(layout="wide", page_title="In vivo Glycine Synthesis Simulator")
st.sidebar.header("🔬 Kinetic Parameters (Literature Fixed)")

# ============================================================================
# 1. ENZYME Vmax (Literature based)
# ============================================================================
st.sidebar.subheader("Enzyme Vmax (mM/s)")
vmax_ftfl = st.sidebar.slider("FtfL (EC 6.3.4.3)", 0.1, 5.0, 1.6, 0.1)
vmax_fch = st.sidebar.slider("Fch (EC 3.5.4.9)", 0.1, 5.0, 1.9, 0.1)
vmax_mtda = st.sidebar.slider("MtdA (EC 1.5.1.5)", 0.1, 10.0, 6.8, 0.1)
vmax_gcs = st.sidebar.slider("GCS System", 0.1, 50.0, 20.0, 0.1)

# ============================================================================
# 2. INITIAL CONCENTRATIONS
# ============================================================================
st.sidebar.subheader("Initial Substrates (mM)")
init_formate = st.sidebar.number_input("Formate", 0.0, 200.0, 20.0)
init_thf = st.sidebar.number_input("THF", 0.0, 2.0, 1.0)
init_nh3 = st.sidebar.number_input("Ammonia", 0.0, 200.0, 100.0)
init_co2 = st.sidebar.number_input("CO2", 0.0, 50.0, 5.0)

st.sidebar.subheader("Cofactor Pools (Constant)")
init_atp = st.sidebar.number_input("ATP", 0.0, 5.0, 5.0)
init_nadph = st.sidebar.number_input("NADPH", 0.0, 5.0, 1.0)
init_nadh = st.sidebar.number_input("NADH", 0.0, 5.0, 0.5)

# ============================================================================
# 3. FIXED KINETIC CONSTANTS (BRENDA)
# ============================================================================
Km_FtfL = {"F": 22.0, "THF": 0.8, "ATP": 0.021}
Km_Fch = {"F10": 0.08, "MH4F": 0.08}
Km_MtdA = {"MH4F": 0.03, "NADPH": 0.01, "mTHF": 0.03, "NADP": 0.01}
Km_GCS = {"mTHF": 0.0677, "NH3": 65.4, "CO2": 3.4, "NADH": 0.02}

# ============================================================================
# 4. RATE EQUATIONS (Unchanged)
# ============================================================================
def rate_ftfl(F, THF, ATP):
    p = Km_FtfL
    num = vmax_ftfl * F * THF * ATP
    den = (p["F"]*p["THF"]*p["ATP"] +
           p["THF"]*p["ATP"]*F +
           p["F"]*p["ATP"]*THF +
           p["F"]*p["THF"]*ATP +
           F*THF*ATP)
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

def rate_gcs(mTHF, NH3, CO2, NADH):
    p = Km_GCS
    return vmax_gcs * (
        mTHF / (p["mTHF"] + mTHF) *
        NH3 / (p["NH3"] + NH3) *
        CO2 / (p["CO2"] + CO2) *
        NADH / (p["NADH"] + NADH)
    )

# ============================================================================
# 5. ODE SYSTEM (Stabilized)
# ============================================================================
def system(t, y):
    F, THF, F10, MH4F, mTHF, Gly, NH3, CO2 = y

    # Constant cofactor pools
    ATP = init_atp
    NADPH = init_nadph
    NADH = init_nadh
    NADP = 0.05

    # Gentle THF/mTHF recycle (prevents overshoot)
    k_shmt = 0.05  # Reduced from 0.5
    v_shmt = k_shmt * (mTHF - THF)

    v1 = rate_ftfl(F, THF, ATP)
    v2 = rate_fch(F10, MH4F)
    v3 = rate_mtda(MH4F, NADPH, mTHF, NADP)
    v4 = rate_gcs(mTHF, NH3, CO2, NADH)

    dF = -v1
    dTHF = -v1 + v4 + v_shmt
    dF10 = v1 - v2
    dMH4F = v2 - v3
    dmTHF = v3 - v4 - v_shmt
    dGly = v4
    dNH3 = -v4
    dCO2 = -v4

    return [dF, dTHF, dF10, dMH4F, dmTHF, dGly, dNH3, dCO2]

# ============================================================================
# 6. SIMULATION
# ============================================================================
st.sidebar.subheader("Simulation Settings")
sim_hours = st.sidebar.slider("Time (hours)", 0.5, 10.0, 3.0)
use_log_scale = st.sidebar.checkbox("Log Scale for Folate", False)

if st.sidebar.button("🚀 Run Simulation", type="primary"):
    # Smooth initialization (tiny non-zero for intermediates)
    y0 = [
        init_formate,
        init_thf,
        1e-6,   # F10 (10-formyl-THF)
        1e-6,   # MH4F (5,10-methenyl-THF)
        1e-6,   # mTHF (5,10-methylene-THF)
        0.0,
        init_nh3,
        init_co2
    ]
    
    t_end = int(sim_hours * 3600)
    t_eval = np.linspace(0, t_end, 2000)

    with st.spinner('Solving ODE system...'):
        sol = solve_ivp(
            system,
            [0, t_end],
            y0,
            method="LSODA",
            t_eval=t_eval,
            rtol=1e-9,
            atol=1e-12
        )

    t = sol.t / 3600
    F, THF, F10, MH4F, mTHF, Gly, NH3, CO2 = sol.y

    # ========================================================================
    # 7. PLOTS
    # ========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Main reaction
    axes[0,0].plot(t, F, 'g-', lw=2, label='Formate')
    axes[0,0].plot(t, NH3, 'c-', lw=2, label='Ammonia')
    axes[0,0].plot(t, Gly, 'r-', lw=3, label='Glycine')
    axes[0,0].set_xlabel('Time (hours)')
    axes[0,0].set_ylabel('Concentration (mM)')
    axes[0,0].set_title('Main Reaction: Formate + NH₃ → Glycine')
    axes[0,0].legend()
    axes[0,0].grid(True, alpha=0.3)
    
    # Folate cycle
    axes[0,1].plot(t, THF, 'c-', lw=2, label='THF')
    axes[0,1].plot(t, F10, 'm-', lw=2, label='10-F-THF')
    axes[0,1].plot(t, MH4F, 'y-', lw=2, label='5,10-M-THF')
    axes[0,1].plot(t, mTHF, 'k-', lw=2, label='5,10-m-THF')
    axes[0,1].set_xlabel('Time (hours)')
    axes[0,1].set_ylabel('Concentration (mM)')
    axes[0,1].set_title('Folate Cycle (Steady-State Behavior)')
    axes[0,1].legend()
    axes[0,1].grid(True, alpha=0.3)
    
    if use_log_scale:
        axes[0,1].set_yscale('log')
        axes[0,1].set_ylim(1e-5, 2)
    
    # Glycine production
    axes[1,0].plot(t, Gly, 'r-', lw=3)
    axes[1,0].set_xlabel('Time (hours)')
    axes[1,0].set_ylabel('Glycine (mM)')
    axes[1,0].set_title('Glycine Accumulation')
    axes[1,0].grid(True, alpha=0.3)
    
    # Formate consumption
    axes[1,1].plot(t, F, 'g-', lw=2)
    axes[1,1].set_xlabel('Time (hours)')
    axes[1,1].set_ylabel('Formate (mM)')
    axes[1,1].set_title('Formate Utilization')
    axes[1,1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    st.pyplot(fig)
    
    # ========================================================================
    # 8. METRICS
    # ========================================================================
    st.subheader("📊 Simulation Results")
    col1, col2, col3 = st.columns(3)
    
    final_gly = Gly[-1]
    formate_used = F[0] - F[-1]
    carbon_yield = (final_gly / formate_used * 100) if formate_used > 0 else 0
    
    col1.metric("Final Glycine", f"{final_gly:.3f} mM")
    col2.metric("Formate Used", f"{formate_used:.3f} mM")
    col3.metric("Carbon Yield", f"{carbon_yield:.1f} %")
    
    # ========================================================================
    # 9. DIAGNOSTICS
    # ========================================================================
    st.subheader("🔍 System Diagnostics")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("**Folate Pool Status**")
        st.write(f"- THF final: {THF[-1]:.4f} mM")
        st.write(f"- 10-F-THF: {F10[-1]:.6f} mM")
        st.write(f"- 5,10-M-THF: {MH4F[-1]:.6f} mM")
        st.write(f"- 5,10-m-THF: {mTHF[-1]:.4f} mM")
        
    with col2:
        st.write("**Physiological Interpretation**")
        if F10[-1] < 1e-4 and MH4F[-1] < 1e-4:
            st.success("✅ Fast equilibrium: 10-F-THF and 5,10-M-THF are transient intermediates (expected)")
        if THF[-1] > 0.1:
            st.success("✅ THF pool maintained")
        if mTHF[-1] > THF[-1]:
            st.info("ℹ️ mTHF accumulates as the primary folate carrier")
    
    # Warning checks
    if init_nh3 < Km_GCS["NH3"]:
        st.warning(f"⚠️ Ammonia ({init_nh3} mM) below Km ({Km_GCS['NH3']} mM)")
    
    if THF[-1] < 0.01:
        st.warning("⚠️ THF pool nearly depleted - consider increasing initial THF")

else:
    st.info("Adjust parameters and click 'Run Simulation'.")