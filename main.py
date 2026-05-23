import streamlit as st
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

st.set_page_config(layout="wide", page_title="Precise Glycine Synthesis Simulator")
st.sidebar.header("🔬 Precise Kinetic Parameters (from Literature)")

# ============================================================================
# 1. ENZYME PARAMETERS (Fixed from literature, but Vmax adjustable)
# ============================================================================
st.sidebar.subheader("Enzyme Vmax (mM/s)")
vmax_ftfl = st.sidebar.slider("FtfL (EC 6.3.4.3)", 0.1, 5.0, 1.6, 0.1)      # Crowther 2008
vmax_fch = st.sidebar.slider("Fch (EC 3.5.4.9)", 0.1, 5.0, 1.9, 0.1)       # Crowther 2008
vmax_mtda = st.sidebar.slider("MtdA (EC 1.5.1.5)", 0.1, 10.0, 6.8, 0.1)    # Crowther 2008
vmax_gcs = st.sidebar.slider("GCS System", 0.1, 50.0, 20.0, 0.1)           # Estimated

# ============================================================================
# 2. INITIAL CONCENTRATIONS (All adjustable)
# ============================================================================
st.sidebar.subheader("Initial Substrates (mM)")
init_formate = st.sidebar.number_input("Formate", 0.0, 200.0, 20.0)
init_thf = st.sidebar.number_input("THF (Taxi)", 0.0, 2.0, 0.2)            # Critical!
init_nh3 = st.sidebar.number_input("Ammonia", 0.0, 200.0, 100.0)           # Must be high
init_co2 = st.sidebar.number_input("CO2", 0.0, 50.0, 5.0)

st.sidebar.subheader("Initial Cofactors (mM)")
init_atp = st.sidebar.number_input("ATP", 0.0, 5.0, 2.0)
init_adp = st.sidebar.number_input("ADP", 0.0, 5.0, 0.1)
init_pi = st.sidebar.number_input("Phosphate", 0.0, 5.0, 1.0)
init_nadph = st.sidebar.number_input("NADPH", 0.0, 5.0, 0.5)
init_nadp = st.sidebar.number_input("NADP+", 0.0, 5.0, 0.05)
init_nadh = st.sidebar.number_input("NADH", 0.0, 5.0, 0.02)
init_nad = st.sidebar.number_input("NAD+", 0.0, 5.0, 0.2)

# ============================================================================
# 3. FIXED KINETIC CONSTANTS (From BRENDA/Literature)
# ============================================================================
Km_FtfL = {"F": 22.0, "THF": 0.8, "ATP": 0.021}
Km_Fch = {"F10": 0.08, "MH4F": 0.08}
Km_MtdA = {"MH4F": 0.03, "NADPH": 0.01, "mTHF": 0.03, "NADP": 0.01}
Km_GCS = {"mTHF": 0.0677, "NH3": 65.4, "CO2": 3.4, "NADH": 0.02}

# ============================================================================
# 4. RATE EQUATIONS (Exact, no simplification)
# ============================================================================
def rate_ftfl(F, THF, ATP):
    p = Km_FtfL
    num = vmax_ftfl * F * THF * ATP
    den = (p["F"]*p["THF"]*p["ATP"] + p["THF"]*p["ATP"]*F + p["F"]*p["ATP"]*THF + p["F"]*p["THF"]*ATP + F*THF*ATP)
    return num / den if den > 0 else 0

def rate_fch(F10, MH4F):
    p = Km_Fch
    v_rev = vmax_fch / 0.54  # Keq = 0.54
    num = (vmax_fch * F10 / p["F10"]) - (v_rev * MH4F / p["MH4F"])
    den = 1 + F10/p["F10"] + MH4F/p["MH4F"]
    return num / den if den > 0 else 0

def rate_mtda(MH4F, NADPH, mTHF, NADP):
    p = Km_MtdA
    v_rev = vmax_mtda / 4.0  # Keq = 4.0
    num = (vmax_mtda * MH4F * NADPH / (p["MH4F"]*p["NADPH"])) - (v_rev * mTHF * NADP / (p["mTHF"]*p["NADP"]))
    den = (1 + MH4F/p["MH4F"] + NADPH/p["NADPH"] + mTHF/p["mTHF"] + NADP/p["NADP"] + MH4F*NADPH/(p["MH4F"]*p["NADPH"]))
    return num / den if den > 0 else 0

def rate_gcs(mTHF, NH3, CO2, NADH):
    p = Km_GCS
    # Four-substrate saturation
    sat_mTHF = mTHF / (p["mTHF"] + mTHF)
    sat_NH3 = NH3 / (p["NH3"] + NH3)
    sat_CO2 = CO2 / (p["CO2"] + CO2)
    sat_NADH = NADH / (p["NADH"] + NADH)
    return vmax_gcs * sat_mTHF * sat_NH3 * sat_CO2 * sat_NADH

# ============================================================================
# 5. ODE SYSTEM
# ============================================================================
def system(t, y):
    F, THF, F10, MH4F, mTHF, Gly, ATP, ADP, Pi, NADPH, NADP, NADH, NAD, NH3, CO2 = y
    
    v1 = rate_ftfl(F, THF, ATP)
    v2 = rate_fch(F10, MH4F)
    v3 = rate_mtda(MH4F, NADPH, mTHF, NADP)
    v4 = rate_gcs(mTHF, NH3, CO2, NADH)
    
    dF = -v1
    dTHF = -v1 + v4   # THF循环：消耗于v1，再生于v4
    dF10 = v1 - v2
    dMH4F = v2 - v3
    dmTHF = v3 - v4
    dGly = v4
    dATP = -v1
    dADP = v1
    dPi = v1
    dNADPH = -v3
    dNADP = v3
    dNADH = -v4
    dNAD = v4
    dNH3 = -v4
    dCO2 = -v4
    
    return [dF, dTHF, dF10, dMH4F, dmTHF, dGly, dATP, dADP, dPi, dNADPH, dNADP, dNADH, dNAD, dNH3, dCO2]

# ============================================================================
# 6. SIMULATION
# ============================================================================
st.sidebar.subheader("Simulation Settings")
sim_hours = st.sidebar.slider("Simulation Time (hours)", 0.5, 10.0, 2.0)

if st.sidebar.button("🚀 Run Precise Simulation", type="primary"):
    # Initial conditions
    y0 = [init_formate, init_thf, 0, 0, 0, 0, init_atp, init_adp, init_pi, 
          init_nadph, init_nadp, init_nadh, init_nad, init_nh3, init_co2]
    
    t_end = int(sim_hours * 3600)
    t_eval = np.linspace(0, t_end, 2000)  # More points = smoother curve
    
    with st.spinner('Solving precise ODE system...'):
        sol = solve_ivp(system, [0, t_end], y0, method='LSODA', t_eval=t_eval, rtol=1e-8, atol=1e-10)
    
    st.success(f'Simulation completed! Time: {sol.t[-1]/3600:.2f} hours')
    
    # Convert time to hours for plotting
    t_hours = sol.t / 3600
    
    # ========================================================================
    # 7. PLOTTING (English labels, clean style)
    # ========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Main Reaction
    ax1 = axes[0, 0]
    ax1.plot(t_hours, sol.y[0], 'g-', lw=2, label='Formate')
    ax1.plot(t_hours, sol.y[1], 'b--', lw=2, label='THF (Taxi)')
    ax1.plot(t_hours, sol.y[13], 'c-', lw=2, label='Ammonia')
    ax1.plot(t_hours, sol.y[5], 'r-', lw=3, label='Glycine')
    ax1.set_xlabel('Time (hours)')
    ax1.set_ylabel('Concentration (mM)')
    ax1.set_title('Main Reaction: Formate + NH3 → Glycine')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Folate Cycle
    ax2 = axes[0, 1]
    ax2.plot(t_hours, sol.y[1], 'c-', lw=2, label='THF')
    ax2.plot(t_hours, sol.y[2], 'm-', lw=2, label='10-formyl-THF')
    ax2.plot(t_hours, sol.y[3], 'y-', lw=2, label='Methenyl-THF')
    ax2.plot(t_hours, sol.y[4], 'k-', lw=2, label='Methylene-THF')
    ax2.set_xlabel('Time (hours)')
    ax2.set_ylabel('Concentration (mM)')
    ax2.set_title('Folate Cycle (Taxi Route)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Cofactors
    ax3 = axes[1, 0]
    ax3.plot(t_hours, sol.y[9], 'g-', lw=2, label='NADPH')
    ax3.plot(t_hours, sol.y[10], 'b-', lw=2, label='NADP+')
    ax3.plot(t_hours, sol.y[11], 'r-', lw=2, label='NADH')
    ax3.plot(t_hours, sol.y[12], 'c-', lw=2, label='NAD+')
    ax3.set_xlabel('Time (hours)')
    ax3.set_ylabel('Concentration (mM)')
    ax3.set_title('Redox Cofactors')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Energy
    ax4 = axes[1, 1]
    ax4.plot(t_hours, sol.y[6], 'g-', lw=2, label='ATP')
    ax4.plot(t_hours, sol.y[7], 'b-', lw=2, label='ADP')
    ax4.plot(t_hours, sol.y[8], 'r-', lw=2, label='Pi')
    ax4.set_xlabel('Time (hours)')
    ax4.set_ylabel('Concentration (mM)')
    ax4.set_title('Energy Status')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    st.pyplot(fig)
    
    # ========================================================================
    # 8. KEY METRICS
    # ========================================================================
    st.subheader("📊 Simulation Results")
    col1, col2, col3, col4 = st.columns(4)
    
    final_gly = sol.y[5, -1]
    formate_used = sol.y[0, 0] - sol.y[0, -1]
    nh3_used = sol.y[13, 0] - sol.y[13, -1]
    
    col1.metric("Final Glycine", f"{final_gly:.4f} mM")
    col2.metric("Formate Consumed", f"{formate_used:.4f} mM")
    col3.metric("Ammonia Consumed", f"{nh3_used:.4f} mM")
    col4.metric("Carbon Yield", f"{(final_gly/formate_used)*100:.1f} %" if formate_used > 0 else "0 %")
    
    # ========================================================================
    # 9. THF STATUS CHECK
    # ========================================================================
    st.subheader("🔍 THF Status Check")
    thf_initial = sol.y[1, 0]
    thf_final = sol.y[1, -1]
    thf_drop = (thf_initial - thf_final) / thf_initial * 100
    
    if thf_drop < 10:
        st.success(f"✅ THF cycle is healthy. THF dropped only {thf_drop:.1f}%.")
    elif thf_drop < 50:
        st.warning(f"⚠️ THF partially depleted ({thf_drop:.1f}%). Consider increasing initial THF.")
    else:
        st.error(f"🚨 THF severely depleted ({thf_drop:.1f}%). First step is bottlenecked!")
        
    # Check if ammonia is limiting
    if init_nh3 < Km_GCS["NH3"]:
        st.warning(f"⚠️ Ammonia ({init_nh3} mM) is below Km ({Km_GCS['NH3']} mM). Increase ammonia!")

else:
    st.info("Adjust parameters in the sidebar and click 'Run Precise Simulation'.")