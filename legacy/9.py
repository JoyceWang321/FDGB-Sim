import streamlit as st
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import time
import pandas as pd

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
    # Clear saved simulation result on reset
    if "simulation_result" in st.session_state:
        del st.session_state.simulation_result

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
# 6. RATE EQUATIONS (PARAMETERIZED)
# ============================================================================
def safe_divide(num, den):
    return num / den if den > 1e-12 else 0

def rate_ftfl(F, THF, ATP, vmax_ftfl_local):
    if F < 1e-8 or THF < 1e-8 or ATP < 1e-8:
        return 0
    p = Km_FtfL
    num = vmax_ftfl_local * F * THF * ATP
    den = (p["F"]*p["THF"]*p["ATP"] + p["THF"]*p["ATP"]*F +
           p["F"]*p["ATP"]*THF + p["F"]*p["THF"]*ATP + F*THF*ATP)
    return safe_divide(num, den)

def rate_fch(F10, MH4F, vmax_fch_local):
    if F10 < 1e-8:
        return 0
    p = Km_Fch
    num = vmax_fch_local * F10 / p["F10"]
    den = 1 + F10/p["F10"] + MH4F/p["MH4F"]
    return safe_divide(num, den)

def rate_mtda(MH4F, NADPH, mTHF, NADP, vmax_mtda_local):
    if MH4F < 1e-8:
        return 0
    p = Km_MtdA
    num = vmax_mtda_local * MH4F * NADPH / (p["MH4F"]*p["NADPH"])
    den = (1 + MH4F/p["MH4F"] + NADPH/p["NADPH"] +
           mTHF/p["mTHF"] + NADP/p["NADP"] +
           MH4F*NADPH/(p["MH4F"]*p["NADPH"]))
    return safe_divide(num, den)

def rate_gcs(mTHF, NH3, NADH, HCO3, Gly, THF, vmax_gcs_local):
    if any(c < 1e-8 for c in [mTHF, NH3, HCO3, NADH, Gly, THF]):
        return 0
    p = Km_GCS
    v_fwd = vmax_gcs_local * (
        mTHF/(p["mTHF"]+mTHF) *
        NH3/(p["NH3"]+NH3) *
        HCO3/(p["HCO3"]+HCO3) *
        NADH/(p["NADH"]+NADH)
    ) * (1/(1+Gly/p["GCS_Gly_inh"]))
    v_rev = vmax_gcs_local * 0.02 * (
        Gly/(p["GCS_Gly_rev"]+Gly) *
        THF/(p["GCS_THF_rev"]+THF)
    )
    return safe_divide(v_fwd - v_rev, 1)

def rate_shmt(Gly, mTHF, Ser, THF, vmax_shmt_local):
    if Gly < 1e-8 and Ser < 1e-8:
        return 0
    p = Km_SHMT
    if THF < 0.01:
        v_fwd = 0
        v_rev = (vmax_shmt_local / Keq_shmt) * (Ser * THF) / (p["Ser"] * p["THF"])
    else:
        v_fwd = vmax_shmt_local * (Gly * mTHF) / (p["Gly"] * p["mTHF"])
        v_rev = (vmax_shmt_local / Keq_shmt) * (Ser * THF) / (p["Ser"] * p["THF"])
    den = (1 + Gly/p["Gly"] + mTHF/p["mTHF"] + 
           Ser/p["Ser"] + THF/p["THF"] +
           (Gly*mTHF)/(p["Gly"]*p["mTHF"]) +
           (Ser*THF)/(p["Ser"]*p["THF"]))
    return safe_divide(v_fwd - v_rev, den)

# ============================================================================
# 7. ODE SYSTEM (PARAMETERIZED)
# ============================================================================
def system(t, y, params):
    y = np.maximum(y, 0)
    F, THF, F10, MH4F, mTHF, Gly, Ser, NH3 = y
    ATP = params['init_atp']
    NADPH = params['init_nadph']
    NADH = params['init_nadh']
    NADP = 0.05
    HCO3 = params['init_hco3']
    if THF < 0.005:
        v1 = 0
    else:
        v1 = rate_ftfl(F, THF, ATP, params['vmax_ftfl'])
    v2 = rate_fch(F10, MH4F, params['vmax_fch'])
    v3 = rate_mtda(MH4F, NADPH, mTHF, NADP, params['vmax_mtda'])
    v4 = rate_gcs(mTHF, NH3, NADH, HCO3, Gly, THF, params['vmax_gcs'])
    v5 = rate_shmt(Gly, mTHF, Ser, THF, params['vmax_shmt'])
    return [-v1, -v1+v4+v5, v1-v2, v2-v3, v3-v4-v5, v4-v5, v5, -v4]

# ============================================================================
# 8. SIMULATION FUNCTIONS (WITH PERSISTENCE)
# ============================================================================
def get_current_params():
    return {
        'vmax_ftfl': st.session_state.vmax_ftfl,
        'vmax_fch': st.session_state.vmax_fch,
        'vmax_mtda': st.session_state.vmax_mtda,
        'vmax_gcs': st.session_state.vmax_gcs,
        'vmax_shmt': st.session_state.vmax_shmt,
        'init_formate': init_formate,
        'init_thf': init_thf,
        'init_gly': init_gly,
        'init_ser': init_ser,
        'init_nh3': init_nh3,
        'init_hco3': init_hco3,
        'init_atp': init_atp,
        'init_nadph': init_nadph,
        'init_nadh': init_nadh,
        'sim_minutes': st.session_state.sim_minutes
    }

def run_simulation_with_params(params):
    mthf_init = max(0.5, params['vmax_gcs'] * 20)
    y0 = [params['init_formate'], params['init_thf'], 1e-6, 1e-6, mthf_init, 
          params['init_gly'], params['init_ser'], params['init_nh3']]
    t_end = params['sim_minutes'] * 60
    t_eval = np.linspace(0, t_end, 3000)
    try:
        sol = solve_ivp(
            lambda t, y: system(t, y, params),
            [0, t_end], y0, method="LSODA", t_eval=t_eval,
            rtol=1e-6, atol=1e-8, max_step=5.0, first_step=0.001, max_num_steps=100000
        )
        if sol.success:
            t = sol.t / 60
            F, THF, F10, MH4F, mTHF, Gly, Ser, NH3 = sol.y
            return {'time': t, 'formate': F, 'thf': THF, 'mthf': mTHF, 'f10': F10, 'mh4f': MH4F,
                    'glycine': Gly, 'serine': Ser, 'ammonia': NH3, 'success': True, 'steps': sol.nfev}
        else:
            return {'success': False, 'message': sol.message}
    except Exception as e:
        return {'success': False, 'message': str(e)}

# ============================================================================
# 9. PLOTTING FUNCTIONS (PERSISTENT)
# ============================================================================
def display_simulation_plots(result):
    t = result['time']
    F = result['formate']
    THF = result['thf']
    mTHF = result['mthf']
    F10 = result['f10']
    MH4F = result['mh4f']
    Gly = result['glycine']
    Ser = result['serine']
    NH3 = result['ammonia']
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes[0,0].plot(t, Gly, 'r-', lw=3, label='Glycine')
    axes[0,0].plot(t, Ser, 'b-', lw=2, label='Serine')
    axes[0,0].fill_between(t, 0, Gly+Ser, alpha=0.1, color='gray', label='Total AA')
    axes[0,0].set_xlabel('Time (min)')
    axes[0,0].set_ylabel('Concentration (mM)')
    axes[0,0].set_title(f'Amino Acid Synthesis (Total: {Gly[-1]+Ser[-1]:.2f} mM)')
    axes[0,0].legend()
    axes[0,0].grid(True, alpha=0.3)
    axes[0,1].plot(t, THF, 'c-', lw=2, label='THF')
    axes[0,1].plot(t, mTHF, 'k-', lw=2, label='mTHF')
    axes[0,1].plot(t, F10, 'm-', lw=2, label='10-F-THF')
    axes[0,1].plot(t, MH4F, 'y-', lw=2, label='5,10-M-THF')
    axes[0,1].axhline(y=0.01, color='r', linestyle='--', alpha=0.5, label='Critical THF')
    axes[0,1].set_xlabel('Time (min)')
    axes[0,1].set_ylabel('Concentration (mM)')
    axes[0,1].set_title('Folate Cycle (THF/mTHF Balance)')
    axes[0,1].legend()
    axes[0,1].grid(True, alpha=0.3)
    axes[1,0].plot(t, F, 'g-', lw=2, label='Formate')
    axes[1,0].plot(t, NH3, 'c-', lw=2, label='Ammonia')
    for substrate, name, color in [(F, 'Formate', 'g'), (NH3, 'NH3', 'c')]:
        idx = np.where(substrate < 0.01)[0]
        if len(idx) > 0:
            axes[1,0].axvline(x=t[idx[0]], color=color, linestyle=':', alpha=0.7, 
                             label=f'{name} depleted @ {t[idx[0]]:.1f}min')
    axes[1,0].set_xlabel('Time (min)')
    axes[1,0].set_ylabel('Concentration (mM)')
    axes[1,0].set_title('Substrate Depletion (Limiting Factors)')
    axes[1,0].legend()
    axes[1,0].grid(True, alpha=0.3)
    v4_vals = [rate_gcs(mTHF[i], NH3[i], init_nadh, init_hco3, Gly[i], THF[i], st.session_state.vmax_gcs) for i in range(len(t))]
    axes[1,1].plot(t, v4_vals, 'r-', lw=2, label='GCS Rate')
    axes[1,1].axhline(y=st.session_state.vmax_gcs, color='gray', linestyle='--', alpha=0.5, label=f'Vmax={st.session_state.vmax_gcs}')
    axes[1,1].set_xlabel('Time (min)')
    axes[1,1].set_ylabel('Rate (mM/min)')
    axes[1,1].set_title('GCS Activity Over Time')
    axes[1,1].legend()
    axes[1,1].grid(True, alpha=0.3)
    plt.tight_layout()
    st.pyplot(fig)

def display_results_summary(result):
    Gly = result['glycine']
    Ser = result['serine']
    THF = result['thf']
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Final Glycine", f"{Gly[-1]:.4f} mM")
    col2.metric("Final Serine", f"{Ser[-1]:.4f} mM")
    col3.metric("Total Amino Acids", f"{Gly[-1]+Ser[-1]:.4f} mM")
    col4.metric("Min THF", f"{np.min(THF):.6f} mM")

# ============================================================================
# 10. SENSITIVITY ANALYSIS (USES PERSISTED SIMULATION DATA)
# ============================================================================
def run_sensitivity_analysis():
    # Check if simulation data exists
    if "simulation_result" not in st.session_state or not st.session_state.simulation_result['success']:
        st.error("Please run a successful simulation first before performing sensitivity analysis.")
        return
    
    st.subheader("📊 Sensitivity Analysis Results")
    base_result = st.session_state.simulation_result
    baseline_gly = base_result['glycine'][-1]
    base_params = get_current_params()

    params_to_test = {
        'vmax_gcs': (0.01, 10.0), 'vmax_ftfl': (0.01, 10.0), 'vmax_mtda': (0.01, 20.0),
        'vmax_shmt': (0.01, 50.0), 'init_formate': (0.0, 200.0), 'init_nh3': (0.0, 200.0),
        'init_hco3': (0.0, 100.0), 'init_nadh': (0.0, 5.0), 'init_thf': (0.0, 5.0)
    }
    
    results = []
    progress_bar = st.progress(0)
    for i, (param_key, (min_val, max_val)) in enumerate(params_to_test.items()):
        current_val = base_params[param_key]
        test_values = [max(min_val, current_val * 0.9), min(max_val, current_val * 1.1)]
        param_results = []
        for test_val in test_values:
            test_params = base_params.copy()
            test_params[param_key] = test_val
            result = run_simulation_with_params(test_params)
            param_results.append(result['glycine'][-1] if result['success'] else baseline_gly)
        delta_param = test_values[1] - test_values[0]
        delta_gly = param_results[1] - param_results[0]
        sensitivity_index = (delta_gly / baseline_gly) / (delta_param / current_val) if delta_param != 0 else 0
        results.append({'Parameter': param_key, 'Current Value': current_val,
                        'Sensitivity Index': sensitivity_index, 'Impact': abs(sensitivity_index)})
        progress_bar.progress((i + 1) / len(params_to_test))
    
    df = pd.DataFrame(results).sort_values('Impact', ascending=False)
    st.dataframe(df.style.format({'Current Value': '{:.3f}', 'Sensitivity Index': '{:.3f}', 'Impact': '{:.3f}'}))
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ['red' if x > 0 else 'blue' for x in df['Sensitivity Index']]
    bars = ax.barh(df['Parameter'], df['Sensitivity Index'], color=colors)
    ax.set_xlabel('Sensitivity Index (∂Gly/∂Param × Param/Gly)')
    ax.set_title(f'Sensitivity Analysis (Baseline Glycine: {baseline_gly:.2f} mM)')
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    for bar, val in zip(bars, df['Sensitivity Index']):
        width = bar.get_width()
        ax.text(width + 0.01 if width >= 0 else width - 0.01, 
                bar.get_y() + bar.get_height()/2, f'{val:.2f}', va='center',
                ha='left' if width >= 0 else 'right')
    plt.tight_layout()
    st.pyplot(fig)
    st.subheader("📝 Interpretation")
    top_pos = df[df['Sensitivity Index'] > 0].head(3)
    top_neg = df[df['Sensitivity Index'] < 0].head(3)
    if not top_pos.empty:
        st.markdown("**🔺 Increasing these parameters increases glycine yield:**")
        for _, row in top_pos.iterrows():
            st.markdown(f"- **{row['Parameter']}**: +1% → +{row['Sensitivity Index']:.2f}% glycine")
    if not top_neg.empty:
        st.markdown("**🔻 Increasing these parameters decreases glycine yield:**")
        for _, row in top_neg.iterrows():
            st.markdown(f"- **{row['Parameter']}**: +1% → -{abs(row['Sensitivity Index']):.2f}% glycine")
    st.subheader("🧪 Experimental Recommendations")
    most_sensitive = df.iloc[0]['Parameter']
    st.markdown(f"""**Highest priority experimental target:** `{most_sensitive}`  
Suggested experimental design:  
1. **Variable control**: Precisely regulate `{most_sensitive}` within physiological range  
2. **Measurements**: Glycine yield, THF/mTHF ratio  
3. **Expected effect**: 1% change in this parameter causes {df.iloc[0]['Sensitivity Index']:.1f}% change in glycine""")

# ============================================================================
# 11. MAIN INTERFACE (ADJUSTED LAYOUT)
# ============================================================================
st.sidebar.subheader("Simulation Settings")
sim_minutes = st.sidebar.slider("Simulation Time (minutes)", 1, 180, 120, key="sim_minutes")

st.subheader("Simulation Control")
control_col1, control_col2 = st.columns([3, 1])  # 左宽右窄

# 左侧列：Run Simulation按钮 + 持久模拟图表
with control_col1:
    if st.button("🚀 Run Simulation", type="primary", key="run_sim_btn"):
        start_time = time.time()
        current_params = get_current_params()
        with st.spinner(f'Simulating {sim_minutes} minutes...'):
            result = run_simulation_with_params(current_params)
        elapsed = time.time() - start_time
        st.write(f"⏱️ Simulation time: {elapsed:.2f} seconds")
        if result['success']:
            st.success(f'✅ Simulation completed! Steps: {result["steps"]}')
            st.session_state.simulation_result = result
        else:
            st.error(f"❌ Simulation failed: {result.get('message', 'Unknown error')}")

    # 持久显示模拟图表（固定在左侧列）
    if "simulation_result" in st.session_state and st.session_state.simulation_result.get("success", False):
        st.subheader("📈 Persistent Simulation Results")
        display_simulation_plots(st.session_state.simulation_result)
        display_results_summary(st.session_state.simulation_result)

# 右侧列：Sensitivity Analysis按钮 + 分析结果
with control_col2:
    st.write("Next, perform sensitivity analysis (tests ±10% perturbations to identify key regulators).")
    if st.button("📊 Sensitivity Analysis", type="secondary", key="sens_btn"):
        run_sensitivity_analysis()