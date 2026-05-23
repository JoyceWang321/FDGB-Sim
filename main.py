import streamlit as st
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'DejaVu Sans'

st.set_page_config(layout="wide", page_title="Glycine Simulator")
st.title("E. coli Glycine Synthesis")

# 侧边栏 - 简化参数
with st.sidebar:
    st.header("Parameters")
    vmax_gcs = st.slider("GCS Vmax", 0.01, 1.0, 0.1, 0.01)
    vmax_dhfr = st.slider("DHFR Vmax", 0.1, 5.0, 1.0, 0.1)
    sim_time = st.slider("Time (min)", 1, 120, 30)
    
    st.header("Initial (mM)")
    formate = st.number_input("Formate", 0.0, 100.0, 20.0)
    thf = st.number_input("THF", 0.0, 5.0, 0.2)
    nadh = st.number_input("NADH", 0.0, 5.0, 0.53)
    nh3 = st.number_input("NH3", 0.0, 100.0, 50.0)

# 简化的 ODE（避免刚性）
def system(t, y):
    F, T, D, G, N = y  # Formate, THF, DHF, Glycine, NADH
    
    # 简化的速率方程
    v_gcs = vmax_gcs * F * T * N / ((F+10)*(T+0.1)*(N+0.01))
    v_dhfr = vmax_dhfr * D * 0.42 / ((D+0.01)*(0.42+0.01))
    v_regen = 0.1 * (1.69 - N)  # NADH 再生
    
    return [
        -v_gcs,           # dF/dt
        -v_gcs + v_dhfr,  # dT/dt
        -v_dhfr,          # dD/dt
        v_gcs,            # dG/dt
        -v_gcs + v_regen  # dN/dt
    ]

if st.button("Run", type="primary"):
    y0 = [formate, thf, 0.05, 0.0, nadh]
    t_span = [0, sim_time*60]  # 分钟转秒
    t_eval = np.linspace(0, sim_time*60, 500)
    
    try:
        sol = solve_ivp(system, t_span, y0, method='RK45', 
                      t_eval=t_eval, rtol=1e-3, atol=1e-6)
        
        fig, ax = plt.subplots(2, 2, figsize=(10, 8))
        
        time = sol.t / 60  # 转回分钟
        
        ax[0,0].plot(time, sol.y[0], 'g-', label='Formate')
        ax[0,0].plot(time, sol.y[3], 'r-', label='Glycine')
        ax[0,0].set_title('Substrates')
        ax[0,0].legend()
        
        ax[0,1].plot(time, sol.y[1], 'b-', label='THF')
        ax[0,1].plot(time, sol.y[2], 'm-', label='DHF')
        ax[0,1].set_title('Folate Cycle')
        ax[0,1].legend()
        
        ax[1,0].plot(time, sol.y[4], 'c-', label='NADH')
        ax[1,0].axhline(y=0.53, color='gray', linestyle='--', label='Initial')
        ax[1,0].set_title('NADH')
        ax[1,0].legend()
        
        ax[1,1].plot(time, sol.y[3], 'r-', linewidth=3)
        ax[1,1].set_title('Glycine Production')
        ax[1,1].set_ylabel('mM')
        
        plt.tight_layout()
        st.pyplot(fig)
        
        st.metric("Final Glycine", f"{sol.y[3,-1]:.3f} mM")
        st.metric("Conversion", f"{(sol.y[3,-1]/formate)*100:.1f}%")
        
    except Exception as e:
        st.error(f"Error: {e}")
        st.info("Try reducing simulation time or adjusting parameters")