import tkinter as tk
from tkinter import messagebox
from psychrolib import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Set SI Units
SetUnitSystem(SI)


def calculate():
    try:
        Tdb = float(entry_temp.get())  # Dry bulb temperature (°C)
        RH = float(entry_rh.get()) / 100.0  # Relative Humidity (0-1)
        P_atm = 101325  # Atmospheric pressure (Pa)

        # Psychrometric calculations
        Psat = GetSatVapPres(Tdb)  # Saturation vapor pressure at Tdb
        Pv = RH * Psat  # Vapor pressure
        W = 0.622 * Pv / (P_atm - Pv)  # Humidity ratio
        Tdp = GetTDewPointFromRelHum(Tdb, RH)  # Dew point temperature (°C)
        Twb = GetTWetBulbFromRelHum(Tdb, RH, P_atm)  # Wet bulb temperature (°C)
        h = GetMoistAirEnthalpy(Tdb, W) / 1000  # Enthalpy (kJ/kg dry air)
        v = GetMoistAirVolume(Tdb, W, P_atm)  # Specific volume (m3/kg dry air)
        rho = 1 / v  # Density (kg/m3)

        # Convert W to g/kg
        W_g_per_kg = W * 1000

        # Saturation temperature is approximately equal to wet bulb temp when RH is low
        Tsat = Twb

        result_text.set(
            f"Atmospheric Pressure: {P_atm:.2f} Pa\n"
            f"Dry Bulb Temperature: {Tdb:.2f} °C\n"
            f"Relative Humidity: {RH*100:.2f} %\n"
            f"Humidity Ratio: {W_g_per_kg:.3f} g/kg dry air\n"
            f"Dew Point Temperature: {Tdp:.2f} °C\n"
            f"Wet Bulb Temperature: {Twb:.2f} °C\n"
            f"Saturation Temperature: {Tsat:.2f} °C (approx)\n"
            f"Enthalpy: {h:.2f} kJ/kg dry air\n"
            f"Vapor Pressure: {Pv:.2f} Pa\n"
            f"Saturation Vapor Pressure: {Psat:.2f} Pa\n"
            f"Specific Volume: {v:.3f} m³/kg dry air\n"
            f"Density: {rho:.3f} kg/m³"
        )

        # Plot the psychrometric chart with the calculated point
        plot_psychrometric_chart(Tdb, RH, P_atm, W, Tdp, Twb, Pv, Psat, h)

    except ValueError as e:
        messagebox.showerror("Error", "Please enter valid numerical values.")
    except Exception as e:
        messagebox.showerror("Error", f"Calculation error: {e}")


def plot_psychrometric_chart(Tdb, RH, P_atm, W, Tdp, Twb, Pv, Psat, enthalpy):
    # Clear previous plot if exists
    if hasattr(plot_psychrometric_chart, "canvas"):
        plot_psychrometric_chart.canvas.get_tk_widget().destroy()

    fig, ax = plt.subplots(figsize=(10, 6))

    # Adjust T_range based on input temperature
    max_T = max(50, Tdb + 20)
    min_T = min(0, Tdb - 20)
    T_range = np.linspace(min_T, max_T, 300)

    # --- Relative humidity curves from 2% to 100% ---
    rh_values = [RH, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for rh in rh_values:
        W_rh = []
        T_valid = []

        for i, T in enumerate(T_range):
            try:
                Psat = GetSatVapPres(T)

                Pv = rh * Psat
                W_val = 0.622 * Pv / (P_atm - Pv)

                # Check if this W value is decreasing (curve starting to bend back)
                if i > 0 and len(W_rh) > 0 and W_val < W_rh[-1]:
                    break

                W_rh.append(W_val)
                T_valid.append(T)
            except:
                # If calculation fails, stop this curve
                break

        # Only plot if we have valid points
        if len(W_rh) > 0:
            W_rh_array = np.array(W_rh) * 1000  # Convert to g/kg
            (line,) = ax.plot(T_valid, W_rh_array, "g", alpha=0.5, lw=0.8)

            # Add label at midpoint of curve
            if len(W_rh) > 2:
                midpoint_idx = len(W_rh) // 2
                ax.text(
                    T_valid[midpoint_idx],
                    W_rh[midpoint_idx] * 1000,
                    f"{rh*100:.1f}%",
                    fontsize=7,
                    color="green",
                    alpha=0.7,
                )

    # --- Constant Enthalpy Lines ---
    h_step = 10  # kJ/kg
    max_h = int(max(50, enthalpy + 20))  # Use the calculated enthalpy here
    for h_val in range(0, max_h, h_step):
        W_h = []
        valid_T = []
        for T in T_range:
            try:
                W_temp = GetHumRatioFromEnthalpyAndTDryBulb(h_val * 1000, T)
                if not np.isnan(W_temp) and W_temp > 0:
                    W_h.append(W_temp)
                    valid_T.append(T)
                else:
                    W_h.append(np.nan)
                    valid_T.append(np.nan)
            except:
                W_h.append(np.nan)
                valid_T.append(np.nan)

        W_h_array = np.array(W_h) * 1000  # Convert to g/kg
        (line,) = ax.plot(valid_T, W_h_array, "r:", alpha=0.4, lw=0.8)

        # Add label at a suitable position
        valid_indices = ~np.isnan(W_h_array)
        if np.any(valid_indices):
            idx = np.where(valid_indices)[0][-1]  # Use last valid point
            ax.text(
                valid_T[idx],
                W_h_array[idx],
                f"{h_val}",
                fontsize=7,
                color="red",
                alpha=0.6,
            )

    # --- Constant wet bulb temperature lines ---
    wb_step = 5  # °C
    max_wb = int(max(30, Twb + 15))
    min_wb = int(min(0, Twb - 15))
    for Tw in range(min_wb, max_wb, wb_step):
        W_wb = []
        valid_T = []
        for T in T_range:
            if Tw <= T:
                try:
                    W_val = GetHumRatioFromTWetBulb(Tw, T, P_atm)
                    W_wb.append(W_val)
                    valid_T.append(T)
                except:
                    W_wb.append(np.nan)
                    valid_T.append(np.nan)
            else:
                W_wb.append(np.nan)
                valid_T.append(np.nan)

        (line,) = ax.plot(
            T_range, np.array(W_wb) * 1000, color="blue", linestyle="--", alpha=0.3
        )

        # Add label at the start of the line
        valid_indices = ~np.isnan(np.array(W_wb))
        if np.any(valid_indices):
            first_idx = np.where(valid_indices)[0][0]
            ax.text(
                T_range[first_idx],
                W_wb[first_idx] * 1000,
                f"Twb={Tw}°C",
                fontsize=7,
                color="blue",
                alpha=0.6,
            )



    # --- Calculated point ---
    ax.plot(Tdb, W * 1000, "ro", markersize=8, label="Current State")
    ax.annotate(
        f"Tdb={Tdb:.1f}°C\nRH={RH*100:.1f}%\nW={W*1000:.1f}g/kg",
        (Tdb, W * 1000),
        textcoords="offset points",
        xytext=(10, 5),
        bbox=dict(boxstyle="round,pad=0.3", fc="yellow", alpha=0.3),
    )

    # --- Chart formatting ---
    ax.set_title(f"Psychrometric Chart (P = {P_atm/1000:.1f} kPa)")
    ax.set_xlabel("Dry Bulb Temperature (°C)")
    ax.set_ylabel("Humidity Ratio (g/kg dry air)")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.5)

    # Set dynamic limits based on the calculated point
    max_y = max(25, W * 1000 * 1.5)
    ax.set_xlim(min_T, max_T)
    ax.set_ylim(0, max_y)

    ax.legend(loc="upper left")

    # Display in Tkinter
    canvas = FigureCanvasTkAgg(fig, master=graph_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    plot_psychrometric_chart.canvas = canvas


# Create the GUI application
def create_gui():
    global root, frame, entry_temp, entry_rh, result_text, graph_frame

    root = tk.Tk()
    root.title("Psychrometric Calculator")
    root.geometry("1000x800")  # Set initial window size

    # Create main frame with padding
    main_frame = tk.Frame(root, padx=10, pady=10)
    main_frame.pack(fill=tk.BOTH, expand=True)

    # Input frame
    input_frame = tk.LabelFrame(main_frame, text="Input Parameters", padx=10, pady=10)
    input_frame.pack(fill=tk.X, padx=5, pady=5)

    # Input fields
    tk.Label(input_frame, text="Dry Bulb Temperature (°C):").grid(
        row=0, column=0, sticky="e", padx=5, pady=5
    )
    entry_temp = tk.Entry(input_frame, width=10)
    entry_temp.grid(row=0, column=1, padx=5, pady=5)
    entry_temp.insert(0, "25")  # Default value

    tk.Label(input_frame, text="Relative Humidity (%):").grid(
        row=0, column=2, sticky="e", padx=5, pady=5
    )
    entry_rh = tk.Entry(input_frame, width=10)
    entry_rh.grid(row=0, column=3, padx=5, pady=5)
    entry_rh.insert(0, "50")  # Default value

    # Calculate button
    calculate_btn = tk.Button(
        input_frame,
        text="Calculate",
        command=calculate,
        bg="#4CAF50",
        fg="white",
        padx=10,
    )
    calculate_btn.grid(row=0, column=4, padx=20, pady=5)

    # Results frame
    results_frame = tk.LabelFrame(
        main_frame, text="Psychrometric Properties", padx=10, pady=10
    )
    results_frame.pack(fill=tk.X, padx=5, pady=5)

    # Results text
    result_text = tk.StringVar()
    result_text.set("Results will appear here after calculation.")
    tk.Label(
        results_frame, textvariable=result_text, justify="left", font=("Courier", 10)
    ).pack(fill=tk.X, padx=5, pady=5)

    # Graph frame
    graph_frame = tk.LabelFrame(
        main_frame, text="Psychrometric Chart", padx=10, pady=10
    )
    graph_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    # Status bar
    status_bar = tk.Label(
        root,
        text="Ready. Enter values and click Calculate.",
        bd=1,
        relief=tk.SUNKEN,
        anchor=tk.W,
    )
    status_bar.pack(side=tk.BOTTOM, fill=tk.X)

    return root


# Run the application
if __name__ == "__main__":
    root = create_gui()
    root.mainloop()
