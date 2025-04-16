import tkinter as tk
from tkinter import messagebox
from psychrolib import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from adjustText import adjust_text
import csv

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
        v = GetMoistAirVolume(Tdb, W, P_atm)  # Specific volume (m³/kg dry air)
        rho = 1 / v  # Density (kg/m³)

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

    fig, ax = plt.subplots(figsize=(8, 8))

    # Adjust T_range based on input temperature
    max_T = max(50, Tdb + 10)
    min_T = min(0, Tdb - 20)
    T_range = np.linspace(min_T, max_T, 300)

    # --- Relative humidity curves from 10% to 100% ---
    rh_values = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
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

            # Make 10%, 50% and 100% RH lines more prominent
            if rh in [0.1, 0.5, 1.0]:
                line_style = "-"
                line_width = 1.5
                line_color = "darkgreen"
                alpha = 0.9
            else:
                line_style = "-"
                line_width = 0.8
                line_color = "g"
                alpha = 0.5

            (line,) = ax.plot(
                T_valid,
                W_rh_array,
                line_color,
                linestyle=line_style,
                alpha=alpha,
                lw=line_width,
            )

            # Add label at midpoint of curve (only for prominent lines)
            if len(W_rh) > 2 and rh in [0.02, 0.05, 0.1, 0.5, 1.0]:
                label_idx = len(W_rh) // 3  # Point at 1/3 of the curve
                x_pos = T_valid[label_idx]
                y_pos = W_rh_array[label_idx]
                y_offset = max(W_rh_array) * 0.00

                # Calculate text rotation angle
                if label_idx > 0 and label_idx < len(W_rh) - 1:
                    dx = T_valid[label_idx + 1] - T_valid[label_idx - 1]
                    dy = W_rh_array[label_idx + 1] - W_rh_array[label_idx - 1]
                    angle = np.degrees(np.arctan2(dy - 0.05, dx))
                else:
                    angle = 0

                ax.text(
                    x_pos,
                    y_pos + y_offset,
                    f"RH{rh*100:.0f}%",
                    fontsize=9,
                    color=line_color,
                    alpha=0.9,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=0.5),
                    rotation=angle,
                    rotation_mode="anchor",
                    ha="center",
                    va="bottom",
                )

    # --- Constant wet bulb temperature lines ---
    wb_step = np.arange(
        max(0, min(T_range)), min(max(T_range), 100), 5
    )  # Lines every 5°C
    for Tw in wb_step:
        W_wb = []
        valid_T = []
        for T in T_range:
            try:
                W_val = GetHumRatioFromTWetBulb(T, Tw, P_atm)
                if W_val > 0:  # Only positive values
                    W_wb.append(W_val)
                    valid_T.append(T)
            except:
                continue

        if len(W_wb) > 3:
            W_wb_array = np.array(W_wb) * 1000  # Convert to g/kg
            (line,) = ax.plot(
                valid_T,
                W_wb_array,
                color="blue",
                linestyle="--",
                alpha=0.6,  # Increased visibility
                lw=1,
            )

            # Better positioned label
            if len(valid_T) > 0:
                label_idx = int(len(valid_T) * 0.8)
                label_idx = min(label_idx, len(valid_T) - 1)

                # Calculate curve angle for text rotation
                if label_idx > 1 and label_idx < len(valid_T) - 1:
                    dx = valid_T[label_idx + 1] - valid_T[label_idx - 1]
                    dy = W_wb_array[label_idx + 1] - W_wb_array[label_idx - 1]
                    angle = np.degrees(np.arctan2(dy-0.3, dx))
                else:
                    angle = 0
                if Tw == 45:
                    ax.text(
                    valid_T[label_idx],
                    W_wb_array[label_idx],
                    f"Wet Bulb Temperature {Tw}°C",
                    fontsize=8,
                    color="blue",
                    alpha=0.8,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=0.3),
                    rotation=angle,
                    rotation_mode="anchor",
                    ha="center",
                    va="center",
                    )
                else:    
                    ax.text(
                        valid_T[label_idx],
                        W_wb_array[label_idx],
                        f"{Tw}°C",
                        fontsize=8,
                        color="blue",
                        alpha=0.8,
                        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=0.3),
                        rotation=angle,
                        rotation_mode="anchor",
                        ha="center",
                        va="center",
                    )

    # --- Constant Enthalpy Lines ---
    h_step = 10  # kJ/kg
    max_h = int(max(50, enthalpy + 100))
    min_h = int(min(0, enthalpy - 30))
    enthalpy_lines = []

    for h_val in range(min_h, max_h, h_step):
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
        enthalpy_lines.append((h_val, line, valid_T, W_h_array))

    # --- Second axis for Enthalpy ---
    ax2 = ax.twinx()
    ax2.set_ylim(min_h, max_h)
    ax2.set_ylabel("Enthalpy [kJ/kg]", color="tab:red", labelpad=20)
    ax2.set_yticks([])

    # Add enthalpy values to the right axis
    for h_val, line, valid_T, W_h_array in enthalpy_lines:
        valid_indices = ~np.isnan(W_h_array)
        if np.any(valid_indices):
            idx = np.where(valid_indices)[0][-1]  # Last valid point
            ax.text(
                valid_T[idx],
                W_h_array[idx],
                f"{h_val}",
                fontsize=7,
                color="red",
                alpha=0.6,
                ha="left",
                va="center",
            )

    # --- Chart formatting ---
    ax.set_title(f"Psychrometric Chart (P = {P_atm/1000:.1f} kPa)")
    ax.set_xlabel("Dry Bulb Temperature (°C)")
    ax.set_ylabel("Humidity Ratio (g/kg dry air)")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.5)

    # Set dynamic limits
    max_y = max(25, W * 1000 * 3.5)
    ax.set_xlim(min_T, max_T)
    ax.set_ylim(0, max_y)

    # Display in Tkinter
    canvas = FigureCanvasTkAgg(fig, master=graph_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    plot_psychrometric_chart.canvas = canvas


def export_plot_as_png():
    if not hasattr(plot_psychrometric_chart, "canvas"):
        messagebox.showerror("Error", "No plot available")
        return

    file_path = tk.filedialog.asksaveasfilename(
        defaultextension=".png",
        filetypes=[("PNG Image", "*.png"), ("All Files", "*.*")],
        title="Save plot as PNG",
    )

    if file_path:
        try:
            plot_psychrometric_chart.canvas.figure.savefig(file_path, dpi=300)
            messagebox.showinfo("Success", f"Plot saved to:\n{file_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save plot:\n{e}")


def export_data_as_csv():
    if not result_text.get().startswith("Atmospheric Pressure:"):
        messagebox.showerror("Error", "No data available")
        return

    file_path = tk.filedialog.asksaveasfilename(
        defaultextension=".csv",
        filetypes=[("CSV File", "*.csv"), ("All Files", "*.*")],
        title="Save results as CSV",
    )

    if file_path:
        try:
            # Extract data from results text
            lines = result_text.get().split("\n")
            data = {
                line.split(":")[0].strip(): line.split(":")[1].strip()
                for line in lines
                if ":" in line
            }

            # Write to CSV
            with open(file_path, mode="w", newline="") as file:
                writer = csv.writer(file)
                writer.writerow(["Parameter", "Value"])  # Header
                for key, value in data.items():
                    writer.writerow([key, value])

            messagebox.showinfo("Success", f"Data saved to:\n{file_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save CSV:\n{e}")


# Create the GUI application
def create_gui():
    global root, frame, entry_temp, entry_rh, result_text, graph_frame

    root = tk.Tk()
    root.title("Psychrometric Calculator")
    root.geometry("800x1000")

    # Main frame
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
    entry_temp.insert(0, "25")

    tk.Label(input_frame, text="Relative Humidity (%):").grid(
        row=0, column=2, sticky="e", padx=5, pady=5
    )
    entry_rh = tk.Entry(input_frame, width=10)
    entry_rh.grid(row=0, column=3, padx=5, pady=5)
    entry_rh.insert(0, "50")

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

    export_frame = tk.Frame(results_frame)
    export_frame.pack(fill=tk.X, pady=5)

    # Export buttons
    btn_export_png = tk.Button(
        export_frame,
        text="Save Plot (PNG)",
        command=export_plot_as_png,
        bg="#2196F3",
        fg="white",
    )
    btn_export_png.pack(side=tk.LEFT, padx=5)

    btn_export_csv = tk.Button(
        export_frame,
        text="Save Results (CSV)",
        command=export_data_as_csv,
        bg="#FF9800",
        fg="white",
    )
    btn_export_csv.pack(side=tk.LEFT, padx=5)

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
