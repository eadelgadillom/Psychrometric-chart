
# Tu código completo con la sección mejorada para evitar superposición de etiquetas en líneas Wet Bulb
# Solo se incluye un ejemplo para la parte relevante (líneas de bulbo húmedo):

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

            # Mejor posicionamiento de etiquetas
            if Tw >= 35 and len(valid_T) > 0:
                label_idx = int(len(valid_T) * 0.8)
                
                # Verificar si la etiqueta se superpondrá
                show_label = True
                x = valid_T[label_idx]
                y = W_wb_array[label_idx]

                # Si el punto está muy cerca del borde inferior, no mostrar etiqueta
                if y < 2:  # 2 g/kg, ajustar según convenga
                    show_label = False

                if show_label:
                    if label_idx > 1 and label_idx < len(valid_T) - 1:
                        dx = valid_T[label_idx + 1] - valid_T[label_idx - 1]
                        dy = W_wb_array[label_idx + 1] - W_wb_array[label_idx - 1]
                        angle = np.degrees(np.arctan2(dy, dx))
                    else:
                        angle = 0
                    if Tw == 45:
                        ax.text(
                            x, y,
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
                            x, y,
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
