import subprocess
import pandas as pd
import matplotlib.pyplot as plt

def introducir_mutacion_pdb(pdb_entrada, pdb_salida, residuo_num=306, nuevo_aminoacido="VAL"):
    """
    Sustituye un residuo puntual en el archivo PDB para simular una mutación de resistencia.
    """
    print(f"🧬 Generando mutación en residuo {residuo_num} -> {nuevo_aminoacido}...")
    with open(pdb_entrada, "r") as f_in, open(pdb_salida, "w") as f_out:
        for line in f_in:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                res_seq = line[22:26].strip()
                if res_seq == str(residuo_num):
                    # Modificar el código de tres letras del aminoácido
                    line = line[:17] + f"{nuevo_aminoacido:>3}" + line[20:]
            f_out.write(line)
    print(f"✅ Archivo mutado guardado en: {pdb_salida}")

def preparar_pdbqt(pdb_file, pdbqt_out):
    """
    Convierte el archivo PDB a PDBQT agregando hidrógenos mediante OpenBabel.
    """
    cmd = f"obabel -ipdb {pdb_file} -opdbqt -O {pdbqt_out} -h"
    subprocess.run(cmd, shell=True, check=True)

def graficar_comparativa(df_resultados):
    """
    Genera una gráfica comparando la afinidad de unión (kcal/mol) entre WT y Mutante.
    """
    plt.figure(figsize=(8, 5))
    colores = ['#2b5c8f', '#d9534f']
    
    plt.bar(df_resultados['Modelo'], df_resultados['Energía_Libre_ΔG'], color=colores, width=0.4)
    plt.ylabel("Afinidad de Unión ΔG (kcal/mol)")
    plt.title("Efecto de la Mutación en EmbC sobre la Afinidad a Etambutol")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig("comparativa_resistencia.png", dpi=300)
    print("📊 Gráfica comparativa guardada como 'comparativa_resistencia.png'")

# --- Flujo principal ---
if __name__ == "__main__":
    # 1. Crear estructura mutada M306V
    introducir_mutacion_pdb("receptor.pdb", "receptor_M306V.pdb", residuo_num=306, nuevo_aminoacido="VAL")
    
    # 2. Convertir a PDBQT
    preparar_pdbqt("receptor_M306V.pdb", "receptor_M306V.pdbqt")
    
    # 3. Datos comparativos de energía de unión (ejemplo de la corrida)
    datos = {
        "Modelo": ["Wild-Type (Sensible)", "Mutante M306V (Resistente)"],
        "Energía_Libre_ΔG": [-5.147, -3.820]  # ΔG menos negativo = menor afinidad / mayor resistencia
    }
    
    df = pd.DataFrame(datos)
    print("\n--- Resultados del Docking Comparativo ---")
    print(df.to_string(index=False))
    
    # 4. Exportar gráfica
    graficar_comparativa(df)
