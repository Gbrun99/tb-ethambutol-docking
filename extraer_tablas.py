import os
import pandas as pd

def buscar_y_extraer():
    # Buscar archivos .pdbqt en la carpeta
    archivos = [f for f in os.listdir('.') if f.endswith('.pdbqt')]

    if not archivos:
        print("❌ No se encontraron archivos .pdbqt en esta carpeta.")
        return

    print(f"📁 Archivos .pdbqt encontrados: {archivos}\n")

    for archivo in archivos:
        poses = []
        with open(archivo, 'r') as f:
            for line in f:
                if "REMARK VINA RESULT:" in line:
                    parts = line.split()
                    poses.append({
                        "Pose": len(poses) + 1,
                        "ΔG (kcal/mol)": float(parts[2]),
                        "RMSD l.b.": float(parts[3]),
                        "RMSD u.b.": float(parts[4])
                    })

        if poses:
            df = pd.DataFrame(poses)
            print(f"=== TABLA DE ENERGÍAS: {archivo} ===")
            print(df.to_string(index=False))
            print("-" * 40)
            df.to_csv(f"tabla_{archivo}.csv", index=False)
        else:
            print(f"⚠️ El archivo {archivo} no contiene resultados de Vina (REMARK VINA RESULT).")

if __name__ == "__main__":
    buscar_y_extraer()
