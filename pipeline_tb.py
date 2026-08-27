import urllib.request
import numpy as np
import matplotlib.pyplot as plt
from openbabel import pybel
from vina import Vina

print("=== 1. DESCARGANDO ESTRUCTURAS (EmbC + ETAMBUTOL) ===")
url_ligand = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/14052/SDF?record_type=3d"
urllib.request.urlretrieve(url_ligand, "ethambutol.sdf")

url_protein = "https://files.rcsb.org/download/3PTY.pdb"
urllib.request.urlretrieve(url_protein, "EmbC_WT.pdb")

print("=== 2. LIMPIANDO PROTEÍNA Y PREPARANDO ARCHIVOS ===")
coords = []
with open("EmbC_WT.pdb") as infile, open("EmbC_clean.pdb", "w") as outfile:
    for line in infile:
        if line.startswith("ATOM  "):
            outfile.write(line)
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])

mol_lig = next(pybel.readfile("sdf", "ethambutol.sdf"))
mol_lig.addh()
mol_lig.write("pdbqt", "ethambutol.pdbqt", overwrite=True)

mol_prot = next(pybel.readfile("pdb", "EmbC_clean.pdb"))
mol_prot.write("pdbqt", "EmbC_temp.pdbqt", overwrite=True)

with open("EmbC_temp.pdbqt") as infile, open("EmbC.pdbqt", "w") as outfile:
    for line in infile:
        if not any(line.startswith(tag) for tag in ["ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF"]):
            outfile.write(line)

cx, cy, cz = np.mean(coords, axis=0)
print(f"Centro del receptor detectado: X={cx:.2f}, Y={cy:.2f}, Z={cz:.2f}")

print("=== 3. EJECUTANDO SIMULACIÓN CON AUTODOCK VINA ===")
v = Vina(sf_name="vina")
v.set_receptor("EmbC.pdbqt")
v.set_ligand_from_file("ethambutol.pdbqt")
v.compute_vina_maps(center=[cx, cy, cz], box_size=[28.0, 28.0, 28.0])
v.dock(exhaustiveness=8, n_poses=9)
v.write_poses("EmbC_Ethambutol_docked.pdbqt", overwrite=True)

energies = v.energies()
affinities = [e[0] for e in energies]

print("\n=== TABLA DE RESULTADOS ===")
for i, aff in enumerate(affinities):
    print(f"Pose {i+1}: {aff:.3f} kcal/mol")

print("\n=== 4. GENERANDO GRÁFICO DE ENERGÍAS ===")
plt.figure(figsize=(8, 5))
plt.bar(range(1, len(affinities)+1), affinities, color="#1b7837", edgecolor="black")
plt.title("Binding Affinity per Pose: Ethambutol vs M. tuberculosis EmbC", fontsize=12, fontweight="bold")
plt.xlabel("Binding Pose (Mode)", fontsize=10)
plt.ylabel("Affinity (kcal/mol)", fontsize=10)
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.savefig("ethambutol_affinities.png", dpi=300, bbox_inches="tight")

print("\n=== 5. GENERANDO VISUALIZACIÓN 3D INTERACTIVA EMBEBIDA ===")
with open("EmbC_clean.pdb") as f:
    pdb_str = f.read().replace("\n", "\\n").replace(", \")

with open("EmbC_Ethambutol_docked.pdbqt") as f:
    pdbqt_str = f.read().replace("\n", "\\n").replace(", \")

html_content = f"""<!DOCTYPE html>
<html>
<head>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
</head>
<body style="margin: 0; padding: 0; background-color: #111;">
    <div id="viewport" style="width: 100vw; height: 100vh;"></div>
    <script>
        let viewer = $3Dmol.createViewer(viewport, {{backgroundColor: #111111}});
        
        let pdbData = {pdb_str};
        let pdbqtData = {pdbqt_str};
        
        viewer.addModel(pdbData, pdb);
        viewer.setStyle({{}}, {{cartoon: {{color: spectrum}}}});
        
        viewer.addModel(pdbqtData, pdbqt);
        viewer.setStyle({{model: 1}}, {{stick: {{colorscheme: greenCarbon, radius: 0.3}}}});
        
        viewer.zoomTo();
        viewer.render();
    </script>
</body>
</html>"""

with open("visualizacion_3d.html", "w") as f:
    f.write(html_content)

print("Visualización 3D embebida generada con éxito.")
print("=== PIPELINE COMPLETADO CON ÉXITO ===")
