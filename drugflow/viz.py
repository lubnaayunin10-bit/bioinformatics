"""
TAHAP 6 — Menampilkan hasil dalam bentuk gambar 3D.

Memakai py3Dmol yang berjalan langsung di dalam halaman web, jadi tidak perlu
memasang PyMOL atau ChimeraX.
"""


def tampilan_3d(protein_pdb: str, pose_pdbqt: str = "",
                center=None, size=None, lebar=760, tinggi=520) -> str:
    """
    Susun potongan HTML berisi tampilan 3D protein + pose hasil docking.
    Kembaliannya HTML siap ditanam ke Streamlit.
    """
    protein_js = _aman(protein_pdb)
    pose_js = _aman(pose_pdbqt) if pose_pdbqt else ""

    kotak = ""
    if center and size:
        kotak = f"""
        viewer.addBox({{
            center: {{x:{center[0]}, y:{center[1]}, z:{center[2]}}},
            dimensions: {{w:{size[0]}, h:{size[1]}, d:{size[2]}}},
            color: 'gold', opacity: 0.18, wireframe: true
        }});"""

    pose = ""
    if pose_js:
        pose = f"""
        viewer.addModel(`{pose_js}`, 'pdbqt');
        viewer.setStyle({{model: 1}}, {{stick: {{colorscheme: 'greenCarbon', radius: 0.22}}}});
        viewer.addStyle({{model: 1}}, {{sphere: {{scale: 0.26}}}});"""

    return f"""
<div style="width:100%;height:{tinggi}px;position:relative;" id="kotak3d"></div>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script>
  let viewer = $3Dmol.createViewer(document.getElementById("kotak3d"),
                                   {{backgroundColor: "#0f172a"}});
  viewer.addModel(`{protein_js}`, 'pdb');
  viewer.setStyle({{model: 0}}, {{cartoon: {{color: 'spectrum', opacity: 0.85}}}});
  {pose}
  {kotak}
  viewer.zoomTo({{model: {1 if pose_js else 0}}});
  viewer.zoom(0.85);
  viewer.render();
</script>
"""


def _aman(teks: str) -> str:
    """Hindari karakter yang merusak template string JavaScript."""
    return teks.replace("\\", "\\\\").replace("`", "'").replace("${", "$ {")


def grafik_skor(skor: list):
    """Diagram batang skor docking tiap pose. Kembaliannya figure matplotlib."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 3.2))
    label = [f"Pose {i+1}" for i in range(len(skor))]
    warna = ["#16a34a" if s == min(skor) else "#94a3b8" for s in skor]

    ax.barh(label, skor, color=warna)
    ax.invert_yaxis()
    ax.set_xlabel("Binding affinity (kcal/mol)")
    ax.axvline(-7, color="#f59e0b", linestyle="--", lw=1.2)
    ax.text(-7, -0.7, " patokan -7", color="#f59e0b", fontsize=8)
    for i, s in enumerate(skor):
        ax.text(s, i, f" {s:.1f}", va="center", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig
