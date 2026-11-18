"""
DrugFlow — alur penemuan obat berbasis fisika komputasi.

Disusun modular: satu modul untuk satu tahap, sesuai Borang 2 nomor 6.

    sources.py   tahap 1 & 2  — ambil data dari RCSB PDB dan PubChem
    prep.py      tahap 3      — siapkan protein dan senyawa
    gridbox.py   tahap 4a     — tentukan kotak pencarian docking
    docking.py   tahap 4b     — jalankan AutoDock Vina
    viz.py       tahap 6      — tampilkan hasil
"""

__version__ = "1.0"
