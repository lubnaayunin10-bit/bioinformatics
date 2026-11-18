"""
TAHAP 3 — Menyiapkan bahan.

Protein : buang molekul air dan ligan lama, (opsional) tambah hidrogen.
Senyawa : SMILES -> struktur 3D -> minimisasi energi, lalu disaring aturan Lipinski.
"""

# Molekul yang bukan bagian protein dan harus dibuang sebelum docking.
BUKAN_PROTEIN = {
    "HOH", "WAT", "DOD",                      # air
    "SO4", "PO4", "CL", "NA", "K", "MG",      # ion & garam
    "CA", "ZN", "MN", "FE", "NI", "CD",
    "GOL", "EDO", "PEG", "DMS", "ACT", "MPD", "TRS",  # sisa kristalisasi
}


# --------------------------------------------------------------------------
# PROTEIN
# --------------------------------------------------------------------------
def bersihkan_protein(teks_pdb: str, rantai: str = None) -> tuple:
    """
    Sisakan hanya atom protein. Kembaliannya: (teks_pdb_bersih, catatan).

    rantai : huruf rantai yang dipertahankan, mis. "A". None = semua rantai.
    """
    baris_simpan, air, hetero = [], 0, 0

    for baris in teks_pdb.splitlines():
        tag = baris[:6].strip()

        if tag == "ATOM":
            resname = baris[17:20].strip()
            ch = baris[21:22]
            altloc = baris[16:17]
            if rantai and ch != rantai:
                continue
            if altloc not in (" ", "A"):     # buang konformasi ganda
                continue
            if resname in BUKAN_PROTEIN:
                continue
            baris_simpan.append(baris)

        elif tag == "HETATM":
            resname = baris[17:20].strip()
            if resname in ("HOH", "WAT", "DOD"):
                air += 1
            else:
                hetero += 1

        elif tag in ("TER", "END"):
            baris_simpan.append(baris)

    if not any(b.startswith("ATOM") for b in baris_simpan):
        raise ValueError(
            "Tidak ada atom protein tersisa setelah pembersihan. "
            "Coba kosongkan pilihan rantai."
        )

    bersih = "\n".join(baris_simpan)
    if not bersih.rstrip().endswith("END"):
        bersih += "\nEND\n"

    catatan = {
        "atom_protein": sum(1 for b in baris_simpan if b.startswith("ATOM")),
        "air_dibuang": air,
        "hetero_dibuang": hetero,
    }
    return bersih, catatan


def tambah_hidrogen(path_pdb_masuk: str, path_pdb_keluar: str) -> str:
    """
    Tambahkan atom hidrogen dan perbaiki residu yang hilang memakai PDBFixer.
    Kalau PDBFixer belum terpasang, berkas disalin apa adanya.
    Kembaliannya: pesan status untuk ditampilkan di layar.
    """
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError:
        import shutil
        shutil.copyfile(path_pdb_masuk, path_pdb_keluar)
        return ("PDBFixer belum terpasang — hidrogen tidak ditambahkan. "
                "Docking tetap bisa jalan, tapi struktur kurang rapi.")

    fixer = PDBFixer(filename=path_pdb_masuk)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)          # pH fisiologis
    with open(path_pdb_keluar, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    return "Hidrogen ditambahkan pada pH 7,4 dan residu yang hilang diperbaiki."


# --------------------------------------------------------------------------
# SENYAWA
# --------------------------------------------------------------------------
def senyawa_3d(smiles: str, seed: int = 42):
    """SMILES -> molekul RDKit 3D yang sudah ditambah hidrogen dan diminimisasi."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        raise ImportError("RDKit belum terpasang. Jalankan: pip install rdkit")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES tidak bisa dibaca: {smiles}")

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        # Sebagian molekul gagal pada percobaan pertama; coba cara longgar.
        if AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=seed) != 0:
            raise ValueError("Gagal membuat struktur 3D dari SMILES ini.")

    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    return mol


def lipinski(mol) -> dict:
    """
    Hitung empat sifat aturan Lipinski.
    Syarat: BM < 500, logP < 5, donor H <= 5, akseptor H <= 10.
    """
    from rdkit.Chem import Descriptors, Lipinski as L

    bm = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    donor = L.NumHDonors(mol)
    akseptor = L.NumHAcceptors(mol)

    langgar = sum([bm >= 500, logp >= 5, donor > 5, akseptor > 10])
    return {
        "Bobot molekul": round(bm, 2),
        "logP": round(logp, 2),
        "Donor H": donor,
        "Akseptor H": akseptor,
        "Pelanggaran": langgar,
        "Lolos": langgar <= 1,     # 1 pelanggaran masih lazim diterima
    }
