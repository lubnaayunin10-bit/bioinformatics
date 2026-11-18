"""
TAHAP 1 & 2 — Mengambil data dari database publik di internet.

Dua database yang dipakai:
  1. RCSB Protein Data Bank  -> struktur 3D protein  (https://www.rcsb.org)
  2. PubChem (NCBI)          -> struktur senyawa obat (https://pubchem.ncbi.nlm.nih.gov)

Keduanya gratis, terbuka, dan tidak perlu API key.
"""

import requests

TIMEOUT = 30
UA = {"User-Agent": "DrugFlow/1.0 (tugas kuliah; python-requests)"}


# --------------------------------------------------------------------------
# 1. RCSB PROTEIN DATA BANK
# --------------------------------------------------------------------------
def ambil_protein(pdb_id: str) -> str:
    """Unduh berkas .pdb mentah dari RCSB. Kembaliannya berupa teks."""
    pdb_id = pdb_id.strip().upper()
    if len(pdb_id) != 4:
        raise ValueError(f"PDB ID harus 4 karakter, bukan '{pdb_id}'.")

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    r = requests.get(url, timeout=TIMEOUT, headers=UA)
    if r.status_code == 404:
        raise ValueError(
            f"PDB ID '{pdb_id}' tidak ditemukan di RCSB. Cek lagi kodenya di rcsb.org."
        )
    r.raise_for_status()
    if "ATOM" not in r.text:
        raise ValueError(f"Berkas {pdb_id}.pdb terunduh tapi tidak berisi koordinat atom.")
    return r.text


def info_protein(pdb_id: str) -> dict:
    """Ambil keterangan ringkas protein: judul, resolusi, metode eksperimen."""
    pdb_id = pdb_id.strip().upper()
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    try:
        r = requests.get(url, timeout=TIMEOUT, headers=UA)
        r.raise_for_status()
        d = r.json()
    except Exception:
        return {}

    resolusi = None
    for kunci in ("rcsb_entry_info", "refine", "em_3d_reconstruction"):
        bagian = d.get(kunci)
        if isinstance(bagian, dict):
            resolusi = bagian.get("resolution_combined") or bagian.get("ls_d_res_high")
        elif isinstance(bagian, list) and bagian:
            resolusi = bagian[0].get("resolution") or bagian[0].get("ls_d_res_high")
        if resolusi:
            break
    if isinstance(resolusi, list) and resolusi:
        resolusi = resolusi[0]

    metode = ""
    exptl = d.get("exptl") or []
    if exptl:
        metode = exptl[0].get("method", "")

    return {
        "id": pdb_id,
        "judul": (d.get("struct") or {}).get("title", ""),
        "resolusi": resolusi,
        "metode": metode,
    }


# --------------------------------------------------------------------------
# 2. PUBCHEM
# --------------------------------------------------------------------------
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# PubChem pernah mengganti nama properti SMILES. Dicoba berurutan supaya
# program tetap jalan di versi API mana pun.
_SMILES_KEYS = ["CanonicalSMILES", "SMILES", "ConnectivitySMILES", "IsomericSMILES"]


def cari_cid(nama: str) -> int:
    """Ubah nama senyawa (mis. 'aspirin') menjadi CID PubChem."""
    url = f"{PUG}/compound/name/{requests.utils.quote(nama.strip())}/cids/JSON"
    r = requests.get(url, timeout=TIMEOUT, headers=UA)
    if r.status_code == 404:
        raise ValueError(f"Senyawa '{nama}' tidak ditemukan di PubChem.")
    r.raise_for_status()
    daftar = r.json().get("IdentifierList", {}).get("CID", [])
    if not daftar:
        raise ValueError(f"Senyawa '{nama}' tidak ditemukan di PubChem.")
    return int(daftar[0])


def ambil_senyawa(cid) -> dict:
    """Ambil SMILES dan keterangan dasar senyawa berdasarkan CID PubChem."""
    cid = str(cid).strip()
    if not cid.isdigit():
        raise ValueError(f"CID harus berupa angka, bukan '{cid}'.")

    props = "Title,MolecularFormula,MolecularWeight," + ",".join(_SMILES_KEYS)
    url = f"{PUG}/compound/cid/{cid}/property/{props}/JSON"
    r = requests.get(url, timeout=TIMEOUT, headers=UA)

    if r.status_code >= 400:
        # Kalau sebagian properti ditolak, coba lagi dengan properti minimal.
        url = f"{PUG}/compound/cid/{cid}/property/Title,MolecularWeight/JSON"
        r = requests.get(url, timeout=TIMEOUT, headers=UA)
        if r.status_code == 404:
            raise ValueError(f"CID {cid} tidak ditemukan di PubChem.")
        r.raise_for_status()

    tabel = r.json().get("PropertyTable", {}).get("Properties", [])
    if not tabel:
        raise ValueError(f"PubChem tidak mengembalikan data untuk CID {cid}.")
    p = tabel[0]

    smiles = next((p[k] for k in _SMILES_KEYS if p.get(k)), None)
    if not smiles:
        smiles = _ambil_smiles_cadangan(cid)

    return {
        "cid": int(cid),
        "nama": p.get("Title", f"CID {cid}"),
        "rumus": p.get("MolecularFormula", ""),
        "bobot_pubchem": p.get("MolecularWeight", ""),
        "smiles": smiles,
    }


def _ambil_smiles_cadangan(cid: str) -> str:
    """Jalan terakhir: unduh SDF lalu baca SMILES-nya dengan RDKit."""
    url = f"{PUG}/compound/cid/{cid}/SDF"
    r = requests.get(url, timeout=TIMEOUT, headers=UA)
    r.raise_for_status()
    try:
        from rdkit import Chem
    except ImportError:
        raise ValueError(
            "PubChem tidak memberi SMILES dan RDKit belum terpasang untuk membacanya."
        )
    mol = Chem.MolFromMolBlock(r.text)
    if mol is None:
        raise ValueError(f"Gagal membaca struktur senyawa CID {cid}.")
    return Chem.MolToSmiles(mol)
