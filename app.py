"""
DrugFlow — aplikasi web alur penemuan obat berbasis fisika komputasi.

Jalankan dengan:   streamlit run app.py
"""

import os
import tempfile
import traceback

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from drugflow import sources, prep, gridbox, docking, viz

st.set_page_config(page_title="DrugFlow", page_icon="🧬", layout="wide")

# --------------------------------------------------------------------------
# Tampilan atas
# --------------------------------------------------------------------------
st.title("DrugFlow")
st.caption(
    "Alur penemuan obat berbasis fisika komputasi — dari ID protein dan senyawa "
    "sampai skor docking, dalam satu klik."
)
st.markdown(
    "<div style='display:inline-block;padding:6px 14px;border-radius:999px;"
    "background:#1e293b;border:1px solid #334155;color:#94a3b8;"
    "font-size:0.85rem;margin-bottom:8px;'>"
    "Kelompok 1 · Biologi Umum · FSK C</div>",
    unsafe_allow_html=True,
)

if "hasil" not in st.session_state:
    st.session_state.hasil = None

# --------------------------------------------------------------------------
# Panel samping
# --------------------------------------------------------------------------
with st.sidebar:
    st.caption("Kelompok 1 · Biologi Umum · FSK C")
    st.header("Masukan")

    pdb_id = st.text_input("PDB ID protein", value="6LU7",
                           help="Contoh: 6LU7 (main protease SARS-CoV-2)")
    rantai = st.text_input("Rantai yang dipakai", value="A",
                           help="Kosongkan untuk memakai semua rantai")

    st.divider()
    cara = st.radio("Senyawa dicari lewat", ["CID PubChem", "Nama senyawa"])
    if cara == "CID PubChem":
        cid_input = st.text_input("CID", value="5280343",
                                  help="5280343 = kuersetin")
        nama_input = ""
    else:
        nama_input = st.text_input("Nama", value="quercetin")
        cid_input = ""

    st.divider()
    st.subheader("Pengaturan docking")
    exhaustiveness = st.slider("Exhaustiveness", 1, 32, 8,
                               help="Makin besar makin teliti, tapi makin lama")
    n_poses = st.slider("Jumlah pose", 1, 20, 5)
    path_vina = st.text_input("Lokasi program Vina (opsional)", value="",
                              placeholder="mis. C:\\vina\\vina.exe")

    st.divider()
    perkakas = docking.periksa_perkakas(path_vina)
    st.subheader("Status perkakas")
    for nama, nilai in [
        ("RDKit", perkakas["rdkit"]),
        ("Meeko", perkakas["meeko"]),
        ("Open Babel", perkakas["openbabel"]),
        ("AutoDock Vina", bool(perkakas["vina_program"] or perkakas["vina_python"])),
        ("PDBFixer", perkakas["pdbfixer"]),
    ]:
        st.write(("✅ " if nilai else "⬜ ") + nama)
    if not (perkakas["vina_program"] or perkakas["vina_python"]):
        st.info("Tanpa Vina, tahap 1–3 tetap jalan. Docking saja yang dilewati.")

    jalan = st.button("Run", type="primary", use_container_width=True)

# --------------------------------------------------------------------------
# Alur utama
# --------------------------------------------------------------------------
if jalan:
    kerja = tempfile.mkdtemp(prefix="drugflow_")
    hasil = {"kerja": kerja}

    try:
        # ---------------- TAHAP 1 & 2 : ambil data ------------------------
        with st.status("Tahap 1–2 · Mengambil data dari PDB dan PubChem...",
                       expanded=True) as kotak:
            st.write(f"Mengunduh {pdb_id.upper()} dari RCSB PDB...")
            teks_pdb = sources.ambil_protein(pdb_id)
            hasil["pdb_mentah"] = teks_pdb
            hasil["info"] = sources.info_protein(pdb_id)

            cid = cid_input or sources.cari_cid(nama_input)
            st.write(f"Mengambil senyawa CID {cid} dari PubChem...")
            hasil["senyawa"] = sources.ambil_senyawa(cid)
            kotak.update(label="Tahap 1–2 · Data berhasil diambil ✅",
                         state="complete", expanded=False)

        # ---------------- TAHAP 3 : siapkan bahan -------------------------
        with st.status("Tahap 3 · Menyiapkan protein dan senyawa...",
                       expanded=True) as kotak:
            bersih, catatan = prep.bersihkan_protein(teks_pdb, rantai.strip() or None)
            hasil["pdb_bersih"], hasil["catatan"] = bersih, catatan
            st.write(f"Dibuang {catatan['air_dibuang']} molekul air "
                     f"dan {catatan['hetero_dibuang']} atom hetero.")

            path_bersih = os.path.join(kerja, "protein_bersih.pdb")
            open(path_bersih, "w").write(bersih)

            path_siap = os.path.join(kerja, "protein_siap.pdb")
            hasil["pesan_h"] = prep.tambah_hidrogen(path_bersih, path_siap)
            st.write(hasil["pesan_h"])

            st.write("Membuat struktur 3D senyawa dari SMILES...")
            mol = prep.senyawa_3d(hasil["senyawa"]["smiles"])
            hasil["lipinski"] = prep.lipinski(mol)
            kotak.update(label="Tahap 3 · Bahan siap ✅",
                         state="complete", expanded=False)

        # ---------------- TAHAP 4 : docking -------------------------------
        hasil["kotak"] = gridbox.hitung_kotak(teks_pdb)

        if perkakas["vina_program"] or perkakas["vina_python"]:
            with st.status("Tahap 4 · Menjalankan docking...", expanded=True) as kotak:
                st.write(f"Grid box dihitung dari {hasil['kotak']['sumber']}.")

                lig_pdbqt = os.path.join(kerja, "ligan.pdbqt")
                cara_lig = docking.ligan_ke_pdbqt(mol, lig_pdbqt)
                st.write(f"Ligan disiapkan dengan {cara_lig}.")

                rec_pdbqt = os.path.join(kerja, "protein.pdbqt")
                cara_rec = docking.protein_ke_pdbqt(path_siap, rec_pdbqt)
                st.write(f"Protein disiapkan dengan {cara_rec}.")

                out_pdbqt = os.path.join(kerja, "hasil_docking.pdbqt")
                hasil["docking"] = docking.jalankan_docking(
                    rec_pdbqt, lig_pdbqt, out_pdbqt,
                    hasil["kotak"]["center"], hasil["kotak"]["size"],
                    exhaustiveness, n_poses, path_vina,
                )
                hasil["pose"] = open(out_pdbqt).read()
                kotak.update(label="Tahap 4 · Docking selesai ✅",
                             state="complete", expanded=False)
        else:
            hasil["docking"] = None

        st.session_state.hasil = hasil

    except Exception as e:
        st.error(f"**Gagal pada salah satu tahap**\n\n{e}")
        with st.expander("Rincian teknis"):
            st.code(traceback.format_exc())

# --------------------------------------------------------------------------
# Menampilkan hasil
# --------------------------------------------------------------------------
h = st.session_state.hasil
if h:
    tab1, tab2, tab3, tab4 = st.tabs(
        ["Ringkasan", "Senyawa & Lipinski", "Hasil docking", "Tampilan 3D"]
    )

    # ---- ringkasan ----
    with tab1:
        info = h.get("info", {})
        a, b, c = st.columns(3)
        a.metric("Protein", info.get("id", "-"))
        b.metric("Resolusi", f"{info.get('resolusi')} Å" if info.get("resolusi") else "-")
        c.metric("Atom protein", h["catatan"]["atom_protein"])
        if info.get("judul"):
            st.write(f"**{info['judul']}**")
        if info.get("metode"):
            st.caption(f"Metode eksperimen: {info['metode']}")

        st.divider()
        k = h["kotak"]
        st.write("**Grid box docking**")
        st.write(f"Titik pusat (x, y, z): `{k['center']}`  ·  "
                 f"Ukuran: `{k['size']}` Å")
        st.caption(f"Ditentukan dari {k['sumber']}, "
                   f"berdasarkan {k['jumlah_atom_acuan']} atom.")

    # ---- senyawa ----
    with tab2:
        s = h["senyawa"]
        st.write(f"### {s['nama']}")
        st.write(f"CID PubChem **{s['cid']}** · Rumus **{s['rumus'] or '-'}**")
        st.code(s["smiles"], language="text")

        lp = h["lipinski"]
        kol = st.columns(4)
        kol[0].metric("Bobot molekul", lp["Bobot molekul"], "< 500")
        kol[1].metric("logP", lp["logP"], "< 5")
        kol[2].metric("Donor H", lp["Donor H"], "≤ 5")
        kol[3].metric("Akseptor H", lp["Akseptor H"], "≤ 10")

        if lp["Lolos"]:
            st.success(f"Lolos aturan Lipinski — {lp['Pelanggaran']} pelanggaran.")
        else:
            st.warning(f"Tidak lolos aturan Lipinski — {lp['Pelanggaran']} pelanggaran.")

    # ---- docking ----
    with tab3:
        d = h.get("docking")
        if not d:
            st.info(
                "Docking dilewati karena AutoDock Vina belum tersedia. "
                "Lihat README untuk cara memasangnya."
            )
        else:
            skor = d["skor"]
            st.metric("Skor terbaik", f"{min(skor):.2f} kcal/mol")
            if min(skor) <= -7:
                st.success("Skor ≤ −7 kcal/mol — tergolong menjanjikan.")
            else:
                st.warning("Skor di atas −7 kcal/mol — ikatannya kurang kuat.")

            st.pyplot(viz.grafik_skor(skor))
            st.dataframe(
                pd.DataFrame({
                    "Pose": [f"Pose {i+1}" for i in range(len(skor))],
                    "Binding affinity (kcal/mol)": skor,
                }),
                use_container_width=True, hide_index=True,
            )
            st.caption(f"Dijalankan memakai {d['mesin']}.")
            st.download_button("⬇️ Unduh pose hasil docking (.pdbqt)",
                               h["pose"], file_name="hasil_docking.pdbqt")

    # ---- 3D ----
    with tab4:
        components.html(
            viz.tampilan_3d(h["pdb_bersih"], h.get("pose", ""),
                            h["kotak"]["center"], h["kotak"]["size"]),
            height=540,
        )
        st.caption("Pita berwarna = protein · batang hijau = senyawa · "
                   "kotak emas = grid box docking. Seret untuk memutar.")
else:
    st.info("Isi panel di sebelah kiri, lalu tekan **Jalankan alur**.")
    