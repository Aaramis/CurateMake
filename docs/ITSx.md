# Pipeline ITSx – Version Factorisée

## Objectif

Identifier et extraire les sous-régions **ITS1**, **5.8S**, et **ITS2** à partir de séquences barcodes (ITS, ITS1, ITS2) grâce à **ITSx**, puis structurer et enregistrer les résultats dans une base de données avec un suivi du statut de traitement.

---

## Commande exécutée

```bash
ITSx --graphical F --preserve T --save_regions all \
     --cpu 30 --multi_thread T --fasta T --summary T \
     -i data/[barcode].fasta \
     -o data/ITS2X
```

---

## Fichiers produits par ITSx

* Séquences complètes : `ITSX.full.fasta`
* Sous-régions : `ITSX.ITS1.fasta`, `ITSX.ITS2.fasta`, `ITSX.5_8S.fasta`
* Fichiers auxiliaires :

  * `ITSX_no_detections.fasta`, `ITSX_no_detections.txt`
  * `ITSX.positions.txt`, `ITSX.problematic.txt`
  * `ITSX.SSU.fasta`, `ITSX.LSU.fasta`
  * `ITSX.summary.txt`

Les détections peuvent être partielles, ambiguës (ITS1 ↔ ITS2), ou inexistantes (notamment pour PLANiTS, CALeDNA).

---

## Headers des séquences

Exemple :

```
>425896|tax_id:3667652|original_id:Ge02721
```

* `425896` → **seq\_id**
* `tax_id`, `original_id` → **métadonnées**

---

# Logique de traitement factorisée

### Étape 1 : Initialisation

* Extraire **toutes les seq\_id** du fichier FASTA d’origine.
* Mettre à jour la base :

  * `barcode_status = "Processed"`

---

### Étape 2 : Traitement générique des sous-régions

* Pour chaque fichier de sous-région (`ITSX.ITS1.fasta`, `ITSX.ITS2.fasta`) :

  1. Charger les séquences.
  2. Valider chaque séquence :

     * `len(seq) > config.min_length`
     * `%N < config.max_n_fraction`
  3. Si valide → insérer une **nouvelle entrée** dans la base :

     * `id_seq = NEW_ID`
     * `source_id = old_seq_id.source_id`
     * `original_id = old_seq_id.original_id`
     * `sequence`, `length`, `%N`
     * `barcode = <region>`
     * `barcode_status = "ITSx"`
     * `parent_seq_id = old_seq_id`

---

### Étape 3 : Reconstruction ITS complet

* Identifier les **id\_seq présents dans ITS1, 5.8S et ITS2**.
* Concaténer les séquences pour former un ITS complet.
* Valider la séquence (mêmes règles).
* Ajouter entrée séquence :

  * `barcode = ITS`
  * `barcode_status = "ITSx"`
  * `parent_seq_id = old_seq_id`

---

### Étape 4 : Séquences non détectées

* Récupérer `ITSX_no_detections.fasta`.
* Mettre à jour la base :

  * `barcode_status = "no_detection"`

---

# Pseudo-code factorisé

```python
def process_ITSx_output(input_fasta, config, db):
    # Étape 1: Initialisation
    seq_ids = get_seq_ids(input_fasta)
    db.update_status(seq_ids, "Processed")

    # Étape 2: Sous-régions génériques
    for region in ["ITS1", "ITS2"]:
        sequences = load_fasta(f"ITSX.{region}.fasta")
        for seq in sequences:
            if validate_sequence(seq, config):
                db.insert_sequence(
                    new_id=generate_id(),
                    source_id=seq.source_id,
                    original_id=seq.original_id,
                    sequence=seq.data,
                    length=len(seq.data),
                    percentage_n=calculate_N_fraction(seq.data),
                    barcode=region,
                    barcode_status="ITSx",
                    parent_seq_id=seq.seq_id
                )

    # Étape 3: Reconstruction ITS complet
    its_full = reconstruct_full_ITS(
        load_fasta("ITSX.ITS1.fasta"),
        load_fasta("ITSX.5_8S.fasta"),
        load_fasta("ITSX.ITS2.fasta")
    )
    for seq in its_full:
        if validate_sequence(seq, config):
            db.insert_sequence(..., barcode="ITS", barcode_status="ITSx")

    # Étape 4: Séquences non détectées
    no_detect = load_fasta("ITSX_no_detections.fasta")
    db.update_status(no_detect.ids, "no_detection")
```